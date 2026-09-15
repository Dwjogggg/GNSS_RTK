#include<iostream>
#include<vector>
#include<cmath>
#include<Eigen/Core>
#include<Eigen/Dense>
#include "SPP_Struct.h"

using namespace std;

int GetSynObs(FILE* FBas, FILE* FRov,RAWDATA* Raw)
{
	double dt;
	int lenR_base, lenR_rover;        /*已接收数据长度，实际接收长度*/
	static int lenD_base = 0, lenD_rover = 0;  /*剩余长度*/
	static unsigned char buff_base[MAXRAWLEN], buff_rover[MAXRAWLEN];      /*最大缓冲区*/
	//1.获取流动站的观测数据，若不成功，继续获取，若文件结束，返回文件结束标记，若成功，得到观测时刻。
	while (!feof(FRov))
	{
		if ((lenR_rover = fread(buff_rover + lenD_rover, sizeof(unsigned char), MAXRAWLEN - lenD_rover, FRov)) < MAXRAWLEN - lenD_rover) return -1;
		lenD_rover += lenR_rover;
		if (DecodeNovOem7Dat(buff_rover, lenD_rover, &Raw->RovEpk, Raw->GpsEph, Raw->BdsEph, &Raw->RovEpk.Res) == 1) break;
	}
	//2.流动站观测时刻与当前基站观测时刻比较，在限差范围内，文件同步成功，返回1；
	dt = (Raw->RovEpk.Time.Week - Raw->BaseEpk.Time.Week) * 604800 + Raw->RovEpk.Time.SecOfWeek- Raw->BaseEpk.Time.SecOfWeek;
	if (fabs(dt) < 0.01) return 1;
	//3.若不在限差范围内，获取基站观测数据，两站的观测时间求差，在限差范围内，同步成功返回1；
	//4.如果不在限差范围内，若基站时间在后，返回0；若基站时间在前，循环获取基站数据，直到成功。
	while (!feof(FBas))
	{
		if ((lenR_base = fread(buff_base + lenD_base, sizeof(unsigned char), MAXRAWLEN - lenD_base, FBas)) < MAXRAWLEN - lenD_base) return -1;
		lenD_base += lenR_base;
		if (DecodeNovOem7Dat(buff_base, lenD_base, &Raw->BaseEpk, Raw->GpsEph, Raw->BdsEph, &Raw->BaseEpk.Res) == 1)
		{
			dt = (Raw->RovEpk.Time.Week - Raw->BaseEpk.Time.Week) * 604800 + Raw->RovEpk.Time.SecOfWeek - Raw->BaseEpk.Time.SecOfWeek;
			if (fabs(dt) < 0.01) return 1;
			else if(dt < 0.5) return 0;
			else;
		}
	}
}

void SDEpochObs(const EPOCHOBSDATA* RovEpk, const EPOCHOBSDATA* BaseEpk, SDEPOCHOBS* SDObs) 
{
	int i, j, SDSatNum = 0;
	// 1、根据流动站卫星数进行循环，对于每颗流动站的卫星
	for(i=0;i<RovEpk->SatNum;i++)
	{
		const SATOBSDATA& rovSat = RovEpk->SatObs[i];
		// 2、检查卫星号和系统号是否正常，卫星观测值是否完整，如果不正常，continue；
		if (rovSat.System != GPS && rovSat.System != BDS) continue;
		if(rovSat.System == GPS&&(rovSat.Prn<1||rovSat.Prn>32)
			|| rovSat.System == BDS && (rovSat.Prn < 1 || rovSat.Prn>63)) continue;
		if (rovSat.P[0] == 0.0 || rovSat.P[1] == 0.0 || rovSat.L[0] == 0.0 || rovSat.L[1] == 0.0) continue;
		// 3、在基站观测值中查找与该流动站相同的卫星，如果未找到，continue；
		for (j = 0;j < BaseEpk->SatNum;j++)
		{
			const SATOBSDATA& baseSat = BaseEpk->SatObs[j];
			if (rovSat.System == baseSat.System && rovSat.Prn == baseSat.Prn)
			{
				if (rovSat.Valid == false || baseSat.Valid == false) break;
				// 4、对同类型和同频率的观测值求差并保存，保存索引号、卫星号和系统号等相关信息，累加单差观测值卫星数量
				SDSATOBS& sdSat = SDObs->SdSatObs[SDSatNum];
				SDSatNum++;
				sdSat.Prn = rovSat.Prn;
				sdSat.System = rovSat.System;
				sdSat.nBas = j;
				sdSat.nRov = i;
				sdSat.dP[0] = rovSat.P[0] - baseSat.P[0];
				sdSat.dP[1] = rovSat.P[1] - baseSat.P[1];
				sdSat.dL[0] = rovSat.L[0] - baseSat.L[0];
				sdSat.dL[1] = rovSat.L[1] - baseSat.L[1];
				sdSat.Valid = true;
				break;
			}
		}
	}
	// 5、循环结束后，对单差观测值的观测时刻和卫星数量赋值。
	SDObs->Time = RovEpk->Time;
	SDObs->SatNum = SDSatNum;
	cout << "SDSatNum: " << SDSatNum << endl;
}

void DetRefSat(const EPOCHOBSDATA* epkA/*Base*/, const EPOCHOBSDATA* epkB/*Rov*/, SDEPOCHOBS* SDObs, DDCOBS* DDObs)
{
	/*每个卫星导航系统各选取一颗卫星作为参考星
	要求：
		伪距和载波相位通过周跳探测，没有粗差和周跳标记
		卫星星历正常，卫星位置计算成功
		高度角大或CN0大
		连续观测时间大于一定时间（6s），没有半周
		上一个历元是否为参考星（KF使用）
	*/
	int i, j, n, BDSNum = 0, GPSNum = 0;
	double Sum[2] = { 0.0 }, MaxSum[2] = { 0.0 };
	int RefPrn[2] = { -1 }, RefIndex[2] = { -1 };

	/*基准星选取*/
	for (int i = 0; i < SDObs->SatNum; i++)
	{
		if (!SDObs->SdSatObs[i].Valid || !epkA->SatPVT[SDObs->SdSatObs[i].nBas].Valid || !epkB->SatPVT[SDObs->SdSatObs[i].nRov].Valid) continue;
		if (epkA->SatObs[SDObs->SdSatObs[i].nBas].LockTime[0] < 6 || epkA->SatObs[SDObs->SdSatObs[i].nBas].LockTime[1] < 6 
		 || epkB->SatObs[SDObs->SdSatObs[i].nRov].LockTime[0] < 6 || epkB->SatObs[SDObs->SdSatObs[i].nRov].LockTime[1] < 6) continue;

		n = SDObs->SdSatObs[i].System == GPS ? 0 : 1;
		if (n == 0) GPSNum++;
		else if (n == 1) BDSNum++;
		else continue;
		Sum[n] = epkB->SatPVT[SDObs->SdSatObs[i].nRov].Elevation * 180 / pi + epkA->SatPVT[SDObs->SdSatObs[i].nBas].Elevation * 180 / pi
		+ epkA->SatObs[SDObs->SdSatObs[i].nBas].cn0[0] + epkA->SatObs[SDObs->SdSatObs[i].nBas].cn0[1]
		+ epkB->SatObs[SDObs->SdSatObs[i].nRov].cn0[0] + epkB->SatObs[SDObs->SdSatObs[i].nRov].cn0[1];

		if (Sum[n] > MaxSum[n]) 
		{
			MaxSum[n] = Sum[n];
			RefPrn[n] = SDObs->SdSatObs[i].Prn;
			RefIndex[n] = i;
		}
	}

	for (i = 0; i < 2; i++) {
		if (MaxSum[i] < 200.0) DDObs->RefPos[i] = -1;
		else DDObs->RefPos[i] = RefIndex[i];
	}
	DDObs->RefPrn[0] = RefPrn[0];   // 参考星卫星号与存储位置，0=GPS; 1=BDS
	DDObs->RefPrn[1] = RefPrn[1];
	// 存储双差卫星数
	DDObs->Sats = GPSNum + BDSNum;
	DDObs->DDSatNum[0] = GPSNum;
	DDObs->DDSatNum[1] = BDSNum;
	cout << "GPSNum: " << GPSNum << " BDSNum: " << BDSNum << endl;
	cout << "RefGPS: " << DDObs->RefPrn[0] << " RefBDS: " << DDObs->RefPrn[1] << endl;
}
// RTK浮点解，使用LS
bool RTKFloat(RAWDATA* Raw, POSRES* Base, POSRES* Rov)
{
	/*1. 设置基站和流动站位置初值
		2. 计算GPS和BDS双差卫星数（单双频混用：各系统每个频率的双差卫
		星数）
		3. 计算基站坐标到所有卫星的几何距离
		4. 计算流动站到参考星的几何距离
		5. 对单差观测值进行循环
		① 参考星不用计算，存在半周不参与计算
		② 线性化双差观测方程得到B矩阵和W向量
		③ 计算权矩阵
		6. 双差观测方程数量大于未知数，则求解
		7. 最小二乘解算
		8. 更新流动站的位置和双差模糊度参数
		9. 流动站位置增量大于阈值或迭代次数小于阈值，返回4迭代计算，否则
		浮点解计算完成。
		10. 保存双差浮点解模糊度及其协因数矩阵，用于LAMBDA模糊度固定
		11. 精度评价*/

	// 1. 设置基站和流动站位置初值
	int GPSnum, BDSnum,i,j;
	double Bxyz[3] = { 0.0 }, Rover[3] = { 0.0 }, dPos[3] = { 0.0 }, Brou[MAXCHANNUM] = {0.0};

	// 2. 计算GPS和BDS双差卫星数（单双频混用：各系统每个频率的双差卫星数）
	for(i=0;i<Raw->SDObs.SatNum;i++)
	{
		if (Raw->SDObs.SdSatObs[i].Valid == true
			&& Raw->SDObs.SdSatObs[i].System == GPS) {
			for (j = 0;j < 3;j++) {
				dPos[j] = Bxyz[j] - Raw->RovEpk.SatPVT[Raw->SDObs.SdSatObs[i].nRov].SatPos[j];
			}
			Brou[i] = sqrt(pow(dPos[0], 2) + pow(dPos[1], 2) + pow(dPos[2], 2));
			GPSnum++;
		}
		if (Raw->SDObs.SdSatObs[i].Valid == true
			&& Raw->SDObs.SdSatObs[i].System == BDS) {
			for (j = 0;j < 3;j++) {
				dPos[j] = Bxyz[j] - Raw->RovEpk.SatPVT[Raw->SDObs.SdSatObs[i].nRov].SatPos[j];
			}
			
			BDSnum++;
		}
		else continue;
	}
	// 3. 计算基站坐标到所有卫星的几何距离

	// 循环计算
	// 4. 计算流动站到参考星的几何距离

	// 5. 对单差观测值进行循环

	// 6. 双差观测方程数量大于未知数，则求解

	// 7. 最小二乘解算

	// 8. 更新流动站的位置和双差模糊度参数

	// 9. 流动站位置增量大于阈值或迭代次数小于阈值，返回4迭代计算，否则浮点解计算完成。

	// 10. 保存双差浮点解模糊度及其协因数矩阵，用于LAMBDA模糊度固定

	// 11. 精度评价
}

