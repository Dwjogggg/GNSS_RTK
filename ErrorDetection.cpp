#include<iostream>
#include<vector>
#include<cmath>
#include<Eigen/Core>
#include<Eigen/Dense>
#include "SPP_Struct.h"

using namespace std;

double Hopfield(const double H, const double Elev)
{
	double T, P, RH, e, hw, hd, Kw, Kd,trop;
	if (H < 0 || H>3e4)return 0;
	RH = Hop_RH0 * exp(-0.0006396 * (H - Hop_H0));
	P = Hop_p0 * pow((1-0.0000226*(H-Hop_H0)),5.225);
	T = Hop_T0 - 0.0065 * (H - Hop_H0);
	e = RH * exp(-37.2465 + 0.213166 * T - 0.000256908 * T * T);

	hw = 11000;
	hd = 40136 + 148.72 * (Hop_T0-273.16);
	Kw = 155.2 * 1e-7 * (4810 / (T * T)) * e * (hw - H);
	Kd = 155.2 * 1e-7 * (P / T)*(hd - H);

	trop = Kd / sin(sqrt(Elev * Elev + 6.25) * pi / 180.0) + Kw / sin(sqrt(Elev * Elev + 2.25) * pi / 180.0);
	return trop;
}

void DetectOutlier(EPOCHOBSDATA* Obs)
{
	/*一个历元中有多颗卫星，通过循环对其中一颗卫星进行粗差探测。
		函数内部可先定义MWGF CurComObs[MAXCHANNUM];用于存放当前历元的计算结果，
		EPOCHOBSDATA结构体中ComObs保存上个历元的平滑结果。*/
	MWGF CurComObs[MAXCHANNUM];
	int i = 0;

	//对于任意一颗卫星的观测数据，其粗差探测的步骤为
	for (i;i < Obs->SatNum;i++)
	{
		//1. 检查该卫星的双频伪距和相位数据是否有效和完整，若不全或为0，将Valid标记为false，continue
		SATOBSDATA& sat = Obs->SatObs[i];
		//sat.Valid = true;
		for (int j = 0; j < 2; j++) {
			if (sat.P[j] == 0.0 || sat.L[j] == 0.0) {
				sat.Valid = false;
				break;
			}
		}
		if (!sat.Valid) { continue; }
		//2. 计算当前历元该卫星的GF和MW组合值
		double lambda1, lambda2, f1, f2;
		if (sat.System == GPS) {
			lambda1 = WL1_GPS;
			lambda2 = WL2_GPS;
			f1 = FG1_GPS, f2 = FG2_GPS;
		}
		else if (sat.System == BDS) {
			lambda1 = WL1_BDS;
			lambda2 = WL3_BDS;
			f1 = FG1_BDS, f2 = FG3_BDS;
		}
		else {
			sat.Valid = false;
			CurComObs[i].PIF = 0.0;
			continue;
		}
		// GF组合（Geometry-Free）：GF = (L1 * λ1 - L2 * λ2)
		double GF = sat.L[0]  - sat.L[1] ;
		// MW组合（Melbourne-Wübbena）
		double MW =  (sat.L[0]*f1 - sat.L[1]*f2) / (f1 - f2) - (sat.P[0] * f1 + sat.P[1] * f2) / (f1 + f2);
		CurComObs[i].Prn = sat.Prn;
		CurComObs[i].Sys = sat.System;
		CurComObs[i].GF = GF;
		CurComObs[i].MW = MW;
		//CurComObs[i].n =]; // 初始化平滑计数
		//3. 从上个历元的MWGF数据中查找该卫星的GF和MW组合值
		double lastGF = 0.0, lastMW = 0.0;
		int lastN = 0;
		bool foundLast = false;
		for (int k = 0; k < MAXCHANNUM; k++) {
			if (Obs->ComObs[k].Prn == sat.Prn && Obs->ComObs[k].Sys == sat.System) {
				lastGF = Obs->ComObs[k].GF;
				lastMW = Obs->ComObs[k].MW;
				lastN = Obs->ComObs[k].n; // 获取上个历元的平滑计数
				foundLast = true;
				break;
			}
		}
		//4. 计算当前历元该卫星GF与上一历元对应GF的差值dGF
		double dGF = 0.0;
		if (foundLast) {
			dGF = GF - lastGF;
		}
		else {
			dGF = 0.0; // 如果没有上一历元数据，可设为0或特殊值
		}
		//5. 计算当前历元该卫星MW与上一历元对应MW平滑值的差值dMW
		double dMW = 0.0;
		if (foundLast) {
			dMW = MW - lastMW;
		}
		else {
			dMW = 0.0; // 如果没有上一历元数据，可设为0或特殊值
		}
		//6. 检查dGF和dMW是否超限，限差阈值建议为5cm和3m。若超限，标记为粗差，将Valid标记为false ，若不超限，标记为可用将Valid标记为true，并计算该卫星的MW平滑值
		if (fabs(dGF) > GF_THRESHOLD || fabs(dMW) > MW_THRESHOLD)
		{
			sat.Valid = false;
		}
		else
		{
			sat.Valid = true;
			if (foundLast) {//找到了上个历元的值，计算平滑值
				CurComObs[i].MW = (lastMW * lastN + MW) / (lastN + 1);
				CurComObs[i].n = lastN + 1;// 更新平滑计数
			}
			else {//没找到使用当前值
				CurComObs[i].MW = MW;
				CurComObs[i].n = 1;// 初始化平滑计数为1
			}
		}
		//7. 对于可用的观测数据，计算伪距的IF组合观测值，用于SPP
		if (sat.Valid) 
		{
			double alpha1, alpha2;
			if (sat.System == GPS) {
				alpha1 = FG1_GPS * FG1_GPS;
				alpha2 = FG2_GPS * FG2_GPS;
			}
			else if (sat.System == BDS) {
				alpha1 = FG1_BDS * FG1_BDS;
				alpha2 = FG3_BDS * FG3_BDS;
			}
			else {
				alpha1 = alpha2 = 0.0;
			}
			if (alpha1 != alpha2 && alpha1 > 0) 
			{
				// IF组合：P_IF = (alpha1 * P1 - alpha2 * P2) / (alpha1 - alpha2)
				double PIF = (alpha1 * sat.P[0] - alpha2 * sat.P[1]) / (alpha1 - alpha2);
				CurComObs[i].PIF = PIF;
			}
			else 
			{
				CurComObs[i].PIF = 0.0;
			}
		}
	    //8. 所有卫星循环计算完成之后，将CurComObs内存拷贝到ComObs，即函数运行结束后， ComObs保存了当前历元的GF和MW平滑值。
		for (int i = 0; i < MAXCHANNUM; i++) {
			Obs->ComObs[i] = CurComObs[i];
		}
	}
}

void SDDetectOutlier(SDEPOCHOBS* Obs)
{
	/*一个历元中有多颗卫星，通过循环对其中一颗卫星进行粗差探测。
		函数内部可先定义MWGF CurComObs[MAXCHANNUM];用于存放当前历元的计算结果，
		EPOCHOBSDATA结构体中ComObs保存上个历元的平滑结果。*/
	MWGF CurComObs[MAXCHANNUM];
	int i = 0;

	//对于任意一颗卫星的观测数据，其粗差探测的步骤为
	for (i;i < Obs->SatNum;i++)
	{
		//1. 检查该卫星的双频伪距和相位数据是否有效和完整，若不全或为0，将Valid标记为false，continue
		SDSATOBS& sat = Obs->SdSatObs[i];
		sat.Valid = true;
		for (int j = 0; j < 2; j++) {
			if (sat.dP[j] == 0.0 || sat.dL[j] == 0.0) {
				sat.Valid = false;
				break;
			}
		}
		if (!sat.Valid) { continue; }
		//2. 计算当前历元该卫星的GF和MW组合值
		double lambda1, lambda2, f1, f2;
		if (sat.System == GPS) {
			lambda1 = WL1_GPS;
			lambda2 = WL2_GPS;
			f1 = FG1_GPS, f2 = FG2_GPS;
		}
		else if (sat.System == BDS) {
			lambda1 = WL1_BDS;
			lambda2 = WL3_BDS;
			f1 = FG1_BDS, f2 = FG3_BDS;
		}
		else {
			sat.Valid = false;
			CurComObs[i].PIF = 0.0;
			continue;
		}
		// GF组合（Geometry-Free）：GF = (L1 * λ1 - L2 * λ2)
		double GF = sat.dL[0] - sat.dL[1];
		// MW组合（Melbourne-Wübbena）
		double MW = (sat.dL[0] * f1 - sat.dL[1] * f2) / (f1 - f2) - (sat.dP[0] * f1 + sat.dP[1] * f2) / (f1 + f2);
		CurComObs[i].Prn = sat.Prn;
		CurComObs[i].Sys = sat.System;
		CurComObs[i].GF = GF;
		CurComObs[i].MW = MW;
		//CurComObs[i].n =]; // 初始化平滑计数
		//3. 从上个历元的MWGF数据中查找该卫星的GF和MW组合值
		double lastGF = 0.0, lastMW = 0.0;
		int lastN = 0;
		bool foundLast = false;
		for (int k = 0; k < MAXCHANNUM; k++) {
			if (Obs->SdCObs[k].Prn == sat.Prn && Obs->SdCObs[k].Sys == sat.System) {
				lastGF = Obs->SdCObs[k].GF;
				lastMW = Obs->SdCObs[k].MW;
				lastN = Obs->SdCObs[k].n; // 获取上个历元的平滑计数
				foundLast = true;
				break;
			}
		}
		//4. 计算当前历元该卫星GF与上一历元对应GF的差值dGF
		double dGF = 0.0;
		if (foundLast) {
			dGF = GF - lastGF;
		}
		else {
			dGF = 0.0; // 如果没有上一历元数据，可设为0或特殊值
		}
		//5. 计算当前历元该卫星MW与上一历元对应MW平滑值的差值dMW
		double dMW = 0.0;
		if (foundLast) {
			dMW = MW - lastMW;
		}
		else {
			dMW = 0.0; // 如果没有上一历元数据，可设为0或特殊值
		}
		//6. 检查dGF和dMW是否超限，限差阈值建议为5cm和3m。若超限，标记为粗差，将Valid标记为false ，若不超限，标记为可用将Valid标记为true，并计算该卫星的MW平滑值
		if (fabs(dGF) > GF_THRESHOLD || fabs(dMW) > MW_THRESHOLD)
		{
			sat.Valid = false;
		}
		else
		{
			sat.Valid = true;
			if (foundLast) {//找到了上个历元的值，计算平滑值
				CurComObs[i].MW = (lastMW * lastN + MW) / (lastN + 1);
				CurComObs[i].n = lastN + 1;// 更新平滑计数
			}
			else {//没找到使用当前值
				CurComObs[i].MW = MW;
				CurComObs[i].n = 1;// 初始化平滑计数为1
			}
		}
		//7. 对于可用的观测数据，计算伪距的IF组合观测值，用于SPP
		if (sat.Valid)
		{
			double alpha1, alpha2;
			if (sat.System == GPS) {
				alpha1 = FG1_GPS * FG1_GPS;
				alpha2 = FG2_GPS * FG2_GPS;
			}
			else if (sat.System == BDS) {
				alpha1 = FG1_BDS * FG1_BDS;
				alpha2 = FG3_BDS * FG3_BDS;
			}
			else {
				alpha1 = alpha2 = 0.0;
			}
			if (alpha1 != alpha2 && alpha1 > 0)
			{
				// IF组合：P_IF = (alpha1 * P1 - alpha2 * P2) / (alpha1 - alpha2)
				double PIF = (alpha1 * sat.dP[0] - alpha2 * sat.dP[1]) / (alpha1 - alpha2);
				CurComObs[i].PIF = PIF;
			}
			else
			{
				CurComObs[i].PIF = 0.0;
			}
		}
		//8. 所有卫星循环计算完成之后，将CurComObs内存拷贝到ComObs，即函数运行结束后， ComObs保存了当前历元的GF和MW平滑值。
		for (int i = 0; i < MAXCHANNUM; i++) {
			Obs->SdCObs[i] = CurComObs[i];
		}
	}
}