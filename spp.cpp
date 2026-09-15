#include<iostream>
#include<vector>
#include<cmath>
#include<Eigen/Core>
#include<Eigen/Dense>
#include "SPP_Struct.h"

using namespace std;

void ComputeSatPVTAtSignalTrans(EPOCHOBSDATA* Epk, GPSEPHREC* Eph, GPSEPHREC* BDSEph,double UserPos[3])
{
	GPSTIME Ts;
	memset(Epk->SatPVT, 0, sizeof(Epk->SatPVT));
	for(int i=0;i<Epk->SatNum;i++)
	{
		Ts.SecOfWeek = Epk->Time.SecOfWeek - Epk->SatObs[i].P[1] / C_Light;
		Ts.Week = Epk->Time.Week;
		//第一次,计算成功返回true
		Epk->SatPVT[i].Valid = CompSatClkOff(Epk->SatObs[i].Prn, Epk->SatObs[i].System, &Ts, Eph, BDSEph, &Epk->SatPVT[i]);
		if (Epk->SatPVT[i].Valid == false)
		{
			continue;
		}
		Ts.SecOfWeek = Epk->Time.SecOfWeek - Epk->SatObs[i].P[1] / C_Light - Epk->SatPVT[i].SatClkOft;
		//第二次，计算成功返回true
		CompSatClkOff(Epk->SatObs[i].Prn, Epk->SatObs[i].System, &Ts, Eph, BDSEph, &Epk->SatPVT[i]);
		Ts.SecOfWeek = Epk->Time.SecOfWeek - Epk->SatObs[i].P[1] / C_Light - Epk->SatPVT[i].SatClkOft;
		//GPST时间容错处理
		while(Ts.SecOfWeek<0)
		{
			Ts.SecOfWeek += 604800.0;
			Ts.Week--;
		}

		double omega = 0;
		if (Epk->SatObs[i].System == GPS)
		{
			int flag = CompGPSSatPVT(Epk->SatObs[i].Prn, &Ts, Eph, &Epk->SatPVT[i]);
			if (flag < 0)
			{
				Epk->SatPVT[i].Valid = false;
				continue;
			}
			omega = GPS_radv;
		}
		else if (Epk->SatObs[i].System == BDS)
		{
			int flag = CompBDSSatPVT(Epk->SatObs[i].Prn, &Ts, BDSEph, &Epk->SatPVT[i]);
			if (flag < 0)
			{
				Epk->SatPVT[i].Valid = false;
				continue;
			}
			if (BDSEph->SVHealth == 1)     //卫星健康状况，0--正常，1--异常
			{
				Epk->SatPVT[i].Valid = false;
				continue;
			}
			omega = BDS_radv;
		}
		
		double p = sqrt(pow(Epk->SatPVT[i].SatPos[0] - UserPos[0], 2) +pow(Epk->SatPVT[i].SatPos[1] - UserPos[1], 2) +pow(Epk->SatPVT[i].SatPos[2] - UserPos[2], 2));
		double dt = p / C_Light;
		Eigen::Matrix<double, 3, 3> Rz;
		Rz << cos(omega*dt),sin(omega * dt),0,
			-sin(omega * dt), cos(omega * dt), 0,
			0, 0, 1;
		Eigen::Vector3d Xs_;
		Xs_ << Epk->SatPVT[i].SatPos[0], Epk->SatPVT[i].SatPos[1], Epk->SatPVT[i].SatPos[2];
		Eigen::Vector3d Xs = Rz * Xs_;      //地球自转改正
		Epk->SatPVT[i].SatPos[0] = Xs(0);
		Epk->SatPVT[i].SatPos[1] = Xs(1);
		Epk->SatPVT[i].SatPos[2] = Xs(2);

        XYZ Xr;  
        Xr.x = UserPos[0];  
        Xr.y = UserPos[1];  
        Xr.z = UserPos[2];
		CompSatElAz(&Xr, Epk->SatPVT[i].SatPos, &Epk->SatPVT[i].Elevation, &Epk->SatPVT[i].Azimuth, &Epk->SatPVT[i].TropCorr);//计算卫星的仰角、方位角、对流层误差

		Epk->SatPVT[i].Valid = true;    //以上计算均成功，标记为有效
	}
}

bool SPP(EPOCHOBSDATA* Epoch, GPSEPHREC* GPSEph, GPSEPHREC* BDSEph, POSRES* Res)
{
	int i,j, n, dim, Iterator, GPSNum, BDSNum;
	double dPos[3],Range,sigma,pdop;
	//计算B和w矩阵
	Eigen::MatrixXd B(MAXCHANNUM, 5), BTB(5,5),BTW(5,1);
	Eigen::VectorXd w(MAXCHANNUM), x(5), X0(5);

	B.setZero();//线性化系数B矩阵初始化

	//X0 << -2267335.9727,5008647.9932,3222372.6377, -2.561, 0.0;//接收机参考值测试用
	X0 << 1000, 1000, 1000, 0.0, 0.0;

	Iterator = 0; //迭代次数
	do {
		x.setZero();
		ComputeSatPVTAtSignalTrans(Epoch, GPSEph, BDSEph, X0.block(0,0,3,1).data());//解算卫星发射时刻坐标
		n = GPSNum = BDSNum = 0;
		for (i = 0;i < Epoch->SatNum;i++)
		{
			if (Epoch->SatPVT[i].Valid == false || Epoch->SatObs[i].Valid == false)  continue;

			for (j = 0;j < 3;j++) dPos[j] = X0(j, 0) - Epoch->SatPVT[i].SatPos[j];
			Range = sqrt(pow(dPos[0], 2) + pow(dPos[1], 2) + pow(dPos[2], 2));

			B(n, 0) = dPos[0] / Range;
			B(n, 1) = dPos[1] / Range;
			B(n, 2) = dPos[2] / Range;
			w(n) = Epoch->ComObs[i].PIF - (Range - Epoch->SatPVT[i].SatClkOft * C_Light + Epoch->SatPVT[i].TropCorr);
			if (Epoch->SatObs[i].System == GPS) {
				B(n, 3) = 1.0; // GPS系统
				B(n, 4) = 0.0; // BDS系统
				w(n) = w(n) - X0(3,0); // 减去GPS钟差
				GPSNum++;
			}
			else if (Epoch->SatObs[i].System == BDS) {
				B(n, 3) = 0.0; // GPS系统
				B(n, 4) = 1.0; // BDS系统
				GPSEPHREC* BDSeph = BDSEph + Epoch->SatObs[i].Prn - 1;
				w(n) = w(n) - X0(4, 0) - C_Light * BDSeph->TGD1 * FG1_BDS * FG1_BDS / (FG1_BDS * FG1_BDS - FG3_BDS * FG3_BDS); // 减去GPS钟差
				BDSNum++;
			}
			else continue;
			n++;
		}
		if (n < 4) {
			return false;
		}
		BTB = B.block(0, 0, n, 5).transpose() * B.block(0, 0, n, 5);
		BTW = B.block(0, 0, n, 5).transpose() * w.head(n);
		dim = 5;
		if (GPSNum == 0) {

			for (i = 0;i < 5;i++) BTB(i, 3) = BTB(i, 4);
			for (i = 0;i < 5;i++)BTB(3, i) = BTB(4, i);
			BTW(3, 0) = BTW(4, 0); // GPS系统无效
			dim--; // GPS系统无效
		}
		else if (BDSNum == 0) {
			dim--; // BDS系统无效
		}
		else;
		x = BTB.block(0, 0, dim, dim).inverse() * BTW.block(0, 0, dim, 1);
		sigma = pdop = 0;
		for (int i = 0;i < n;i++) {
			sigma += w(i) * w(i); // 残差平方和
		}
		sigma = sqrt(sigma / (n - dim)); // 残差标准差
		pdop = sqrt(BTB.block(0, 0, dim, dim).inverse()(0, 0) + BTB.block(0, 0, dim, dim).inverse()(1, 1) + BTB.block(0, 0, dim, dim).inverse()(2, 2));
		for (i = 0;i < 3;i++) X0(i, 0) += x(i, 0); // 更新用户位置
		if (dim == 5) {
			X0(3, 0) += x(3, 0); // 更新GPS钟差
			X0(4, 0) += x(4, 0); // 更新BDS钟差
		}
		else if (dim == 4) {
			if (GPSNum == 0) X0(4, 0) += x(3, 0); // 更新BDS钟差
			if (BDSNum == 0) X0(3, 0) += x(3, 0); // 更新GPS钟差
		}
		BTB.setZero();
		BTW.setZero();
		w.setZero();
		Iterator++;
	} while (Iterator < 10 && x.block(0,0,3,1).norm() > 1e-4); // 收敛条件
	Res->Pos[0] = X0(0);
	Res->Pos[1] = X0(1);
	Res->Pos[2] = X0(2);
	Res->Sigma = sigma;
	Res->PDOP = pdop;

	return true;
}

bool SPV(EPOCHOBSDATA* Epoch,POSRES* Res)
{
	int i, j, n;
	double dPos[3], Range,Rangedot,SigmaVel;
	//计算B和w矩阵
	Eigen::MatrixXd B(MAXCHANNUM, 4), BTB(4, 4), BTW(4, 1);
	Eigen::VectorXd w(MAXCHANNUM), x(4);
	n = 0;
	for (i = 0;i < Epoch->SatNum;i++)
	{
		if (Epoch->SatPVT[i].Valid == false || Epoch->SatObs->Valid == false)  continue;

		for (j = 0;j < 3;j++) dPos[j] = Epoch->SatPVT[i].SatPos[j] - Res->Pos[j];
		Range = sqrt(pow(dPos[0], 2) + pow(dPos[1], 2) + pow(dPos[2], 2));
		Rangedot = (dPos[0] * Epoch->SatPVT[i].SatVel[0] + dPos[1] * Epoch->SatPVT[i].SatVel[1] + dPos[2] * Epoch->SatPVT[i].SatVel[2]) / Range;

		B(n, 0) = dPos[0] / Range;
		B(n, 1) = dPos[1] / Range;
		B(n, 2) = dPos[2] / Range;
		B(n, 3) = 1;
		w(n) = Epoch->SatObs->D[0] - (Rangedot - C_Light * Epoch->SatPVT[i].SatClkSft);
		n++;
	}
	if (n < 4) {
		return false;
	}
	BTB = B.block(0, 0, n, 4).transpose() * B.block(0, 0, n, 4);
	BTW = B.block(0, 0, n, 4).transpose() * w.head(n);

	SigmaVel = 0;
	for (int i = 0;i < n;i++) {
		SigmaVel += w(i) * w(i); // 残差平方和
	}
	SigmaVel = sqrt(SigmaVel / (n - 4));

	x = BTB.inverse() * BTW;

	Res->Vel[0] = x(0);
	Res->Vel[1] = x(1);
	Res->Vel[2] = x(2);
	Res->SigmaVel = SigmaVel;

	return true;
}

//RTK代码部分
//一份SPP代码，完成基站和流动站定位
//➢ 保存基站和流动站所有卫星的卫星位置、钟差等中间计算结果，用于RTK定位
//➢ 保存基站和流动站的定位结果
//bool RTKSPP(RAWDATA* RawEpoch,POSRES* Base, POSRES* Rov)
//{
//	EPOCHOBSDATA* Epoch;
//	GPSEPHREC* GPSEph = RawEpoch->GpsEph;
//	GPSEPHREC* BDSEph = RawEpoch->BdsEph;
//	POSRES* Res;
//	for(int k=0;k<2;k++)
//	{
//		int i, j, n, dim, Iterator, GPSNum, BDSNum;
//		double dPos[3], Range, sigma, pdop;
//		//计算B和w矩阵
//		Eigen::MatrixXd B(MAXCHANNUM, 5), BTB(5, 5), BTW(5, 1);
//		Eigen::VectorXd w(MAXCHANNUM), x(5), X0(5);
//
//		B.setZero();//线性化系数B矩阵初始化
//
//		X0 << 1000, 1000, 1000, 0.0, 0.0;
//		//选择数据
//		if(k==0){
//			*Epoch = RawEpoch->BaseEpk;//基准站
//			Res = Base;
//		}
//		else if (k == 1) {
//			*Epoch = RawEpoch->RovEpk;//流动站
//			Res = Rov;
//		}
//		else return false;
//
//		Iterator = 0; //迭代次数
//		do {
//			x.setZero();
//			ComputeSatPVTAtSignalTrans(Epoch, GPSEph, BDSEph, X0.block(0, 0, 3, 1).data());//解算卫星发射时刻坐标
//			n = GPSNum = BDSNum = 0;
//			for (i = 0;i < Epoch->SatNum;i++)
//			{
//				if (Epoch->SatPVT[i].Valid == false || Epoch->SatObs[i].Valid == false)  continue;
//
//				for (j = 0;j < 3;j++) dPos[j] = X0(j, 0) - Epoch->SatPVT[i].SatPos[j];
//				Range = sqrt(pow(dPos[0], 2) + pow(dPos[1], 2) + pow(dPos[2], 2));
//
//				B(n, 0) = dPos[0] / Range;
//				B(n, 1) = dPos[1] / Range;
//				B(n, 2) = dPos[2] / Range;
//				w(n) = Epoch->ComObs[i].PIF - (Range - Epoch->SatPVT[i].SatClkOft * C_Light + Epoch->SatPVT[i].TropCorr);
//				if (Epoch->SatObs[i].System == GPS) {
//					B(n, 3) = 1.0; // GPS系统
//					B(n, 4) = 0.0; // BDS系统
//					w(n) = w(n) - X0(3, 0); // 减去GPS钟差
//					GPSNum++;
//				}
//				else if (Epoch->SatObs[i].System == BDS) {
//					B(n, 3) = 0.0; // GPS系统
//					B(n, 4) = 1.0; // BDS系统
//					GPSEPHREC* BDSeph = BDSEph + Epoch->SatObs[i].Prn - 1;
//					w(n) = w(n) - X0(4, 0) - C_Light * BDSeph->TGD1 * FG1_BDS * FG1_BDS / (FG1_BDS * FG1_BDS - FG3_BDS * FG3_BDS); // 减去GPS钟差
//					BDSNum++;
//				}
//				else continue;
//				n++;
//			}
//			if (n < 4) {
//				return false;
//			}
//			BTB = B.block(0, 0, n, 5).transpose() * B.block(0, 0, n, 5);
//			BTW = B.block(0, 0, n, 5).transpose() * w.head(n);
//			dim = 5;
//			if (GPSNum == 0) {
//
//				for (i = 0;i < 5;i++) BTB(i, 3) = BTB(i, 4);
//				for (i = 0;i < 5;i++)BTB(3, i) = BTB(4, i);
//				BTW(3, 0) = BTW(4, 0); // GPS系统无效
//				dim--; // GPS系统无效
//			}
//			else if (BDSNum == 0) {
//				dim--; // BDS系统无效
//			}
//			else;
//			x = BTB.block(0, 0, dim, dim).inverse() * BTW.block(0, 0, dim, 1);
//			sigma = pdop = 0;
//			for (int i = 0;i < n;i++) {
//				sigma += w(i) * w(i); // 残差平方和
//			}
//			sigma = sqrt(sigma / (n - dim)); // 残差标准差
//			pdop = sqrt(BTB.block(0, 0, dim, dim).inverse()(0, 0) + BTB.block(0, 0, dim, dim).inverse()(1, 1) + BTB.block(0, 0, dim, dim).inverse()(2, 2));
//			for (i = 0;i < 3;i++) X0(i, 0) += x(i, 0); // 更新用户位置
//			if (dim == 5) {
//				X0(3, 0) += x(3, 0); // 更新GPS钟差
//				X0(4, 0) += x(4, 0); // 更新BDS钟差
//			}
//			else if (dim == 4) {
//				if (GPSNum == 0) X0(4, 0) += x(3, 0); // 更新BDS钟差
//				if (BDSNum == 0) X0(3, 0) += x(3, 0); // 更新GPS钟差
//			}
//			BTB.setZero();
//			BTW.setZero();
//			w.setZero();
//			Iterator++;
//		} while (Iterator < 10 && x.block(0, 0, 3, 1).norm() > 1e-4); // 收敛条件
//		Res->Pos[0] = X0(0);
//		Res->Pos[1] = X0(1);
//		Res->Pos[2] = X0(2);
//		Res->Sigma = sigma;
//		Res->PDOP = pdop;
//
//		Res->SatNum = n;
//		Res->GPSSatNum = GPSNum;
//		Res->BDSSatNum = BDSNum;
//		Res->IsSuccess = true;
//	}

	//int i, j, n, dim, Iterator, GPSNum, BDSNum;
	//double dPos[3], Range, sigma, pdop;
	////计算B和w矩阵
	//Eigen::MatrixXd B(MAXCHANNUM, 5), BTB(5, 5), BTW(5, 1);
	//Eigen::VectorXd w(MAXCHANNUM), x(5), X0(5);

	//B.setZero();//线性化系数B矩阵初始化

	////X0 << -2267335.9727,5008647.9932,3222372.6377, -2.561, 0.0;//接收机参考值测试用
	//X0 << 1000, 1000, 1000, 0.0, 0.0;

	//Iterator = 0; //迭代次数
	//do {
	//	x.setZero();
	//	ComputeSatPVTAtSignalTrans(Epoch, GPSEph, BDSEph, X0.block(0, 0, 3, 1).data());//解算卫星发射时刻坐标
	//	n = GPSNum = BDSNum = 0;
	//	for (i = 0;i < Epoch->SatNum;i++)
	//	{
	//		if (Epoch->SatPVT[i].Valid == false || Epoch->SatObs[i].Valid == false)  continue;

	//		for (j = 0;j < 3;j++) dPos[j] = X0(j, 0) - Epoch->SatPVT[i].SatPos[j];
	//		Range = sqrt(pow(dPos[0], 2) + pow(dPos[1], 2) + pow(dPos[2], 2));

	//		B(n, 0) = dPos[0] / Range;
	//		B(n, 1) = dPos[1] / Range;
	//		B(n, 2) = dPos[2] / Range;
	//		w(n) = Epoch->ComObs[i].PIF - (Range - Epoch->SatPVT[i].SatClkOft * C_Light + Epoch->SatPVT[i].TropCorr);
	//		if (Epoch->SatObs[i].System == GPS) {
	//			B(n, 3) = 1.0; // GPS系统
	//			B(n, 4) = 0.0; // BDS系统
	//			w(n) = w(n) - X0(3, 0); // 减去GPS钟差
	//			GPSNum++;
	//		}
	//		else if (Epoch->SatObs[i].System == BDS) {
	//			B(n, 3) = 0.0; // GPS系统
	//			B(n, 4) = 1.0; // BDS系统
	//			GPSEPHREC* BDSeph = BDSEph + Epoch->SatObs[i].Prn - 1;
	//			w(n) = w(n) - X0(4, 0) - C_Light * BDSeph->TGD1 * FG1_BDS * FG1_BDS / (FG1_BDS * FG1_BDS - FG3_BDS * FG3_BDS); // 减去GPS钟差
	//			BDSNum++;
	//		}
	//		else continue;
	//		n++;
	//	}
	//	if (n < 4) {
	//		return false;
	//	}
	//	BTB = B.block(0, 0, n, 5).transpose() * B.block(0, 0, n, 5);
	//	BTW = B.block(0, 0, n, 5).transpose() * w.head(n);
	//	dim = 5;
	//	if (GPSNum == 0) {

	//		for (i = 0;i < 5;i++) BTB(i, 3) = BTB(i, 4);
	//		for (i = 0;i < 5;i++)BTB(3, i) = BTB(4, i);
	//		BTW(3, 0) = BTW(4, 0); // GPS系统无效
	//		dim--; // GPS系统无效
	//	}
	//	else if (BDSNum == 0) {
	//		dim--; // BDS系统无效
	//	}
	//	else;
	//	x = BTB.block(0, 0, dim, dim).inverse() * BTW.block(0, 0, dim, 1);
	//	sigma = pdop = 0;
	//	for (int i = 0;i < n;i++) {
	//		sigma += w(i) * w(i); // 残差平方和
	//	}
	//	sigma = sqrt(sigma / (n - dim)); // 残差标准差
	//	pdop = sqrt(BTB.block(0, 0, dim, dim).inverse()(0, 0) + BTB.block(0, 0, dim, dim).inverse()(1, 1) + BTB.block(0, 0, dim, dim).inverse()(2, 2));
	//	for (i = 0;i < 3;i++) X0(i, 0) += x(i, 0); // 更新用户位置
	//	if (dim == 5) {
	//		X0(3, 0) += x(3, 0); // 更新GPS钟差
	//		X0(4, 0) += x(4, 0); // 更新BDS钟差
	//	}
	//	else if (dim == 4) {
	//		if (GPSNum == 0) X0(4, 0) += x(3, 0); // 更新BDS钟差
	//		if (BDSNum == 0) X0(3, 0) += x(3, 0); // 更新GPS钟差
	//	}
	//	BTB.setZero();
	//	BTW.setZero();
	//	w.setZero();
	//	Iterator++;
	//} while (Iterator < 10 && x.block(0, 0, 3, 1).norm() > 1e-4); // 收敛条件
	//Res->Pos[0] = X0(0);
	//Res->Pos[1] = X0(1);
	//Res->Pos[2] = X0(2);
	//Res->Sigma = sigma;
	//Res->PDOP = pdop;

	//return true;
//}
