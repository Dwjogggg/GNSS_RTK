#include<iostream>
#include<vector>
#include<cmath>
#include<Eigen/Core>
#include<Eigen/Dense>
#include "SPP_Struct.h"

using namespace std;

bool CompSatClkOff(const int Prn, const GNSSys Sys, const GPSTIME* t, GPSEPHREC* GPSEph, GPSEPHREC* BDSEph, SATPVT* Mid)
{
	double dt, LimT = 7500.0;
	GPSTIME CurT;
	GPSEPHREC* EPH;

	CurT = *t;
	if (Sys == GPS) EPH = GPSEph + Prn - 1;
	else if (Sys == BDS) {
		EPH = BDSEph + Prn - 1;
		CurT.Week -= 1356;
		CurT.SecOfWeek -= 14;
		LimT = 3900.0;
	}
	else return false;

	if (Prn != EPH->PRN || EPH->System != Sys)
	{
		return false;
	}

	dt = (CurT.Week - EPH->TOC.Week) * 604800.0 + CurT.SecOfWeek - EPH->TOC.SecOfWeek;
	if (fabs(dt) > LimT || EPH->SVHealth != 0) return false;

	//计算卫星钟差
	Mid->SatClkOft = EPH->ClkBias + EPH->ClkDrift * dt + EPH->ClkDriftRate * pow(dt, 2);
	Mid->SatClkSft = EPH->ClkDrift + 2.0 * EPH->ClkDriftRate * dt;
	Mid->Valid = true;
	return true;
}

int CompGPSSatPVT( int Prn, const GPSTIME* t, const GPSEPHREC* Eph, SATPVT* Mid)
{
	Prn--;
	if (Eph[Prn].PRN < 1 || Eph[Prn].PRN > 32 || Eph[Prn].PRN-1 != Prn || Eph[Prn].System != GPS)
	{
		return -1; //卫星号不合法
	}
	//卫星位置与钟差计算
	double t0 = (t->Week- Eph[Prn].TOE.Week)*604800.0 + (t->SecOfWeek - Eph[Prn].TOE.SecOfWeek);
	double A = Eph[Prn].SqrtA * Eph[Prn].SqrtA;
	double n0 = sqrt(GPS_GM / (A * A * A));
	double n = n0 + Eph[Prn].DetlaN;
	double M = Eph[Prn].M0 + n * t0;
	
	double E = M;
	bool flag = true;
	while (flag)
	{
		double E1 = M + Eph[Prn].e * sin(E);
		if (abs(E1 - E) < 1e-12)
		{
			flag = false;
		}
		E = E1;
	}

	double v = atan2(sqrt(1.0 - Eph[Prn].e * Eph[Prn].e) * sin(E), cos(E) - Eph[Prn].e);
	double phi = v + Eph[Prn].omega;

	double u = phi + Eph[Prn].Cuc * cos(2.0 * phi) + Eph[Prn].Cus * sin(2.0 * phi);
	double r = A * (1.0 - Eph[Prn].e * cos(E)) + Eph[Prn].Crc * cos(2.0 * phi) + Eph[Prn].Crs * sin(2.0 * phi);
	double i = Eph[Prn].i0 + Eph[Prn].iDOT * t0 + Eph[Prn].Cic * cos(2.0 * phi) + Eph[Prn].Cis * sin(2.0 * phi);
	double OMEGA = Eph[Prn].OMEGA0 + (Eph[Prn].OMEGADot - GPS_radv) * t0 - GPS_radv * Eph[Prn].TOE.SecOfWeek;

	//卫星位置
	Mid->SatPos[0] = r*cos(u)*cos(OMEGA) - r*sin(u)*cos(i)*sin(OMEGA);
	Mid->SatPos[1] = r*cos(u)*sin(OMEGA) + r*sin(u)*cos(i)*cos(OMEGA);
	Mid->SatPos[2] = r*sin(u)*sin(i);
	
	//计算卫星速度
	double dE = n / (1.0 - Eph[Prn].e * cos(E));
	double dphi = (sqrt(1.0 - Eph[Prn].e * Eph[Prn].e) * dE) / (1.0 - Eph[Prn].e * cos(E));

	double du = 2 * (Eph[Prn].Cus*cos(2*phi)- Eph[Prn].Cuc*sin(2*phi))*dphi+dphi;
	double dr = Eph[Prn].e * A * sin(E) * dE + 2 * (Eph[Prn].Crs * cos(2 * phi) - Eph[Prn].Crc * sin(2 * phi)) * dphi;
	double di = Eph[Prn].iDOT + 2 * (Eph[Prn].Cis * cos(2 * phi) - Eph[Prn].Cic * sin(2 * phi)) * dphi;
	double dOMEGA = Eph[Prn].OMEGADot - GPS_radv;

	double xk_ = r * cos(u);
	double yk_ = r * sin(u);
	double xkdot = dr * cos(u) - r * du * sin(u);
	double ykdot = dr * sin(u) + r * du * cos(u);

	Eigen::Matrix<double,3,4> R;
	R << cos(OMEGA),-sin(OMEGA)*cos(i),-(xk_*sin(OMEGA)+yk_*cos(OMEGA)*cos(i)),yk_*sin(OMEGA)*sin(i),
		 sin(OMEGA),cos(OMEGA)*cos(i),xk_*cos(OMEGA)-yk_*sin(OMEGA)*cos(i),-yk_*cos(OMEGA)*sin(i),
		 0,sin(i),0,yk_*cos(i);
	Eigen::Vector4d V;
	V << xkdot,ykdot,dOMEGA,di;

	//卫星速度
	Eigen::Vector3d vel = R * V;
	Mid->SatVel[0] = vel(0);
	Mid->SatVel[1] = vel(1);
	Mid->SatVel[2] = vel(2);
	//相对论改正
	double dtr = F_RELATIV * Eph[Prn].e * Eph[Prn].SqrtA * sin(E);
	Mid->SatClkOft = Eph[Prn].ClkBias + Eph[Prn].ClkDrift * t0 + Eph[Prn].ClkDriftRate * pow(t0, 2) + dtr;
	double ddtr = F_RELATIV * Eph[Prn].e * Eph[Prn].SqrtA * cos(E) * dE;
	Mid->SatClkSft = Eph[Prn].ClkDrift + 2.0 * Eph[Prn].ClkDriftRate * t0 + ddtr;

	return 0;
}

int CompBDSSatPVT( int Prn, const GPSTIME* t, const GPSEPHREC* Eph, SATPVT* Mid)
{
	Prn--;
	if (Eph[Prn].PRN < 1 || Eph[Prn].PRN > 63 || Eph[Prn].PRN-1 != Prn || Eph[Prn].System != BDS)
	{
		return -1; //卫星号不合法
	}
	//卫星位置与钟差计算
	double t0 = (t->Week - 1356.0 - Eph[Prn].TOE.Week) * 604800.0 + (t->SecOfWeek - 14.0 - Eph[Prn].TOE.SecOfWeek);
	double A = Eph[Prn].SqrtA * Eph[Prn].SqrtA;
	double n0 = sqrt(BDS_GM / (A * A * A));
	double n = n0 + Eph[Prn].DetlaN;
	double M = Eph[Prn].M0 + n * t0;

	double E = M;
	bool flag = true;
	while (flag)
	{
		double E1 = M + Eph[Prn].e * sin(E);
		if (abs(E1 - E) < 1e-12)
		{
			flag = false;
		}
		E = E1;
	}

	double v = atan2(sqrt(1.0 - Eph[Prn].e * Eph[Prn].e) * sin(E), cos(E) - Eph[Prn].e);
	double phi = v + Eph[Prn].omega;

	double u = phi + Eph[Prn].Cuc * cos(2.0 * phi) + Eph[Prn].Cus * sin(2.0 * phi);
	double r = A * (1.0 - Eph[Prn].e * cos(E)) + Eph[Prn].Crc * cos(2.0 * phi) + Eph[Prn].Crs * sin(2.0 * phi);
	double i = Eph[Prn].i0 + Eph[Prn].iDOT * t0 + Eph[Prn].Cic * cos(2.0 * phi) + Eph[Prn].Cis * sin(2.0 * phi);

	double xk = r * cos(u);
	double yk = r * sin(u);

	//GEO卫星
	if((1<= Eph[Prn].PRN && Eph[Prn].PRN <=5)||(59<= Eph[Prn].PRN && Eph[Prn].PRN <=63))
	{
		double OMEGA = Eph[Prn].OMEGA0 + Eph[Prn].OMEGADot * t0 - BDS_radv * Eph[Prn].TOE.SecOfWeek;
		double Xgk = xk * cos(OMEGA) - yk * sin(OMEGA) * cos(i);
		double Ygk = xk * sin(OMEGA) + yk * cos(OMEGA) * cos(i);
		double Zgk = yk * sin(i);

		double sita = (-5 / 180.0) * pi;
		Eigen::Matrix3d Rx;
		Rx << 1,0,0,
			0, cos(sita), sin(sita),
			0, -sin(sita),cos(sita);
		double pika = BDS_radv * t0;
		Eigen::Matrix3d Rz;
		Rz << cos(pika),sin(pika), 0,
			-sin(pika), cos(pika), 0,
			0, 0, 1;
		Eigen::Vector3d Gk;
		Gk << Xgk, Ygk, Zgk;
		Eigen::Vector3d Gk1 = Rz * Rx * Gk;

		Mid->SatPos[0] = Gk1(0);
		Mid->SatPos[1] = Gk1(1);
		Mid->SatPos[2] = Gk1(2);

		double dE = n / (1.0 - Eph[Prn].e * cos(E));
		double dphi = (sqrt(1.0 - Eph[Prn].e * Eph[Prn].e) * dE) / (1.0 - Eph[Prn].e * cos(E));

		double du = 2 * (Eph[Prn].Cus * cos(2 * phi) - Eph[Prn].Cuc * sin(2 * phi)) * dphi + dphi;
		double dr = Eph[Prn].e * A * sin(E) * dE + 2 * (Eph[Prn].Crs * cos(2 * phi) - Eph[Prn].Crc * sin(2 * phi)) * dphi;
		double di = Eph[Prn].iDOT + 2 * (Eph[Prn].Cis * cos(2 * phi) - Eph[Prn].Cic * sin(2 * phi)) * dphi;
		double dOMEGA = Eph[Prn].OMEGADot;

		double xkdot = dr * cos(u) - r * du * sin(u);
		double ykdot = dr * sin(u) + r * du * cos(u);

		Eigen::Vector3d Vel;
		Vel << xkdot * cos(OMEGA) - xk * sin(OMEGA) * dOMEGA - ykdot * cos(i) * sin(OMEGA) + yk * sin(i) * di * sin(OMEGA),
			xkdot* sin(OMEGA) + xk * cos(OMEGA) * dOMEGA + ykdot * cos(i) * cos(OMEGA) - yk * sin(i) * di * cos(OMEGA) - yk * cos(i) * sin(OMEGA) * dOMEGA,
			ykdot* sin(i) + yk * cos(i) * di;
		Eigen::Matrix3d dRz;
		dRz << cos(BDS_radv), sin(BDS_radv), 0,
			-sin(BDS_radv), cos(BDS_radv), 0,
			0, 0, 1;
		Eigen::Vector3d GEO_Vel = Rz * Rx * Vel + dRz * Rx * Gk;
		Mid->SatVel[0] = GEO_Vel(0);
		Mid->SatVel[1] = GEO_Vel(1);
		Mid->SatVel[2] = GEO_Vel(2);

		double dtr = F_RELATIV * Eph[Prn].e * Eph[Prn].SqrtA * sin(E);
		Mid->SatClkOft = Eph[Prn].ClkBias + Eph[Prn].ClkDrift * t0 + Eph[Prn].ClkDriftRate * pow(t0, 2) + dtr;
		double ddtr = F_RELATIV * Eph[Prn].e * Eph[Prn].SqrtA * cos(E) * dE;
		Mid->SatClkSft = Eph[Prn].ClkDrift + 2.0 * Eph[Prn].ClkDriftRate * t0 + ddtr;

		return 0;
	}
	//MEO/IGSO卫星
	else
	{
		double OMEGA = Eph[Prn].OMEGA0 + (Eph[Prn].OMEGADot - BDS_radv) * t0 - BDS_radv * Eph[Prn].TOE.SecOfWeek;

		//卫星位置
		Mid->SatPos[0] = r * cos(u) * cos(OMEGA) - r * sin(u) * cos(i) * sin(OMEGA);
		Mid->SatPos[1] = r * cos(u) * sin(OMEGA) + r * sin(u) * cos(i) * cos(OMEGA);
		Mid->SatPos[2] = r * sin(u) * sin(i);

		//计算卫星速度
		double dE = n / (1.0 - Eph[Prn].e * cos(E));
		double dphi = (sqrt(1.0 - Eph[Prn].e * Eph[Prn].e) * dE) / (1.0 - Eph[Prn].e * cos(E));

		double du = 2 * (Eph[Prn].Cus * cos(2 * phi) - Eph[Prn].Cuc * sin(2 * phi)) * dphi + dphi;
		double dr = Eph[Prn].e * A * sin(E) * dE + 2 * (Eph[Prn].Crs * cos(2 * phi) - Eph[Prn].Crc * sin(2 * phi)) * dphi;
		double di = Eph[Prn].iDOT + 2 * (Eph[Prn].Cis * cos(2 * phi) - Eph[Prn].Cic * sin(2 * phi)) * dphi;
		double dOMEGA = Eph[Prn].OMEGADot - BDS_radv;

		double xk_ = r * cos(u);
		double yk_ = r * sin(u);
		double xkdot = dr * cos(u) - r * du * sin(u);
		double ykdot = dr * sin(u) + r * du * cos(u);

		Eigen::Matrix<double, 3, 4> R;
		R << cos(OMEGA), -sin(OMEGA) * cos(i), -(xk_ * sin(OMEGA) + yk_ * cos(OMEGA) * cos(i)), yk_* sin(OMEGA)* sin(i),
			sin(OMEGA), cos(OMEGA)* cos(i), xk_* cos(OMEGA) - yk_ * sin(OMEGA) * cos(i), -yk_ * cos(OMEGA) * sin(i),
			0, sin(i), 0, yk_* cos(i);
		Eigen::Vector4d V;
		V << xkdot, ykdot, dOMEGA, di;

		//卫星速度
		Eigen::Vector3d vel = R * V;
		Mid->SatVel[0] = vel(0);
		Mid->SatVel[1] = vel(1);
		Mid->SatVel[2] = vel(2);
		//相对论改正
		double dtr = F_RELATIV * Eph[Prn].e * Eph[Prn].SqrtA * sin(E);
		Mid->SatClkOft = Eph[Prn].ClkBias + Eph[Prn].ClkDrift * t0 + Eph[Prn].ClkDriftRate * pow(t0, 2) + dtr;
		double ddtr = F_RELATIV * Eph[Prn].e * Eph[Prn].SqrtA * cos(E) * dE;
		Mid->SatClkSft = Eph[Prn].ClkDrift + 2.0 * Eph[Prn].ClkDriftRate * t0 + ddtr;

		return 0;
	}
}