#include<iostream>
#include<vector>
#include<cmath>
#include "SPP_Struct.h"

using namespace std;

void BLHToXYZ(const GEOCOOR* blh, XYZ* xyz, const double R, const double F)
{
	if (blh == nullptr || xyz == nullptr)
	{
		cout << "Error: BLHToXYZ: Null pointer passed." << endl;
		return;
	}

	double e2 = 2 * F - F * F;
	double N = R / sqrt(1 - e2 * sin(blh -> latitude) * sin(blh->latitude));
	xyz->x = (N + blh->height) * cos(blh->latitude) * cos(blh->longtitude);
	xyz->y = (N + blh->height) * cos(blh->latitude) * sin(blh->longtitude);
	xyz->z = (N * (1 - e2) + blh->height) * sin(blh->latitude);
}

void XYZToBLH(const XYZ* xyz, GEOCOOR* blh, const double R, const double F)
{
	if (xyz == nullptr || blh == nullptr)
	{
		cout << "Error: XYZToBLH: Null pointer passed." << endl;
		return;
	}

	double e2 = 2 * F - F * F;

	int iter = 0;
	double dz,delta;
	delta = 0;
	dz = e2 * xyz->z;
	do {
		double sinB = (xyz->z + dz) / sqrt(xyz->x * xyz->x + xyz->y * xyz->y + (xyz->z + dz) * (xyz->z + dz));
		double N = R / sqrt(1 - e2 * sinB * sinB);
		delta = dz - e2 * N * sinB;
		dz = e2 * N * sinB;
		iter += 1;
	} while (iter < 15 && delta > 1e-6);
	blh->longtitude = atan2(xyz->y, xyz->x);
	blh->latitude = atan2(xyz->z + dz, sqrt(xyz->x * xyz->x + xyz->y * xyz->y));
	blh->height = sqrt(xyz->x * xyz->x + xyz->y * xyz->y + (xyz->z + dz) * (xyz->z + dz)) - R / sqrt(1 - e2 * sin(blh->latitude) * sin(blh->latitude));
}

void BLHToNEUMat(const GEOCOOR* Blh, double Mat[])
{
	Mat[0] = -sin(Blh->longtitude);
	Mat[1] = cos(Blh->longtitude);
	Mat[2] = 0.0;
	Mat[3] = -sin(Blh->latitude) * cos(Blh->longtitude);
	Mat[4] = -sin(Blh->latitude) * sin(Blh->longtitude);
	Mat[5] = cos(Blh->latitude);
	Mat[6] = cos(Blh->latitude) * cos(Blh->longtitude);
	Mat[7] = cos(Blh->latitude) * sin(Blh->longtitude);
	Mat[8] = sin(Blh->latitude);
}

void CompSatElAz(const XYZ* Xr, const double Xs[], double* Elev, double* Azim,double* trop)
{
	GEOCOOR blh;
	double dE, dN, dH, dx, dy, dz;
	XYZToBLH(Xr, &blh, R_WGS84, F_WGS84);
	double Mat[9];
	BLHToNEUMat(&blh, Mat);
	dx = Xs[0] - Xr->x;
	dy = Xs[1] - Xr->y;
	dz = Xs[2] - Xr->z;
	dE = Mat[0] * dx + Mat[1] * dy + Mat[2] * dz;
	dN = Mat[3] * dx + Mat[4] * dy + Mat[5] * dz;
	dH = Mat[6] * dx + Mat[7] * dy + Mat[8] * dz;

	*Elev = atan(dH / sqrt(dE * dE + dN * dN));
	*Azim = atan2(dE, dN);

	*trop = Hopfield(blh.height, *Elev*180/pi);
}

void CompEnuPos(const XYZ* X0, const double Xr[], double dNeu[])
{
	GEOCOOR blh; 
	double dE, dN, dH, dx, dy, dz;
	XYZToBLH(X0, &blh, R_WGS84, F_WGS84);
	double Mat[9];
	BLHToNEUMat(&blh, Mat);
	dx = Xr[0] - X0->x;
	dy = Xr[1] - X0->y;
	dz = Xr[2] - X0->z;
	dE = Mat[0] * dx + Mat[1] * dy + Mat[2] * dz;
	dN = Mat[3] * dx + Mat[4] * dy + Mat[5] * dz;
	dH = Mat[6] * dx + Mat[7] * dy + Mat[8] * dz;
	dNeu[0] = dN;
	dNeu[1] = dE;
	dNeu[2] = dH;
}