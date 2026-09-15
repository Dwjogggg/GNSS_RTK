#include<iostream>
#include<cmath>
#include<vector>
#include "SPP_Struct.h"

using namespace std;

void CommonTimeToMjdTime(const COMMONTIME& ct, MJDTIME& mjd)
{
	int y, m;
	if (ct.Month <= 2) {
		y = ct.Year - 1;
		m = ct.Month + 12;
	}
	else {
		y = ct.Year;
		m = ct.Month;
	}
	mjd.Days = (int)(365.25 * y) + (int)(30.6001 * (m + 1)) + ct.Day - 679019;
	mjd.FracDay = (ct.Hour + ct.Minute / 60.0 + ct.Second / 3600.0) / 24.0;
}

void MjdTimeToCommonTime(const MJDTIME& mjd, COMMONTIME& ct)
{
	int a, b, c, d, e;
	a = mjd.Days + 2400001;
	b = a + 1537;
	c = (int)((b - 122.1) / 365.25);
	d = (int)(365.25 * c);
	e = (int)((b - d) / 30.6001);
	ct.Day = b - d - (int)(30.6001 * e);
	ct.Month = e - 1 - 12 * (int)(e / 14.0);
	ct.Year = c - 4715 - (int)((7 + ct.Month) / 10.0);
	ct.Hour = (int)(24 * mjd.FracDay);
	ct.Minute = (int)(60 * (24 * mjd.FracDay - ct.Hour));
	ct.Second = 3600 * (24 * mjd.FracDay - ct.Hour - ct.Minute / 60.0);
}

void MjdTimeToGPSTime(const MJDTIME& mjd, GPSTIME& gps) 
{
	gps.Week = (int)((mjd.Days - 44244) / 7.0);
	gps.SecOfWeek = (mjd.Days - 44244 - gps.Week * 7) * 86400 + mjd.FracDay * 86400;
}

void GPSTimeToMjdTime(const GPSTIME& gps, MJDTIME& mjd)
{
	mjd.Days = 44244 + gps.Week * 7 + (int)(gps.SecOfWeek / 86400.0);
	mjd.FracDay = (gps.SecOfWeek - (int)(gps.SecOfWeek / 86400.0) * 86400)/86400.0;
}

void CommonTimeToGPSTime(const COMMONTIME& ct, GPSTIME& gps)
{
	MJDTIME mjd;
	CommonTimeToMjdTime(ct, mjd);
	MjdTimeToGPSTime(mjd, gps);
}

void GPSTimeToCommonTime(const GPSTIME& gps, COMMONTIME& ct)
{
	MJDTIME mjd;
	GPSTimeToMjdTime(gps, mjd);
	MjdTimeToCommonTime(mjd, ct);
}

