#pragma once
#include<iostream>
#include<vector>
#include<cmath>
#include<Eigen/Core>
#include<Eigen/Dense>

using namespace std;

#define pi 3.141592653589793
#define R_WGS84 6378137.0
#define F_WGS84 (1.0 / 298.257223563)
#define POLYCRC32   0xEDB88320u         /* CRC32 polynomial */
#define FILEMODE 1
#define MAXBDSNUM 63                
#define MAXGPSNUM  32
#define MAXCHANNUM 36
#define MAXBUFFNUM 60000
#define MAXRAWLEN 40960

#define C_Light 299792458.0      
#define  FG1_GPS  1575.42E6             /* L1信号频率 */
#define  FG2_GPS  1227.60E6             /* L2信号频率 */
#define  WL1_GPS  (C_Light/FG1_GPS)
#define  WL2_GPS  (C_Light/FG2_GPS)

#define  FG1_BDS  1561.098E6            /* B1信号的基准频率 */
#define  FG3_BDS  1268.520E6            /* B3信号的基准频率 */
#define  WL1_BDS  (C_Light/FG1_BDS)
#define  WL3_BDS  (C_Light/FG3_BDS)       // 波长

#define  GPS_GM 3.986005E14
#define GPS_radv 7.2921151467E-5

#define  BDS_GM 3.986004418E14
#define BDS_radv 7.2921150E-5

#define F_RELATIV -4.442807633E-10       //相对论改正

#define Hop_H0 0.000                    /*Hopfield海平面（m）*/
#define Hop_T0 15+273.16                /*Hopfield温度（K）*/
#define Hop_p0 1013.25                  /*Hopfield气压（mbar）*/
#define Hop_RH0 0.500                   /*Hopfield相对湿度*/

#define GF_THRESHOLD 0.05
#define MW_THRESHOLD 3.0


enum GNSSys {BDS,GPS,UNKS};             

/*通用时间定义*/
struct COMMONTIME
{
	short          Year;
	unsigned short Month;
	unsigned short Day;
	unsigned short Hour;
	unsigned short Minute;
	double         Second;
};

/*简化儒略日*/
struct MJDTIME 
{
	int    Days;
	double FracDay;
	MJDTIME()
	{
		Days = 0;
		FracDay = 0.0;
	}
};

/*GPS时间定义*/
struct GPSTIME
{
	unsigned short    Week;
	double            SecOfWeek;
	GPSTIME()
	{
		Week = 0;
		SecOfWeek = 0.0;
	}
};
/*地心地固坐标系*/
union XYZ 
{
	struct
	{
		double x;
		double y;
		double z;
	};
	double xyz[3];
};
/*大地坐标系*/
union GEOCOOR
{
	struct
	{
		double latitude;
		double longtitude;
		double height;
	};
	double Blh[3];
};
/*站心坐标系*/
union NEU {
	struct
	{
		double dN;
		double dE;
		double dU;
	};
	double Neu[3];
};
/*卫星位置与钟差结构体定义*/
struct SATPVT
{
	double SatPos[3], SatVel[3];
	double SatClkOft, SatClkSft;
	double Elevation, Azimuth;
	double TropCorr;
	double Tgd1, Tgd2;
	bool Valid;//false=没有星历或星历过期，true=计算成功

	SATPVT()
	{
		SatPos[0] = SatPos[1] = SatPos[2] = 0.0;
		SatVel[0] = SatVel[1] = SatVel[2] = 0.0;
		SatClkOft = SatClkSft = 0.0;
		Elevation = pi / 2.0;
		Azimuth = 0.0;

		double TropCorr = Tgd1 = Tgd2 = 0.0;
		Valid = false;
	}
};

/* 每颗卫星的观测数据定义 */
struct SATOBSDATA 
{
	short Prn;
	GNSSys System;
	double P[2], L[2], D[2];
	double cn0[2], LockTime[2];// C/N0, 载波锁定时间
	unsigned char half[2];     // ParityFlag 半周跳
	bool Valid;

	SATOBSDATA()
	{
		Prn = 0;
		System = UNKS;
		for (int i = 0; i < 2; i++)
			P[i] = L[i] = D[i] = 0.0;
		Valid = false;
	}
};
/* 每颗卫星的组合观测数据定义 */
struct MWGF
{
	short Prn;//卫星号
	GNSSys Sys;
	double MW;
	double GF;
	double PIF;
	int n;    //平滑计数

	MWGF()
	{
		Prn = n = 0;
		Sys = UNKS;
		MW = GF = PIF = 0.0;
	}
};

/* 每个历元的定位结果结构体定义 */
struct POSRES
{
	GPSTIME Time;
	double Pos[3];
	double Vel[3];
	double PDOP, Sigma, SigmaVel;
	double RcvClkOft[2];               /* 0 为GPS钟差; 1=BDS钟差 */
	double RcvClkSft;
	int GPSSatNum;
	int BDSSatNum;      /* 定位使用的GPS卫星数 */
	int SatNum, SVs, solnSVs;
	bool IsSuccess;                /* 定位是否成功, 1为成功, 0为失败 */

	POSRES()
	{
		for (int i = 0; i < 3; i++)	Pos[i] = Vel[i] = 0.0;
		RcvClkOft[0] = RcvClkOft[1] = RcvClkSft = 0.0;
		PDOP = Sigma = SigmaVel = 999.9;
		GPSSatNum = BDSSatNum = SatNum = 0;
		IsSuccess = false;
	}
};

/* 每个历元的观测数据定义 */
struct EPOCHOBSDATA
{
	GPSTIME Time;
	short SatNum;
	SATOBSDATA SatObs[MAXCHANNUM];
	MWGF ComObs[MAXCHANNUM];
	SATPVT SatPVT[MAXCHANNUM]; //卫星位置与钟差
	POSRES Res; //保存基站或NovAtel接收机定位结果
	EPOCHOBSDATA()
	{
		SatNum = 0;
		Res.Pos[0] = Res.Pos[1] = Res.Pos[2] = 0.0;
	}
};

/*  每颗卫星的单差观测数据定义  */
struct SDSATOBS
{
	short    Prn;
	GNSSys  System;
	bool    Valid;
	double   dP[2], dL[2];   // 伪距单差、相位单差（m）
	short    nBas, nRov;   // 存储单差观测值对应的基准和流动站的数值索引号

	SDSATOBS()
	{
		Prn = nBas = nRov = 0;
		System = UNKS;
		dP[0] = dL[0] = dP[1] = dL[1] = 0.0;
		Valid = false;
	}
};

/*  每个历元的单差观测数据定义  */
struct SDEPOCHOBS
{
	GPSTIME    Time;
	short      SatNum;
	SDSATOBS   SdSatObs[MAXCHANNUM];
	MWGF       SdCObs[MAXCHANNUM];

	SDEPOCHOBS()
	{
		SatNum = 0;
	}
};

/* 双差相关的数据定义 */
struct DDCOBS
{
	int RefPrn[2], RefPos[2];         // 参考星卫星号与存储位置，0=GPS; 1=BDS
	int Sats, DDSatNum[2];            // 待估的双差模糊度数量，0=GPS; 1=BDS
	double FixedAmb[MAXCHANNUM * 4];  // 包括双频最优解[0,AmbNum]和次优解[AmbNum,2*AmbNum]
	double ResAmb[2], Ratio;          // LAMBDA浮点解中的模糊度残差
	float  FixRMS[2];                 // 固定解定位中rms误差
	double dPos[3];                   // 基线向量
	bool bFixed;                      // true为固定，false为未固定

	DDCOBS()
	{
		int i;
		for (i = 0; i < 2; i++) {
			DDSatNum[i] = 0;    // 各卫星系统的双差数量
			RefPos[i] = RefPrn[i] = -1;
		}
		Sats = 0;              // 双差卫星总数
		dPos[0] = dPos[1] = dPos[2] = 0.0;
		ResAmb[0] = ResAmb[1] = FixRMS[0] = FixRMS[1] = Ratio = 0.0;
		bFixed = false;
		for (i = 0; i < MAXCHANNUM * 2; i++)
		{
			FixedAmb[2 * i + 0] = FixedAmb[2 * i + 1] = 0.0;
		}
	}
};

/* 每颗卫星的历元数据定义 */
struct GPSEPHREC
{
	short PRN;
	GNSSys System;
	GPSTIME TOC, TOE;
	double ClkBias, ClkDrift, ClkDriftRate;
	double IODE, IODC;
	double SqrtA, M0, e, OMEGA0, i0, omega;
	double Crs, Cuc, Cus, Cic, Cis, Crc;
	double DetlaN, OMEGADot, iDOT;
	int SVHealth;
	double TGD1, TGD2;
};

/*  RTK定位的数据定义  */
struct RAWDATA 
{
	EPOCHOBSDATA BaseEpk;
	EPOCHOBSDATA RovEpk;
	SDEPOCHOBS SDObs;
	DDCOBS DDObs;
	GPSEPHREC GpsEph[MAXGPSNUM], BdsEph[MAXBDSNUM];
};

// 定义结构体用来存储上一时刻locktime和当前时刻观测数据
struct LockTimeRecord 
{
	EPOCHOBSDATA preObs;
	EPOCHOBSDATA curObs;
	//初始化
	LockTimeRecord() {
		preObs.SatNum = 0;
		curObs.SatNum = 0;//以卫星数量为0初始化
	}
};

//时间转换函数声明
void CommonTimeToMjdTime(const COMMONTIME& ct, MJDTIME& mjd);
void MjdTimeToCommonTime(const MJDTIME& mjd, COMMONTIME& ct);
void MjdTimeToGPSTime(const MJDTIME& mjd, GPSTIME& gps);
void GPSTimeToMjdTime(const GPSTIME& gps, MJDTIME& mjd);
void CommonTimeToGPSTime(const COMMONTIME& ct, GPSTIME& gps);
void GPSTimeToCommonTime(const GPSTIME& gps, COMMONTIME& ct);

//坐标转换函数声明
void BLHToXYZ(const GEOCOOR* blh, XYZ* xyz, const double R, const double F);
void XYZToBLH(const XYZ* xyz, GEOCOOR* blh, const double R, const double F);
void BLHToNEUMat(const GEOCOOR* Blh, double Mat[]);
void CompSatElAz(const XYZ* Xr, const double Xs[], double* Elev, double* Azim, double* trop);
void CompEnuPos(const XYZ* X0, const double Xr[], double dNeu[]);

//小端解码函数
double R8(unsigned char* p);
float R4(unsigned char* p);
int I4(unsigned char* p);
unsigned int UI4(unsigned char* p);
short I2(unsigned char* p);
unsigned short UI2(unsigned char* p);
unsigned char UI1(unsigned char* p);

// CRC32校验函数
unsigned int crc32(const unsigned char* buff, int len);

//解码主函数
int DecodeNovOem7Dat(unsigned char Buff[], int& Len, EPOCHOBSDATA*obs, GPSEPHREC geph[], GPSEPHREC beph[], POSRES* pos);
//解码子函数
void decode_rangeb_oem7(unsigned char* data, EPOCHOBSDATA* obs);//解码OEM7格式的观测数据
//GPS星历
int decode_gpsephem(unsigned char* buff, GPSEPHREC* geph);
//BDS星历
int decode_bdsephem(unsigned char* buff, GPSEPHREC* beph);
//接收机解算结果
void decode_psrpos(unsigned char* buff, POSRES* pos);

//卫星钟差与卫星钟速计算
bool CompSatClkOff(const int Prn, const GNSSys Sys, const GPSTIME* t, GPSEPHREC* GPSEph, GPSEPHREC* BDSEph, SATPVT* Mid);
//卫星位置与钟差
int CompGPSSatPVT(/*const*/int Prn, const GPSTIME* t, const GPSEPHREC* Eph, SATPVT* Mid);
int CompBDSSatPVT(/*const*/ int Prn, const GPSTIME* t, const GPSEPHREC* Eph, SATPVT* Mid);

//对流层延迟改正
double Hopfield(const double H,const double Elev);
//粗差探测
void DetectOutlier(EPOCHOBSDATA* Obs);

//信号发射时刻卫星位置计算
void ComputeSatPVTAtSignalTrans(EPOCHOBSDATA* Epk, GPSEPHREC* Eph, GPSEPHREC* BDSEph, double UserPos[3]);
//SPP函数
bool SPP(EPOCHOBSDATA* Epoch, GPSEPHREC* GPSEph, GPSEPHREC* BDSEph, POSRES* Res);
//SPV函数
bool SPV(EPOCHOBSDATA* Epoch, POSRES* Res);

// RTK
bool RTKSPP(EPOCHOBSDATA* Epoch, GPSEPHREC* GPSEph, GPSEPHREC* BDSEph, POSRES* Res);
int GetSynObs(FILE* FBas, FILE* FRov, RAWDATA* Raw);
void SDDetectOutlier(SDEPOCHOBS* Obs);
void SDEpochObs(const EPOCHOBSDATA* RovEpk, const EPOCHOBSDATA* BaseEpk, SDEPOCHOBS* SDObs);
void DetRefSat(const EPOCHOBSDATA* epkA, const EPOCHOBSDATA* epkB, SDEPOCHOBS* SDObs, DDCOBS* DDObs);
void DDEpochObs(const EPOCHOBSDATA* RovEpk, const EPOCHOBSDATA* BaseEpk, SDEPOCHOBS* SDObs, DDCOBS* DDObs);
bool RTKFloat(RAWDATA* Raw, POSRES* Base, POSRES* Rov);