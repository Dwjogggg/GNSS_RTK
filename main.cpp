#include<iostream>
#include<vector>
#include<cmath>
#include<iomanip>
#include <fstream>
#include"sockets.h"
#include "SPP_Struct.h"

using namespace std;

int mode;//文件读取模式或网络接收模式

int main()
{
	FILE *Rovfobs, *Basefobs;
	unsigned char buff[MAXBUFFNUM];
	int LenD, LenR,state;
	RAWDATA Raw;
	EPOCHOBSDATA obs;
	GPSEPHREC gephs[MAXGPSNUM], bephs[MAXBDSNUM];
	SDEPOCHOBS SDObs;
	DDCOBS DDObs;
	POSRES pos, bestPos;
	string sys;

	std::cout << "请选择模式：\n1. 读取文件模式\n0. 网络接收模式\n";
	std::cin >> mode;

	if (mode == 1) 
	{
		if (fopen_s(&Rovfobs, "D:\\大三\\卫导程序设计2\\zero-baseline-data\\oem719-202202021500-rover.bin", "rb") != 0)
		{
			printf("文件打开失败");
			return 0;
		}
		if (fopen_s(&Basefobs, "D:\\大三\\卫导程序设计2\\zero-baseline-data\\oem719-202202021500-base.bin", "rb") != 0)
		{
			printf("文件打开失败");
			return 0;
		}

		LenD = 0;
		while (!feof(Rovfobs)&& !feof(Basefobs))
		{
		/*if ((LenR = fread(buff + LenD, 1, MAXBUFFNUM - LenD, fobs)) != MAXBUFFNUM - LenD)    return 0;
		LenD += LenR;
		DecodeNovOem7Dat(buff, LenD, &obs, gephs, bephs, &bestPos);
		DetectOutlier(&obs);*/
			state=GetSynObs(Basefobs, Rovfobs, &Raw);
			//cout << "state: " << state << endl;
			cout << Raw.RovEpk.Time.Week << " " << Raw.RovEpk.Time.SecOfWeek << endl;
			DetectOutlier(&Raw.RovEpk);
			DetectOutlier(&Raw.BaseEpk);
			SDEpochObs(&Raw.RovEpk, &Raw.BaseEpk, &SDObs);
			SDDetectOutlier(&SDObs);
			SPP(&Raw.RovEpk, Raw.GpsEph, Raw.BdsEph, &pos);
			SPP(&Raw.BaseEpk, Raw.GpsEph, Raw.BdsEph, &bestPos);
			DetRefSat(&Raw.BaseEpk, &Raw.RovEpk, &SDObs, &DDObs);
		}
		fclose(Rovfobs);
		fclose(Basefobs);
	}

	else if (mode == 0)
	{
		SOCKET NetGps;
		if (OpenSocket(NetGps, "47.114.134.129", 7190) == false) {
			printf("Cannot connect to the server.\n");
		}
		LenD = 0;
		do {
			Sleep(980);
			if ((LenR = recv(NetGps, (char*)buff, MAXBUFFNUM, 0)) > 0)
			{
				memcpy(buff + LenD, buff, LenR);
				LenD += LenR;
				DecodeNovOem7Dat(buff, LenD, &obs, gephs, bephs, &bestPos);
				DetectOutlier(&obs);
			}
		} while (true);
	}
	else
	{
		printf("Invalid mode selected.\n");
		return 1;
	}

	return 0;
}