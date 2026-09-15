#include<iostream>
#include<vector>
#include<cmath>
#include "SPP_Struct.h"

using namespace std;

extern int mode;

double R8(unsigned char* p)
{
	double r;
	memcpy(&r, p, 8);
	return r;
}

float R4(unsigned char* p)
{
	float r;
	memcpy(&r, p, 4);
	return r;
}

int I4(unsigned char* p)
{
	int r;
	memcpy(&r, p, 4);
	return r;
}

unsigned int UI4(unsigned char* p)
{
	unsigned int r;
	memcpy(&r, p, 4);
	return r;
}

short I2(unsigned char* p)
{
	short r;
	memcpy(&r, p, 2);
	return r;
}

unsigned short UI2(unsigned char* p)
{
	unsigned short r;
	memcpy(&r, p, 2);
	return r;
}

unsigned char UI1(unsigned char* p)
{
	unsigned char r;
	memcpy(&r, p, 1);
	return r;
}

unsigned int crc32(const unsigned char* buff, int len)
{
	int i, j;
	unsigned int crc = 0;
	for (i = 0; i < len; i++)
	{
		crc ^= buff[i];
		for (j = 0; j < 8; j++)
		{
			if (crc & 1)
				crc = (crc >> 1) ^ POLYCRC32;
			else
				crc >>= 1;
		}
	}
	return crc;
}

void decode_rangeb_oem7(unsigned char* data, EPOCHOBSDATA* Epkobs)
{
	EPOCHOBSDATA obs;
	int i, j, n, k, ObsNum, Freq, Prn;
	GNSSys sys;
	double wl;
	int PhaseLockFlag, CodeLockedFlag, ParityFlag, SatSystem, SigType;
	unsigned int ChanStatus;
	unsigned char* p = data + 28;

	LockTimeRecord Epk;
	
	//1. 从消息头中解码得到观测时刻，该时刻为接收机钟表面时，用GPSTIME结构体表示。
	obs.Time.Week = UI2(data + 14);
	obs.Time.SecOfWeek = UI4(data + 16) * 1E-3;
	//    2. 解码得到观测值数量，为所有卫星所有信号观测值的总数。
	ObsNum = UI4(p);  // 观测卫星数量
	memset(obs.SatObs, 0, MAXCHANNUM * sizeof(SATOBSDATA));
	//    3. 对所有信号观测值进行循环解码
	for (i = 0, p += 4; i < ObsNum; i++, p += 44) {
		//    ① 解码得到跟踪状态标记，从中取出Phase lock flag / Code locked
		//    flag / Parity known flag / Satellite system / signal type等数据
		ChanStatus = UI4(p + 40);
		ParityFlag = (ChanStatus >> 11) & 0x01;
		PhaseLockFlag = (ChanStatus >> 10) & 0x01;
		CodeLockedFlag = (ChanStatus >> 12) & 0x01;
		SatSystem = (ChanStatus >> 16) & 0x07;
		SigType = (ChanStatus >> 21) & 0x1F;
		//    ② 如果卫星系统不是GPS或BDS， continue至①

	//    ③ 如果GPS卫星的信号类型不是L1 C / A或者L2P（Y）， BDS卫星不是B1I
	//    或B3I， continue至①，并记录信号频率类型，第一频率s = 0，第二频
	//    率s = 1；
		if (SatSystem == 0) {
			sys = GPS;
			if (SigType == 0) {
				Freq = 0; wl = WL1_GPS;
			}
			else if (SigType == 9) { Freq = 1; wl = WL2_GPS; }
			else continue;
		}
		else if (SatSystem == 4) {
			sys = BDS;
			if (SigType == 0 || SigType == 4) {
				Freq = 0; wl = WL1_BDS;
			}
			else if (SigType == 2 || SigType == 6) { Freq = 1; wl = WL3_BDS; }
			else continue;
		}
		else continue;

		//    ④ 解码得到卫星号Prn以及卫星系统号，在当前观测值结构体中进行搜索，
		//    如果找到相同的卫星，将解码的观测值填充到该卫星对应的数组中；如
		//    果在当前已解码的卫星数据中没有发现，则填充到现有数据的末尾。
		Prn = UI2(p);
		for (j = 0; j < MAXCHANNUM; j++) {
			if (obs.SatObs[j].System == sys && obs.SatObs[j].Prn == Prn) {
				n = j; break;
			}
			if (obs.SatObs[j].Prn == 0) {
				k = n = j; break;
			}
		}
		obs.SatObs[n].Prn = Prn;
		obs.SatObs[n].System = sys;
		obs.SatObs[n].P[Freq] = CodeLockedFlag == 1 ? R8(p + 4) : 0.0;
		obs.SatObs[n].L[Freq] = -wl * (PhaseLockFlag == 1 ? R8(p + 16) : 0.0);
		obs.SatObs[n].D[Freq] = -wl * R4(p + 28);
		obs.SatObs[n].cn0[Freq] = R4(p + 32);     // C/N0
		obs.SatObs[n].LockTime[Freq] = R4(p + 36);// 载波锁定时间
		obs.SatObs[n].half[Freq] = ParityFlag;    // 半周跳
		obs.SatObs[n].Valid = true;               // 解到数据认为有效，后续逐步排除
	}
	obs.SatNum = k + 1;
	// Locktime：判断历元间相位数据的连续性，没有发生中断时，Locktime(k + 1) >= Locktime(k)，数据正常；反之有周跳
	Epk.curObs = obs;
	Epk.preObs = *Epkobs;
	// 判断是否为第一个历元，是则设置为ture进行下个历元检查
	if (Epk.preObs.SatNum == 0) 
	{
		Epk.preObs = Epk.curObs;
		*Epkobs = Epk.curObs;
		return;
	}
	// 若不是则与上个历元每个卫星进行比对
	else 
	{
		for (i = 0;i < Epk.curObs.SatNum;i++)
		{
			SATOBSDATA& curSat = Epk.curObs.SatObs[i];
			// Parity：判断相位数据是否存在半周，Parity = 0，可能存在半周，valid = false
			if (curSat.System == BDS && curSat.Prn == 16)
			{
				int bb = 0;
			}

			if (curSat.half[0] == 0 || curSat.half[1] == 0)
			{
				curSat.Valid = false;
				continue;
			}
			// 在上个历元中寻找相同卫星，比对Locktime
			bool findComSat = false;
			for(j=0;j<Epk.preObs.SatNum;j++)
			{
				if(curSat.System == Epk.preObs.SatObs[j].System && curSat.Prn == Epk.preObs.SatObs[j].Prn)
				{
					findComSat = true;
					// 找到相同卫星，比较Locktime，Locktime(k + 1) >= Locktime(k)，数据正常；反之有周跳
					if (curSat.LockTime[0] < Epk.preObs.SatObs[j].LockTime[0] || curSat.LockTime[1] < Epk.preObs.SatObs[j].LockTime[1])
					{
						// 去掉这颗卫星的观测数据
						curSat.Valid = false;
						break;
					}
					else
					{
						break;
					}
				}
			}
			// 若找不到，说明该卫星是新进卫星，跳过
			if (!findComSat) continue;
		}
	}
	// 保存当前历元数据为上个历元数据
	Epk.preObs = Epk.curObs;
	*Epkobs = Epk.curObs;
}

int decode_gpsephem(unsigned char* buff, GPSEPHREC* geph)
{
	int prn;
	unsigned char* p = buff + 28;
	GPSEPHREC* eph;

	prn = UI4(p);
	if (prn < 1 || prn > MAXGPSNUM) return 0;

	eph = geph + prn - 1;
	eph->PRN = prn;
	eph->System = GPS;

	eph->SVHealth = UI4(p + 12);

	eph->TOC.Week = eph->TOE.Week = UI4(p + 24);
	eph->TOC.SecOfWeek = R8(p + 164);
	eph->TOE.SecOfWeek = R8(p + 32);

	eph->ClkBias = R8(p + 180);
	eph->ClkDrift = R8(p + 188);
	eph->ClkDriftRate = R8(p + 196);

	eph->IODE = UI4(p + 16);
	eph->IODC = UI4(p + 160);

	eph->SqrtA = sqrt(R8(p + 40));
	eph->M0 = R8(p + 56);
	eph->e = R8(p + 64);
	eph->omega = R8(p + 72);//近地幅角
	eph->OMEGA0 = R8(p + 144);//升交点经度
	eph->i0 = R8(p + 128);

	eph->Cuc = R8(p + 80);
	eph->Cus = R8(p + 88);
	eph->Crc = R8(p + 96);
	eph->Crs = R8(p + 104);
	eph->Cic = R8(p + 112);
	eph->Cis = R8(p + 120);
	eph->iDOT = R8(p + 136);
	eph->OMEGADot = R8(p + 152);
	eph->DetlaN = R8(p + 48);

	return 1;
}

int decode_bdsephem(unsigned char* buff, GPSEPHREC* beph)
{
	int prn;
	unsigned char* p = buff + 28;
	GPSEPHREC* eph;

	prn = UI4(p);//
	if (prn < 1 || prn > MAXBDSNUM) return 0;

	eph = beph + prn - 1;
	eph->PRN = prn;
	eph->System = BDS;

	eph->TOC.Week = eph->TOE.Week = UI4(p + 4);

	eph->SVHealth = UI4(p + 16);

	eph->TGD1 = R8(p + 20);
	eph->TGD2 = R8(p + 28);

	eph->TOC.SecOfWeek = UI4(p + 40);
	eph->TOE.SecOfWeek = UI4(p + 72);

	eph->ClkBias = R8(p + 44);
	eph->ClkDrift = R8(p + 52);
	eph->ClkDriftRate = R8(p + 60);

	eph->IODC = UI4(p + 36);  //对于BDS，这实际上是AODC
	eph->IODE = UI4(p + 68);  //对于BDS，这实际上是AODE

	eph->SqrtA = R8(p + 76);
	eph->M0 = R8(p + 108);
	eph->e = R8(p + 84);
	eph->omega = R8(p + 92);  //近地幅角
	eph->OMEGA0 = R8(p + 116);//升交点经度
	eph->i0 = R8(p + 132);

	eph->Crs = R8(p + 172);
	eph->Cuc = R8(p + 148);
	eph->Cus = R8(p + 156);
	eph->Cic = R8(p + 180);
	eph->Cis = R8(p + 188);
	eph->Crc = R8(p + 164);
	eph->DetlaN = R8(p + 100);
	eph->OMEGADot = R8(p + 124);
	eph->iDOT = R8(p + 140);

	return 1;
}

void decode_psrpos(unsigned char* buff, POSRES* pos)
{
	GEOCOOR blh;
	XYZ xyz;
	unsigned char* p = buff + 28;
	pos->Time.Week = UI2(buff + 14);
	pos->Time.SecOfWeek = UI4(buff + 16);
	blh.latitude = R8(p + 8)*pi/180;
	blh.longtitude = R8(p + 16) * pi / 180;
	blh.height = R8(p + 24)+R4(p+32);

	BLHToXYZ(&blh, &xyz, R_WGS84, F_WGS84);

	pos->Pos[0] = xyz.x;
	pos->Pos[1] = xyz.y;
	pos->Pos[2] = xyz.z;
	/*pos->SigmaPos[0] = R4(p + 40);
	pos->SigmaPos[1] = R4(p + 44);
	pos->SigmaPos[2] = R4(p + 48);
	pos->SVs = UI1(p + 64);
	pos->solnSVs = UI1(p + 65);*/
}

int DecodeNovOem7Dat(unsigned char Buff[], int& Len, EPOCHOBSDATA* obs, GPSEPHREC geph[], GPSEPHREC beph[], POSRES* pos)
{
	if (Buff == nullptr || obs == nullptr || geph == nullptr || beph == nullptr)
	{
		cout << "Error: DecodeNovOem7Dat: Null pointer passed." << endl;
		return -1;
	}

	int i,MsgLen,MsgID,Status;
	i = 0;
	while (1)
	{
		//    1.设置循环变量i = 0，开始查找AA 44 12同步字符
		for (;i < Len - 2;i++) 
		{
			if (Buff[i] == 0xAA && Buff[i + 1] == 0x44 && Buff[i + 2] == 0x12)break;
		}
		//    2.找到同步字符后，获取消息头长度的字符28字节， 若字节数量不足即i + 28 > len，跳出循环（break）至第6步
		if (i + 28 > Len) break;
		//    3. 从消息头中解码消息长度MsgLen和消息类型MsgID，获得整条消息buff[i, i + 28 + MsgLen + 4]，若字节数量不足即i + 28 + MsgLen + 4 > len，跳出循环（break） 到第6步
		MsgID = UI2(Buff + i + 4);
		MsgLen = UI2(Buff + i + 8);
		if (i + 28 + MsgLen + 4 > Len) break;

		//    4. CRC检验，若不通过，跳过同步字符3个字节，即i = i + 3，返回到第1步
		if (crc32(Buff + i, 28 + MsgLen) != UI4(Buff + i + 28 + MsgLen))
		{
			i += 3; //跳过当前AA 44 12
			continue;
		}
		//    5. 根据消息ID，调用对应的解码函数，若解码得到观测值，跳出循环至第6步；否则，跳过整条消息，即i = i + 28 + MsgLen + 4，返回第1步
		Status = 0;
		switch(MsgID)
		{
		case 43:
			decode_rangeb_oem7(Buff + i, obs);
			Status = 1;//历元只读一次
			break;
		case 7:
			decode_gpsephem(Buff + i, geph);
			break;
		case 1696:
			decode_bdsephem(Buff + i, beph);
			break;
		case 42:
			decode_psrpos(Buff + i, pos);
			break;
		default:
			printf("Unknown message ID: 0x%04X\n", MsgID);
			break;
		}
		i = i + 28 + MsgLen + 4;
		if (Status == 1 && mode == 1) break;

		//    6. 循环结束后，将不足一条消息的剩余字节，拷贝至buff缓冲区的开始处，将len设置为剩余字节数量，并返回给主函数。
	}
	memcpy(Buff + 0, Buff + i, Len - i);
	Len = Len - i;
	return Status;
}
