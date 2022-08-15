#include "myutils.h"
#include "seqdb.h"
#include "xdpmem.h"
#include "objmgr.h"
#include "seq.h"
#include "alndata.h"
#include "pathinfo.h"

#define LOG_ALNS	0

void RevCompInPlace(string &Seq);

/***
Similar approach to:
A simple algorithm for detecting circular permutations in proteins
S. Uliel, A. Fliess, A. Amir, R. Unger
Bioinformatics 1999
https://doi.org/10.1093/bioinformatics/15.11.930
***/

float ViterbiFastMem(XDPMem &Mem, const string &A, const string &B, PathInfo &PI);
void RevCompSeq(const byte *Seq, unsigned L, byte *RCSeq);
void WriteAln(FILE *f, const byte *A, const byte *B,
  const char *Path, unsigned ColCount);

static float FixPath(const string &Q, const string &T,
  const string &DoublePath, uint &StartPosQ, string &QRow, string &TRow)
	{
	const uint QL = SIZE(Q);
	const uint TL = SIZE(T);

	QRow.clear();
	TRow.clear();

	const uint ColCount = SIZE(DoublePath);

// For validation only
	{
	uint NM = 0;
	uint ND = 0;
	uint NI = 0;
	for (uint i = 0; i < ColCount; ++i)
		{
		char c = DoublePath[i];
		if (c == 'M')
			++NM;
		else if (c == 'D')
			++ND;
		else if (c == 'I')
			++NI;
		else
			asserta(false);
		}

	asserta(NM + ND == 2*QL);
	asserta(NM + NI == TL);
	}

	uint StartCol = UINT_MAX;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = DoublePath[Col];
		if (StartCol == UINT_MAX)
			{
			if (c == 'M')
				{
				StartCol = Col;
				break;
				}
			}
		}
	asserta(StartCol != UINT_MAX);

	uint TPos = 0;
	StartPosQ = 0;
	for (uint Col = 0; Col < StartCol; ++Col)
		{
		char c = DoublePath[Col];
		switch (c)
			{
		case 'M':
			asserta(false);
			break;
		case 'D':
			++StartPosQ;
			break;
		case 'I':
			++TPos;
			break;
		default:
			asserta(false);
			}
		}

	uint IdCount = 0;
	uint QPos = 0;
	for (uint Col = StartCol; Col < ColCount; ++Col)
		{
		asserta(TPos <= TL && QPos <= QL);
		uint CirclePosQ = (StartPosQ + QPos)%QL;
		char c = DoublePath[Col];
		switch (c)
			{
		case 'M':
			{
			byte q = Q[CirclePosQ];
			byte t = T[TPos++];
			++QPos;
			QRow += q;
			TRow += t;
			if (toupper(q) == toupper(t))
				++IdCount;
			asserta(TPos <= TL && QPos <= QL);
			break;
			}

		case 'D':
			{
			byte q = Q[CirclePosQ];
			++QPos;
			QRow += q;
			TRow += '-';
			asserta(TPos <= TL && QPos <= QL);
			break;
			}

		case 'I':
			{
			byte t = T[TPos++];
			TRow += t;
			QRow += '-';
			asserta(TPos <= TL && QPos <= QL);
			break;
			}

		default:
			asserta(false);
			}

		if (TPos == TL || QPos == QL)
			break;
		}

	while (TPos < TL)
		{
		QRow += '-';
		TRow += T[TPos++];
		}
	while (QPos < QL)
		{
		QRow += Q[QPos++];
		TRow += '-';
		}

	uint SingleColCount = SIZE(QRow);
	asserta(SIZE(TRow) == SingleColCount);
	float FractId = float(IdCount)/SingleColCount;
	return FractId;
	}

void AlignPairPlusOnly(const string &LabelQ, const string &Q,
  const string &LabelT, const string &T, AlnData &AD)
	{
	AD.Clear();

	XDPMem Mem;

	const uint QL = SIZE(Q);
	const uint TL = SIZE(T);
	const uint QL2 = QL*2;

// Q2 is tandem duplication of Q
	string Q2 = Q + Q;

	PathInfo *PI = ObjMgr::GetPathInfo();

	ViterbiFastMem(Mem, Q2, T, *PI);

	PI->GetPathStr(AD.m_DoublePath);

	uint StartPosQ = UINT_MAX;

	string QRow;
	string TRow;
	float FractId = FixPath(Q, T, AD.m_DoublePath, StartPosQ, QRow, TRow);

#if LOG_ALNS
	Log("______________________________________________________________________\n");
	Log("Q>%s\n", LabelQ.c_str());
	Log("T>%s\n", LabelT.c_str());
	Log("\nDouble:\n");
	WriteAln(g_fLog, Q2, T, DoublePath.c_str(), SIZE(DoublePath));
	Log("\nSingle (Start=%u, QL=%d):\n", StartPosQ, QL);
	Log("Q:  %s\n", QRow.c_str());
	Log("T:  %s\n", TRow.c_str());
#endif

	AD.m_LabelQ = LabelQ;
	AD.m_LabelT = LabelT;
	AD.m_SeqQ = Q;
	AD.m_SeqT = T;
	AD.m_Plus = true;
	AD.m_StartPosQ = StartPosQ;
	AD.m_FractId = FractId;
	AD.m_QRow = QRow;
	AD.m_TRow = TRow;
	}

void AlignPair(
  const string &LabelQ, const string &SeqQ,
  const string &LabelT, const string &SeqT,
  AlnData &AD)
	{
	string SeqQRC = SeqQ;
	RevCompInPlace(SeqQRC);

	AlnData ADPlus;
	AlnData ADMinus;

	AlignPairPlusOnly(LabelQ, SeqQ, LabelT, SeqT, ADPlus);

	if (opt(plusonly))
		ADMinus.m_FractId = -1.0;
	else
		AlignPairPlusOnly(LabelQ, SeqQRC, LabelT, SeqT, ADMinus);

	if (ADPlus.m_FractId >= ADMinus.m_FractId)
		{
		AD = ADPlus;
		AD.m_Plus = true;
		}
	else
		{
		AD = ADMinus;
		AD.m_Plus = false;
		}
	}
