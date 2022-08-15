#include "myutils.h"
#include "alndata.h"

bool ReportOverhangs(FILE *f, const AlnData &AD)
	{
	if (f == 0)
		return false;

	const string &Q = AD.m_SeqQ;
	const string &T = AD.m_SeqT;
	const string &DoublePath = AD.m_DoublePath;
	const uint QL = SIZE(Q);
	const uint TL = SIZE(T);

	string Q2 = Q + Q;

	const uint ColCount = SIZE(DoublePath);

	uint FirstMCol = UINT_MAX;
	uint LastMCol = UINT_MAX;
	uint QCount = 0;
	uint NM = 0;
	uint ND = 0;
	uint NI = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = DoublePath[Col];
		switch (c)
			{
		case 'M':
			{
			++NM;
			if (FirstMCol == UINT_MAX)
				FirstMCol = Col;
			LastMCol = Col;
			break;
			}

		case 'D':
			{
			++ND;
			break;
			}

		case 'I':
			{
			++NI;
			break;
			}

		default:
			asserta(false);
			continue;
			}
		}
	asserta(FirstMCol != UINT_MAX);
	asserta(LastMCol != UINT_MAX);
	asserta(NM + ND == 2*QL);
	asserta(NM + NI == TL);

	uint NQ = 0;
	for (uint Col = FirstMCol; Col <= LastMCol; ++Col)
		{
		char c = DoublePath[Col];
		if (c == 'M' || c == 'D')
			++NQ;
		}

	if (NQ == QL)
		return false;

	string QRow;
	string TRow;
	string AnnotRow;
	string PathRow;
	uint QPos = 0;
	uint TPos = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = DoublePath[Col];
		PathRow += c;
		switch (c)
			{
		case 'M':
			{
			char q = Q2[QPos++];
			char t = T[TPos++];
			if (q == t)
				AnnotRow += '|';
			else
				AnnotRow += ' ';
			QRow += q;
			TRow += t;
			break;
			}
		case 'D':
			{
			QRow += Q2[QPos++];
			TRow += '-';
			AnnotRow += ' ';
			break;
			}
		case 'I':
			{
			QRow += '-';
			TRow += T[TPos++];
			AnnotRow += ' ';
			break;
			}
		default:
			asserta(false);
			}
		}

	asserta(TPos == TL);
	asserta(f != 0);

	fprintf(f, "\n");
	fprintf(f, "================================================\n");
	fprintf(f, "Q>%s\n", AD.m_LabelQ.c_str());
	fprintf(f, "T>%s\n", AD.m_LabelT.c_str());
	fprintf(f, "QL %u, QL2 %u, TL %u, NQ %u, overhang %+d\n",
	  QL, 2*QL, TL, NQ, (int) NQ - (int) QL);
	fprintf(f, "\n");
	fprintf(f, "q  %s\n", AD.m_SeqQ.c_str());
	fprintf(f, "Q  %s\n", QRow.c_str());
	fprintf(f, "   %s\n", AnnotRow.c_str());
	fprintf(f, "T  %s\n", TRow.c_str());
	fprintf(f, "   %s\n", PathRow.c_str());

	return true;
	}
