#include "myutils.h"
#include "seq.h"
#include "alndata.h"

void SeqToFasta(FILE *f, const char *Label, const byte *Seq, unsigned L);
void SeqToFasta(FILE *f, const string &Label, const string &Seq);

void AlnData::WriteFasta(FILE *f) const
	{
	if (f == 0)
		return;

	string RotatedQ;
	StripGapsFromQRow(RotatedQ);
	SeqToFasta(f, m_LabelQ, RotatedQ);
	}

void AlnData::WriteTsv(FILE *f) const
	{
	if (f == 0)
		return;
	fprintf(f, "%s", m_LabelQ.c_str());
	fprintf(f, "\t%s", m_LabelT.c_str());
	fprintf(f, "\t%u", m_StartPosQ);
	fprintf(f, "\t%c", pom(m_Plus));
	fprintf(f, "\t%.1f", 100.0*m_FractId);
	fprintf(f, "\n");
	}

void AlnData::StripGapsFromQRow(string &Q) const
	{
	Q.clear();
	const uint ColCount = SIZE(m_QRow);
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char q = toupper(m_QRow[Col]);
		if (q != '-')
			Q += q;
		}
	}

void AlnData::WriteAln(FILE *f) const
	{
	if (f == 0)
		return;

	const uint ColCount = SIZE(m_QRow);
	asserta(SIZE(m_TRow) == ColCount);

	fprintf(f, "\n");
	fprintf(f, "Q  >%s (%u nt) start %u strand %c\n",
	 m_LabelQ.c_str(), SIZE(m_SeqQ), m_StartPosQ + 1, pom(m_Plus));
	fprintf(f, "T  >%s (%u nt)\n", m_LabelT.c_str(), SIZE(m_SeqT));
	uint IdCount = 0;
	uint GapCount = 0;
	string AnnotRow;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char q = toupper(m_QRow[Col]);
		char t = toupper(m_TRow[Col]);
		if (q == '-' || t == '-')
			++GapCount;
		else if (q == t)
			++IdCount;

		if (q == t)
			AnnotRow += '|';
		else
			AnnotRow += ' ';
		}

	const uint BLOCK_SIZE = 80;

	const char *QR = m_QRow.c_str();
	const char *TR = m_TRow.c_str();
	const char *AR = AnnotRow.c_str();

	uint ColFrom = 0;
	for (;;)
		{
		if (ColFrom >= ColCount)
			break;
		uint n = ColCount - ColFrom;
		if (n > BLOCK_SIZE)
			n = BLOCK_SIZE;
		fprintf(f, "\n");
		fprintf(f, "Q  %*.*s\n", n, n, QR + ColFrom);
		fprintf(f, "   %*.*s\n", n, n, AR + ColFrom);
		fprintf(f, "T  %*.*s\n", n, n, TR + ColFrom);

		ColFrom += BLOCK_SIZE;
		}

	double PctId = double(IdCount*100)/ColCount;
	double GapPct = double(GapCount*100)/ColCount;

	fprintf(f, "%u/%u ids (%.1f%%), %u gaps (%.1f%%)\n",
	  IdCount, ColCount, PctId, GapCount, GapPct);
	}
