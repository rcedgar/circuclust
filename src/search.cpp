#include "myutils.h"
#include "seqdb.h"
#include "alndata.h"
#include "objmgr.h"
#include "seq.h"

void SetSubstMx(const SeqDB &DB);
void USortS(const string &Query, uint MinTargetLength, uint MaxTargetLength, 
  const SeqDB &DB, vector<unsigned> &WordCounts, vector<unsigned> &Order);
void AlignPair(const string &LabelQ, const string &SeqQ,
  const string &LabelT, const string &SeqT, AlnData &AD);
void SeqToFasta(FILE *f, const char *Label, const byte *Seq, unsigned L);
void SeqToFasta(FILE *f, const string &Label, const string &Seq);
bool ReportOverhangs(FILE *f, const AlnData &AD);

static FILE *g_fTsv;
static FILE *g_fRotatedFa;
static FILE *g_fOv;

static const uint MINL = 32;

static void OnHit(const AlnData &AD, uint TargetSeqIndex)
	{
	AD.WriteFasta(g_fRotatedFa);

	if (g_fTsv == 0)
		return;

	fprintf(g_fTsv, "H");
	fprintf(g_fTsv, "\t%s", AD.m_LabelQ.c_str());
	fprintf(g_fTsv, "\t%s", AD.m_LabelT.c_str());
	fprintf(g_fTsv, "\t%.1f", 100.0*AD.m_FractId);
	fprintf(g_fTsv, "\n");
	}

static void OnNoHit(const string &QueryLabel)
	{
	if (g_fTsv == 0)
		return;

	fprintf(g_fTsv, "N");
	fprintf(g_fTsv, "\t%s", QueryLabel.c_str());
	fprintf(g_fTsv, "\t.");
	fprintf(g_fTsv, "\t.");
	fprintf(g_fTsv, "\n");
	}

void cmd_search()
	{
	SeqDB Query;
	SeqDB DB;

	if (!optset_id)
		Die("-id required");

	double MinFractId = opt(id);
	if (MinFractId < 0 || MinFractId > 1)
		Die("-id must be in range 0.0 to 1.0");

	Progress("Reading query...");
	Query.FromFasta(opt(search));
	Progress(" done.\n");

	Progress("Reading db...");
	DB.FromFasta(opt(db));
	Progress(" done.\n");

	asserta(!optset_explode);
	asserta(!optset_overhangs);

	g_fTsv = CreateStdioFile(opt(tsvout));
	g_fRotatedFa = CreateStdioFile(opt(rotated));
	const uint MaxOv = (optset_maxoverhangs ? opt(maxoverhangs) : 50);

	SetSubstMx(Query);
	const unsigned QuerySeqCount = Query.GetSeqCount();

	uint ShortCount = 0;
	uint ClusterCount = 0;
	uint HitCount = 0;
	uint OverhangCount = 0;
	for (unsigned QuerySeqIndex = 0; QuerySeqIndex < QuerySeqCount; ++QuerySeqIndex)
		{
		double Pct = GetPct(HitCount, QuerySeqIndex+1);
		ProgressStep(QuerySeqIndex, QuerySeqCount, "Searching, %u hits (%.1f%%)",
		  HitCount, Pct);

		string QueryLabel;
		string QuerySeq;
		Query.GetStrs(QuerySeqIndex, QueryLabel, QuerySeq);

		const uint QL = SIZE(QuerySeq);
		if (QL < MINL)
			{
			++ShortCount;
			continue;
			}
		uint MaxDeltaL = uint(QL*(1 - MinFractId));
		if (MaxDeltaL < 8)
			MaxDeltaL = 8;
		uint MinTL = QL - MaxDeltaL;
		uint MaxTL = QL + MaxDeltaL;

		vector<uint> WordCounts;
		vector<uint> Order;
		USortS(QuerySeq, MinTL, MaxTL, DB, WordCounts, Order);

		uint TopWordCount = UINT_MAX;
		const uint N = SIZE(Order);
		if (N > 0)
			TopWordCount = WordCounts[Order[0]];

		bool HitFound = false;
		for (uint i = 0; i < N; ++i)
			{
			uint TargetSeqIndex = Order[i];
			uint WordCount = WordCounts[TargetSeqIndex];
			if (WordCount < TopWordCount/2)
				break;

			string TargetLabel;
			string TargetSeq;
			DB.GetStrs(TargetSeqIndex, TargetLabel, TargetSeq);

			AlnData AD;
			AlignPair(QueryLabel, QuerySeq, TargetLabel, TargetSeq, AD);

			if (AD.m_FractId >= MinFractId)
				{
				++HitCount;
				HitFound = true;
				OnHit(AD, TargetSeqIndex);
				break;
				}
			}
		if (!HitFound)
			OnNoHit(QueryLabel);
		}
	if (ShortCount > 0)
		Warning("%u short sequences < %unt discarded", ShortCount, MINL);

	CloseStdioFile(g_fTsv);
	CloseStdioFile(g_fRotatedFa);
	}
