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
void RotateLongorf(string &QuerySeq);

SeqDB *g_RefDB;

static SeqDB g_DB;
static FILE *g_fTsv;
static FILE *g_fFa;
static FILE *g_fAln;
static FILE *g_fOv;
static bool g_Explode;
static string g_ExplodeDir;

static const uint MINL = 32;

static void OnHit(const AlnData &AD, uint TargetSeqIndex)
	{
	if (g_fTsv != 0)
		{
		fprintf(g_fTsv, "H");
		fprintf(g_fTsv, "\t%s", AD.m_LabelQ.c_str());
		fprintf(g_fTsv, "\t%s", AD.m_LabelT.c_str());
		fprintf(g_fTsv, "\t%.1f", 100.0*AD.m_FractId);
		fprintf(g_fTsv, "\n");
		}

	if (g_fAln != 0)
		AD.WriteAln(g_fAln);

	if (g_Explode)
		{
		string FileName = g_ExplodeDir;
		Psa(FileName, "/cluster%u", TargetSeqIndex);
		FILE *f = fopen(FileName.c_str(), "a");
		asserta(f != 0);
		AD.WriteFasta(f);
		CloseStdioFile(f);
		}
	}

static void OnNewCluster(const string &QueryLabel, 
  const string &QuerySeq, uint TargetSeqIndex)
	{
	if (g_fTsv != 0)
		{
		fprintf(g_fTsv, "C");
		fprintf(g_fTsv, "\t%s", QueryLabel.c_str());
		fprintf(g_fTsv, "\t.");
		fprintf(g_fTsv, "\t.");
		fprintf(g_fTsv, "\n");
		}

	if (g_Explode)
		{
		string FileName = g_ExplodeDir;
		Psa(FileName, "/cluster%u", TargetSeqIndex);
		FILE *f = CreateStdioFile(FileName);
		SeqToFasta(f, QueryLabel, QuerySeq);
		CloseStdioFile(f);
		}
	}

void cmd_cluster()
	{
	SeqDB Input;
	Input.FromFasta(opt(cluster));

	if (!optset_id)
		Die("-id required");

	asserta(!optset_output);
	double MinFractId = opt(id);
	if (MinFractId < 0 || MinFractId > 1)
		Die("-id must be in range 0.0 to 1.0");

	g_Explode = optset_explode;
	if (g_Explode)
		g_ExplodeDir = string(opt(explode));

	g_fTsv = CreateStdioFile(opt(tsvout));
	g_fFa = CreateStdioFile(opt(fastaout));
	g_fAln = CreateStdioFile(opt(alnout));

	if (optset_overhangs)
		{
		g_fOv = CreateStdioFile(opt(overhangs));
		setbuf(g_fOv, 0);
		}
	const uint MaxOv = (optset_maxoverhangs ? opt(maxoverhangs) : 50);

	g_RefDB = &g_DB;

	SetSubstMx(Input);
	const unsigned QuerySeqCount = Input.GetSeqCount();

	uint ShortCount = 0;
	uint ClusterCount = 0;
	uint HitCount = 0;
	uint OverhangCount = 0;
	for (unsigned QuerySeqIndex = 0; QuerySeqIndex < QuerySeqCount; ++QuerySeqIndex)
		{
		if (optset_overhangs)
			ProgressStep(QuerySeqIndex, QuerySeqCount, "Clustering, %u overhangs", OverhangCount);
		else
			ProgressStep(QuerySeqIndex, QuerySeqCount, "Clustering %u clusters %u hits",
			  ClusterCount, HitCount);

		string QueryLabel;
		string QuerySeq;
		Input.GetStrs(QuerySeqIndex, QueryLabel, QuerySeq);

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
		USortS(QuerySeq, MinTL, MaxTL, g_DB, WordCounts, Order);

		uint TopWordCount = UINT_MAX;
		const uint N = SIZE(Order);
		if (N > 0)
			TopWordCount = WordCounts[Order[0]];
		bool AssignedToCluster = false;
		for (uint i = 0; i < N; ++i)
			{
			uint TargetSeqIndex = Order[i];
			uint WordCount = WordCounts[TargetSeqIndex];
			if (WordCount < TopWordCount/2)
				break;

			string TargetLabel;
			string TargetSeq;
			g_DB.GetStrs(TargetSeqIndex, TargetLabel, TargetSeq);

			AlnData AD;
			AlignPair(QueryLabel, QuerySeq, TargetLabel, TargetSeq, AD);

			if (AD.m_FractId >= MinFractId)
				{
				++HitCount;
				AssignedToCluster = true;
				OnHit(AD, TargetSeqIndex);
				if (optset_overhangs)
					{
					bool Ov = ReportOverhangs(g_fOv, AD);
					if (Ov)
						{
						++OverhangCount;
						if (OverhangCount > MaxOv)
							{
							ProgressLog("\n\n>%u overhangs\n", MaxOv);
							exit(0);
							}
						}
					}
				break;
				}
			}
		if (!AssignedToCluster)
			{
			if (opt(rotate_longorf))
				RotateLongorf(QuerySeq);
			OnNewCluster(QueryLabel, QuerySeq, ClusterCount);
			++ClusterCount;
			g_DB.AddSeq(QueryLabel.c_str(), (const byte *) QuerySeq.c_str(), SIZE(QuerySeq));
			SeqToFasta(g_fFa, QueryLabel, QuerySeq.c_str());
			}
		}
	if (ShortCount > 0)
		Warning("%u short sequences < %unt discarded", ShortCount, MINL);

	CloseStdioFile(g_fTsv);
	CloseStdioFile(g_fFa);
	CloseStdioFile(g_fAln);
	}
