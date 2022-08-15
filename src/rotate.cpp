#include "myutils.h"
#include "seqdb.h"
#include "xdpmem.h"
#include "objmgr.h"
#include "seq.h"
#include "alndata.h"
#include "translator.h"

void AlignPair(const string &LabelQ, const string &SeqQ,
  const string &LabelT, const string &SeqT, AlnData &AD);
void SetSubstMx(const SeqDB &DB);
void RotateLongorf(string &Seq);
void StrToFasta(FILE *f, const string &Label, const string &Seq);

void cmd_rotate()
	{
	SeqDB Input;
	Input.FromFasta(opt(rotate));
	FILE *fOut = CreateStdioFile(opt(output));

	SetSubstMx(Input);

	const unsigned SeqCount = Input.GetSeqCount();

	string CentroidLabel;
	string CentroidSeq;

	Translator T;
	for (unsigned i = 0; i < SeqCount; ++i)
		{
		ProgressStep(i, SeqCount, "Rotating");

		if (i == 0)
			{
			Input.GetStrs(0, CentroidLabel, CentroidSeq);
			RotateLongorf(CentroidSeq);
			StrToFasta(fOut, CentroidLabel, CentroidSeq);
			continue;
			}

		string LabelQ, SeqQ;
		Input.GetStrs(i, LabelQ, SeqQ);

		AlnData AD;
		AlignPair(LabelQ, SeqQ, CentroidLabel, CentroidSeq, AD);

		string RotatedQ;
		AD.StripGapsFromQRow(RotatedQ);

		int Frame;
		uint OrfIndex;
		T.SetStr(RotatedQ, true);

		T.GetLongestOrf(Frame, OrfIndex);
		if (Frame == 1)
			{
			StrToFasta(fOut, LabelQ, RotatedQ);
			continue;
			}

		if (Frame == 2 || Frame == 3)
			{
			T.SetFrame(Frame);
			uint CodonCount = SIZE(SeqQ)/3;
			string NtStr;
			T.GetNtSubStr(0, CodonCount, NtStr);
			StrToFasta(fOut, LabelQ, NtStr);
			continue;
			}

		Log(">%s\n", LabelQ.c_str());
		Warning("Longest ORF in minus frame");
		}

	CloseStdioFile(fOut);
	}
