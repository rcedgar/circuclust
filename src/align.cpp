#include "myutils.h"
#include "seqdb.h"
#include "xdpmem.h"
#include "objmgr.h"
#include "seq.h"
#include "alndata.h"

void AlignPair(const string &LabelQ, const string &SeqQ,
  const string &LabelT, const string &SeqT, AlnData &AD);
void SetSubstMx(const SeqDB &DB);

void cmd_align()
	{
	SeqDB Query;
	SeqDB DB;
	Query.FromFasta(opt(align));
	DB.FromFasta(opt(db));
	FILE *fOut = CreateStdioFile(opt(output));
	FILE *fTsv = CreateStdioFile(opt(tsvout));
	FILE *fFa = CreateStdioFile(opt(rotated));
	FILE *fAln = CreateStdioFile(opt(alnout));

	SetSubstMx(DB);

	const unsigned QuerySeqCount = Query.GetSeqCount();
	const unsigned DBSeqCount = DB.GetSeqCount();

	for (unsigned i = 0; i < QuerySeqCount; ++i)
		{
		ProgressStep(i, QuerySeqCount, "Aligning");

		string LabelQ, SeqQ;
		Query.GetStrs(i, LabelQ, SeqQ);

		for (unsigned j = 0; j < DBSeqCount; ++j)
			{
			string LabelT, SeqT;
			DB.GetStrs(j, LabelT, SeqT);

			AlnData AD;
			AlignPair(LabelQ, SeqQ, LabelT, SeqT, AD);
			AD.WriteAln(fOut);
			AD.WriteTsv(fTsv);
			AD.WriteFasta(fFa);
			AD.WriteAln(fAln);
			}
		}

	CloseStdioFile(fOut);
	CloseStdioFile(fTsv);
	CloseStdioFile(fAln);
	CloseStdioFile(fFa);
	}
