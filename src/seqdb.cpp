/* (C) Copyright 2016 Robert C. Edgar, all rights reserved. */

#include "myutils.h"
#include "seqdb.h"
#include "alpha.h"
#include "sfasta.h"
#include "seq.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "fastaseqsource.h"

void SeqToFasta(FILE *f, const char *Label, const byte *Seq, unsigned L)
	{
	if (f == 0)
		return;

	const unsigned ROWLEN = 80;
	if (Label != 0)
		fprintf(f, ">%s\n", Label);
	unsigned BlockCount = (L + ROWLEN - 1)/ROWLEN;
	for (unsigned BlockIndex = 0; BlockIndex < BlockCount; ++BlockIndex)
		{
		unsigned From = BlockIndex*ROWLEN;
		unsigned To = From + ROWLEN;
		if (To >= L)
			To = L;
		for (unsigned Pos = From; Pos < To; ++Pos)
			fputc(Seq[Pos], f);
		fputc('\n', f);
		}
	}

void StrToFasta(FILE *f, const string &Label, const string &Seq)
	{
	SeqToFasta(f, Label.c_str(), (const byte *) Seq.c_str(), SIZE(Seq));
	}

void SeqToFasta(FILE *f, const string &Label, const string &Seq)
	{
	SeqToFasta(f, Label.c_str(), (const byte *) Seq.c_str(), SIZE(Seq));
	}

SeqDB::~SeqDB()
	{
	Clear();
	}

SeqDB::SeqDB()
	{
	Clear(true);
	}

void SeqDB::Clear(bool ctor)
	{
	if (!ctor)
		{
		for (unsigned i = 0; i < m_SeqCount; ++i)
			{
			unsigned n = strlen(m_Labels[i]);
			MYFREE(m_Labels[i], n, SeqDB);
			MYFREE(m_Seqs[i], m_SeqLengths[i], SeqDB);
			}
		MYFREE(m_Labels, m_Size, SeqDB);
		MYFREE(m_Seqs, m_Size, SeqDB);
		MYFREE(m_SeqLengths, m_Size, SeqDB);
		}

	m_FileName.clear();
	m_SeqCount = 0;
	m_Size = 0;

	m_Labels = 0;
	m_Seqs = 0;
	m_SeqLengths = 0;

	m_Aligned = false;
	m_IsNucleo = false;
	m_IsNucleoSet = false;
	}

void SeqDB::InitEmpty(bool Nucleo)
	{
	Clear();
	m_IsNucleo = Nucleo;
	m_IsNucleoSet = true;
	}

//static void OnSeq(const string &Label, const string &Seq, void *vpDB)
//	{
//	SeqDB *DB = (SeqDB *) vpDB;
//	uint L = SIZE(Seq);
//	uint N = DB->GetSeqCount();
//	for (uint i = 0; i < N; ++i)
//		{
//		string Label;
//		DB->GetLabelStr(i, Label);
//		}
//	DB->AddSeq(Label.c_str(), (const byte *) Seq.c_str(), L);
//	for (uint i = 0; i < N+1; ++i)
//		{
//		string Label;
//		DB->GetLabelStr(i, Label);
//		}
//	}

void SeqDB::FromFasta(const string &FileName, bool AllowGaps)
	{
	Clear();
	m_FileName = FileName;

	//ReadFasta(FileName, OnSeq, (void *) this);
	FASTASeqSource SS;
	SS.Open(FileName);

	SeqInfo *SI = ObjMgr::GetSeqInfo();
	uint64 FileSize = GetStdioFileSize64(SS.m_LR.m_f);
	uint SeqCount = 0;
	ProgressStep(0, 1002, "Reading %s", FileName.c_str());
	uint LastMills = 0;
	while (SS.GetNext(SI))
		{
		if (SeqCount++%1000 == 0)
			{
			uint64 Pos = GetStdioFilePos64(SS.m_LR.m_f);
			uint Mills = uint(Pos*1000.0/double(FileSize));
			if (Mills > LastMills)
				{
				ProgressStep(Mills, 1002, "Reading %s", FileName.c_str());
				LastMills = Mills;
				}
			}

		AddSeq(SI->m_Label, SI->m_Seq, SI->m_L);
		}
	ProgressStep(1001, 1002, "Reading %s", FileName.c_str());

	SS.Close();
	}

void SeqDB::ToFasta(const string &FileName) const
	{
	FILE *f = CreateStdioFile(FileName);
	for (unsigned SeqIndex = 0; SeqIndex < GetSeqCount(); ++SeqIndex)
		ToFasta(f, SeqIndex);
	CloseStdioFile(f);
	}

void SeqDB::SeqToFasta(FILE *f, unsigned SeqIndex, bool WithLabel) const
	{
	if (WithLabel)
		fprintf(f, ">%s\n", GetLabel(SeqIndex));

	const unsigned ROWLEN = 80;

	unsigned L = GetSeqLength(SeqIndex);
	const byte *Seq = GetSeq(SeqIndex);
	unsigned BlockCount = (L + ROWLEN - 1)/ROWLEN;
	for (unsigned BlockIndex = 0; BlockIndex < BlockCount; ++BlockIndex)
		{
		unsigned From = BlockIndex*ROWLEN;
		unsigned To = From + ROWLEN;
		if (To >= L)
			To = L;
		for (unsigned Pos = From; Pos < To; ++Pos)
			fputc(Seq[Pos], f);
		fputc('\n', f);
		}
	}

void SeqDB::ToFasta(FILE *f, unsigned SeqIndex) const
	{
	asserta(SeqIndex < m_SeqCount);
	fprintf(f, ">%s\n", GetLabel(SeqIndex));
	SeqToFasta(f, SeqIndex);
	}

unsigned SeqDB::GetMaxLabelLength() const
	{
	const unsigned SeqCount = GetSeqCount();
	unsigned MaxL = 0;
	for (unsigned Index = 0; Index < SeqCount; ++Index)
		{
		unsigned L = (unsigned) strlen(m_Labels[Index]);
		if (L > MaxL)
			MaxL = L;
		}
	return MaxL;
	}

unsigned SeqDB::GetMaxSeqLength() const
	{
	const unsigned SeqCount = GetSeqCount();
	unsigned MaxL = 0;
	for (unsigned Index = 0; Index < SeqCount; ++Index)
		{
		unsigned L = m_SeqLengths[Index];
		if (L > MaxL)
			MaxL = L;
		}
	return MaxL;
	}

void SeqDB::LogMe() const
	{
	Log("\n");
	const unsigned SeqCount = GetSeqCount();
	Log("SeqDB %u seqs, aligned=%c\n", SeqCount, tof(m_Aligned));
	if (SeqCount == 0)
		return;

	Log("Index             Label  Length  Seq\n");
	Log("-----  ----------------  ------  ---\n");
	for (unsigned Index = 0; Index < SeqCount; ++Index)
		{
		Log("%5u", Index);
		Log("  %16.16s", m_Labels[Index]);
		unsigned L = m_SeqLengths[Index];
		Log("  %6u", L);
		Log("  %*.*s", L, L, m_Seqs[Index]);
		Log("\n");
		}
	}

bool SeqDB::GuessIsNucleo() const
	{
	const unsigned SeqCount = GetSeqCount();
	unsigned N = 0;
	for (unsigned i = 0; i < 100; ++i)
		{
		unsigned SeqIndex = unsigned(rand()%SeqCount);
		const byte *Seq = GetSeq(SeqIndex);
		unsigned L = GetSeqLength(SeqIndex);
		const unsigned Pos = unsigned(rand()%L);
		byte c = Seq[Pos];

		if (g_IsNucleoChar[c])
			++N;
		}
	return (N > 80);
	}

void SeqDB::SetIsNucleo()
	{
	const unsigned SeqCount = GetSeqCount();
	unsigned N = 0;
	for (unsigned i = 0; i < 100; ++i)
		{
		unsigned SeqIndex = unsigned(rand()%SeqCount);
		const byte *Seq = GetSeq(SeqIndex);
		unsigned L = GetSeqLength(SeqIndex);
		const unsigned Pos = unsigned(rand()%L);
		byte c = Seq[Pos];

		if (g_IsNucleoChar[c])
			++N;
		}
	m_IsNucleo = (N > 80);
	m_IsNucleoSet = true;
	}

unsigned SeqDB::GetTotalLength() const
	{
	const unsigned SeqCount = GetSeqCount();
	unsigned TotalLength = 0;
	for (unsigned Id = 0; Id < SeqCount; ++Id)
		TotalLength += GetSeqLength(Id);
	return TotalLength;
	}

void SeqDB::GetStrs(uint SeqIndex, string &Label, string &Seq) const
	{
	GetLabelStr(SeqIndex, Label);
	GetSeqStr(SeqIndex, Seq);
	}

void SeqDB::GetSeqStr(uint SeqIndex, string &Seq) const
	{
	Seq.clear();
	const byte *ByteSeq = GetSeq(SeqIndex);
	const uint L = GetSeqLength(SeqIndex);
	for (uint i = 0; i < L; ++i)
		Seq += ByteSeq[i];
	}

void SeqDB::GetLabelStr(uint SeqIndex, string &Label) const
	{
	Label.clear();
	const char *CharLabel = GetLabel(SeqIndex);
	const uint L = ustrlen(CharLabel);
	for (uint i = 0; i < L; ++i)
		Label += CharLabel[i];
	}

unsigned SeqDB::AddSeq(const char *Label, const byte *Seq, unsigned L)
	{
	StartTimer(AddSeq);
	if (m_SeqCount >= m_Size)
		{
		unsigned NewSize = unsigned(m_Size*1.5) + 1024;
		char **NewLabels = MYALLOC(char *, NewSize, SeqDB);
		byte **NewSeqs = MYALLOC(byte *, NewSize, SeqDB);
		unsigned *NewSeqLengths = MYALLOC(unsigned, NewSize, SeqDB);

		for (unsigned i = 0; i < m_SeqCount; ++i)
			{
			NewLabels[i] = m_Labels[i];
			NewSeqs[i] = m_Seqs[i];
			NewSeqLengths[i] = m_SeqLengths[i];
			}

		MYFREE(m_Labels, m_SeqCount, SeqDB);
		MYFREE(m_Seqs, m_SeqCount, SeqDB);
		MYFREE(m_SeqLengths, m_SeqCount, SeqDB);

		m_Labels = NewLabels;
		m_Seqs = NewSeqs;
		m_SeqLengths = NewSeqLengths;
		m_Size = NewSize;
		}

	unsigned Index = m_SeqCount++;
	m_Seqs[Index] = MYALLOC(byte, L, SeqDB);
	memcpy(m_Seqs[Index], Seq, L);

	unsigned n = strlen(Label) + 1;
	m_Labels[Index] = MYALLOC(char, n, SeqDB);
	memcpy(m_Labels[Index], Label, n);

	if (Index == 0)
		m_Aligned = true;
	else
		m_Aligned = (m_Aligned && L == m_SeqLengths[0]);

	m_SeqLengths[Index] = L;
	EndTimer(AddSeq);
	return Index;
	}

unsigned SeqDB::GetIndex(const char *Label) const
	{
	for (unsigned i = 0; i < m_SeqCount; ++i)
		if (strcmp(Label, m_Labels[i]) == 0)
			return i;
	Die("SeqDB::GetIndex(%s), not found", Label);
	return UINT_MAX;
	}

void SeqDB::MakeLabelToIndex(map<string, unsigned> &LabelToIndex)
	{
	LabelToIndex.clear();
	for (unsigned i = 0; i < m_SeqCount; ++i)
		{
		const string &Label = string(GetLabel(i));
		if (LabelToIndex.find(Label) != LabelToIndex.end())
			Die("Duplicate label: %s", Label.c_str());
		LabelToIndex[Label] = i;
		}
	}

void SeqDB::FromSeqDBSubset(const SeqDB &DB, const unsigned *SeqIndexes, unsigned N)
	{
	m_SeqCount = 0;
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned SeqIndex = SeqIndexes[i];
		const char *Label = DB.GetLabel(SeqIndex);
		const byte *Seq = DB.GetSeq(SeqIndex);
		unsigned L = DB.GetSeqLength(SeqIndex);

		AddSeq(Label, Seq, L);
		}
	}
