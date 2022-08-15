/* (C) Copyright 2016 Robert C. Edgar, all rights reserved. */

#include "myutils.h"
#include "seqdb.h"
#include "seq.h"
#include "alpha.h"
#include "sort.h"
#include "partition.h"
#include <algorithm>

extern unsigned opt_w;
extern bool g_IsNucleo;
extern const byte *g_CharToLetter;
extern uint g_AlphaSize;

const vector<float> *g_SortVecFloat;
static byte *g_QueryHasWord;
static unsigned g_DictSize;
unsigned g_QueryUniqueWordCount;
static PartMem *g_Mem;

unsigned GetWord(const byte *Seq)
	{
	unsigned Word = 0;
	const byte *Front = Seq;
	for (unsigned i = 0; i < opt_w; ++i)
		{
		unsigned Letter = g_CharToLetter[*Front++];
		if (Letter >= g_AlphaSize)
			return UINT_MAX;
		Word = (Word*g_AlphaSize) + Letter;
		assert(Word < g_DictSize);
		}
	return Word;
	}

static bool CmpDescVecFloat(unsigned i, unsigned j)
	{
	return (*g_SortVecFloat)[i] > (*g_SortVecFloat)[j];
	}

void GetQueryUniqueWords(const SeqData &Query, vector<unsigned> &Words)
	{
	Words.clear();

	if (g_QueryHasWord == 0)
		{
		g_DictSize = g_AlphaSize;
		for (unsigned i = 1; i < opt_w; ++i)
			g_DictSize *= g_AlphaSize;

		g_QueryHasWord = myalloc(byte, g_DictSize);
		}

	memset(g_QueryHasWord, 0, g_DictSize);

	if (SIZE(Query.Seq) <= opt_w)
		return;

	const unsigned L = SIZE(Query.Seq) - opt_w + 1;
	const byte *Seq = (const byte *) Query.Seq.c_str();
	for (unsigned i = 0; i < L; ++i)
		{
		unsigned Word = GetWord(Seq++);
		if (Word == UINT_MAX)
			continue;
		if (g_QueryHasWord[Word] == 0)
			++g_QueryUniqueWordCount;
		assert(Word < g_DictSize);
		g_QueryHasWord[Word] = 1;
		}

	for (unsigned i = 0; i < g_DictSize; ++i)
		{
		if (g_QueryHasWord[i])
			Words.push_back(i);
		}
	}

void SortDescending(const vector<float> &Values, vector<unsigned> &Order)
	{
	const unsigned N = SIZE(Values);
	Range(Order, N);
	g_SortVecFloat = &Values;
	sort(Order.begin(), Order.end(), CmpDescVecFloat);
	}

static void SetQuery(const string &Query)
	{
	if (g_QueryHasWord == 0)
		{
		g_DictSize = g_AlphaSize;
		for (unsigned i = 1; i < opt_w; ++i)
			g_DictSize *= g_AlphaSize;

		g_QueryHasWord = myalloc(byte, g_DictSize);
		}

	memset(g_QueryHasWord, 0, g_DictSize);

	if (SIZE(Query) <= opt_w)
		return;

	const unsigned L = SIZE(Query) - opt_w + 1;
	const byte *Seq = (const byte *) Query.c_str();
	for (unsigned i = 0; i < L; ++i)
		{
		unsigned Word = GetWord(Seq++);
		if (Word == UINT_MAX)
			continue;
		assert(Word < g_DictSize);
		if (g_QueryHasWord[Word] == 0)
			++g_QueryUniqueWordCount;
		assert(Word < g_DictSize);
		g_QueryHasWord[Word] = 1;
		}
	}

static unsigned GetUniqueWordsInCommon(const string &Target)
	{
	if (SIZE(Target) <= opt_w)
		return 0;

	unsigned Count = 0;
	const unsigned L = SIZE(Target) - opt_w + 1;
	const byte *Seq = (const byte *) Target.c_str();
	for (unsigned i = 0; i < L; ++i)
		{
		unsigned Word = GetWord(Seq++);
		if (Word == UINT_MAX)
			continue;
		assert(Word < g_DictSize);
		if (g_QueryHasWord[Word])
			++Count;
		}
	return Count;
	}

void USortS(const string &Query, uint MinTargetLength, uint MaxTargetLength, 
  const SeqDB &DB, vector<unsigned> &WordCounts, vector<unsigned> &Order)
	{
	//CheckDB();

	WordCounts.clear();
	Order.clear();

	//CheckDB();
	SetQuery(Query);
	//CheckDB();

	const unsigned SeqCount = DB.GetSeqCount();
	float MaxU = 0;
	float NextU = 0;
	vector<float> WordCountsf;
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		//CheckDB();
		uint TargetLength = DB.GetSeqLength(SeqIndex);
		if (TargetLength < MinTargetLength || TargetLength > MaxTargetLength)
			{
			WordCountsf.push_back(0);
			continue;
			}
		string Target;
		DB.GetSeqStr(SeqIndex, Target);
		//CheckDB();
		unsigned WordCount = GetUniqueWordsInCommon(Target);
		//CheckDB();
		WordCountsf.push_back(float(WordCount));
		//CheckDB();
		}
	//CheckDB();
	SortDescending(WordCountsf, Order);
	for (unsigned i = 0; i < SeqCount; ++i)
		WordCounts.push_back(unsigned(WordCountsf[i]));
	//CheckDB();
	}

void USort(const string &Query, const SeqDB &DB, vector<unsigned> &WordCounts, 
  vector<unsigned> &Order)
	{
	WordCounts.clear();
	Order.clear();

	SetQuery(Query);

	const unsigned SeqCount = DB.GetSeqCount();
	float MaxU = 0;
	float NextU = 0;
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		string Target;
		DB.GetSeqStr(SeqIndex, Target);
		unsigned WordCount = GetUniqueWordsInCommon(Target);
		WordCounts.push_back(WordCount);
		}
//	SortDescending(WordCounts, Order);

	unsigned N = SIZE(WordCounts);
	unsigned *v = myalloc(unsigned, N);
	if (g_Mem == 0)
		g_Mem = new PartMem;
	unsigned m = PartOrderDesc(WordCounts.data(), SIZE(WordCounts), *g_Mem, v);
	for (unsigned i = 0; i < m; ++i)
		Order.push_back(v[i]);
	myfree(v);
	}

//unsigned GetUS(const string &Query, unsigned *ptrTargetIndexes, unsigned *ptrWordCounts)
//	{
//	extern SeqDB *g_RefDB;
//
//	vector<unsigned> WordCounts;
//	vector<unsigned> Order;
//
//	g_QueryUniqueWordCount = 0;
//	USortS(Query, *g_RefDB, WordCounts, Order);
//	unsigned n = SIZE(Order);
//	for (unsigned i = 0; i < n; ++i)
//		{
//		unsigned k = Order[i];
//		ptrTargetIndexes[i] = k;
//		unsigned WordCount = (unsigned) WordCounts[k];
//		ptrWordCounts[i] = WordCount;
//		}
//#if	0
//	{
//	extern SeqDB *g_RefDB;
//	Log("\n");
//	Log("Q>%s\n", Query->Label);
//	Log("%u hot\n", n);
//	for (unsigned k = 0; k < n; ++k)
//		{
//		unsigned TargetIndex = ptrTargetIndexes[k];
//		unsigned WordCount = ptrWordCounts[k];
//		const char *Label = g_RefDB->GetLabel(TargetIndex);
//		Log("%7u  %7u  %s\n", TargetIndex, WordCount, Label);
//		}
//	}
//#endif
//	return n;
//	}

unsigned GetU(const string &Query, unsigned *ptrTargetIndexes, unsigned *ptrWordCounts)
	{
	extern SeqDB *g_RefDB;

	vector<unsigned> WordCounts;
	vector<unsigned> Order;

	g_QueryUniqueWordCount = 0;
	USort(Query, *g_RefDB, WordCounts, Order);
	unsigned n = SIZE(Order);
	for (unsigned i = 0; i < n; ++i)
		{
		unsigned k = Order[i];
		ptrTargetIndexes[i] = k;
		unsigned WordCount = (unsigned) WordCounts[k];
		ptrWordCounts[i] = WordCount;
		}
#if	0
	{
	extern SeqDB *g_RefDB;
	Log("\n");
	Log("Q>%s\n", Query->Label);
	Log("%u hot\n", n);
	for (unsigned k = 0; k < n; ++k)
		{
		unsigned TargetIndex = ptrTargetIndexes[k];
		unsigned WordCount = ptrWordCounts[k];
		const char *Label = g_RefDB->GetLabel(TargetIndex);
		Log("%7u  %7u  %s\n", TargetIndex, WordCount, Label);
		}
	}
#endif
	return n;
	}

unsigned GetQueryUniqueWordCount(SeqData *Query)
	{
	unsigned L = SIZE(Query->Seq);
	if (L <= opt_w)
		return 1;
	return L - opt_w + 1;
	}

void ZeroDict()
	{
	memset(g_QueryHasWord, 0, g_DictSize);
	}

unsigned GetTopTargetIndex(const vector<unsigned> &Words, unsigned n)
	{
	extern SeqDB *g_RefDB;
//	memset(g_QueryHasWord, 0, g_DictSize);
	for (unsigned i = 0; i < n; ++i)
		g_QueryHasWord[Words[i]] = true;

	const unsigned SeqCount = g_RefDB->GetSeqCount();
	unsigned MaxU = 0;
	unsigned TopTargetIndex = UINT_MAX;
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		string Target;
		g_RefDB->GetSeqStr(SeqIndex, Target);
		unsigned U = GetUniqueWordsInCommon(Target);
		if (U > MaxU)
			{
			MaxU = U;
			TopTargetIndex = SeqIndex;
			}
		}
	for (unsigned i = 0; i < n; ++i)
		g_QueryHasWord[Words[i]] = false;
	return TopTargetIndex;
	}
