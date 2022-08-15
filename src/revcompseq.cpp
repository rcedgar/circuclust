#include "myutils.h"
#include "alpha.h"
#include "seq.h"

void RevCompSeq(const byte *Seq, unsigned L, byte *RCSeq)
	{
	for (unsigned i = 0; i < L/2; ++i)
		{
		byte c1 = Seq[i];
		byte c2 = Seq[L-i-1];

		byte rc1 = g_CharToCompChar[c1];
		byte rc2 = g_CharToCompChar[c2];

		RCSeq[i] = rc2;
		RCSeq[L-i-1] = rc1;
		}
	if (L%2 != 0)
		{
		byte c = Seq[L/2];
		RCSeq[L/2] = g_CharToCompChar[c];
		}
	}

void RevCompSeqInPlace(byte *Seq, unsigned L)
	{
	RevCompSeq((const byte *) Seq, L, Seq);
	}

void SeqData::RevCompInPlace()
	{
	RevCompSeqInPlace((byte *) Seq.c_str(), SIZE(Seq));
	}

void RevCompInPlace(string &Seq)
	{
	RevCompSeqInPlace((byte *) Seq.c_str(), SIZE(Seq));
	}
