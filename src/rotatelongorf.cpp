#include "myutils.h"
#include "translator.h"

void RotateLongorf(string &Seq)
	{
	Translator T;
	T.SetStr(Seq, true);

	vector<int> EndlessFrames;
	vector<bool> HasMets;
	T.GetEndlessFrames(EndlessFrames, HasMets);
	if (SIZE(EndlessFrames) > 0)
		{
		int Frame = EndlessFrames[0];
		T.SetFrame(Frame);

		string NtSeq_EndlessFrame;
		T.GetNtStr(NtSeq_EndlessFrame);

		Seq = NtSeq_EndlessFrame;
		return;
		}

	int Frame;
	uint OrfIndex;
	T.GetLongestOrf(Frame, OrfIndex);

	if (Frame == 0)
		return;

	T.SetFrame(Frame);
	asserta(OrfIndex  < SIZE(T.m_OrfAaLengths));
	uint CodonCount = T.m_OrfAaLengths[OrfIndex];
	uint AaPos = T.m_OrfAaStarts[OrfIndex];
	T.GetNtSubStr(AaPos, CodonCount, Seq);
	}
