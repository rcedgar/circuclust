/* (C) Copyright 2016 Robert C. Edgar, all rights reserved. */

#include "sfasta.h"
#include "alpha.h"

void ReadFasta(const string &FileName, fn_OnSeq OnSeq, void *UserData)
	{
	string Label;
	string Seq;

	FILE *f = OpenStdioFile(FileName);
	const uint64 PROGSIZE = 10*1024*1024;
	uint64 FileSize = GetStdioFileSize64(f);
	if (FileSize > PROGSIZE)
		ProgressStep(0, 1001, "Reading %s", FileName.c_str());

	uint64 SumL = 0;
	string Line;
	uint64 LastPos = 0;
	uint SeqIndex = 0;
	while (ReadLineStdioFile(f, Line))
		{
		if (Line[0] == '>')
			{
 			if (!Seq.empty())
				{
				SumL += SIZE(Seq);
				OnSeq(Label, Seq, UserData);
				}
			Label = Line.substr(1, string::npos);
			Seq.clear();
			if (FileSize > PROGSIZE && SeqIndex%100 == 0)
				{
				uint64 Pos = GetStdioFilePos64(f);
				if (Pos != LastPos)
					{
					uint64 n = (Pos*1000)/FileSize;
					if (n > 1001)
						n = 1001;
					ProgressStep((uint) n, 1002, "Reading %s (%s)", FileName.c_str(), MemBytesToStr(SumL));
					LastPos = Pos;
					}
				}
			++SeqIndex;
			}
		else
			{
			for (uint i = 0; i < SIZE(Line); ++i)
				{
				char c = Line[i];
				if (isalpha(c) || c == '*')
					Seq += c;
				else if (isspace(c))
					continue;
				else
					Warning("Invalid byte 0x%02x in FASTA sequence data", c);
				}
			}
		}
	if (FileSize > PROGSIZE)
		ProgressStep(1001, 1002, "Reading %s (%s)", FileName.c_str(), MemBytesToStr(SumL));

	if (!Seq.empty())
		OnSeq(Label, Seq, UserData);

	CloseStdioFile(f);
	}
