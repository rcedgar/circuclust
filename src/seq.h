/* (C) Copyright 2016 Robert C. Edgar, all rights reserved. */

#ifndef seq_h
#define seq_h

struct ORFData;

struct SeqData
	{
	string Label;
	string Seq;

	void RevCompInPlace();

	const byte *GetSeq() const
		{
		return (const byte *) Seq.c_str();
		}

	const uint GetL() const
		{
		return SIZE(Seq);
		}
	};

#endif // seq_h
