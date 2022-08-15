#include "myutils.h"
#include "translator.h"
#include "alpha.h"

bool Translator::ValidFrame(int Frame) const
	{
	if (m_Circular && m_NtL%3 != 0)
		return Frame == +1 || Frame == -1;

	switch (Frame)
		{
	case -3: case -2: case -1: case +1: case +2: case +3:
		return true;
		}
	return false;
	}

void Translator::SetSeq(const byte *NtSeq, uint NtL, bool Circular)
	{
	m_NtSeq = NtSeq;
	m_NtL = NtL;
	m_Circular = Circular;
	SetFrame(1);
	}

void Translator::SetStr(const string &NtStr, bool Circular)
	{
	uint NtL = SIZE(NtStr);
	SetSeq((const byte *) NtStr.c_str(), NtL, Circular);
	}

// Nt position of first base in frame.
// Range is 0 .. (L-1)
uint Translator::CalcFirstNtPos(int Frame) const
	{
	if (Frame > 0)
		return uint(Frame-1);
	else
		{
		uint uFrame = uint(-Frame);
		return uint(m_NtL - uFrame);
		}
	}

uint Translator::CalcLastNtPos(int Frame) const
	{
	uint FirstNtPos = CalcFirstNtPos(Frame);
	uint CodonCount = CalcCodonCount(Frame);
	if (Frame > 0)
		{
		uint LastNtPos = FirstNtPos + 3*CodonCount - 1;
		asserta(LastNtPos < m_NtL);
		return LastNtPos;
		}
	else
		{
		uint n = 3*CodonCount - 1;
		asserta(FirstNtPos >= n);
		uint LastNtPos = FirstNtPos - n;
		return LastNtPos;
		}
	}

uint Translator::CalcCodonCount(int Frame) const
	{
	if (m_Circular)
		{
		if (m_NtL%3 == 0)
			return m_NtL/3;
		else
			return m_NtL;
		}

	uint n = uint(abs(Frame) - 1);
	uint N = (m_NtL - n)/3;
	return N;
	}

uint Translator::GetNtPos(uint AAPos, uint BaseNr) const
	{
	asserta(BaseNr < 3);
	uint n = 3*AAPos + BaseNr;
	if (m_Frame > 0)
		{
		uint NtPos = m_FirstNtPos + n;
		if (m_Circular)
			NtPos = NtPos%m_NtL;
		else
			asserta(NtPos < m_NtL);
		return NtPos;
		}
	else
		{
		if (m_Circular)
			{
			uint FirstNtPos = m_FirstNtPos;
			while (FirstNtPos < n)
				FirstNtPos += m_NtL;
			uint NtPos = FirstNtPos - n;
			return NtPos;
			}
		else
			{
			asserta(m_FirstNtPos >= n);
			uint NtPos = m_FirstNtPos - n;
			return NtPos;
			}
		}
	}

void Translator::GetCodonStr(uint AAPos, string &Codon) const
	{
	Codon.clear();
	Codon += (char) GetNt(AAPos, 0);
	Codon += (char) GetNt(AAPos, 1);
	Codon += (char) GetNt(AAPos, 2);
	}

byte Translator::GetAa(uint AAPos) const
	{
	byte b1 = GetNt(AAPos, 0);
	byte b2 = GetNt(AAPos, 1);
	byte b3 = GetNt(AAPos, 2);
	byte aa = TranslateCodon(b1, b2, b3);
	return aa;
	}

byte Translator::GetNt(uint AAPos, uint BaseNr) const
	{
	asserta(BaseNr < 3);
	uint NtPos = GetNtPos(AAPos, BaseNr);
	if (m_Frame > 0)
		return m_NtSeq[NtPos];
	else
		{
		byte c = m_NtSeq[NtPos];
		byte cc = g_CharToCompChar[c];
		return cc;
		}
	}

void Translator::SetFrame(int Frame, bool DoSetOrfs)
	{
	asserta(ValidFrame(Frame));
	m_Frame = Frame;
	m_FirstNtPos = CalcFirstNtPos(Frame);
	m_CodonCount = CalcCodonCount(Frame);
	if (DoSetOrfs)
		{
		uint MinCodonCount = 30;
		if (optset_mincodons)
			MinCodonCount = opt(mincodons);
		SetOrfs(MinCodonCount);
		}
	}

void Translator::GetAaStr(string &AaStr) const
	{
	AaStr.clear();
	asserta(m_Frame != 0 && m_CodonCount != UINT_MAX);
	for (uint AAPos = 0; AAPos < m_CodonCount; ++AAPos)
		{
		char aa = GetAa(AAPos);
		AaStr += aa;
		}
	}

void Translator::GetNtStr(string &NtStr) const
	{
	NtStr.clear();
	asserta(m_Frame != 0 && m_CodonCount != UINT_MAX);
	for (uint AAPos = 0; AAPos < m_CodonCount; ++AAPos)
		{
		for (uint k = 0; k < 3; ++k)
			{
			char nt = GetNt(AAPos, k);
			NtStr += nt;
			}
		}
	}

bool Translator::IsStop(uint AaPos) const
	{
	char aa = GetAa(AaPos);
	return aa == '*';
	}

bool Translator::IsMet(uint AaPos) const
	{
	char aa = GetAa(AaPos);
	return aa == 'M';
	}

void Translator::GetNtSubStr(uint AaStart, uint AaLength, string &NtStr) const
	{
	NtStr.clear();
	for (uint Offset = 0; Offset < AaLength; ++Offset)
		{
		for (uint BaseNr = 0; BaseNr < 3; ++BaseNr)
			NtStr += GetNt(AaStart + Offset, BaseNr);
		}
	}

void Translator::GetAaSubStr(uint AaStart, uint AaLength, string &AaStr) const
	{
	AaStr.clear();
	for (uint Offset = 0; Offset < AaLength; ++Offset)
		{
		char aa = GetAa(AaStart + Offset);
		AaStr += aa;
		}
	}

void Translator::SetOrfs(uint MinCodonCount)
	{
	m_OrfMinCodonCount = MinCodonCount;

	m_OrfAaStarts.clear();
	m_OrfAaLengths.clear();
	m_OrfHasStops.clear();
	m_OrfHasMets.clear();
	m_OrfEndless.clear();

	vector<uint> PosVec;
	vector<uint> IsMetVec;
	uint StopCount = 0;
	uint MetCount = 0;
	for (uint AaPos = 0; AaPos < m_CodonCount; ++AaPos)
		{
		if (IsMet(AaPos))
			{
			++MetCount;
			PosVec.push_back(AaPos);
			IsMetVec.push_back(true);
			}
		else if (IsStop(AaPos))
			{
			++StopCount;
			PosVec.push_back(AaPos);
			IsMetVec.push_back(false);
			}
		}

	if (StopCount == 0)
		{
	// Add ORF terminated by end of sequence
		uint StartPos = UINT_MAX;
		if (MetCount > 0)
			{
			m_OrfHasMets.push_back(true);
			StartPos = PosVec[0];
			m_OrfAaStarts.push_back(StartPos);
			}
		else
			{
			StartPos = 0;
			m_OrfHasMets.push_back(false);
			m_OrfAaStarts.push_back(0);
			}
		m_OrfAaLengths.push_back(m_CodonCount - StartPos);
		m_OrfHasStops.push_back(false);
		if (m_Circular)
			m_OrfEndless.push_back(true);
		else
			m_OrfEndless.push_back(false);
		return;
		}

	uint MetPos = UINT_MAX;
	uint LastMetPos = UINT_MAX;
	uint FirstStopPos = UINT_MAX;
	const uint N = SIZE(PosVec);
	for (uint i = 0; i < N; ++i)
		{
		uint Pos = PosVec[i];
		bool IsMet = IsMetVec[i];
		if (IsMet)
			{
			if (MetPos == UINT_MAX)
				MetPos = Pos;
			LastMetPos = Pos;
			}
		else
			{
			if (FirstStopPos == UINT_MAX)
				FirstStopPos = Pos;
			if (MetPos != UINT_MAX)
				{
				uint AaLength = Pos - MetPos;
				if (AaLength >= MinCodonCount)
					{
					m_OrfAaStarts.push_back(MetPos);
					m_OrfAaLengths.push_back(AaLength);
					m_OrfHasStops.push_back(true);
					m_OrfHasMets.push_back(true);
					m_OrfEndless.push_back(false);
					}
				}
			MetPos = UINT_MAX;
			LastMetPos = UINT_MAX;
			}
		}

	if (m_Circular)
		{
	// Circular, special case at start/end of sequence
		if (FirstStopPos != UINT_MAX && LastMetPos != UINT_MAX)
			{
			uint AaLength_Start = FirstStopPos;
			uint AaLength_End = m_CodonCount - LastMetPos;
			uint AaLength = AaLength_Start + AaLength_End;
			if (AaLength >= MinCodonCount)
				{
				m_OrfAaStarts.push_back(LastMetPos);
				m_OrfAaLengths.push_back(AaLength);
				m_OrfHasStops.push_back(true);
				m_OrfHasMets.push_back(true);
				m_OrfEndless.push_back(false);
				}
			}
		}

	ValidateOrfs();
	}

void Translator::ValidateOrf(uint OrfIndex) const
	{
	uint AaStart = m_OrfAaStarts[OrfIndex];
	uint AaLength = m_OrfAaLengths[OrfIndex];
	bool HasStop = m_OrfHasStops[OrfIndex];
	bool HasMet = m_OrfHasMets[OrfIndex];
	bool Endless = m_OrfEndless[OrfIndex];

	asserta(AaLength >= m_OrfMinCodonCount);

	if (HasMet)
		asserta(IsMet(AaStart));

	if (HasStop)
		asserta(IsStop(AaStart+AaLength));
	else
		{
		asserta(m_Circular);
		asserta(Endless);
		}

	uint NtStart = GetNtPos(AaStart, 0);
	uint NtEnd = GetNtPos(AaStart + AaLength - 1, 2);
	if (!m_Circular)
		{
		if (m_Frame > 0)
			{
			asserta(NtEnd == NtStart + 3*AaLength - 1);
			asserta(NtStart < NtEnd);
			}
		else
			{
			asserta(NtStart == NtEnd + 3*AaLength - 1);
			assert(NtStart > NtEnd);
			}
		}
	}

void Translator::ValidateOrfs() const
	{
	const uint OrfCount = SIZE(m_OrfAaStarts);
	asserta(SIZE(m_OrfAaLengths) == OrfCount);
	asserta(SIZE(m_OrfHasStops) == OrfCount);
	asserta(SIZE(m_OrfHasMets) == OrfCount);
	asserta(SIZE(m_OrfEndless) == OrfCount);

	for (uint i = 0; i < OrfCount; ++i)
		ValidateOrf(i);
	}

void Translator::GetLongestOrf(int &Frame, uint &OrfIndex)
	{
	Frame = 0;
	OrfIndex = UINT_MAX;

	vector<int> Frames;
	GetFrames(Frames);
	const uint N = SIZE(Frames);
	uint MaxLen = 0;
	for (uint i = 0; i < N; ++i)
		{
		SetFrame(Frames[i], true);
		const uint N = SIZE(m_OrfAaLengths);
		for (uint j = 0; j < N; ++j)
			{
			uint Len = m_OrfAaLengths[j];
			if (Len >  MaxLen)
				{
				MaxLen = Len;
				OrfIndex = j;
				Frame = Frames[i];
				}
			}
		}
	}

uint Translator::GetMaxOrfLength() const
	{
	uint MaxLen = 0;
	const uint N = SIZE(m_OrfAaLengths);
	for (uint i = 0; i < N; ++i)
		MaxLen = max(m_OrfAaLengths[i], MaxLen);
	return MaxLen;
	}

bool Translator::CurrentFrameHasMet() const
	{
	bool HasMet = false;
	for (uint AaPos = 0; AaPos < m_CodonCount; ++AaPos)
		{
		char aa = GetAa(AaPos);
		if (aa == 'M')
			return true;
		}
	return false;
	}

bool Translator::CurrentFrameIsEndless() const
	{
	bool StopFound = false;
	bool HasMet = false;
	for (uint AaPos = 0; AaPos < m_CodonCount; ++AaPos)
		{
		char aa = GetAa(AaPos);
		if (aa == 'M')
			HasMet = true;
		if (IsStop(AaPos))
			return false;
		}
	return true;
	}

void Translator::GetFrames(vector<int> &Frames) const
	{
	Frames.clear();
	if (ValidFrame(1)) Frames.push_back(1);
	if (ValidFrame(2)) Frames.push_back(2);
	if (ValidFrame(3)) Frames.push_back(3);
	if (ValidFrame(-1)) Frames.push_back(-1);
	if (ValidFrame(-2)) Frames.push_back(-2);
	if (ValidFrame(-3)) Frames.push_back(-3);
	}

void Translator::GetEndlessFrames2(int &PlusFrame, int &MinusFrame)
	{
	PlusFrame = 0;
	MinusFrame = 0;
	vector<int> Frames;
	GetFrames(Frames);
	for (uint k = 0; k < SIZE(Frames); ++k)
		{
		int Frame = Frames[k];
		SetFrame(Frame);
		if (CurrentFrameIsEndless() && CurrentFrameHasMet())
			{
			if (Frame > 0)
				PlusFrame = Frame;
			else
				MinusFrame = Frame;
			}
		}
	}

void Translator::GetEndlessFrames(vector<int> &EndlessFrames,
  vector<bool> &HasMets)
	{
	EndlessFrames.clear();
	vector<int> Frames;
	GetFrames(Frames);
	for (uint k = 0; k < SIZE(Frames); ++k)
		{
		int Frame = Frames[k];
		SetFrame(Frame);
		if (CurrentFrameIsEndless())
			{
			bool HasMet = CurrentFrameHasMet();
			EndlessFrames.push_back(Frame);
			HasMets.push_back(HasMet);
			}
		}
	}

void Translator::HasEndlessFrame(bool &Plus, bool &Minus)
	{
	Plus = false;
	Minus = false;

	vector<int> Frames;
	GetFrames(Frames);
	for (uint k = 0; k < SIZE(Frames); ++k)
		{
		int Frame = Frames[k];
		SetFrame(Frame);
		if (CurrentFrameIsEndless())
			{
			if (Frame > 0)
				Plus = true;
			else
				Minus = true;
			}
		}
	}

void Translator::PrettyPrint(FILE *f) const
	{
	if (f == 0)
		return;

	Pf(f, "\n");
	Pf(f, "Frame %+d\n", m_Frame);

	string NtStr, AaStr;
	GetAaStr(AaStr);

	const uint CodonCount = SIZE(AaStr);

	string AaRow;
	string NtRow;
	string CoordRow;

	for (uint AaPos = 0; AaPos < CodonCount; ++AaPos)
		{
		uint NtPos = GetNtPos(AaPos, 0);
		string s;
		Ps(s, "%4u", NtPos);
		CoordRow += s;

		AaRow += "   ";
		AaRow += AaStr[AaPos];

		NtRow += ' ';

		string CodonStr;
		GetCodonStr(AaPos, CodonStr);

		NtRow += CodonStr;
		}
	Pf(f, "%s\n", CoordRow.c_str());
	Pf(f, "%s\n", AaRow.c_str());
	Pf(f, "%s\n", NtRow.c_str());
	}
