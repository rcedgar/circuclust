#pragma once

struct OrfSpec
	{
	uint Frame = 0;
	uint NtStart = UINT_MAX;
	uint NtLength = UINT_MAX;
	bool HasStop = false;
	};

class Translator
	{
public:
	const byte *m_NtSeq;
	uint m_NtL = 0;
	bool m_Circular = false;

// Per-frame
	int m_Frame = 0;
	uint m_FirstNtPos = UINT_MAX;
	uint m_CodonCount = UINT_MAX;

	uint m_OrfMinCodonCount = UINT_MAX;
	vector<uint> m_OrfAaStarts;
	vector<uint> m_OrfAaLengths;
	vector<bool> m_OrfHasStops;
	vector<bool> m_OrfHasMets;
	vector<bool> m_OrfEndless;

public:
	void SetSeq(const byte *NtSeq, uint NtL, bool Circular);
	void SetStr(const string &NtStr, bool Circular);
	void SetFrame(int Frame, bool DoSetOrfs = false);
	bool ValidFrame(int Frame) const;

public:
	byte GetNt(uint AAPos, uint BaseNr) const;
	byte GetAa(uint AAPos) const;
	void GetCodonStr(uint AAPos, string &Codon) const;
	uint GetNtPos(uint AAPos, uint BaseNr) const;
	bool IsStop(uint AaPos) const;
	bool IsMet(uint AaPos) const;
	void GetFrames(vector<int> &Frames) const;

	uint CalcFirstNtPos(int Frame) const;
	uint CalcLastNtPos(int Frame) const;
	uint CalcCodonCount(int Frame) const;
	void GetAaStr(string &AaStr) const;
	void GetNtStr(string &NtStr) const;
	void GetAaSubStr(uint AaStart, uint AaLength, string &AaStr) const;
	void GetNtSubStr(uint AaStart, uint AaLength, string &NtStr) const;
	void SetOrfs(uint MinCodonCount);
	void ValidateOrfs() const;
	void ValidateOrf(uint i) const;
	bool CurrentFrameIsEndless() const;
	bool CurrentFrameHasMet() const;
	uint GetMaxOrfLength() const;
	void HasEndlessFrame(bool &Plus, bool &Minus);
	void GetEndlessFrames(vector<int> &Frames, vector<bool> &HasMets);
	void GetEndlessFrames2(int &PlusFrame, int &MinusFrame);
	void GetLongestOrf(int &Frame, uint &OrfIndex);
	void PrettyPrint(FILE *f) const;
	};

byte TranslateCodon(byte c1, byte c2, byte c3);
byte TranslateCodon3(const char Codon[3]);
uint GetNCBICodeNr();
void SetCode(uint NCBICodeNr);
uint GetPhase(int Endless_Plus, int Endless_Minus);
void RevCompStr(const string &Seq, string &RCSeq);
