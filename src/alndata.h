#pragma once

struct SeqData;

class AlnData
	{
public:
	string m_LabelQ;
	string m_LabelT;
	string m_SeqQ;
	string m_SeqT;

	string m_QRow;
	string m_TRow;
	uint m_StartPosQ;
	bool m_Plus;
	float m_FractId;

	string m_DoublePath;

public:
	AlnData() {}
	void Clear() {}

	void StripGapsFromQRow(string &Q) const;
	void WriteAln(FILE *f) const;
	void WriteTsv(FILE *f) const;
	void WriteFasta(FILE *f) const;
	};
