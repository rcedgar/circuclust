#include "myutils.h"
#include "alpha.h"

static const char *g_TableStrings[] =
	{
#include "table.inc"
	};

static const uint LINES_PER_CODE = 6;
static const uint g_TableStringSize = sizeof(g_TableStrings)/sizeof(g_TableStrings[0]);
static const uint g_CodeCount = g_TableStringSize/LINES_PER_CODE;

static char **g_CodeNames;
static byte **g_CodonIntToAAChar;
static uint *g_NCBICodonInts;

static char *g_ActiveCodeName = 0;
static uint g_ActiveCode = UINT_MAX;
static uint g_ActiveCode_NCBI = UINT_MAX;
static byte *g_ActiveCodonIntToAAChar = 0;
static byte *g_ActiveCodonIntToAACharUpper = 0;

static uint GetCodonInt(byte b1, byte b2, byte b3)
	{
	uint i1 = g_CharToLetterNucleo[b1];
	if (i1 == INVALID_LETTER)
		return UINT_MAX;
	uint i2 = g_CharToLetterNucleo[b2];
	if (i2 == INVALID_LETTER)
		return UINT_MAX;
	uint i3 = g_CharToLetterNucleo[b3];
	if (i3 == INVALID_LETTER)
		return UINT_MAX;
	assert(i1 < 4 && i2 < 4 && i3 < 4);
	return i1*16 + i2*4 + i3;
	}

byte TranslateCodon_Ints(uint i1, uint i2, uint i3)
	{
	uint CodonInt = i1*16 + i2*4 + i3;
	byte aa = g_ActiveCodonIntToAACharUpper[CodonInt];
	return aa;
	}

byte TranslateCodon(byte b1, byte b2, byte b3)
	{
	if (b1 == '-' && b2 == '-' && b3 == '-')
		return '-';
	uint CodonInt = GetCodonInt(b1, b2, b3);
	if (CodonInt == UINT_MAX)
		return 'X';
	asserta(CodonInt < 64);
	byte aa = g_ActiveCodonIntToAACharUpper[CodonInt];
	return aa;
	}

byte TranslateCodon3(const char Codon[3])
	{
	byte b1 = byte(Codon[0]);
	byte b2 = byte(Codon[1]);
	byte b3 = byte(Codon[2]);

	return TranslateCodon(b1, b2, b3);
	}

void CodonIntToNtInts(uint CodonInt, uint &i1, uint &i2, uint &i3)
	{
	asserta(CodonInt < 64);

	i1 = (CodonInt/16)%4;
	i2 = (CodonInt/4)%4;
	i3 = CodonInt%4;
	}

void LogCode(uint Code, bool NCBIOrder = false)
	{
	asserta(Code < g_CodeCount);

	Log("\n");
	Log("[%u] %s\n", Code, g_CodeNames[Code]);

	string AAs;
	string Base1s;
	string Base2s;
	string Base3s;
	string Starts;

	for (uint k = 0; k < 64; ++k)
		{
		uint CodonInt = (NCBIOrder ? g_NCBICodonInts[k] : k);
		uint i1, i2, i3;
		CodonIntToNtInts(CodonInt, i1, i2, i3);

		byte b1 = g_LetterToCharNucleo[i1];
		byte b2 = g_LetterToCharNucleo[i2];
		byte b3 = g_LetterToCharNucleo[i3];

		byte aa = g_CodonIntToAAChar[Code][CodonInt];
		if (aa == '*')
			Starts += '*';
		else if (islower(aa))
			Starts += 'M';
		else
			Starts += '.';

		AAs += aa;
		Base1s += b1;
		Base2s += b2;
		Base3s += b3;
		}
	Log("    AA  %s\n", AAs.c_str());
	Log("Starts  %s\n", Starts.c_str());
	Log(" Base1  %s\n", Base1s.c_str());
	Log(" Base2  %s\n", Base2s.c_str());
	Log(" Base3  %s\n", Base3s.c_str());
	}

static void InitCode(uint Code)
	{
	uint Base = Code*LINES_PER_CODE;
	const char *Name = g_TableStrings[Base];
	const byte *AA = (const byte *) g_TableStrings[Base+1];
	const byte *Starts = (const byte *) g_TableStrings[Base+2];
	const byte *Base1 = (const byte *) g_TableStrings[Base+3];
	const byte *Base2 = (const byte *) g_TableStrings[Base+4];
	const byte *Base3 = (const byte *) g_TableStrings[Base+5];

	g_CodeNames[Code] = mystrsave(Name);

	vector<bool> Done(64, false);
	for (uint i = 0; i < 64; ++i)
		{
		byte aa = AA[i];
		byte start = Starts[i];
		byte b1 = Base1[i];
		byte b2 = Base2[i];
		byte b3 = Base3[i];

		byte AminoInt = (aa == '*' ? 20 : g_CharToLetterAmino[aa]);
		asserta(AminoInt <= 20);

		uint CodonInt = GetCodonInt(b1, b2, b3);
		g_NCBICodonInts[i] = CodonInt;
		asserta(CodonInt < 64);
		asserta(!Done[CodonInt]);
		Done[CodonInt] = true;

		if (start == 'M')
			aa = tolower(aa);
		else if (start == '*')
			aa = '*';
		else
			asserta(start == '-');

		g_CodonIntToAAChar[Code][CodonInt] = aa;
		}
	}

uint GetNCBICodeNr()
	{
	return g_ActiveCode_NCBI;
	}

void SetCode(uint NCBICodeNr)
	{
	g_ActiveCode = UINT_MAX;
	g_ActiveCode_NCBI = UINT_MAX;
	for (uint i = 0; i < g_CodeCount; ++i)
		{
		if (atoi(g_CodeNames[i]) == NCBICodeNr)
			{
			g_ActiveCode = i;
			break;
			}
		}
	if (g_ActiveCode == UINT_MAX)
		Die("Code %u not found", NCBICodeNr);

	g_ActiveCode_NCBI = NCBICodeNr;
	g_ActiveCodonIntToAAChar = g_CodonIntToAAChar[g_ActiveCode];
	g_ActiveCodeName = g_CodeNames[g_ActiveCode];

	if (g_ActiveCodonIntToAACharUpper == 0)
		g_ActiveCodonIntToAACharUpper = myalloc(byte, 64);
	for (uint i = 0; i < 64; ++i)
		{
		byte c = g_ActiveCodonIntToAAChar[i];
		g_ActiveCodonIntToAACharUpper[i] = toupper(c);
		}
	}

void InitCodes()
	{
	asserta(g_CodeCount*LINES_PER_CODE == g_TableStringSize);

	g_NCBICodonInts = myalloc(uint, 64);
	g_CodonIntToAAChar = myalloc(byte *, g_CodeCount);
	g_CodeNames = myalloc(char *, g_CodeCount);

	for (uint Code = 0; Code < g_CodeCount; ++Code)
		{
		g_CodonIntToAAChar[Code] = myalloc(byte, 64);
		InitCode(Code);
		}

	if (optset_code)
		SetCode(opt(code));
	else
		SetCode(1);
	}

void cmd_log_tables()
	{
	opt(log_tables);

	InitCodes();
	for (uint Code = 0; Code < g_CodeCount; ++Code)
		{
		LogCode(Code, true);
		LogCode(Code, false);
		}
	}
