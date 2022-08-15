#include "myutils.h"
#include "seqdb.h"
#include "alpha.h"

float **g_SubstMx;
bool g_IsNucleo;
unsigned opt_w;
const byte *g_CharToLetter;
uint g_AlphaSize;

const float MATCH_SCORE = 1.0;
const float MISMATCH_SCORE = -2.0;

void SetSubstMx(const SeqDB &DB)
	{
	g_IsNucleo = DB.GuessIsNucleo();
	if (g_IsNucleo)
		{
		ProgressLog("\nNucleotide sequences\n");
		opt_w = 8;
		g_CharToLetter = g_CharToLetterNucleo;
		g_AlphaSize = 4;
		}
	else
		{
		ProgressLog("\nAmino acid sequences\n");
		opt_w = 3;
		g_CharToLetter = g_CharToLetterAmino;
		g_AlphaSize = 20;
		}
	if (g_SubstMx == 0)
		{
		g_SubstMx = myalloc(float *, 256);
		for (uint i = 0; i < 256; ++i)
			{
			g_SubstMx[i] = myalloc(float, 256);
			for (uint j = 0; j < 256; ++j)
				g_SubstMx[i][j] = 0;
			}
		}

	for (uint i = 0; i < 256; ++i)
		{
		byte Ci = toupper(i);
		for (uint j = 0; j < 256; ++j)
			{
			byte Cj = toupper(j);
			g_SubstMx[Ci][Cj] = (Ci == Cj ? MATCH_SCORE : MISMATCH_SCORE);
			}
		}
	}
