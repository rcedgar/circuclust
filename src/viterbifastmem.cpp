#include "myutils.h"
#include "xdpmem.h"
#include "tracebit.h"
#include "pathinfo.h"

//const float MATCH_SCORE = 1.0;
//const float MISMATCH_SCORE = -2.0;
const float TERM_GAP_OPEN = -2.0;
const float TERM_GAP_EXT = -0.5;
const float GAP_OPEN = -10.0;
const float GAP_EXT = -1.0;

extern float **g_SubstMx;

float ViterbiFastMem(XDPMem &Mem, const string &A, const string &B, PathInfo &PI)
	{
	const uint LA = SIZE(A);
	const uint LB = SIZE(B);
	if (LA*LB > 100*1000*1000)
		Die("ViterbiFastMem, seqs too long LA=%u, LB=%u", LA, LB);

	Mem.Alloc(LA, LB);
	PI.Alloc2(LA, LB);
	
	StartTimer(ViterbiFastMem);

	float OpenA = TERM_GAP_OPEN;
	float ExtA = TERM_GAP_EXT;

	float *Mrow = Mem.GetDPRow1();
	float *Drow = Mem.GetDPRow2();
	byte **TB = Mem.GetTBBit();

// Use Mrow[-1], so...
	Mrow[-1] = MINUS_INFINITY;
	for (unsigned j = 0; j <= LB; ++j)
		{
		Mrow[j] = MINUS_INFINITY;
		Drow[j] = MINUS_INFINITY;
		}
	
// Main loop
	float M0 = float (0);
	for (unsigned i = 0; i < LA; ++i)
		{
		byte a = A[i];
		const float *MxRow = g_SubstMx[a];
		float OpenB = TERM_GAP_OPEN;
		float ExtB = TERM_GAP_EXT;
		float I0 = MINUS_INFINITY;

		byte *TBrow = TB[i];
		for (unsigned j = 0; j < LB; ++j)
			{
			byte b = B[j];
			byte TraceBits = 0;
			float SavedM0 = M0;

		// MATCH
			{
			float xM = M0;
			if (Drow[j] > xM)
				{
				xM = Drow[j];
				TraceBits = TRACEBITS_DM;
				}
			if (I0 > xM)
				{
				xM = I0;
				TraceBits = TRACEBITS_IM;
				}
			M0 = Mrow[j];

		// Mrow[j] = xM + (toupper(a) == toupper(b) ? MATCH_SCORE : MISMATCH_SCORE);
			Mrow[j] = xM + MxRow[b];
			}
			
		// DELETE
			{
			float md = SavedM0 + OpenB;
			Drow[j] += ExtB;
			if (md >= Drow[j])
				{
				Drow[j] = md;
				TraceBits |= TRACEBITS_MD;
				}
			}
			
		// INSERT
			{
			float mi = SavedM0 + OpenA;
			I0 += ExtA;
			if (mi >= I0)
				{
				I0 = mi;
				TraceBits |= TRACEBITS_MI;
				}
		// I0 = DPI[i][j+1]
			}
			
			OpenB = GAP_OPEN;
			ExtB = GAP_EXT;
			
			TBrow[j] = TraceBits;
			}
		
	// Special case for end of Drow[]
		{
		TBrow[LB] = 0;
		float md = M0 + TERM_GAP_OPEN;
		Drow[LB] += TERM_GAP_EXT;
		if (md >= Drow[LB])
			{
			Drow[LB] = md;
			TBrow[LB] = TRACEBITS_MD;
			}
		}
		
		M0 = MINUS_INFINITY;

		OpenA = GAP_OPEN;
		ExtA = GAP_EXT;
		}
	
// Special case for last row of DPI
	byte *TBrow = TB[LA];
	float I1 = MINUS_INFINITY;
	for (unsigned j = 1; j < LB; ++j)
		{
	// Mrow[j-1] = DPM[LA][j]
	// I1 = DPI[LA][j]
		
		TBrow[j] = 0;
		float mi = Mrow[int(j)-1] + TERM_GAP_OPEN;
		I1 += TERM_GAP_EXT;
		if (mi > I1)
			{
			I1 = mi;
			TBrow[j] = TRACEBITS_MI;
			}
		}
	
	float FinalM = Mrow[LB-1];
	float FinalD = Drow[LB];
	float FinalI = I1;
// FinalM = DPM[LA][LB]
// FinalD = DPD[LA][LB]
// FinalI = DPI[LA][LB]
	
	float Score = FinalM;
	byte State = 'M';
	if (FinalD > Score)
		{
		Score = FinalD;
		State = 'D';
		}
	if (FinalI > Score)
		{
		Score = FinalI;
		State = 'I';
		}

	EndTimer(ViterbiFastMem);

	TraceBackBitMem(Mem, LA, LB, State, PI);

	return Score;
	}
