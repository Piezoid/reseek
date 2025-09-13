#include "myutils.h"
#include "arrays.h"
#include "xdpmem.h"

float GetBlosum62Score(char a, char b);

int SWFastGapless_Int(XDPMem &Mem, const Matrix<int8_t> &SMx, uint LA, uint LB,
  uint &Besti, uint &Bestj)
	{
	Mem.Alloc(LA+1, LB+1);
	asserta(SMx.Rows() == LA);
	asserta(SMx.Cols() == LB);
	// Use Matrix directly - need to change function signature

	Besti = UINT_MAX;
	Bestj = UINT_MAX;

	int *Mrow = Mem.GetDPRow1Int();

// Use Mrow[-1], so...
	Mrow[-1] = MINUS_INFINITY_INT;
	for (uint j = 0; j <= LB; ++j)
		Mrow[j] = MINUS_INFINITY_INT;
	
	int BestScore = 0;

// Main loop
	int M0 = 0;
	for (uint i = 0; i < LA; ++i)
		{
		span<const int8_t> SMxRow = SMx[i];
		for (uint j = 0; j < LB; ++j)
			{
			int SavedM0 = M0;

		// MATCH
			{
			int xM = M0;
			if (xM < 0)
				xM = 0;

			M0 = Mrow[j];
			xM += SMxRow[j];
			if (xM > BestScore)
				{
				BestScore = xM;
				Besti = i;
				Bestj = j;
				}

			Mrow[j] = xM;
			}
			}
		
		M0 = MINUS_INFINITY_INT;
		}
	
	return BestScore;
	}
