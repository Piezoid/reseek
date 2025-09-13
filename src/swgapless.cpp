#include "myutils.h"
#include "arrays.h"
#include "xdpmem.h"

float GetBlosum62Score(char a, char b);

static float SWGaplessNoTB(XDPMem &Mem, const Matrix<float> &SMx, uint LA, uint LB)
	{
#if TRACE && !DOTONLY
	SMx.LogMe();
#endif

	Mem.Alloc(LA+32, LB+32);

	span<float> Mrow = Mem.GetDPRow1();

	*(Mrow.data() - 1) = MINUS_INFINITY;
	for (uint j = 0; j <= LB; ++j)
		Mrow[j] = MINUS_INFINITY;
	
	float BestScore = 0.0f;
	float M0 = float(0);
	for (uint i = 0; i < LA; ++i)
		{
		span<const float> SMxRow = SMx[i];
		float I0 = MINUS_INFINITY;
		for (uint j = 0; j < LB; ++j)
			{
			float xM = M0;
			if (0.0f >= xM)
				xM = 0.0f;

			M0 = Mrow[j];
			xM += SMxRow[j];
			if (xM > BestScore)
				BestScore = xM;

			Mrow[j] = xM;
			}
		M0 = MINUS_INFINITY;
		}
	return BestScore;
	}

float SWFastGapless(XDPMem &Mem, const Matrix<float> &SMx, uint LA, uint LB,
  uint &Besti, uint &Bestj)
	{
	Mem.Alloc(LA+1, LB+1);
	asserta(SMx.Rows() == LA);
	asserta(SMx.Cols() == LB);

	Besti = UINT_MAX;
	Bestj = UINT_MAX;

	span<float> Mrow = Mem.GetDPRow1();

// Use Mrow[-1], so...
	*(Mrow.data() - 1) = MINUS_INFINITY;
	for (uint j = 0; j <= LB; ++j)
		Mrow[j] = MINUS_INFINITY;
	
	float BestScore = 0.0f;

// Main loop
	float M0 = float (0);
	for (uint i = 0; i < LA; ++i)
		{
		span<const float> SMxRow = SMx[i];
		for (uint j = 0; j < LB; ++j)
			{
			float SavedM0 = M0;

		// MATCH
			{
			float xM = M0;
			if (xM < 0.0f)
				xM = 0.0f;

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
		
		M0 = MINUS_INFINITY;
		}
	
	return BestScore;
	}

// Recursion DP[i+1][j+1] = max { 0, DP[i][j] + S[i][j] }
float SWGapless(Matrix<float> &DPMx, const Matrix<float> &SMx, uint LA, uint LB,
  uint &Loi, uint &Loj, uint &ColCount)
	{
	DPMx = Matrix<float>::Allocate(LA+1, LB+1);
#if DEBUG
	for (uint i = 0; i <= LA; ++i) {
		for (uint j = 0; j <= LB; ++j) {
			DPMx[i][j] = FLT_MAX;
		}
	}
#endif
	// Use Matrix directly

	Loi = UINT_MAX;
	Loj = UINT_MAX;
	ColCount = UINT_MAX;
	for (uint i = 0; i <= LA; ++i)
		DPMx[i][0] = 0;
	for (uint j = 0; j <= LB; ++j)
		DPMx[0][j] = 0;

	float BestScore = 0;
	uint Besti = UINT_MAX;
	uint Bestj = UINT_MAX;
	for (uint i = 0; i < LA; ++i)
		{
		for (uint j = 0; j < LB; ++j)
			{
			float sij = SMx[i][j];
			float m = DPMx[i][j] + sij;
			if (m > 0)
				{
				DPMx[i+1][j+1] = m;
				if (m > BestScore)
					{
					BestScore = m;
					Besti = i;
					Bestj = j;
					}
				}
			else
				DPMx[i+1][j+1] = 0;
			}
		}

	if (BestScore == 0)
		return 0;

	assert(Besti < LA && Bestj < LB);
	uint Min = min(Besti, Bestj);
	float Sum = 0;
	float MaxSum = 0;
	Loi = Besti;
	Loj = Bestj;
	uint i = Besti;
	uint j = Bestj;
	//DPMx.LogMe();
	for (;;)
		{
		float s = SMx[i][j];
		Sum += s;
		if (Sum >= BestScore - 0.01 || i == 0 || j == 0)
			{
			Loi = i;
			Loj = j;
			ColCount = Besti - i + 1;
			MaxSum = Sum;
			break;
			}
		--i;
		--j;
		}

	return BestScore;
	}

#if 0
static void MakeBlosumS(const string &A, const string &B,
  Matrix<float> &MxS)
	{
	uint LA = SIZE(A);
	uint LB = SIZE(B);
	MxS = Matrix<float>::Allocate(LA, LB);
#if DEBUG
	for (uint i = 0; i < LA; ++i) {
		for (uint j = 0; j < LB; ++j) {
			MxS[i][j] = FLT_MAX;
		}
	}
#endif
	// Use Matrix directly
	for (uint i = 0; i < LA; ++i)
		{
		char a = A[i];
		for (uint j = 0; j < LB; ++j)
			{
			char b = B[j];
			MxS[i][j] = GetBlosum62Score(a, b);
			}
		}
	}

static void LogAln(const string &A, const string &B,
  uint Loi, uint Loj, uint ColCount)
	{
	const uint LA = SIZE(A);
	const uint LB = SIZE(B);
	asserta(Loi + ColCount <= LA);
	asserta(Loj + ColCount <= LB);
	Log("A ");
	for (uint Col = 0; Col < ColCount; ++Col)
		Log("%c", A[Loi+Col]);
	Log(" %u .. %u\n", Loi + 1, Loi + ColCount);

	Log("B ");
	for (uint Col = 0; Col < ColCount; ++Col)
		Log("%c", B[Loj+Col]);
	Log(" %u .. %u\n", Loj + 1, Loj + ColCount);
	}

static void Test2(const string &A, const string &B)
	{
	uint LA = SIZE(A);
	uint LB = SIZE(B);

	Matrix<float> SMx;
	Matrix<float> DPMx;
	MakeBlosumS(A, B, SMx);

	uint Loi, Loj, ColCount;
	float Score = SWGapless(DPMx, SMx, LA, LB, Loi, Loj, ColCount);
	Log("\n");
	Log("A=%s\n", A.c_str());
	Log("B=%s\n", B.c_str());
	Log("Score %.1f, %u, %u, %u\n", Score, Loi, Loj, ColCount);
	if (Score > 0)
		LogAln(A, B, Loi, Loj, ColCount);

	XDPMem Mem;
	uint Besti, Bestj;
	float Score2 = SWFastGapless(Mem, SMx, LA, LB, Besti, Bestj);
	Log("Score2 %.1f, %u, %u\n", Score2, Besti, Bestj);
	asserta(feq(Score, Score2));
	}

void cmd_testsw()
	{
	Test2("LQNGSEQVENCE", "LQNGSEQVENCE");
	Test2("QNGSEQVENCE", "LQNGSEQVENCE");
	Test2("LQNGSEQVENCE", "QNGSEQVENCE");
	Test2("LQNGSEQVENC", "QNGSEQVENCE");
	Test2("SEQVENCE", "QVE");
	Test2("QVE", "SEQVENCE");
	}
#endif // 0
