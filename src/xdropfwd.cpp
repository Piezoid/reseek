#include "myutils.h"
#include "mx.h"
#include "xdpmem.h"
#include "tracebit.h"
#include "pathinfo.h"
#include "swtrace.h"

void GetPathCounts(const string &Path, uint &M, uint &D, uint &I);

static void TraceBack(XDPMem &Mem, uint Besti, uint Bestj, char StartState, string &Path)
	{
	Path.clear();
	byte **TB_M = Mem.GetTBM();
	byte **TB_D = Mem.GetTBD();
	byte **TB_I = Mem.GetTBI();
	uint i = Besti;
	uint j = Bestj;
	char State = StartState;
	uint step = 0;
	
	Log("TRACEBACK_STEP: Starting at (%u,%u) state=%c\n", i, j, State);
	
	for (;;)
		{
		Path += State;
		byte t; // The tracebit for the current state and cell
		
		char NextState = '?';
		switch (State)
			{
		case 'M':
			t = TB_M[i][j];
			Log("TRACEBACK_STEP: step=%u at (%u,%u) state=%c tb=0x%02x path=\"%s\"\n", 
				step, i, j, State, t, Path.c_str());
			// Check for the start-of-alignment marker first.
			if (t & TRACEBITS_SM)
				{
				Log("TRACEBACK_STEP: Found start-of-alignment marker, ending traceback\n");
				goto end_traceback;
				}
			NextState = GetTBBitM(t);
			Log("TRACEBACK_STEP: M->%c transition from (%u,%u)\n", NextState, i, j);
			asserta(i > 0 && j > 0);
			--i;
			--j;
			break;

		case 'D':
			t = TB_D[i][j];
			Log("TRACEBACK_STEP: step=%u at (%u,%u) state=%c tb=0x%02x path=\"%s\"\n", 
				step, i, j, State, t, Path.c_str());
			NextState = GetTBBitD(t);
			Log("TRACEBACK_STEP: D->%c transition from (%u,%u)\n", NextState, i, j);
			asserta(i > 0);
			--i;
			break;

		case 'I':
			t = TB_I[i][j];
			Log("TRACEBACK_STEP: step=%u at (%u,%u) state=%c tb=0x%02x path=\"%s\"\n", 
				step, i, j, State, t, Path.c_str());
			NextState = GetTBBitI(t);
			Log("TRACEBACK_STEP: I->%c transition from (%u,%u)\n", NextState, i, j);
			asserta(j > 0);
			--j;
			break;

		default:
			Die("TraceBackBit, invalid state %c", State);
			}
		State = NextState;
		++step;
		
		// Safety check to prevent infinite loops
		if (step > 10000)
			{
			Log("ERROR: Traceback exceeded 10000 steps, possible infinite loop\n");
			Die("Traceback infinite loop detected");
			}
		}
end_traceback:; // Label for goto.
	std::reverse(Path.begin(), Path.end());
	Log("TRACEBACK_COMPLETE: Final path=\"%s\" (length=%u)\n", Path.c_str(), (uint)Path.length());
	}

// After traceback start of path is *ptrSegLoA,*ptrSegLoB
// may be different from LoA,LoB
float XDropFwd(XDPMem &Mem,
  float X, float Open, float Ext, 
  fn_SubstScore SubFn, void *UserData,
  uint LoA, uint aLA, uint LoB, uint aLB,
  uint *ptrSegLoA, uint *ptrSegLoB,
  string &Path)
	{
	*ptrSegLoA = UINT_MAX;
	*ptrSegLoB = UINT_MAX;
	asserta(LoA < aLA);
	asserta(LoB < aLB);
	uint LA = aLA - LoA;
	uint LB = aLB - LoB;

	Path.clear();

	if (LA == 1 || LB == 1)
		{
		float Score = SubFn(UserData, LoA, LoB);
		if (Score > 0)
			Path.push_back('M');
		return Score;
		}

	if (Open > 0.0f || Ext > 0.0f)
		Warning("XDropFwdFast(): non-negative open %.1f, ext %.1f", Open, Ext);
	const float AbsOpen = -Open;
	const float AbsExt = -Ext;

	Mem.Alloc(LA+1, LB+1);

	byte **TB_M = Mem.GetTBM();
	byte **TB_D = Mem.GetTBD();
	byte **TB_I = Mem.GetTBI();
	INIT_TRACE(LA, LB, TB_M);
	INIT_TRACE(LA, LB, TB_D);
	INIT_TRACE(LA, LB, TB_I);

	float *Mrow = Mem.GetDPRow1();
	float *Drow = Mem.GetDPRow2();
	float *Irow = Mem.GetIrow();

	Mrow[-1] = MINUS_INFINITY;
	TRACE_M(0, -1, MINUS_INFINITY);

	// Explicitly initialize all of row 0 (i=0) for local alignment
	// This ensures proper boundary conditions regardless of band expansion
	for (uint j = 0; j <= LB; ++j)
		{
		Mrow[j] = 0;           // M(0,j) = 0 for all j (local alignment boundary)
		Drow[j] = MINUS_INFINITY; // D(0,j) = -∞ for all j (no deletions from boundary)
		Irow[j] = MINUS_INFINITY; // I(0,j) = -∞ for all j (no insertions from boundary)
		TRACE_M(0, j, Mrow[j]);
		TRACE_D(0, j, Drow[j]);
		}

// Main loop
	float BestScore = 0;
	TRACE_M(1, 1, 0);
	uint Besti = 0;
	uint Bestj = 0;
	char BestState = 'M';

	uint prev_jlo = 0;
	uint prev_jhi = 0;
	uint jlo = 1;
	uint jhi = 1;

// Inner loop does this:
//	Mrow[j] = DPM[i][j+1] -> DPM[i+1][j+1]
//	Drow[j] = DPD[i][j]   -> DPD[i+1][j]

	float M0 = BestScore;
	for (uint i = 1; i <= LA; ++i)
		{
#if TRACE && !DOTONLY
		Log("XDrop i=%u j=%u .. %u\n", i, jlo, jhi);
#endif
		if (jlo == prev_jlo)
			{
			asserta(jlo>0);
			Mrow[jlo-1] = MINUS_INFINITY;
			Drow[jlo] = MINUS_INFINITY;
			TRACE_M(i, jlo-1, MINUS_INFINITY);
			TRACE_D(i, jlo, MINUS_INFINITY);
			}

		uint endj = min(prev_jhi+1,LB);
		for (uint j = endj+1; j <= min(jhi+1, LB); ++j)
			{
			Mrow[j-1] = MINUS_INFINITY;
			Drow[j] = MINUS_INFINITY;
			Irow[j] = MINUS_INFINITY;
			TRACE_M(i, j-1, MINUS_INFINITY);
			TRACE_D(i, j, MINUS_INFINITY);
			}

		uint next_jlo = UINT_MAX;
		uint next_jhi = UINT_MAX;

		float I0 = MINUS_INFINITY;

		byte *TB_M_row = TB_M[i];
		byte *TB_D_row = TB_D[i];
		byte *TB_I_row = TB_I[i];
		asserta(jlo>0);
		asserta(jlo<=jhi);
		float SavedM0 = UNINIT;
		
		// Canonical state variables for proper Gotoh implementation
		// Diagonal predecessors from (i-1,j-1)
		float m_diag = MINUS_INFINITY; // M(i-1,j-1)
		float d_diag = MINUS_INFINITY; // D(i-1,j-1)  
		float i_diag = MINUS_INFINITY; // I(i-1,j-1)
		
		// Current row I-scores (i,j-1) -> (i,j)
		float i_curr = MINUS_INFINITY; // I(i,j-1) for current j
		
		// Pipeline M(i,j-1) scores within the inner loop
		float m_left = MINUS_INFINITY; // Represents M(i, j-1)
		
		// Initialize diagonal predecessors for the first j in the band
		if (jlo > 0)
			{
			m_diag = Mrow[jlo-1]; // M(i-1,j-1) for j=jlo
			d_diag = (jlo > 1) ? Drow[jlo-1] : MINUS_INFINITY; // D(i-1,j-1) for j=jlo
			i_diag = (jlo > 1) ? Irow[jlo-1] : MINUS_INFINITY; // I(i-1,j-1) for j=jlo
			}

		for (uint j = jlo; j <= jhi; ++j)
			{
			// Separate trace bits for each state
			byte tb_m = 0; // Trace bits for M(i,j)
			byte tb_d = 0; // Trace bits for D(i,j)  
			byte tb_i = 0; // Trace bits for I(i,j)

			// Load pristine scores from previous row (i-1,j) before they get overwritten
			float m_up = Mrow[j]; // M(i-1,j)
			float d_up_pristine = Drow[j]; // D(i-1,j) - preserve before modification
			float i_up_pristine = Irow[j]; // I(i-1,j) - preserve before modification
			float m_new = MINUS_INFINITY; // M(i,j) - will be computed
			
#if DEBUG
			// Only log cells with meaningful scores or interesting traceback patterns
			bool has_meaningful_input = (m_diag > MINUS_INFINITY) || (d_diag > MINUS_INFINITY) || 
										(i_diag > MINUS_INFINITY) || (m_up > MINUS_INFINITY) || 
										(d_up_pristine > MINUS_INFINITY) || (i_up_pristine > MINUS_INFINITY) || 
										(i_curr > MINUS_INFINITY);
#endif

		// MATCH: M(i,j) = S(A_i, B_j) + max(M(i-1,j-1), D(i-1,j-1), I(i-1,j-1), 0)
			{
			float xM = m_diag;
			tb_m = 0; // Default to M->M transition
			if (d_diag > xM)
				{
				xM = d_diag;
				tb_m = TRACEBITS_DM;
				}
			if (i_diag > xM)
				{
				xM = i_diag;
				tb_m = TRACEBITS_IM;
				}

			// Apply local alignment choice: max(0, predecessor_score)
			// If predecessor score is negative, start a new alignment
			if (xM < 0)
				{
				xM = 0;
				tb_m = TRACEBITS_SM; // Mark as start of new alignment
				}

			float s = SubFn(UserData, LoA + i-1, LoB + j-1);
			TRACE_Sub(i-1, j-1, s);
			s += xM;
			
			m_new = s; // M(i,j) - assign to outer scope variable
			TRACE_M(i, j, s);

			float h = s - BestScore + X;
		// Match-Match
			if (h > 0)
				{
				next_jlo = min(next_jlo, j+1);
				next_jhi = j+1;
				}

		// Match-Delete
			if (h > AbsOpen)
				next_jlo = min(next_jlo, j);

		// Match-Insert potentially extends current row
			if (h > AbsExt && j == jhi && jhi + 1 < LB)
				{
				++jhi;
				uint new_endj = min(jhi+1, LB);
				new_endj = max(new_endj, endj);
				for (uint j2 = endj+1; j2 <= new_endj; ++j2)
					{
				// Nasty special case for j=j2-1, Mrow[j] has already
				// been updated for current i.
					if (j2-1 > j)
						{
						Mrow[j2-1] = MINUS_INFINITY;
						TRACE_M(i, j2-1, MINUS_INFINITY);
						}

					Drow[j2] = MINUS_INFINITY;
					TRACE_M(i, j2, MINUS_INFINITY);
					}
				endj = new_endj;
				}
			}

		// DELETE: D(i,j) = max(M(i-1,j) + Open, D(i-1,j) + Ext)
			{
			float d_new = d_up_pristine + Ext;
			tb_d = TRACEBITS_DD; // Explicitly set D->D (extend)

			// A gap can only be opened from a valid M score.
			if (m_up > MINUS_INFINITY)
				{
				float md = m_up + Open;
				if (md >= d_new)
					{
					d_new = md;
					tb_d = TRACEBITS_MD; // M->D (open)
					}
				}
			Drow[j] = d_new;
			TRACE_D(i, j, d_new);

			float h = d_new - BestScore + X;
		// Delete-Match
			if (h > 0)
				{
				next_jlo = min(next_jlo, j);
				next_jhi = max(next_jhi, j+1);
				}
			}
			
		// INSERT: I(i,j) = max(M(i,j-1) + Open, I(i,j-1) + Ext)
			{
			float i_new = i_curr + Ext;
			tb_i = TRACEBITS_II; // Explicitly set I->I (extend)

			// A gap can only be opened from a valid M score, M(i, j-1), held in m_left.
			if (m_left > MINUS_INFINITY)
				{
				float mi = m_left + Open;
				if (mi >= i_new)
					{
					i_new = mi;
					tb_i = TRACEBITS_MI; // M->I (open)
					}
				}
			Irow[j] = i_new;
			i_curr = i_new;

			float h = i_curr - BestScore + X;
		// Insert-Match
			if (h > 0)
				{
				next_jlo = min(next_jlo, j);
				next_jhi = max(next_jhi, j+1);
				}

		// Insert-Insert potentially extends current row
			if (h > AbsExt && j == jhi && jhi + 1 < LB)
				{
				++jhi;
				uint new_endj = min(jhi+1, LB);
				new_endj = max(new_endj, endj);
				for (uint j2 = endj+1; j2 <= new_endj; ++j2)
					{
					Mrow[j2-1] = MINUS_INFINITY;
					Drow[j2] = MINUS_INFINITY;
					TRACE_M(i, j2-1, MINUS_INFINITY);
					TRACE_D(i, j, MINUS_INFINITY);
					}
				endj = new_endj;
				}
			}
		
			// Update BestScore by checking all three states for cell (i,j)
			if (m_new >= BestScore)
				{
				float old_best = BestScore;
				BestScore = m_new;
				Besti = i;
				Bestj = j;
				BestState = 'M';
				Log("BEST_UPDATE: M(%u,%u)=%.3f (prev=%.3f, state=M) [tb=0x%02x]\n", 
					i, j, m_new, old_best, tb_m);
				}
			if (Drow[j] >= BestScore)
				{
				float old_best = BestScore;
				BestScore = Drow[j];
				Besti = i;
				Bestj = j;
				BestState = 'D';
				Log("BEST_UPDATE: D(%u,%u)=%.3f (prev=%.3f, state=D) [tb=0x%02x]\n", 
					i, j, Drow[j], old_best, tb_d);
				}
			if (i_curr >= BestScore)
				{
				float old_best = BestScore;
				BestScore = i_curr;
				Besti = i;
				Bestj = j;
				BestState = 'I';
				Log("BEST_UPDATE: I(%u,%u)=%.3f (prev=%.3f, state=I) [tb=0x%02x]\n", 
					i, j, i_curr, old_best, tb_i);
				}
		
			// Store M(i,j) in the matrix
			Mrow[j] = m_new;
			
			// Pass M(i,j) to the next iteration, where it will be M(i,j-1)
			m_left = m_new;
			
			// Store trace bits separately for each state
			TB_M_row[j] = tb_m;
			TB_D_row[j] = tb_d;
			TB_I_row[j] = tb_i;
			
#if DEBUG
			// Only log cells with meaningful scores, interesting traceback bits, or potential issues
			bool has_meaningful_output = (m_new > MINUS_INFINITY) || (Drow[j] > MINUS_INFINITY) || (i_curr > MINUS_INFINITY);
			bool has_interesting_traceback = (tb_m & TRACEBITS_SM) || (tb_d & TRACEBITS_MD) || (tb_i & TRACEBITS_MI);
			bool is_potential_issue = (m_new > 0 && m_new < 1.0) || (Drow[j] > 0 && Drow[j] < 1.0) || (i_curr > 0 && i_curr < 1.0);
			
			if (has_meaningful_input || has_meaningful_output || has_interesting_traceback || is_potential_issue)
				{
				// Compact format: Cell(i,j) [inputs] -> [outputs] [traceback]
				Log("(%u,%u) [m_d=%.1f d_d=%.1f i_d=%.1f m_u=%.1f d_u=%.1f i_u=%.1f i_c=%.1f] -> [M=%.2f D=%.2f I=%.2f] [tb_m=0x%02x tb_d=0x%02x tb_i=0x%02x]\n",
					i, j, 
					m_diag > MINUS_INFINITY ? m_diag : 0, d_diag > MINUS_INFINITY ? d_diag : 0, i_diag > MINUS_INFINITY ? i_diag : 0,
					m_up > MINUS_INFINITY ? m_up : 0, d_up_pristine > MINUS_INFINITY ? d_up_pristine : 0, 
					i_up_pristine > MINUS_INFINITY ? i_up_pristine : 0, i_curr > MINUS_INFINITY ? i_curr : 0,
					m_new > MINUS_INFINITY ? m_new : 0, Drow[j] > MINUS_INFINITY ? Drow[j] : 0, i_curr > MINUS_INFINITY ? i_curr : 0,
					tb_m, tb_d, tb_i);
				}
#endif
			
			// Update state for next iteration (j+1)
			// The values from (i-1,j) become diagonal predecessors for (i-1,j+1)
			m_diag = m_up; // M(i-1,j) becomes M(i-1,j-1) for next j
			d_diag = d_up_pristine; // D(i-1,j) becomes D(i-1,j-1) for next j
			i_diag = i_up_pristine; // I(i-1,j) becomes I(i-1,j-1) for next j
			// i_curr will be updated to I(i,j+1) in the next iteration
			}

	// Special case for end of Drow[]
		if (jhi < LB)
			{
			const uint jhi1 = jhi+1;
			TB_D[i][jhi1] = 0;
			float md = Mrow[jhi1] + Open; // M(i-1,jhi1) + Open
			Drow[jhi1] += Ext; // D(i-1,jhi1) + Ext
			if (md >= Drow[jhi1])
				{
				Drow[jhi1] = md;
				TB_D[i][jhi1] = TRACEBITS_MD;
				TRACE_D(i, jhi1, md);
				}
			}

		if (next_jlo == UINT_MAX)
			break;

		prev_jlo = jlo;
		prev_jhi = jhi;
		jlo = next_jlo;
		jhi = next_jhi;
		if (jlo > LB)
			jlo = LB;
		if (jhi > LB)
			jhi = LB;
		asserta(jlo <= jhi);
		asserta(jlo >= prev_jlo);

		if (jlo == prev_jlo)
			{
			// No change in jlo, initialize for next row
			// M0 will be set from Mrow[jlo-1] in the loop
			Drow[jlo] = MINUS_INFINITY;
			Irow[jlo] = MINUS_INFINITY;
			TRACE_D(i, jlo, MINUS_INFINITY);
			}
		else
			{
			assert(jlo > prev_jlo);
			// jlo increased, M0 will be set from Mrow[jlo-1] in the loop
			}
		}

	DONE_TRACE(BestScore, Besti, Bestj, TB_M);
	if (BestScore <= 0.0f)
		return 0.0f;

	Log("TRACEBACK_START: BestScore=%.3f, Besti=%u, Bestj=%u, BestState=%c\n", 
		BestScore, Besti, Bestj, BestState);
	Log("TRACEBACK_START: About to call TraceBack(Mem, %u, %u, %c, Path)\n", 
		Besti, Bestj, BestState);

	TraceBack(Mem, Besti, Bestj, BestState, Path);
	uint nM, nD, nI;
	GetPathCounts(Path, nM, nD, nI);
	uint Loi = LoA + Besti - nM - nD;
	uint Loj = LoB + Bestj - nM - nI;
	*ptrSegLoA = Loi;
	*ptrSegLoB = Loj;

	{
	const uint ColCount = SIZE(Path);
	uint PosA = Loi;
	uint PosB = Loj;
	float Score2WithGaps = 0;
	char prevState = '\0'; // Previous state to track gap transitions (null for first character)
	
	Log("SCORE_VERIFICATION: Starting recalculation\n");
	Log("SCORE_VERIFICATION: Path=\"%s\" (length=%u)\n", Path.c_str(), ColCount);
	Log("SCORE_VERIFICATION: StartPos=(%u,%u), BestPos=(%u,%u), BestState=%c\n", 
		Loi, Loj, Besti, Bestj, BestState);
	Log("SCORE_VERIFICATION: Path counts: M=%u D=%u I=%u\n", nM, nD, nI);
	
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		float step_score = 0;
		switch (c)
			{
		case 'M': 
			step_score = SubFn(UserData, PosA, PosB);
			Score2WithGaps += step_score;
			Log("SCORE_VERIFICATION: step=%u M(%u,%u) score=%.3f total=%.3f\n", 
				Col, PosA, PosB, step_score, Score2WithGaps);
			++PosA;
			++PosB;
			break;

		case 'D':
			// Add gap penalty: Open for first gap, Ext for extension
			if (prevState == 'D')
				{
				step_score = Ext; // Gap extension
				Log("SCORE_VERIFICATION: step=%u D(%u,%u) gap_ext=%.3f total=%.3f\n", 
					Col, PosA, PosB, step_score, Score2WithGaps + step_score);
				}
			else
				{
				step_score = Open; // Gap opening (including first gap)
				Log("SCORE_VERIFICATION: step=%u D(%u,%u) gap_open=%.3f total=%.3f\n", 
					Col, PosA, PosB, step_score, Score2WithGaps + step_score);
				}
			Score2WithGaps += step_score;
			++PosA;
			break;

		case 'I':
			// Add gap penalty: Open for first gap, Ext for extension
			if (prevState == 'I')
				{
				step_score = Ext; // Gap extension
				Log("SCORE_VERIFICATION: step=%u I(%u,%u) gap_ext=%.3f total=%.3f\n", 
					Col, PosA, PosB, step_score, Score2WithGaps + step_score);
				}
			else
				{
				step_score = Open; // Gap opening (including first gap)
				Log("SCORE_VERIFICATION: step=%u I(%u,%u) gap_open=%.3f total=%.3f\n", 
					Col, PosA, PosB, step_score, Score2WithGaps + step_score);
				}
			Score2WithGaps += step_score;
			++PosB;
			break;
			}
		prevState = c;
		}
	// Compare full score (including gap penalties) with DP result
	float score_diff = fabs(Score2WithGaps - BestScore);
	Log("SCORE_VERIFICATION: BestScore=%.3f, RecalculatedScore=%.3f, Diff=%.3f\n", 
		BestScore, Score2WithGaps, score_diff);
	if (score_diff >= 0.1)
		{
		Log("ERROR: Score mismatch detected! This will trigger assertion failure.\n");
		Log("ERROR: BestScore=%.3f vs RecalculatedScore=%.3f (diff=%.3f)\n", 
			BestScore, Score2WithGaps, score_diff);
		Log("ERROR: Path=\"%s\"\n", Path.c_str());
		Log("ERROR: BestPos=(%u,%u) state=%c, StartPos=(%u,%u)\n", 
			Besti, Bestj, BestState, Loi, Loj);
		}
	asserta(score_diff < 0.1);
	}
	return BestScore;
	}
