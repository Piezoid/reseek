#include "myutils.h"
#include "pdbchain.h"
#include <cmath>
#include <vector>

// Helper struct for 3D coordinates to improve data locality.
struct Vec3d {
  float x, y, z;
};

inline float dist_sq(const Vec3d &a, const Vec3d &b) {
  const float dx = a.x - b.x;
  const float dy = a.y - b.y;
  const float dz = a.z - b.z;
  return dx * dx + dy * dy + dz * dz;
}

struct LocCounter {
  uint nr_preserved;
  uint nr_considered;
  Vec3d query;
  Vec3d target;
};

// Pre-calculated constants to avoid magic numbers and repeated calculations.
static constexpr float g_LDDT_R0 = 15.0;
static constexpr float g_LDDT_R0_sq = g_LDDT_R0 * g_LDDT_R0; // 225.0
static constexpr float g_LDDT_thresholds[4] = {0.5, 1.0, 2.0, 4.0};
static constexpr uint g_nr_thresholds = 4;

double GetLDDT_mu(const PDBChain &Q, const PDBChain &T,
                  const std::vector<uint> &PosQs,
                  const std::vector<uint> &PosTs, bool DaliScorerCompatible) {
  const uint nr_aligned = SIZE(PosQs);
  if (nr_aligned < 2) {
    return 0.0;
  }
  asserta(SIZE(PosTs) == nr_aligned);

  // --- Optimization 1: Improve Data Locality ---
  // Copy the coordinates of only the aligned residues into contiguous arrays.
  // This makes subsequent memory access much faster (cache-friendly) and avoids
  // indirect lookups via PDBChain::GetCoord() inside the hot loop.
  std::vector<LocCounter> locs;
  locs.reserve(nr_aligned);

  for (uint i = 0; i < nr_aligned; ++i) {
    const uint posQ = PosQs[i];
    const uint posT = PosTs[i];
    if (posQ != UINT_MAX && posT != UINT_MAX) {
      float xq, yq, zq, xt, yt, zt;
      Q.GetXYZ(posQ, xq, yq, zq);
      T.GetXYZ(posT, xt, yt, zt);
      locs.push_back({0, 0, {xq, yq, zq}, {xt, yt, zt}});
    }
  }

  const uint nr_cols_considered = SIZE(locs);
  if (nr_cols_considered < 2) {
    return 0.0;
  }

  for (uint i = 0; i < nr_cols_considered; ++i) {
    for (uint j = i + 1; j < nr_cols_considered; ++j) {
      // --- Optimization 3: Avoid sqrt() with Squared Distances ---
      // Calculate squared distances first. This is much cheaper than sqrt().
      const float d1_sq = dist_sq(locs[i].query, locs[j].query);
      const float d2_sq = dist_sq(locs[i].target, locs[j].target);

      // Apply the 15A cutoff using the cheap squared distance check.
      // This filter avoids the expensive sqrt() for most pairs.
      if (DaliScorerCompatible) {
        if (d1_sq > g_LDDT_R0_sq)
          continue;
      } else {
        if (d1_sq > g_LDDT_R0_sq && d2_sq > g_LDDT_R0_sq)
          continue;
      }

      // This pair of residues (i, j) contributes to the local neighborhood
      // score of *both* residue i and residue j.
      locs[i].nr_considered += g_nr_thresholds;
      locs[j].nr_considered += g_nr_thresholds;

      // Only now, for the few pairs that pass the cutoff, do we compute the
      // sqrt().
      const float d1 = std::sqrt(static_cast<float>(d1_sq));
      const float d2 = std::sqrt(static_cast<float>(d2_sq));
      const float diff = std::abs(d1 - d2);

      // FIXME: gcc currently use inefficient comiss + setbe chain for this.
      // Although this path is not that hot
      const uint nr_preserved =
          (diff <= g_LDDT_thresholds[0]) + (diff <= g_LDDT_thresholds[1]) +
          (diff <= g_LDDT_thresholds[2]) + (diff <= g_LDDT_thresholds[3]);

      locs[i].nr_preserved += nr_preserved;
      locs[j].nr_preserved += nr_preserved;
    }
  }

  // --- Final Calculation ---
  // Calculate the final score by averaging the individual residue scores.
  double total_score = 0.0;
  for (uint i = 0; i < nr_cols_considered; ++i) {
    const LocCounter &loc_count = locs[i];
    if (loc_count.nr_considered > 0) {
      total_score += (double)loc_count.nr_preserved / loc_count.nr_considered;
    }
  }

  return total_score / nr_cols_considered;
}

#if 0
#include "daliscorer.h"
#include "seqdb.h"

void cmd_test()
	{
	asserta(optset_input);

	string Name;
	GetStemName(g_Arg1, Name);

	SeqDB MSA;
	MSA.FromFasta(g_Arg1, true);

	FILE* fOut = CreateStdioFile(opt(output));
	const bool MissingSeqOk = opt(missingtestseqok);

	DALIScorer DS;
	DS.LoadChains(opt(input));
	bool Ok = DS.SetMSA(Name, MSA, false, MissingSeqOk);
	if (!Ok)
		Die("SetMSA failed");

	const uint SeqCount = MSA.GetSeqCount();
	double Sum_Z = 0;
	double Sum_Z15 = 0;
	double Sum_LDDT_mu = 0;
	double Sum_LDDT_fm = 0;

	uint PairCount = 0;
	for (uint SeqIdx1 = 0; SeqIdx1 < SeqCount; ++SeqIdx1)
		{
		const uint ChainIdx1 = DS.m_SeqIdxToChainIdx[SeqIdx1];
		const char *Label1 = MSA.GetLabel(SeqIdx1).c_str();
		const vector<uint> &ColToPos1 = DS.m_ColToPosVec[SeqIdx1];
		for (uint SeqIdx2 = SeqIdx1 + 1; SeqIdx2 < SeqCount; ++SeqIdx2)
			{
			const uint ChainIdx2 = DS.m_SeqIdxToChainIdx[SeqIdx2];
			const vector<uint> &ColToPos2 = DS.m_ColToPosVec[SeqIdx2];
			const char *Label2 = MSA.GetLabel(SeqIdx2).c_str();
			if (ChainIdx1 == UINT_MAX || ChainIdx2 == UINT_MAX)
				{
				Log("%s\t%s\tERROR_structure_not_found\n", Label1, Label2);
				continue;
				}
			++PairCount;

			double LDDT_mu = DS.GetLDDTChainPair_muscle(
			  ChainIdx1, ChainIdx2, ColToPos1, ColToPos2);

			const PDBChain &Chain1 = *DS.m_Chains[ChainIdx1];
			const PDBChain &Chain2 = *DS.m_Chains[ChainIdx2];
			double LDDT_mu2 = GetLDDT_mu(
			  Chain1, Chain2, ColToPos1, ColToPos2, true);
			double LDDT_mu3 = GetLDDT_mu(
			  Chain1, Chain2, ColToPos1, ColToPos2, false);

			Log("label1=%s\tlabel2=%s\tLDDT_mu=%.4f\tLDDT_mu2=%.4f\tLDDT_mu3=%.4f\n",
			  Label1, Label2, LDDT_mu, LDDT_mu2, LDDT_mu3);
			
			asserta(feq(LDDT_mu, LDDT_mu2));
			}
		}
	}
#endif // 0
