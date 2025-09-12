#include "myutils.h"
#include "pdbchain.h"
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
static constexpr float g_LDDT_R0_squared = g_LDDT_R0 * g_LDDT_R0; // 225.0
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
        if (d1_sq > g_LDDT_R0_squared)
          continue;
      } else {
        if (d1_sq > g_LDDT_R0_squared && d2_sq > g_LDDT_R0_squared)
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

double GetLDDT_mu_fast(const PDBChain &Q, const PDBChain &T,
  const vector<uint> &PosQs, const vector<uint> &PosTs)
	{
	const uint nr_cols = SIZE(PosQs);
	if (nr_cols == 0)
		return 0;
	asserta(SIZE(PosTs) == nr_cols);
	uint *nr_considered_vec = myalloc(uint, nr_cols);
	uint *nr_preserved_vec = myalloc(uint, nr_cols);
	zero_array(nr_considered_vec, nr_cols);
	zero_array(nr_preserved_vec, nr_cols);
	for (uint coli = 0; coli < nr_cols; ++coli)
		{
		uint pos1i = PosQs[coli];
		uint pos2i = PosTs[coli];
		assert(pos1i != UINT_MAX);
		assert(pos2i != UINT_MAX);

		for (uint colj = coli + 1; colj < nr_cols; ++colj)
			{
			uint pos1j = PosQs[colj];
			uint pos2j = PosTs[colj];
			assert(pos1j != UINT_MAX);
			assert(pos2j != UINT_MAX);

			float d1_squared = Q.GetDist2(pos1i, pos1j);
			float d2_squared = T.GetDist2(pos2i, pos2j);
			if (d1_squared > g_LDDT_R0_squared && d2_squared > g_LDDT_R0_squared)
				continue;

			float d1 = sqrtf(d1_squared);
			float d2 = sqrtf(d2_squared);
			for (uint k = 0; k < g_nr_thresholds; ++k)
				{
				float t = g_LDDT_thresholds[k];
				float diff = abs(d1 - d2);
				if (diff <= t)
					{
					nr_preserved_vec[coli] += 1;
					nr_preserved_vec[colj] += 1;
					}
				}
			nr_considered_vec[coli] += g_nr_thresholds;
			nr_considered_vec[colj] += g_nr_thresholds;
			}
		}

	float total = 0;
	for (uint col = 0; col < nr_cols; ++col)
		{
		float score = 0;
		uint nr_preserved = nr_preserved_vec[col];
		uint nr_considered = nr_considered_vec[col];
		if (nr_considered > 0)
			score = float(nr_preserved)/nr_considered;
		total += score;
		}
	myfree(nr_considered_vec);
	myfree(nr_preserved_vec);
	float avg = total/nr_cols;
	return avg;
	}
