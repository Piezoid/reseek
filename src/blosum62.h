#pragma once

#include "arrays.h"

// Type for dense storage of BLOSUM62 scores
using blosum62_t = int8_t;
// Type for summing BLOSUM62 scores
using blosum_sum_t = int;

extern SquareMatrix<blosum62_t> g_SubstMx;


void SetBLOSUM62();
blosum62_t GetBlosum62Score(char a, char b);


blosum_sum_t GetBlosum62PathScore(
  const string &A, uint LoA,
  const string &B, uint LoB,
  const string &Path);