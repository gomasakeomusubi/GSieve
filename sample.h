#ifndef INCLUDE_SAMPLE
#define INCLUDE_SAMPLE
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/matrix.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_RR.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <random>
#include <list>

#include <chrono>
#include <ctime>

#include "tool.h"

namespace sample{
  using NTL::ZZ;
  using NTL::vec_ZZ;
  using NTL::mat_ZZ;
  using std::cout;
  using std::endl;
  using std::default_random_engine;
  using std::uniform_real_distribution;
  using NTL::power2_ZZ;
  using NTL::to_ZZ;
  using NTL::SetSeed;
  using NTL::RandomBnd;
  using std::max;


  void gen_rref_mat(mat_ZZ& A, int seed);

  ZZ max(mat_ZZ& A);

  ZZ max(vec_ZZ& v);

  void gen_random_unimodular2(mat_ZZ &L,int dim,int seed,int bits,int vl=0);

  void sample(mat_ZZ B, LatticeVector *ret);
  void sample_set(mat_ZZ B, vector<LatticeVector*>& V, long num);
}
#endif
