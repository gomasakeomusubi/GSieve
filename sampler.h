#ifndef __SAMPLE__
#define __SAMPLE__

#include "tool.h"
#include <NTL/ZZ.h>
#include <NTL/vec_RR.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>

NTL_CLIENT

class KleinSampler {
    public:
        void Init(const mat_ZZ &B);
        LatticeVector* Sample();
    private:
        long n_;
        long m_;
        double t_;
        mat_ZZ B_;
        mat_RR mu_;
        vec_RR coef_;
        vec_RR s_prime_square_;
        long SampleZ(RR c, RR s_square);
};

#endif