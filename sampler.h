#ifndef __SAMPLE__
#define __SAMPLE__

#include "tool.h"
#include <NTL/ZZ.h>
#include <NTL/vec_RR.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/mat_ZZ.h>

NTL_CLIENT

class KleinSampler { 
    public:
        void Init(const mat_ZZ &B);
        ListPoint* Sample();
    private:
        long n_;
        long m_;
        double t_;
        mat_int64 B_;
        mat_double mu_;
        vec_double coef_;
        vec_double s_prime_square_;
        long SampleZ(double c, double s_square);
};

#endif