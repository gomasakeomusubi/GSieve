#ifndef INCLUDE_GSIEVE
#define INCLUDE_GSIEVE

#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/vec_ZZ.h>
#include <NTL/vec_RR.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_RR.h>
#include <vector>
#include <queue>
#include <chrono>
#include <thread>
#include "tool.h"

NTL_CLIENT

class GSieve{
    private:
        mat_ZZ B;
        vector<LatticeVector*> L, V;
        queue<LatticeVector*> S;
        LatticeVector *min_vector;
        ZZ thresh;
        int simu_samp;      // 一度のループで初めにサンプルする数。固定。
        int concurrency;
    public:
        GSieve(mat_ZZ B, ZZ thresh, int simu_samp, int concurrency) :B(B), thresh(sqr(thresh)), simu_samp(simu_samp), concurrency(concurrency){}
        ~GSieve(){
            deleteList(L);
            deleteList(V);
            deleteQueue(S);
        }
        void Setup();

        void SampleReduce(LatticeVector *p);
        void SampleReduce_Parallel();

        void ListReduce(LatticeVector *p);
        void ListReduce_Parallel();

        void VectorReduce_Parallel();

        void GaussReduce(LatticeVector *p);
        void GaussReduce_Parallel();

        LatticeVector *GaussSieve(vector<double> &chk_time, long &num_sample, long &cnt);
        LatticeVector *GaussSieve_Parallel(vector<double> &chk_time, long &num_sample, long &cnt);
        
        mat_ZZ getBasis(){ return B; }
        long getLsize(){ return L.size(); }
        long getSsize(){ return S.size(); }
};

#endif