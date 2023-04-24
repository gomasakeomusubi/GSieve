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
        long n_;
        long m_;
        vector<LatticeVector*> L, V;
        queue<LatticeVector*> S;
        LatticeVector *min_vector_;
        ZZ goal_norm2_;
        int simu_samp_;      // 一度のループで初めにサンプルする数。固定。
        int concurrency_;
        void SampleReduce(LatticeVector *p);
        void SampleReduce_Parallel();
        void ListReduce(LatticeVector *p);
        void ListReduce_Parallel();
        void VectorReduce_Parallel();
        void GaussReduce(LatticeVector *p);
        void GaussReduce_Parallel();
        long max_list_size_;
        long collisions_;
        long iterations_;
        long sample_vectors_;
    public:
        GSieve(mat_ZZ B) :B(B){}
        ~GSieve(){
            CleanUp();
        }
        void CleanUp(){
            deleteList(L);
            deleteList(V);
            deleteQueue(S);
        }
        void Setup();
        void SetConcurrency(long num){ concurrency_ = num; }
        void SetSimultaneousSamples(long num){ simu_samp_ = num; }

        LatticeVector *GaussSieve(vector<double> &chk_time);
        LatticeVector *GaussSieve_Parallel(vector<double> &chk_time);

        long getIterations(){ return iterations_; }
        long getSampleVectors(){ return sample_vectors_; }
        ZZ getMinNorm2(){ return min_vector_->norm2; }
};

#endif