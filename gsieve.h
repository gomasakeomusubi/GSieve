#ifndef __GSIEVE__
#define __GSIEVE__

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
#include "sampler.h"

NTL_CLIENT

class GSieve{
    private:
        long n_;
        long m_;
        vector<ListPoint*> L;
        vector<ListPoint*> V, V_, V__;
        queue<ListPoint*> S;
        KleinSampler* sampler_;
        int64 goal_norm_;
        int64 min_norm_;
        int simu_samp_;
        int concurrency_;
        // void VectorReduce_Parallel();
        int64 GaussReduce(ListPoint* p);
        int64 GaussReduce_Parallel();
        // statistics
        long max_list_size_;
        long collisions_;
        long iterations_;
        long sample_vectors_;
        // vector<double> chk_time_;
    public:
        ~GSieve(){
            CleanUp();
        }
        void CleanUp();
        void Init(const mat_ZZ &B, KleinSampler* sampler);
        void SetGoalNorm(long norm){ goal_norm_ = norm; }
        void SetConcurrency(long num){ concurrency_ = num; }
        void SetSimultaneousSamples(long num){ simu_samp_ = num; }

        void GaussSieve();
        void GaussSieve_Parallel();

        void printL();
        void printV();
        long getIterations(){ return iterations_; }
        long getCollisions(){ return collisions_; }
        long getListSize(){ return L.size(); }
        long getSampleVectors(){ return sample_vectors_; }
        ListPoint* getMinVec(){ return *L.begin(); }
        // vector<double> getChkTime(){ return chk_time_; }
};

#endif