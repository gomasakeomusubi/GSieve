#ifndef __GSIEVE__
#define __GSIEVE__

#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/vec_ZZ.h>
#include <NTL/vec_RR.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_RR.h>
#include <vector>
#include <algorithm>
#include <stack>
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
        vector<ListPoint*> V;
        stack<ListPoint*> S;
        KleinSampler* sampler_;
        int64 goal_norm_;
        int64 min_norm_;
        int simu_samp_;
        int concurrency_;
        // void VectorReduce_Parallel();
        int64 GaussReduce(ListPoint* p);
        int64 TripleReduce_2red(ListPoint* p, int &p_pos);
        bool check_2red2(const ListPoint *p1, const ListPoint *p2);
        int64 TripleReduce(ListPoint* p);
        // int64 Reduce_3red(ListPoint* p);
        int64 GaussReduce_Parallel();
        // statistics
        long max_list_size_;
        long collisions_;
        long iterations_;
        long sample_vectors_;
        vector<double> chk_time_;
        double timeL2V;
        double timeV2V;
        double timeV2L;
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
        long getListSize(){ return max_list_size_; }
        long getSampleVectors(){ return sample_vectors_; }
        ListPoint* getMinVec();
        vector<double> getChkTime(){ return chk_time_; }
        double getTimeL2V(){ return timeL2V; }
        double getTimeV2V(){ return timeV2V; }
        double getTimeV2L(){ return timeV2L; }
};

#endif