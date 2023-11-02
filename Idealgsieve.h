#ifndef __IDEALGSIEVE__
#define __IDEALGSIEVE__

#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/vec_ZZ.h>
#include <NTL/vec_RR.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/LLL.h>
#include <vector>
#include <algorithm>
#include <stack>
#include <chrono>
#include <thread>
#include "tool.h"
#include "sampler.h"

NTL_CLIENT

class IdealGSieve{
    private:
        long n_;
        long m_;
        vector<ListPoint*> L;
        vector<ListPoint*> V;
        stack<ListPoint*> S;
        KleinSampler* sampler_;
        int64 goal_norm_;
        int64 min_norm_;
        long num_rots_;
        vec_int64 modf_;
        size_t simu_samp_;
        int concurrency_;
        int64 IdealGaussReduce_rot1(ListPoint* p);
        int64 IdealGaussReduce_rot2(ListPoint* p);
        int64 IdealGaussReduce_rot2_modified(ListPoint* p);
        int64 IdealGaussReduce_rot2_parallel();
        int64 IdealGaussReduce_anti_cyclic(ListPoint* p);
        int64 IdealGaussReduce_anti_cyclic_parallel();
        int64 IdealTripleReduce_rot2_2red(ListPoint* p);
        int64 IdealTripleReduce_rot2_2red_parallel();
        int64 IdealTripleReduce_anti_cyclic_2red(ListPoint* p, size_t &p_pos);
        bool Ideal_check_red3(const ListPoint *p1, const ListPoint *p2, const ListPoint *p3, ListPoint *p_new);
        int adjust_order_Ideal_check_red3(const ListPoint *p1, const ListPoint *p2, const ListPoint *p3, ListPoint *p_new);
        bool Ideal_check_red3_anti_cyclic(const ListPoint *p1, const ListPoint *p2, const ListPoint *p3, ListPoint *p_new);
        int64 IdealTripleReduce_rot2_rot1(ListPoint* p);
        int64 IdealTripleReduce_rot2_rot2(ListPoint* p);
        int64 IdealTripleReduce_rot2_rot3(ListPoint* p);
        int64 IdealTripleReduce_rot2_rot3_parallel();
        int64 IdealTripleReduce_anti_cyclic_rot1(ListPoint* p);
        int64 IdealTripleReduce_anti_cyclic_rot2(ListPoint* p);
        int64 IdealTripleReduce_anti_cyclic_rot3(ListPoint* p);
        // statistics
        long max_list_size_;
        long collisions_;
        long iterations_;
        long sample_vectors_;
        long reductions_;
        long updates_;
        vector<double> chk_time_;
        double timeL2V;
        double timeV2V;
        double timeV2L;
    public:
        ~IdealGSieve(){
            CleanUp();
        }
        void CleanUp();
        void Init(const mat_ZZ &B, KleinSampler* sampler, long num_rots, const vec_ZZ &modf);
        void SetGoalNorm(long norm){ goal_norm_ = norm; }
        void SetConcurrency(long num){ concurrency_ = num; }
        void SetSimultaneousSamples(long num){ simu_samp_ = num; }

        void IdealGaussSieve();
        void IdealGaussSieve_parallel();

        long getIterations(){ return iterations_; }
        long getCollisions(){ return collisions_; }
        long getListSize(){ return max_list_size_; }
        long getSampleVectors(){ return sample_vectors_; }
        long getReductions(){ return reductions_; }
        long getUpdates(){ return updates_; }
        ListPoint* getMinVec();
        vector<double> getChkTime(){ return chk_time_; }
        double getTimeL2V(){ return timeL2V; }
        double getTimeV2V(){ return timeV2V; }
        double getTimeV2L(){ return timeV2L; }

        void TestRotation(ListPoint *p1, long rep);
        void genRotations(const ListPoint *p1, vector<ListPoint*> &res);
        void genRotations_anti_cyclic(const ListPoint *p1, vector<ListPoint*> &res);
};

#endif