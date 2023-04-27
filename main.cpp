#include <sstream>
#include <fstream>
#include <iostream>
#include <unistd.h>
#include <getopt.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>
#include <chrono>

#include "sampler.h"
#include "gsieve.h"
#include "math.h"
#include "tool.h"

NTL_CLIENT

int main(int argc, char** argv){
    char *filename = NULL;
    long goal_norm = 0;
    int concurrency = 1;
    int simu_samp = 1;

    int option;
    while((option = getopt(argc, argv, "f:g:c:s:")) != -1){
        switch(option){
            case 'f': filename = optarg; break;
            case 'g': goal_norm = atol(optarg); break;
            case 'c': concurrency = atoi(optarg); break;
            case 's': simu_samp = atoi(optarg); break;
        }
    };

    ostringstream oss;
    ifstream ifs;
    mat_ZZ B;
    oss << filename;
    ifs.open(oss.str(), ios::in);
    if(ifs.is_open()){
        ifs >> B;
        ifs.close();
    }
    else cin >> B;

    // G_BKZ_FP(B, 0.99, 20);
    ZZ det2;
    LLL(det2, B, 99, 100, 0);

    int exp_time = 1;   // 繰り返し回数
    vector<string> index = {"time(ms)", "list_size", "sample_vectors", "collisions", "iterations", "norm"};        // 評価項目
    vector<int> list_concurrency = {12};
    vector<int> list_simu_samp = {10, 20, 40, 80, 160};
    int size_index = (int)index.size();
    chrono::system_clock::time_point start, end;

    #if 1
    // GS
    vector<double> rec1[size_index];
    for(int i = 0; i < exp_time; i++){
        KleinSampler sampler;
        GSieve gs;

        gs.Init(B, &sampler);

        gs.SetGoalNorm(goal_norm);

        start = chrono::system_clock::now();
        gs.GaussSieve();
        end = chrono::system_clock::now();
        double elapsed = chrono::duration_cast<chrono::microseconds>(end-start).count()/1000; 

        ListPoint* lp = gs.getMinVec();
        cout << "vec: " << lp->v << endl << "norm: " << lp->norm << "/" << sqrt(lp->norm) << endl;

        // vector<double> chk_time = gs.getChkTime();
        rec1[0].emplace_back(elapsed);
        rec1[1].emplace_back(gs.getListSize());
        rec1[2].emplace_back(gs.getSampleVectors());
        rec1[3].emplace_back(gs.getCollisions());
        rec1[4].emplace_back(gs.getIterations());
        rec1[5].emplace_back(gs.getMinVec()->norm);
    }
    string denotes = "dim," + to_string(B.NumRows());
    out2csv("GS_test", rec1, index, denotes);
    
    #endif

    #if 0
    // GS_P
    for(int i = 0; i < exp_time; i++){
        KleinSampler sampler;
        GSieve gs;

        gs.Init(B, &sampler);

        gs.SetGoalNorm(goal_norm);
        gs.SetConcurrency(concurrency);
        gs.SetSimultaneousSamples(simu_samp);

        start = chrono::system_clock::now();
        gs.GaussSieve_Parallel();
        end = chrono::system_clock::now();
        double elapsed = chrono::duration_cast<chrono::microseconds>(end-start).count()/1000;

        rec2[0].emplace_back(elapsed);
        rec2[1].emplace_back(chk_time[0]);
        rec2[2].emplace_back(chk_time[1]);
        rec2[3].emplace_back(chk_time[2]);
        rec2[4].emplace_back(gs->getSampleVectors());
        rec2[5].emplace_back(gs->getIterations());
        rec2[6].emplace_back(sqrt(to_double(gs->getMinNorm2())));
    }
    string denotes = "dim," + to_string(B.NumRows()) + ",sampling," + to_string(simu_samp) + ",concurrency," + to_string(list_concurrency[j]);
    out2csv("GS_P_test", rec2, index, denotes);

    #endif

    return 0;
}