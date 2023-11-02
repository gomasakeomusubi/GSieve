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

    G_BKZ_FP(B, 0.99, 20);
    // ZZ det2;
    // LLL(det2, B, 99, 100, 0);

    int exp_time = 1;   // 繰り返し回数
    // 評価項目
    // vector<string> index = {"time(ms)", "time1", "time2", "time3", "L->V",
    //                          "V->V", "V->L", "list_size", "sample_vectors",
    //                           "collisions", "iterations", "norm"};
    vector<int> list_concurrency = {1, 2, 4, 6, 8, 12, 16, 20};
    vector<int> list_simu_samp = {120, 160, 200, 240};
    // int size_index = (int)index.size();
    chrono::system_clock::time_point start, end;

    #if 1
    // GS
    vector<string> index_GS = {"time(ms)", "list_size", "sample_vectors",
                              "collisions", "iterations", "norm"};
    int size_index_GS = (int)index_GS.size();
    cout << "GS:" << endl;
    vector<double> rec1[size_index_GS];
    for(int i = 0; i < exp_time; i++){
        KleinSampler sampler;
        GSieve gs;

        start = chrono::system_clock::now();
        gs.Init(B, &sampler);
        gs.SetGoalNorm(goal_norm);
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

        cout << "Innerpros : " << gs.getInnerpros() << endl;
        cout << "Reductions: " << gs.getReductions() << endl;
    }
    string denotes1 = "dim," + to_string(B.NumRows());
    out2csv("GS_test", rec1, index_GS, denotes1);
    
    #endif

    #if 0
    // GS_P
    cout << "GS_P:" << endl;
    vector<double> rec2[size_index];
    for(int i = 0; i < 3; i++){
        KleinSampler sampler;
        GSieve gs;
        // simu_samp = list_simu_samp[i];

        gs.Init(B, &sampler);

        gs.SetGoalNorm(goal_norm);
        gs.SetConcurrency(concurrency);
        gs.SetSimultaneousSamples(simu_samp);

        start = chrono::system_clock::now();
        gs.GaussSieve_Parallel();
        end = chrono::system_clock::now();
        double elapsed = chrono::duration_cast<chrono::microseconds>(end-start).count()/1000;
        
        ListPoint* lp = gs.getMinVec();
        cout << "vec: " << lp->v << endl << "norm: " << lp->norm << "/" << sqrt(lp->norm) << endl;

        rec2[0].emplace_back(elapsed);
        rec2[1].emplace_back(gs.getChkTime()[0]);
        rec2[2].emplace_back(gs.getChkTime()[1]);
        rec2[3].emplace_back(gs.getChkTime()[2]);
        rec2[4].emplace_back(gs.getTimeL2V());
        rec2[5].emplace_back(gs.getTimeV2V());
        rec2[6].emplace_back(gs.getTimeV2L());
        rec2[7].emplace_back(gs.getListSize());
        rec2[8].emplace_back(gs.getSampleVectors());
        rec2[9].emplace_back(gs.getCollisions());
        rec2[10].emplace_back(gs.getIterations());
        rec2[11].emplace_back(gs.getMinVec()->norm);
    }
    string denotes2 = "dim," + to_string(B.NumRows()) + ",sampling," + to_string(simu_samp) + ",concurrency," + to_string(concurrency);
    out2csv("GS_P_test", rec2, index, denotes2);

    #endif

    return 0;
}