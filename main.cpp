#include <sstream>
#include <fstream>
#include <iostream>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/mat_ZZ.h>
#include <chrono>

#include "sample.h"
#include "gsieve.h"
#include "tool.h"

NTL_CLIENT

int main(void){
    ostringstream oss;
    ifstream ifs;
    mat_ZZ B;

    oss.str("");
    oss << "input.txt";
    ifs.open(oss.str(), ios::in);
    ifs >> B;
    ifs.close();

    ZZ thresh;
    thresh = 1515;      // 閾値　BKZの出力を目安に設定
    // thresh = 1703;
    int simu_samp = 80;
    int concurrency = 12;
    LatticeVector *v = newLatticeVector(B.NumRows());


    int exp_time = 1;   // 繰り返し回数
    vector<string> index = {"time(ms)", "time1(ms)", "time2(ms)", "time3(ms)", "num_sample", "loop", "norm"};        // 評価項目
    vector<int> list_concurrency = {12};
    vector<int> list_simu_samp = {10, 20, 40, 80, 160};
    int size_index = (int)index.size();
    chrono::system_clock::time_point start, end;

    #if 1
    // GS
    vector<double> rec1[size_index];
    for(int i = 0; i < exp_time; i++){
        GSieve *gs = new GSieve(B, thresh, simu_samp, concurrency);
        vector<double> chk_time;
        long num_sample = 0, cnt = 0;

        start = chrono::system_clock::now();
        v = gs->GaussSieve(chk_time, num_sample, cnt);
        end = chrono::system_clock::now();
        double elapsed = chrono::duration_cast<chrono::microseconds>(end-start).count()/1000;

        rec1[0].emplace_back(elapsed);
        rec1[1].emplace_back(chk_time[0]);
        rec1[2].emplace_back(chk_time[1]);
        rec1[3].emplace_back(chk_time[2]);
        rec1[4].emplace_back(num_sample);
        rec1[5].emplace_back(cnt);
        rec1[6].emplace_back(sqrt(to_double(v->norm2)));

        delete gs;
    }
    string denotes = "dim," + to_string(B.NumRows());
    out2csv("GS_test", rec1, index, denotes);

    #endif

    #if 1
    // GS_P
    for(size_t j = 0; j < list_concurrency.size(); j++){
        vector<double> rec2[size_index];
        for(int i = 0; i < exp_time; i++){
            GSieve *gs = new GSieve(B, thresh, simu_samp, list_concurrency[j]);
            vector<double> chk_time;
            long num_sample = 0, cnt = 0;

            start = chrono::system_clock::now();
            v = gs->GaussSieve_Parallel(chk_time, num_sample, cnt);
            end = chrono::system_clock::now();
            double elapsed = chrono::duration_cast<chrono::microseconds>(end-start).count()/1000;

            rec2[0].emplace_back(elapsed);
            rec2[1].emplace_back(chk_time[0]);
            rec2[2].emplace_back(chk_time[1]);
            rec2[3].emplace_back(chk_time[2]);
            rec2[4].emplace_back(num_sample);
            rec2[5].emplace_back(cnt);
            rec2[6].emplace_back(sqrt(to_double(v->norm2)));

            delete gs;
        }
        string denotes = "dim," + to_string(B.NumRows()) + ",sampling," + to_string(simu_samp) + ",concurrency," + to_string(list_concurrency[j]);
        out2csv("GS_P_test", rec2, index, denotes);
    }
    #endif

    return 0;
}