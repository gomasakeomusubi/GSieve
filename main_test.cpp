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
#include <time.h>

#include "sample.h"
#include "sampler.h"
#include "gsieve.h"
#include "tool.h"

NTL_CLIENT

int main(int argc, char** argv){
    char *filename = NULL;
    ZZ goal_norm2; goal_norm2 = 0;
    int concurrency = 1;
    int simu_samp = 1;

    int option;
    while((option = getopt(argc, argv, "f:g:c:s:")) != -1){
        switch(option){
            case 'f': filename = optarg; break;
            case 'g': goal_norm2 = to_ZZ(atol(optarg)); break;
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
    
    KleinSampler sampler;
    sampler.Init(B);

    LatticeVector* lp;
    clock_t start = clock();
    for(int i = 0; i < 10000; i++){
        lp = sampler.Sample();
    }
    clock_t end = clock();

    cout << (double)(end-start)*1000/CLOCKS_PER_SEC << "ms" << endl;

    return 0;
}