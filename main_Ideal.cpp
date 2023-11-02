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
#include "Idealgsieve.h"
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
    vec_ZZ modf;
    long dim, index, seed;
    oss << filename;
    ifs.open(oss.str(), ios::in);
    if(ifs.is_open()){
        ifs >> B >> dim >> index >> seed;
        modf.SetLength(index);
        ifs >> modf;
        ifs.close();
    }
    else{
        cin >> B >> index;
        cin >> modf;
        dim = B.NumRows();
        seed = -1;
    }
    cout << "dim: " << dim << ", index: " << index << endl;
    long num_rots = dim + 1;
    cout << num_rots << endl;
    // anti-cyclicならdim
    long chk_dim = dim;
    while(chk_dim % 2 == 0) chk_dim /= 2;
    if(chk_dim == 1) num_rots--;
    
    // latticeの種類をチェック
    vector<string> lattice_name = {"cyclic", "anti-cyclic", "prime-cyclic", "common"};
    int lattice_type = check_lattice(modf);
    if(lattice_type == -1){
        cout << "error: can't check lattice type." << endl;
        return -1;
    }
    cout << lattice_name[lattice_type]  << " lattice" << endl;

    cout << "sv bound:" << 1.05 * svbound(B, to_RR(B[0][0])) << endl;;
    G_BKZ_FP(B, 0.99, 20);
    // ZZ det2;
    // LLL(det2, B, 99, 100, 0);

    chrono::system_clock::time_point start, end, ss, se;
    KleinSampler sampler;
    IdealGSieve Igs;
    start = chrono::system_clock::now();
    ss = chrono::system_clock::now();
    Igs.Init(B, &sampler, num_rots, modf);
    se = chrono::system_clock::now();
    Igs.SetGoalNorm(goal_norm);
    Igs.SetConcurrency(concurrency);
    Igs.SetSimultaneousSamples(simu_samp);
    // Igs.IdealGaussSieve_parallel();
    Igs.IdealGaussSieve();
    end = chrono::system_clock::now();
    double elapsed = chrono::duration_cast<chrono::microseconds>(end-start).count()/1000;
    double time_init = chrono::duration_cast<chrono::microseconds>(se-ss).count()/1000;

    ListPoint* lp = Igs.getMinVec();
    cout << "vector    : " << lp->v << endl;
    cout << "norm      : " << lp->norm << "/" << sqrt(lp->norm) << endl;
    cout << "times(ms) : " << elapsed << endl;
    cout << "init time : " << time_init << endl;
    cout << "size of L : " << Igs.getListSize() << endl;
    cout << "samples   : " << Igs.getSampleVectors() << endl;
    cout << "collisions: " << Igs.getCollisions() << endl;
    cout << "Iterations: " << Igs.getIterations() << endl;
    cout << "Reductions: " << Igs.getReductions() << endl;
    cout << "Updatess  : " << Igs.getUpdates() << endl;

    vector<string> item = {"norm", "time", "|L| ", "samp", "coll", "iter"};
    vector<double> record(item.size());
    record[0] = sqrt(lp->norm);
    record[1] = elapsed;
    record[2] = Igs.getListSize();
    record[3] = Igs.getSampleVectors();
    record[4] = Igs.getCollisions();
    record[4] = Igs.getIterations();
    string denotes = "dim," + to_string(dim) + ",index," + to_string(index);
    out2csv("IGS_test", record, item, denotes);

    // Igs.TestRotation(Igs.getMinVec(), 64);

    return 0;
}