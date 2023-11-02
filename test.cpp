#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <queue>
#include <ctime>
#include <chrono>
#include <cstdlib>
#include <unistd.h>
#include <algorithm>
#include <omp.h>
#include "tool.h"

NTL_CLIENT

template<typename T>
void print(vector<T> &v, string name){
    cout << name << ":" << endl;
    for(auto n: v) cout << n << " ";
    cout << endl;
}

int main(){
    // cout << "使用可能な最大スレッド数：" << omp_get_max_threads() << endl;

    // int i, a[100];
    // #pragma omp parallel for ordered num_threads(4)
    // for(i = 0; i < 100; i++){
    //     a[i] = 0;
    //     #pragma omp ordered
    //     printf("i=%d thread_num=%d\n", i, omp_get_thread_num());
    // }

    vector<ListPoint*> A;
    sort(A.begin(), A.end());

    return 0;
}