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
#include "tool.h"

NTL_CLIENT

template<typename T>
void print(vector<T> &v, string name){
    cout << name << ":" << endl;
    for(auto n: v) cout << n << " ";
    cout << endl;
}

int main(){
    int exp_time = 100;
    double all_time = 0;
    for(int d = 0; d < exp_time; d++){
        chrono::system_clock::time_point start, end;
        int N = 10000;
        vector<int> A(N);
        vector<int*> A_ad(N), B;
        for(int i = 0; i < N; i++) A[i] = i;
        for(int i = 0; i < N; i++) A_ad[i] = &i;
        // print(A, "A");

        srand(time(NULL));
        vector<bool> c(N, false);
        for(int i = 0; i < N; i++){
            if(rand() % 2) c[i] = true;
        }
        // print(c, "c");

        start = chrono::system_clock::now();
#ifdef DEBUG
        for(int i = 0; i < A_ad.size(); i++){
            if(c[i]){
                B.push_back(A_ad[i]);
                A_ad.erase(A_ad.begin() + i);
                c.erase(c.begin() + i);
                i--;
            }
        }
#endif
#ifndef DEBUG
        vector<int*> A_ad2;
        for(int i = 0; i < A_ad.size(); i++){
            if(c[i]){
                B.push_back(A_ad[i]);
            }
            else{
                A_ad2.push_back(A_ad[i]);
            }
        }
        A_ad.swap(A_ad2);
#endif
        end = chrono::system_clock::now();
        double elapsed = chrono::duration_cast<chrono::microseconds>(end-start).count()/1000;
        all_time += elapsed;
        // print(A, "A");
        // print(B, "B");
    }

    cout << all_time << endl;

    return 0;
}