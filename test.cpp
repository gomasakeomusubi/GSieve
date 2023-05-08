#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <NTL/RR.h>
#include <NTL/mat_ZZ.h>
#include <queue>

NTL_CLIENT

int main(){
    vector<int> A = {1, 2};
    vector<int> B = {1, 3, 4};

    cout << "before" << endl;
    cout << "A: ";
    for(auto v: A) cout << v << " ";
    cout << endl;
    cout << "B: ";
    for(auto v: B) cout << v << " ";
    cout << endl;

    vector<int> C;
    int i = 0, j = 0;
    while(i < A.size() && j < B.size()){
        if(A[i] < B[j]) {
            C.emplace_back(A[i]);
            i++;
        }
        else {
            C.emplace_back(B[j]);
            j++;
        }
    }

    while(i < A.size()){
        C.emplace_back(A[i]);
        i++;
    }
    while(j < B.size()){
        C.emplace_back(B[j]);
        j++;
    }

    queue<int> D;
    for(int i = 0; i < 10; i++) D.push(i);

    while(!D.empty()){
        cout << D.front() << " ";
        D.pop();
    }
    cout << endl;

    return 0;
}