#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <queue>
#include "tool.h"

NTL_CLIENT

int main(){
    vec_ZZ v;
    cin >> v;

    int dim = v.length();
    dim--;

    mat_ZZ rot;
    rot.SetDims(dim, dim);
    clear(rot);
    for(int i = 0; i < dim-1; i++) rot[i][i+1] = 1;
    for(int i = 0; i < dim; i++) rot[dim-1][i] = -v[i];

    cout << rot << endl;

    mat_ZZ rot_inv;
    rot_inv.SetDims(dim, dim);
    clear(rot_inv);
    for(int i = 0; i < dim-1; i++) rot_inv[i+1][i] = 1;
    for(int i = 0; i < dim; i++) rot_inv[0][i] = -v[dim-1-i];

    cout << rot_inv << endl;

    cout << rot * rot_inv << endl;

    return 0;
}