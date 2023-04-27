#ifndef __TOOL__
#define __TOOL__

#include <NTL/ZZ.h>
#include <NTL/vec_double.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/LLL.h>
#include <vector>
#include <queue>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

typedef int64_t int64;

NTL_CLIENT
NTL_vector_decl(int64, vec_int64)
NTL_io_vector_decl(int64, vec_int64)
NTL_vector_decl(vec_int64, vec_vec_int64)
NTL_matrix_decl(int64, vec_int64, vec_vec_int64, mat_int64)
NTL_vector_decl(vec_double, vec_vec_double)
NTL_matrix_decl(double, vec_double, vec_vec_double, mat_double)

struct ListPoint {
    vec_int64 v;
    int64 norm;
};

ListPoint* NewListPoint(long dims);
void DeleteListPoint(ListPoint* p);

void VecZZToListPoint(const vec_ZZ &v, ListPoint* p);
void MatInt64FromMatZZ(const mat_ZZ B, mat_int64 &A);
void MatDoubleFromMatRR(const mat_RR B, mat_double &A);

bool reduceVector(ListPoint *p1, const ListPoint *p2);

void out2csv(string filename, vector<double> rec[], vector<string> index, string denotes);

#endif