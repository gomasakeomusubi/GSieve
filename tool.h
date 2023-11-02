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
void CopyListPoint(ListPoint* new_p, const ListPoint* old_p);
void DeleteListPoint(ListPoint* p);
void printList(vector<ListPoint*> List, string listname);

void VecZZToListPoint(const vec_ZZ &v, ListPoint* p);
void MatInt64FromMatZZ(const mat_ZZ B, mat_int64 &A);
void MatDoubleFromMatRR(const mat_RR B, mat_double &A);

int check_lattice(const vec_ZZ &modf);
double svbound(const mat_ZZ &B, const RR p);

bool reduceVector(ListPoint *p1, const ListPoint *p2);
bool reduceVectorDot(ListPoint *p1, const ListPoint *p2, int64 dot_p1p2);
bool check_red2(const ListPoint *p1, const ListPoint *p2);
bool check_red3(const ListPoint *p1, const ListPoint *p2, const ListPoint *p3, ListPoint *p_new);

void rotation_anti_cyclic(ListPoint* p1);
void rotation(ListPoint *p1, const vec_int64 &modf);
void rotation(ListPoint* p1,  const ListPoint* p2, const vec_int64 &modf);
void direct_rotation(ListPoint *p1, const vec_int64 &modf, int num_of_rot);
void rotation_inv(ListPoint *p1, const vec_int64 &modf);
int calc_rot_num(const ListPoint* p, const vec_int64 &modf);
bool Ideal_check_red2(const ListPoint *p1, const ListPoint *p2, int num_rots, const vec_int64 &modf);
int IdealreduceVector(ListPoint *p1, ListPoint *p2, const vec_int64 &modf, int num_rots);
bool IdealreduceVector_anti_cyclic(ListPoint *p1, const ListPoint *p2, int num_rots);

void out2csv(string filename, vector<double> rec[], vector<string> index, string denotes);
void out2csv(string filename, vector<double> rec, vector<string> index, string denotes);

#endif