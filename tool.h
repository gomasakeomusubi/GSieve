#ifndef INCLUDE_TOOL
#define INCLUDE_TOOL

#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <vector>
#include <queue>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

NTL_CLIENT

struct LatticeVector{
    vec_ZZ vec;
    ZZ norm2;
};

LatticeVector *newLatticeVector(int dim);
void deleteLatticeVector(vector<LatticeVector*> List, int index);
void deleteList(vector<LatticeVector*> List);
void deleteQueue(queue<LatticeVector*> Que);

ZZ norm2(vec_ZZ v);

bool reduceVector(LatticeVector *a, LatticeVector *b);

void printList(vector<LatticeVector*> L, string name);

void out2csv(string filename, vector<double> rec[], vector<string> index, string denotes);

#endif