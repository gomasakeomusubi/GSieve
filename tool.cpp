#include "tool.h"

NTL_CLIENT

ZZ norm2(vec_ZZ v){
    int n = v.length();
    ZZ sum; sum = 0;
    
    for(int i = 0; i < n; i++){
        sum += sqr(v[i]);
    }

    return sum;
}

LatticeVector *newLatticeVector(int dim){
    LatticeVector *plv = new LatticeVector;
    plv->vec.SetLength(dim);
    for(int i = 0; i < dim; i++){
        plv->vec[i] = 0;
    }
    plv->norm2 = 0;

    return plv;
}

void deleteLatticeVector(vector<LatticeVector*> List, int index){
    delete List[index];
    List.erase(List.begin()+index);
}

void deleteList(vector<LatticeVector*> List){
    int size_list = List.size();
    while(!List.empty()){
        delete List.back();
        List.pop_back();
    }
}

void deleteQueue(queue<LatticeVector*> Que){
    while(!Que.empty()){
        delete Que.front();
        Que.pop();
    }
}

bool reduceVector(LatticeVector *a, LatticeVector *b){
    ZZ q, dot;

    if(b->norm2 == 0){
        return false;
    }

    dot = a->vec*b->vec;
    if(2*dot <= b->norm2){
        return false;
    }
    q = round(to_double(dot)/to_double(b->norm2));
    a->vec = a->vec - q * b->vec; 
    a->norm2 = a->norm2 - 2*q*dot + q*q*b->norm2;
    return true;
}

void printList(vector<LatticeVector*> L, string name){
    cout << name << ": " << flush;
    for(int i = 0; i < L.size(); i++){
        cout << L[i]->norm2 << " " << flush;
    }
    cout << endl;
}

void out2csv(string filename, vector<double> rec[], vector<string> index, string denotes){
    ostringstream oss;
    ofstream ofs;
    oss.str("");
    oss << "output/" + filename + ".csv";
    ofs.open(oss.str(), ios::out | ios::app);

    int size_index = (int)index.size();
    int exp_time = rec[0].size();
    ofs << denotes << endl;
    ofs << "count,";
    for(int i = 1; i <= exp_time; i++) ofs << i << ",";
    ofs << endl;

    for(int i = 0; i < size_index; i++){
        ofs << index[i] << ",";
        for(int j = 0; j < exp_time; j++){
            ofs << rec[i][j] << ",";
        }
        ofs << endl;
    }
    
    ofs.close();
}