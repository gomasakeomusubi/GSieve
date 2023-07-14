#include "tool.h"

NTL_vector_impl(int64, vec_int64)
NTL_io_vector_impl(int64, vec_int64)
NTL_vector_impl(vec_int64, vec_vec_int64)
NTL_matrix_impl(int64, vec_int64, vec_vec_int64, mat_int64)
NTL_vector_impl(vec_double, vec_vec_double)
NTL_matrix_impl(double, vec_double, vec_vec_double, mat_double)

ListPoint* NewListPoint(long dims){
    ListPoint* p = new ListPoint;
    p->norm = 0;
    p->v.SetLength(dims);
    for(int i = 0; i < dims; i++){
        p->v[i] = 0;
    }

    return p;
}

void DeleteListPoint(ListPoint* p){
    delete p;
}

void VecZZToListPoint(const vec_ZZ &v, ListPoint* p){
    long dims = v.length();
    p->v.SetLength(dims);
    p->norm = 0;
    for(int i = 0; i < dims; i++){
        p->v[i] = to_long(v[i]);
        p->norm += p->v[i] * p->v[i];
    }
}

void MatInt64FromMatZZ(const mat_ZZ B, mat_int64 &A){
    long n_ = B.NumRows();
    long m_ = B.NumCols();
    A.SetDims(n_, m_);
    for(int i = 0; i < n_; i++){
        for(int j = 0; j < m_; j++){
            A[i][j] = to_long(B[i][j]);
        }
    }
}

void MatDoubleFromMatRR(const mat_RR B, mat_double &A){
    long n_ = B.NumRows();
    long m_ = B.NumCols();
    A.SetDims(n_, m_);
    for(int i = 0; i < n_; i++){
        for(int j = 0; j < m_; j++){
            A[i][j] = to_double(B[i][j]);
        }
    }
}

bool reduceVector(ListPoint *p1, const ListPoint *p2){
    long dims = p1->v.length();
    int64 dot = 0;

    if(p2->norm == 0){
        return false;
    }

    for(int i = 0; i < dims; i++){
        dot += p1->v[i] * p2->v[i];
    }
    if(2 * abs(dot) <= p2->norm){
        return false;
    }
    long q = round((double)dot / p2->norm);
    for(int i = 0; i < dims; i++){
        p1->v[i] -= q * p2->v[i];
    }
    p1->norm = p1->norm - 2 * q * dot + q * q * p2->norm;
    
    return true;
}

void rotation_anti_cyclic(ListPoint* p1){
    long dims = p1->v.length();
    int64 tmp = p1->v[dims-1];

    for(int i = dims - 1; i > 0; i--){
        p1->v[i] = p1->v[i - 1];
    }
    p1->v[0] = -tmp;
}

void rotation(ListPoint* p1, const vec_int64 &modf){
    long dims = p1->v.length();
    int64 norm = 0;
    int64 tmp = p1->v[dims-1];
    for(int i = dims - 1; i > 0; i--){
        p1->v[i] = p1->v[i-1] - tmp * modf[i];
        norm += p1->v[i] * p1->v[i];
    }
    p1->v[0] = -tmp * modf[0];
    norm += p1->v[0] * p1->v[0];
    p1->norm = norm;
}

void rotation_inv(ListPoint* p1, const vec_int64 &modf){
    long dims = p1->v.length();
    int64 norm = 0;
    int64 tmp = p1->v[0];
    for(int i = 0; i < dims - 1; i++){
        p1->v[i] = p1->v[i+1] - tmp * modf[dims-i-1];
        norm += p1->v[i] * p1->v[i];
    }
    p1->v[dims-1] = -tmp * modf[0];
    norm += p1->v[dims-1] * p1->v[dims-1];
    p1->norm = norm;
}

bool IdealreduceVector(ListPoint *p1, const ListPoint *p2, const vec_int64 &modf, int index){
    long n = p2->v.length();
    ListPoint *lp = NewListPoint(n);
    for(int i = 0; i < n; i++) lp->v[i] = p2->v[i];
    lp->norm = p2->norm;

    bool vec_change = false;

    for(int i = 0; i < index; i++, rotation(lp, modf)){
        if(lp->norm > p1->norm){
            continue;
            // break;
        }
        if(reduceVector(p1, lp)){
            vec_change = true;
        }
    }

    DeleteListPoint(lp);

    return vec_change;
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