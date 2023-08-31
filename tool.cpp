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

void CopyListPoint(ListPoint* new_p, const ListPoint* old_p){
    long dims = old_p->v.length();
    // new_p->v.SetLength(dims);
    for(int i = 0; i < dims; i++) new_p->v[i] = old_p->v[i];
    new_p->norm = old_p->norm;
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

bool reduceVectorDot(ListPoint *p1, const ListPoint *p2, int64 &dot_p1p2){
    long dims = p1->v.length();
    int64 dot = 0;

    if(p2->norm == 0){
        return false;
    }

    for(int i = 0; i < dims; i++){
        dot += p1->v[i] * p2->v[i];
    }
    dot_p1p2 = dot;
    if(2 * abs(dot) <= p2->norm){
        return false;
    }
    long q = round((double)dot / p2->norm);
    for(int i = 0; i < dims; i++){
        p1->v[i] -= q * p2->v[i];
    }
    p1->norm = p1->norm - 2 * q * dot + q * q * p2->norm;

    dot_p1p2 -= q * p2->norm;
    
    return true;
}

bool check_red2(const ListPoint *p1, const ListPoint *p2){
    // input: p1, p2 (p1 <= p2)
    // output: 2-reduced p1, p2
    //         if reduced 1 else 0.

    long dims = p1->v.length();
    int64 dot = 0;

    for(int i = 0; i < dims; i++){
        dot += p1->v[i] * p2->v[i];
    }
    if(2 * abs(dot) <= p2->norm){       // p1? p2?...
        return true;    // 2-reduced.
    }
    else return false;  // not 2-reduced.
}

bool check_red3(const ListPoint *p1, const ListPoint *p2, const ListPoint *p3, ListPoint *p_new){
    // input : p1, p2, p3 (p1 <= p2 <= p3)
    // output: 3-reduced p_new,
    //         if reduced 1 else 0.

    if(!check_red2(p1, p2)) return false;
    if(!check_red2(p1, p3)) return false;
    if(!check_red2(p2, p3)) return false;

    int64 dot12 = 0;
    int64 dot13 = 0;
    int64 dot23 = 0;
    long dims = p1->v.length();
    for(int i = 0; i < dims; i++){
        dot12 += p1->v[i] * p2->v[i];
        dot13 += p1->v[i] * p3->v[i];
        dot23 += p2->v[i] * p3->v[i];
    }

    int sign12 = 1, sign13 = 1, sign23 = 1;
    if(dot12 <= 0) sign12 = -1;
    if(dot13 <= 0) sign13 = -1;
    if(dot23 <= 0) sign23 = -1;

    if(sign12 * sign13 * sign23 == 1) return false;
    
    p_new->norm = 0;
    // p_newのノルムが最小になるように
    // dot の値により最大2パターン？
    for(int i = 0; i < dims; i++){
        p_new->v[i] = p1->v[i] - sign12 * p2->v[i] - sign13 * p3->v[i];
        p_new->norm += p_new->v[i] * p_new->v[i];
    }

    if(p_new->norm < p3->norm) return true;
    else return false;
}

void rotation_anti_cyclic(ListPoint* p1){
    long dims = p1->v.length();
    int64 tmp = p1->v[dims-1];

    for(int i = dims - 1; i > 0; i--){
        p1->v[i] = p1->v[i - 1];
    }
    p1->v[0] = -tmp;
}

// anti-cyclicはn回転で元のベクトル*-1に戻る
// prime cyclotonicはn+1回転で元のベクトルに戻る
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

void rotation(ListPoint* p1,  const ListPoint* p2, const vec_int64 &modf){
    long dims = p2->v.length();
    p1->v.SetLength(dims);
    p1->norm = 0;
    for(int i = dims - 1; i > 0; i--){
        p1->v[i] = p2->v[i-1] - p2->v[dims-1] * modf[i];
        p1->norm += p1->v[i] * p1->v[i];
    }
    p1->v[0] = -p2->v[dims-1] * modf[0];
    p1->norm += p1->v[0] * p1->v[0];
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

int calc_rot_num(const ListPoint* p, const vec_int64 &modf){
    int dim = p->v.length();
    ListPoint *lp = NewListPoint(dim);
    for(int i = 0; i < dim; i++) lp->v[i] = p->v[i];
    lp->norm = p->norm;

    int itr = 0;
    while(1){
        itr++;
        rotation(lp, modf);
        // rot(v) == -v
        bool same = true;
        for(int i = 0; i < dim; i++){
            if(lp->v[i] != -p->v[i]){
                same = false;
                break;
            }
        }
        if(same) break;
        // rot(v) == v
        same = true;
        for(int i = 0; i < dim; i++){
            if(lp->v[i] != p->v[i]){
                same = false;
                break;
            }
        }
        if(same) break;
    }

    return itr - 1;
}

// rotationでp以下のノルムをもつvectorが出たときのみreduce
bool IdealreduceVector(ListPoint *p1, const ListPoint *p2, const vec_int64 &modf, int num_rots){
    ListPoint *lp = NewListPoint(p2->v.length());
    CopyListPoint(lp, p2);

    bool vec_change = false;
    bool loop = true;
    while(loop){
        loop = false;
        for(int i = 0; i < num_rots; i++, rotation(lp, modf)){
            if(lp->norm > p1->norm){
                continue;
            }
            if(reduceVector(p1, lp)){
                vec_change = true;
                loop = true;
            }
        }
    }

    DeleteListPoint(lp);

    return vec_change;
}

// rotationで出たベクトルとpを比較してnormが長い方をreduce
// output: if p1(p2/no one) was reduced, 1(2/0).
int IdealreduceVector2(ListPoint *p1, ListPoint *p2, const vec_int64 &modf, int num_rots){
    long dims = p2->v.length();
    ListPoint* lp = NewListPoint(dims);
    CopyListPoint(lp, p2);
    int64 init_p1 = p1->norm;

    ListPoint* tmp = NewListPoint(dims);
    for(int i = 0; i < num_rots; i++){
        CopyListPoint(tmp, lp);
        if(tmp->norm > p1->norm){
            if(reduceVector(tmp, p1) && tmp->norm < p2->norm && tmp->norm > 0){
                CopyListPoint(p2, tmp);
                DeleteListPoint(tmp);
                DeleteListPoint(lp);
                return 2;
            }
        }
        else{
            if(reduceVector(p1, tmp)){
                if(p1->norm == 0 && i > 0){
                    CopyListPoint(p1, tmp);
                }
                else{
                    DeleteListPoint(tmp);
                    DeleteListPoint(lp);
                    return 1;
                }
            }
        }
        rotation(lp, modf);
    }

    DeleteListPoint(tmp);
    DeleteListPoint(lp);

    return 0;
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

void out2csv(string filename, vector<double> rec, vector<string> index, string denotes){
    ostringstream oss;
    ofstream ofs;
    oss.str("");
    oss << "output/" + filename + ".csv";
    ofs.open(oss.str(), ios::out | ios::app);

    int size_index = (int)index.size();
    ofs << denotes << endl;

    for(int i = 0; i < size_index; i++){
        ofs << index[i] << "," << rec[i] << endl;
    }
    ofs << endl;
    
    ofs.close();
}