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

void printList(vector<ListPoint*> list, string listname){
    cout << listname << endl;
    for(auto v: list) cout << v->norm << " ";
    cout << endl;
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

int check_lattice(const vec_ZZ &modf){
    vector<bool> Lattice_type(4, true); // cyclic, anti-cyclic, prime-cyclotomic, others

    long dim = modf.length();
    if(modf[0] == 1) Lattice_type[0] = false;
    else if(modf[0] == -1) Lattice_type[1] = false, Lattice_type[2] = false; 
    for(int i = 1; i < dim-1; i++){
        if(modf[i] != 1) Lattice_type[2] = false;
        if(modf[i] != 0) Lattice_type[0] = false, Lattice_type[1] = false;
    }
    if(modf[dim-1] != 1) Lattice_type[0] = false, Lattice_type[1] = false, Lattice_type[2] = false;

    for(int i = 0; i < 4; i++){
        if(Lattice_type[i]) return i;
    }

    return -1;
}

double svbound(const mat_ZZ &B, const RR p){
    double alpha = sqrt(B.NumRows() / (double)(2 * M_PI * exp(1)));
    double beta = to_double(pow(p, 1 / to_RR(B.NumRows())));
    return alpha * beta;
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

bool reduceVectorDot(ListPoint *p1, const ListPoint *p2, int64 dot_p1p2){
    long dims = p1->v.length();
    int64 dot = 0;

    if(p2->norm == 0){
        return false;
    }

    if(2 * abs(dot_p1p2) <= p2->norm){
        return false;
    }
    long q = round((double)dot_p1p2 / p2->norm);
    for(int i = 0; i < dims; i++){
        p1->v[i] -= q * p2->v[i];
    }
    p1->norm = p1->norm - 2 * q * dot + q * q * p2->norm;
    
    return true;
}

bool check_red2(const ListPoint *p1, const ListPoint *p2){
    // input: p1, p2 (p1 <= p2)
    // output: 2-reduced p1, p2
    //         if unchanged 1 else 0.

    long dims = p1->v.length();
    int64 dot = 0;

    for(int i = 0; i < dims; i++){
        dot += p1->v[i] * p2->v[i];
    }
    if(2 * abs(dot) <= p1->norm){       // p1? p2?...
        return true;    // 2-reduced.
    }
    else{
        return false;   // not 2-reduced.
    }
}

// O(n)
bool check_red3(const ListPoint *p1, const ListPoint *p2, const ListPoint *p3, ListPoint *p_new){
    // input : p1, p2, p3 (p1 <= p2 <= p3)
    // output: 3-reduced p_new,
    //         if reduced 1 else 0.

    // if(!check_red2(p1, p2)) return false;
    // if(!check_red2(p1, p3)) return false;
    // if(!check_red2(p2, p3)) return false;

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

void direct_rotation(ListPoint* p1, const vec_int64 &modf, int num_of_rot){
    long dims = p1->v.length();
    int64 e = p1->v[dims - num_of_rot];
    ListPoint* tmp = NewListPoint(dims);
    for(int i = 0; i < dims; i++) tmp->v[i] = p1->v[i] - e;
    for(int i = 0; i < dims; i++) p1->v[i] = tmp->v[(i - num_of_rot + dims) % dims];
    p1->v[dims - num_of_rot] = - e;
    DeleteListPoint(tmp);
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

// rotationで出たベクトルとpを比較してnormが長い方をreduce
// output: if p1(p2/no one) was reduced, 1(2/0).
// O(n^2)
int IdealreduceVector(ListPoint *p1, ListPoint *p2, const vec_int64 &modf, int num_rots){
    long dims = p2->v.length();
    ListPoint* lp = NewListPoint(dims);
    CopyListPoint(lp, p2);

    ListPoint* tmp = NewListPoint(dims);
    for(int i = 0; i < num_rots; i++){
        CopyListPoint(tmp, lp);
        if(tmp->norm > p1->norm){
            if(reduceVector(tmp, p1) && (tmp->norm < p2->norm)){
                CopyListPoint(p2, tmp);
                DeleteListPoint(tmp);
                DeleteListPoint(lp);
                return 2;
            }
        }
        else{
            if(reduceVector(p1, tmp)){
                // ここを有効にするとcollisionが減るため実行時間は増加、リストサイズは増減？
                // 一応rot(p2)同士で確認するためリストサイズは減る...はず...
                // 53dim, prime では i > 1 が一番効率いい（謎）
                // if(p1->norm == 0 && (i > 0)){
                //     CopyListPoint(p1, tmp);
                // }
                // else{
                    DeleteListPoint(tmp);
                    DeleteListPoint(lp);
                    return 1;
                // }
            }
        }
        
        rotation(lp, modf);
        if(p2->norm > lp->norm){
            CopyListPoint(p2, lp);
            DeleteListPoint(tmp);
            DeleteListPoint(lp);
            return 2;
        }
    }

    DeleteListPoint(tmp);
    DeleteListPoint(lp);

    return 0;
}

// p1 > p2
bool IdealreduceVector_anti_cyclic(ListPoint *p1, const ListPoint *p2, int num_rots){
    ListPoint *lp = NewListPoint(p2->v.length());
    CopyListPoint(lp, p2);

    for(int i = 0; i < num_rots; i++){
        if(reduceVector(p1, lp)){
            DeleteListPoint(lp);
            return true;
        }
        rotation_anti_cyclic(lp);
    }

    DeleteListPoint(lp);

    return false;
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