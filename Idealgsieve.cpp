#include "Idealgsieve.h"

ListPoint* IdealGSieve::getMinVec(){
    for(ListPoint* v: L){
        if(v->norm == min_norm_) return v;
    }
    for(ListPoint* v: V){
        if(v->norm == min_norm_) return v;
    }
    while(!S.empty()){
        if(S.top()->norm == min_norm_) return S.top();
        if(S.top()->norm < min_norm_) cout << S.top()->norm << endl;
        DeleteListPoint(S.top());
        S.pop();
    }

    cout << "error:min vec " << min_norm_ << " isn't here..." << endl;
    return L[0];
}

void IdealGSieve::CleanUp(){
    for(size_t i = 0; i < L.size(); i++) DeleteListPoint(L[i]);
    L.clear();
    for(size_t i = 0; i < V.size(); i++) DeleteListPoint(V[i]);
    V.clear();
    while(!S.empty()){
        DeleteListPoint(S.top());
        S.pop();
    }
}

void IdealGSieve::Init(const mat_ZZ &B, KleinSampler* sampler, long num_rots, const vec_ZZ &modf){
    n_ = B.NumRows();
    m_ = B.NumCols();
    sampler_ = sampler;
    num_rots_ = num_rots;
    iterations_ = 0;
    collisions_ = 0;
    sample_vectors_ = 0;
    reductions_ = 0;
    updates_ = 0;
    CleanUp();

    sampler_->Init(B);
    min_norm_ = to_long(B[0] * B[0]);

    modf_.SetLength(n_);
    for(int i = 0; i < n_; i++) modf_[i] = to_long(modf[i]);

    ListPoint *ltmp = NewListPoint(n_);
    VecZZToListPoint(B[0], ltmp);
    num_rots_ = calc_rot_num(ltmp, modf_);
    DeleteListPoint(ltmp);
    cout << "num_rots: " << num_rots_ << endl;
    
    ListPoint* p;
    int64 current_norm;
    for(int i = 0; i < n_; i++){
        p = NewListPoint(m_);
        VecZZToListPoint(B[i], p);
        // current_norm = IdealGaussReduce_rot1(p);
        current_norm = IdealGaussReduce_rot2(p);
        // current_norm = IdealGaussReduce_rot2_modified(p);
        // current_norm = IdealGaussReduce_anti_cyclic(p);
        // current_norm = IdealTripleReduce(p);
        // current_norm = IdealTripleReduce_allrot(p);
        // current_norm = IdealTripleReduce_rot2_rot1(p);
        // current_norm = IdealTripleReduce_rot2_rot2(p);
        // current_norm = IdealTripleReduce_rot2_rot3(p);
        // current_norm = IdealTripleReduce_anti_cyclic_rot0(p);
        // current_norm = IdealTripleReduce_anti_cyclic_rot1(p);
        // current_norm = IdealTripleReduce_anti_cyclic_rot2(p);
        // current_norm = IdealTripleReduce_anti_cyclic(p);
        if(current_norm == 0){
            collisions_++;
        }
        else if(current_norm < min_norm_){
            min_norm_ = current_norm;
        }
    }
    max_list_size_ = L.size();
    concurrency_ = 1;
    simu_samp_ = 1;
    goal_norm_ = 0;
    timeL2V = 0;
    timeV2V = 0;
    timeV2L = 0;
}

void IdealGSieve::TestRotation(ListPoint *p1, long rep){
    int dim = p1->v.length();
    ListPoint *lp = NewListPoint(dim);
    CopyListPoint(lp, p1);
    for(int i = 0; i < rep; i++){
        cout << i << " : " << lp->norm << " / " << lp->v << endl;
        rotation(lp, modf_);
    }
    cout << endl;

    DeleteListPoint(lp);
}

void IdealGSieve::genRotations(const ListPoint* p, vector<ListPoint*> &res){
    ListPoint* tmp = NewListPoint(m_);
    CopyListPoint(tmp, p);
    for(size_t i = 0; i < res.size(); i++){
        CopyListPoint(res[i], tmp);
        rotation(tmp, modf_);
    }
    DeleteListPoint(tmp);
}

void IdealGSieve::genRotations_anti_cyclic(const ListPoint* p, vector<ListPoint*> &res){
    ListPoint* tmp = NewListPoint(m_);
    CopyListPoint(tmp, p);
    for(size_t i = 0; i < res.size(); i++){
        CopyListPoint(res[i], tmp);
        rotation_anti_cyclic(tmp);
    }
    DeleteListPoint(tmp);
}

int64 IdealGSieve::IdealGaussReduce_rot1(ListPoint* p){
    vector<bool> vec_change_L(L.size(), false);
    vector<ListPoint*> rots(num_rots_+1);
    for(int i = 0; i < num_rots_+1; i++) rots[i] = NewListPoint(m_);

    bool loop = true;
    while(loop){
        loop = false;    
        // p <-> rot(L)
        bool vec_change = true;
        while(vec_change){
            vec_change = false;
            for(size_t id_L = 0; id_L < L.size(); id_L++){
                genRotations(L[id_L], rots);
                for(int rot = 0; rot <= num_rots_; rot++){
                    if(rots[rot]->norm < L[id_L]->norm) CopyListPoint(L[id_L], rots[rot]);
                    if(p->norm >= rots[rot]->norm){
                        if(reduceVector(p, rots[rot])){
                            vec_change = true;
                            updates_++;
                        }
                    }
                    else{
                        if(reduceVector(rots[rot], p)){
                            if(rots[rot]->norm < L[id_L]->norm){
                                CopyListPoint(L[id_L], rots[rot]);
                                vec_change_L[id_L] = true;
                                vec_change = true;
                                updates_++;
                            }
                        }
                    }
                    reductions_++;
                }
            }
        }

        // rot(p) -> L
        vec_change = true;
        while(vec_change){
            vec_change = false;
            genRotations(p, rots);
            for(size_t id_L = 0; id_L < L.size(); id_L++){
                for(int rot = 0; rot <= num_rots_; rot++){
                    if(rots[rot]->norm < p->norm) CopyListPoint(p, rots[rot]);
                    if(L[id_L]->norm >= rots[rot]->norm){
                        if(reduceVector(L[id_L], rots[rot])){
                            vec_change_L[id_L] = true;
                            vec_change = true;
                            loop = true;
                            updates_++;
                        }
                    }
                    else{
                        if(reduceVector(rots[rot], L[id_L])){
                            if(rots[rot]->norm < p->norm){
                                CopyListPoint(p, rots[rot]);
                                vec_change = true;
                                loop = true;
                                updates_++;
                            }
                        }
                    }
                    reductions_++;
                }
            }
        }
    }

    for(int i = 0; i < num_rots_+1; i++) DeleteListPoint(rots[i]);

    // L -> S
    vector<ListPoint*> L_tmp;
    for(size_t i = 0; i < L.size(); i++){
        if(vec_change_L[i]){
            if(L[i]->norm == 0){
                DeleteListPoint(L[i]);
                collisions_++;
            }
            else S.push(L[i]);
        }
        else L_tmp.emplace_back(L[i]);
    }
    L.swap(L_tmp);
    if(p->norm == 0){
        DeleteListPoint(p);
        return 0;
    }
    else{
        auto itr = lower_bound(L.begin(), L.end(), p, [](const ListPoint* i, const ListPoint* j){
            return i->norm < j->norm;
        });
        if(itr == L.end()) L.emplace_back(p);
        else L.insert(itr, p);
        return p->norm;
    }
}

int64 IdealGSieve::IdealGaussReduce_rot2(ListPoint* p){
    vector<bool> vec_change_L(L.size(), false);
        
    // rot(p) <-> rot(L)
    vector<ListPoint*> prots(num_rots_+1), lrots(num_rots_+1);
    for(int i = 0; i < num_rots_+1; i++){
        prots[i] = NewListPoint(m_);
        lrots[i] = NewListPoint(m_);
    }
    bool vec_change = true;
    while(vec_change){
        vec_change = false;
        if(p->norm == 0) break;
        genRotations(p, prots);
        for(size_t id_L = 0; id_L < L.size(); id_L++){
            if(p->norm == 0) break;
            if(L[id_L]->norm == 0) continue;
            genRotations(L[id_L], lrots);
            for(int rot1 = 0; rot1 <= num_rots_; rot1++){
                if(prots[rot1]->norm < p->norm) {CopyListPoint(p, prots[rot1]); updates_++;}
                for(int rot2 = 0; rot2 <= num_rots_; rot2++){
                    if(lrots[rot2]->norm < L[id_L]->norm) {CopyListPoint(L[id_L], lrots[rot2]); updates_++;}
                    if(prots[rot1]->norm >= lrots[rot2]->norm){
                        if(reduceVector(prots[rot1], lrots[rot2])){
                            if(prots[rot1]->norm == 0){
                                if(p->norm >= L[id_L]->norm) CopyListPoint(p, prots[rot1]);
                                else CopyListPoint(L[id_L], prots[rot1]);
                                rot1 = num_rots_; rot2 = num_rots_;     // これがないとp, l共に0になる可能性がある
                                updates_++;
                            }
                            else if(prots[rot1]->norm < p->norm){
                                CopyListPoint(p, prots[rot1]);
                                vec_change = true;
                                updates_++;
                            }
                        }
                    }
                    else{
                        if(reduceVector(lrots[rot2], prots[rot1])){
                            if(lrots[rot2]->norm == 0){
                                if(p->norm >= L[id_L]->norm) CopyListPoint(p, lrots[rot2]);
                                else CopyListPoint(L[id_L], lrots[rot2]);
                                rot1 = num_rots_; rot2 = num_rots_;
                                updates_++;
                            }
                            else if(lrots[rot2]->norm < L[id_L]->norm){
                                CopyListPoint(L[id_L], lrots[rot2]);
                                vec_change_L[id_L] = true;
                                vec_change = true;
                                updates_++;
                            }
                        }
                    }
                    reductions_++;
                }
            }
        }
    }

    for(int i = 0; i < num_rots_+1; i++){
        DeleteListPoint(prots[i]);
        DeleteListPoint(lrots[i]);
    }

    // L -> S
    vector<ListPoint*> L_tmp;
    for(size_t i = 0; i < L.size(); i++){
        if(L[i]->norm == 0){
            DeleteListPoint(L[i]);
            collisions_++;
        }
        else{
            if(vec_change_L[i]) S.push(L[i]);
            else L_tmp.emplace_back(L[i]);
        }
    }
    L.swap(L_tmp);
    if(p->norm == 0){
        DeleteListPoint(p);
        return 0;
    }
    else{
        auto itr = lower_bound(L.begin(), L.end(), p, [](const ListPoint* i, const ListPoint* j){
            return i->norm < j->norm;
        });
        if(itr == L.end()) L.emplace_back(p);
        else L.insert(itr, p);
        return p->norm;
    }
}

int64 IdealGSieve::IdealGaussReduce_rot2_modified(ListPoint* p){
    vector<bool> vec_change_L(L.size(), false);
        
    // rot(p) <-> rot(L)
    vector<ListPoint*> prots(num_rots_+1), lrots(num_rots_+1);
    for(int i = 0; i < num_rots_+1; i++){
        prots[i] = NewListPoint(m_);
        lrots[i] = NewListPoint(m_);
    }
    int64 dot, sum_p, sum_l, last_p, last_l, norm_p, norm_l;
    bool vec_change_p;
    bool vec_change = true;
    while(vec_change){
        if(p->norm == 0) break;
        vec_change = false;
        vec_change_p = false;
        genRotations(p, prots);
        // set p
        sum_p = 0;
        for(int d = 0; d < m_; d++) sum_p += p->v[d];
        norm_p = p->norm;
        last_p = p->v[m_-1];
        for(size_t id_L = 0; id_L < L.size(); id_L++){
            if(vec_change_p) break;
            if(p->norm == 0) break;
            if(L[id_L]->norm == 0) continue;
            genRotations(L[id_L], lrots);
            for(int rot1 = 0; rot1 <= num_rots_; rot1++){
                // set l
                sum_l = 0;
                for(int d = 0; d < m_; d++) sum_l += lrots[rot1]->v[d];
                norm_l = lrots[rot1]->norm;
                last_l = lrots[rot1]->v[m_-1];
                // calc dot
                dot = 0;
                for(int d = 0; d < m_; d++) dot += p->v[d] * lrots[rot1]->v[d];
                bool once = true;
                // if(prots[rot1]->norm < p->norm) {CopyListPoint(p, prots[rot1]); updates_++;}
                for(int rot2 = 0; rot2 <= num_rots_; rot2++){
                    // if(lrots[rot2]->norm < L[id_L]->norm) {CopyListPoint(L[id_L], lrots[rot2]); updates_++;}
                    if(norm_p >= norm_l){
                        if(reduceVectorDot(prots[rot2], lrots[(rot1+rot2) % (num_rots_+1)], dot)){
                            if(prots[rot2]->norm == 0){
                                if(p->norm >= L[id_L]->norm) {CopyListPoint(p, prots[rot2]); vec_change_p = true;}
                                else CopyListPoint(L[id_L], prots[rot2]);
                                rot1 = num_rots_; rot2 = num_rots_;
                                updates_++;
                            }
                            else if(prots[rot2]->norm < p->norm){
                                CopyListPoint(p, prots[rot2]);
                                vec_change = true;
                                rot1 = num_rots_; rot2 = num_rots_;
                                updates_++;
                                vec_change_p = true;
                            }
                        }
                    }
                    else{
                        if(reduceVectorDot(lrots[(rot1+rot2) % (num_rots_+1)], prots[rot2], dot)){
                            if(lrots[(rot1+rot2) % (num_rots_+1)]->norm == 0){
                                if(p->norm >= L[id_L]->norm) {CopyListPoint(p, lrots[(rot1+rot2) % (num_rots_+1)]); vec_change_p = true;}
                                else CopyListPoint(L[id_L], lrots[(rot1+rot2) % (num_rots_+1)]);
                                rot1 = num_rots_; rot2 = num_rots_;
                                updates_++;
                            }
                            else if(lrots[(rot1+rot2) % (num_rots_+1)]->norm < L[id_L]->norm){
                                CopyListPoint(L[id_L], lrots[(rot1+rot2) % (num_rots_+1)]);
                                vec_change_L[id_L] = true;
                                vec_change = true;
                                rot1 = num_rots_; rot2 = num_rots_;
                                updates_++;
                            }
                        }
                    }
                    if(once){
                        // cout << p->norm << "/" << p->v << endl;
                        // cout << L[id_L]->norm << "/" << L[id_L]->v << endl;
                        // int64 tmp = 0;
                        // for(int _ = 0; _ < m_; _++) tmp += p->v[_] * L[id_L]->v[_];
                        // cout << "pl: " << tmp << endl;
                        // cout << dot << "/" << norm_p << "/" << norm_l << "/" \
                        // << sum_p << "/" << sum_l << "/" << last_p << "/" << last_l << endl;
                        once = false;
                    }
                    if(rot2 < num_rots_){
                        // dot = <x*v, x*l>
                        dot = dot - last_l * sum_p - last_p * sum_l + (m_ + 1) * last_p * last_l;
                        int64 tmp = 0;
                        for(int d = 0; d < m_; d++) tmp += prots[rot2+1]->v[d] * lrots[(rot1+rot2+1) % (num_rots_+1)]->v[d];
                        cout << "the: " << tmp << ", pre: " << dot << endl;
                        // |v| = |x*v|, |l| = |x*l|
                        norm_p = prots[rot2+1]->norm;
                        norm_l = lrots[(rot1+rot2+1) % (num_rots_+1)]->norm;
                        // sum(v) = sum(x*v), sum(l) = sum(x*l)
                        sum_p -= (m_ + 1) * last_p;
                        sum_l -= (m_ + 1) * last_l;
                        // last(v) = last(x*v), last(l) = last(x*l)
                        last_p = prots[rot2+1]->v[m_-1];
                        last_l = lrots[(rot1+rot2+1) % (num_rots_+1)]->v[m_-1];
                        // cout << dot << "/" << norm_p << "/" << norm_l << "/" \
                        // << sum_p << "/" << sum_l << "/" << last_p << "/" << last_l << endl;
                    }
                    
                    reductions_++;
                    if(reductions_ == 200) exit(0);
                }
            }
        }
    }

    for(int i = 0; i < num_rots_+1; i++){
        DeleteListPoint(prots[i]);
        DeleteListPoint(lrots[i]);
    }

    // L -> S
    vector<ListPoint*> L_tmp;
    for(size_t i = 0; i < L.size(); i++){
        if(L[i]->norm == 0){
            DeleteListPoint(L[i]);
            collisions_++;
        }
        else{
            if(vec_change_L[i]) S.push(L[i]);
            else L_tmp.emplace_back(L[i]);
        }
    }
    L.swap(L_tmp);
    if(p->norm == 0){
        DeleteListPoint(p);
        return 0;
    }
    else{
        auto itr = lower_bound(L.begin(), L.end(), p, [](const ListPoint* i, const ListPoint* j){
            return i->norm < j->norm;
        });
        if(itr == L.end()) L.emplace_back(p);
        else L.insert(itr, p);
        return p->norm;
    }
}

int64 IdealGSieve::IdealGaussReduce_rot2_parallel(){   
    vector<bool> vec_change_L(L.size(), false);
    vector<bool> vec_change_V(V.size(), false);
    
    // V <- rot(L)
    #pragma omp parallel for
    for(size_t id_V = 0; id_V < V.size(); id_V++){
        vector<ListPoint*> prots(num_rots_+1), lrots(num_rots_+1);
        for(int i = 0; i < num_rots_+1; i++){
            prots[i] = NewListPoint(m_);
            lrots[i] = NewListPoint(m_);
        }
        bool vec_change = true;
        while(vec_change){
            vec_change = false;
            genRotations(V[id_V], prots);
            for(size_t id_L = 0; id_L < L.size(); id_L++){
                genRotations(L[id_L], lrots);
                for(int rot1 = 0; rot1 <= num_rots_; rot1++){
                    if(prots[rot1]->norm < V[id_V]->norm) CopyListPoint(V[id_V], prots[rot1]);
                    for(int rot2 = 0; rot2 <= num_rots_; rot2++){
                        if(prots[rot1]->norm < lrots[rot2]->norm) continue;
                        if(reduceVector(prots[rot1], lrots[rot2])){
                            if(prots[rot1]->norm < V[id_V]->norm){
                                CopyListPoint(V[id_V], prots[rot1]);
                                vec_change = true;
                            }
                        }
                    }
                }
            }
        }
        for(int i = 0; i < num_rots_+1; i++){
            DeleteListPoint(prots[i]);
            DeleteListPoint(lrots[i]);
        }
    }

    // V <- rot(V)
    vector<ListPoint*> prots1(num_rots_+1), prots2(num_rots_+1);
    for(int i = 0; i < num_rots_+1; i++){
        prots1[i] = NewListPoint(m_);
        prots2[i] = NewListPoint(m_);
    }
    for(size_t i = 0; i < V.size(); i++){
        if(vec_change_V[i]) continue;
        genRotations(V[i], prots1);
        for(size_t j = i + 1; j < V.size(); j++){
            if(vec_change_V[j]) continue;
            genRotations(V[j], prots2);
            for(int rot1 = 0; rot1 <= num_rots_; rot1++){
                if(prots1[rot1]->norm < V[i]->norm) CopyListPoint(V[i], prots1[rot1]);
                for(int rot2 = 0; rot2 <= num_rots_; rot2++){
                    if(prots2[rot2]->norm < V[j]->norm) CopyListPoint(V[j], prots2[rot2]);
                    if(prots1[rot1]->norm >= prots2[rot2]->norm){
                        if(reduceVector(prots1[rot1], prots2[rot2])){
                            if(prots1[rot1]->norm < V[i]->norm){
                                CopyListPoint(V[i], prots1[rot1]);
                                vec_change_V[i] = true;
                            }
                        }
                    }
                    else{
                        if(reduceVector(prots2[rot2], prots1[rot1])){
                            if(prots2[rot2]->norm < V[j]->norm){
                                CopyListPoint(V[j], prots2[rot2]);
                                vec_change_V[j] = true;
                            }
                        }
                    }
                }
            }
        }
    }
    for(int i = 0; i < num_rots_+1; i++){
        DeleteListPoint(prots1[i]);
        DeleteListPoint(prots2[i]);
    }

    // V -> S
    vector<ListPoint*> tmp_V;
    for(size_t i = 0; i < V.size(); i++){
        if(V[i]->norm == 0){
            DeleteListPoint(V[i]);
            collisions_++;
        }
        else{
            if(vec_change_V[i]) S.push(V[i]);
            else tmp_V.emplace_back(V[i]);
        }
    }
    V.swap(tmp_V);
    sort(V.begin(), V.end(), [](const ListPoint* i, const ListPoint* j){
        return i->norm < j->norm;
    });

    // rot(V) -> L
    #pragma omp parallel for
    for(size_t id_L = 0; id_L < L.size(); id_L++){
        vector<ListPoint*> prots(num_rots_+1), lrots(num_rots_+1);
        for(int i = 0; i < num_rots_+1; i++){
            prots[i] = NewListPoint(m_);
            lrots[i] = NewListPoint(m_);
        }
        genRotations(L[id_L], lrots);
        for(size_t id_V = 0; id_V < V.size(); id_V++){
            genRotations(V[id_V], prots);
            for(int rot1 = 0; rot1 <= num_rots_; rot1++){
                if(lrots[rot1]->norm < L[id_L]->norm) CopyListPoint(L[id_L], lrots[rot1]);
                for(int rot2 = 0; rot2 <= num_rots_; rot2++){
                    if(lrots[rot1]->norm < prots[rot2]->norm) continue;
                    if(reduceVector(lrots[rot1], prots[rot2])){
                        if(lrots[rot1]->norm < L[id_L]->norm){
                            CopyListPoint(L[id_L], lrots[rot1]);
                            vec_change_L[id_L] = true;
                        }
                    }
                }
            }
        }
        for(int i = 0; i < num_rots_+1; i++){
            DeleteListPoint(prots[i]);
            DeleteListPoint(lrots[i]);
        }
    }
    
    // L -> S
    vector<ListPoint*> L_tmp;
    for(size_t i = 0; i < L.size(); i++){
        if(L[i]->norm == 0){
            DeleteListPoint(L[i]);
            collisions_++;
        }
        else{
            if(vec_change_L[i]) S.push(L[i]);
            else L_tmp.emplace_back(L[i]);
        }
    }
    L.swap(L_tmp);
    
    // L += V
    if(V.empty()) return min_norm_;
    vector<ListPoint*> tmp_L;
    size_t i = 0, j = 0;
    int64 current_norm = V[0]->norm;
    while(!(i == L.size() && j == V.size())){
        while(i < L.size() && (j == V.size() || L[i]->norm <= V[j]->norm)){
            tmp_L.emplace_back(L[i]);
            i++;
        }
        while(j < V.size() && (i == L.size() || L[i]->norm > V[j]->norm)){
            tmp_L.emplace_back(V[j]);
            j++;
        }
    }
    V.clear();
    L.swap(tmp_L);

    return current_norm;    
}

int64 IdealGSieve::IdealGaussReduce_anti_cyclic(ListPoint* p){   
    // <x^i*p,l> = <p,x^(num_rots-i)*l>より各<p,l>のペアで片方のrotationだけ試せばよい
    vector<ListPoint*> rots(num_rots_+1);
    for(int i = 0; i < num_rots_+1; i++) rots[i] = NewListPoint(m_);

    // p <- rot(L)
    bool vec_change = true;
    size_t id_L;
    while(vec_change){
        vec_change = false;
        for(id_L = 0; id_L < L.size(); id_L++){
            if(p->norm < L[id_L]->norm) break;
            genRotations_anti_cyclic(L[id_L], rots);
            for(int rot = 0; rot <= num_rots_; rot++){
                if(reduceVector(p, rots[rot])){
                    vec_change = true;
                }
            }
        }
    }

    if(p->norm == 0){
        for(int i = 0; i < num_rots_+1; i++) DeleteListPoint(rots[i]);
        DeleteListPoint(p);
        return 0;
    }

    // insert p into L
    L.insert(L.begin()+id_L, p);
    id_L++;

    // rot(p) -> L
    vector<bool> vec_change_L(L.size(), false);
    genRotations(p, rots);
    for(; id_L < L.size(); id_L++){
        for(int rot = 0; rot <= num_rots_; rot++){
            if(reduceVector(L[id_L], rots[rot])){
                vec_change_L[id_L] = true;
            }
        }
    }

    for(int i = 0; i < num_rots_+1; i++) DeleteListPoint(rots[i]);

    // L -> S
    vector<ListPoint*> tmp_L;
    for(size_t i = 0; i < L.size(); i++){
        if(L[i]->norm == 0){
            DeleteListPoint(L[i]);
            collisions_++;
        }
        else{
            if(vec_change_L[i]) S.push(L[i]);
            else tmp_L.emplace_back(L[i]);
        }
    }
    L.swap(tmp_L);

    return p->norm;
}

int64 IdealGSieve::IdealGaussReduce_anti_cyclic_parallel(){   
    // V <- rot(L)
    vector<int> thresh_V(V.size(), 0);
    #pragma omp parallel for
    for(size_t id_V = 0; id_V < V.size(); id_V++){
        vector<ListPoint*> rots(num_rots_+1);
        for(int i = 0; i < num_rots_+1; i++) rots[i] = NewListPoint(m_);
        bool vec_change = true;
        size_t id_L;
        while(vec_change){
            vec_change = false;
            for(id_L = 0; id_L < L.size(); id_L++){
                if(V[id_V]->norm < L[id_L]->norm) break;
                genRotations_anti_cyclic(L[id_L], rots);
                for(int rot = 0; rot <= num_rots_; rot++){
                    if(reduceVector(V[id_V], rots[rot])){
                        vec_change = true;
                    }
                }
            }
        }
        for(int i = 0; i < num_rots_+1; i++) DeleteListPoint(rots[i]);
        thresh_V[id_V] = id_L;
    }

    // rot(V) -> L
    vector<bool> vec_change_L(L.size(), false);
    #pragma omp parallel for
    for(size_t id_L = 0; id_L < L.size(); id_L++){
        vector<ListPoint*> rots(num_rots_+1);
        for(int i = 0; i < num_rots_+1; i++) rots[i] = NewListPoint(m_);
        for(size_t id_V = 0; id_V < V.size(); id_V++){
            if(L[id_L]->norm < V[id_V]->norm) continue;
            genRotations(V[id_V], rots);
            for(int rot = 0; rot <= num_rots_; rot++){
                if(reduceVector(L[id_L], rots[rot])){
                    vec_change_L[id_L] = true;
                }
            }
        }
        for(int i = 0; i < num_rots_+1; i++) DeleteListPoint(rots[i]);
    }

    // V <- rot(V)
    vector<bool> vec_change_V(V.size(), false);
    vector<ListPoint*> rots(num_rots_+1);
    for(int i = 0; i < num_rots_+1; i++) rots[i] = NewListPoint(m_);
    bool vec_change = true;
    while(vec_change){
        vec_change = false;
        for(size_t i = 0; i < V.size(); i++){
            if(V[i]->norm == 0) continue;
            genRotations_anti_cyclic(V[i], rots);
            for(size_t j = i + 1; j < V.size(); j++){
                if(V[i]->norm == 0 || V[j]->norm == 0) continue;
                for(int rot = 0; rot <= num_rots_; rot++){
                    if(V[i]->norm >= V[j]->norm){
                        if(reduceVector(rots[rot], V[j])){
                            if(rots[rot]->norm == 0){
                                V[i]->norm = 0;
                                break;
                            }
                            CopyListPoint(V[i], rots[rot]);
                            vec_change_V[i] = true;
                            vec_change = true;
                        }
                    }
                    else{
                        if(reduceVector(V[j], rots[rot])){
                            if(V[j]->norm == 0) break;
                            vec_change_V[j] = true;
                            vec_change = true;
                        }
                    }
                }
            }
        }
    }
    for(int i = 0; i < num_rots_+1; i++) DeleteListPoint(rots[i]);

    // V -> S
    vector<ListPoint*> tmp_V;
    for(size_t i = 0; i < V.size(); i++){
        if(V[i]->norm == 0){
            DeleteListPoint(V[i]);
            collisions_++;
        }
        else{
            if(vec_change_V[i]) S.push(V[i]);
            else tmp_V.emplace_back(V[i]);
        }
    }
    V.swap(tmp_V);

    // L -> S
    vector<ListPoint*> tmp_L;
    for(size_t i = 0; i < L.size(); i++){
        if(L[i]->norm == 0){
            DeleteListPoint(L[i]);
            collisions_++;
        }
        else{
            if(vec_change_L[i]) S.push(L[i]);
            else tmp_L.emplace_back(L[i]);
        }
    }
    L.swap(tmp_L);

    // L += V
    // cout << "L += V" << endl;
    if(V.empty()) return min_norm_;
    vector<ListPoint*> tmp_LV;
    size_t i = 0, j = 0;
    int64 current_norm = V[0]->norm;
    while(!(i == L.size() && j == V.size())){
        while(i < L.size() && (j == V.size() || L[i]->norm <= V[j]->norm)){
            tmp_LV.emplace_back(L[i]);
            i++;
        }
        while(j < V.size() && (i == L.size() || L[i]->norm > V[j]->norm)){
            tmp_LV.emplace_back(V[j]);
            j++;
        }
    }
    V.clear();
    L.swap(tmp_LV);

    return current_norm;
}

int64 IdealGSieve::IdealTripleReduce_rot2_2red(ListPoint* p){
    vector<bool> vec_change_L(L.size(), false);
        
    // rot(p) <-> rot(L)
    vector<ListPoint*> prots(num_rots_+1), lrots(num_rots_+1);
    for(int i = 0; i < num_rots_+1; i++){
        prots[i] = NewListPoint(m_);
        lrots[i] = NewListPoint(m_);
    }
    bool vec_change = true;
    while(vec_change){
        vec_change = false;
        if(p->norm == 0) break;
        genRotations(p, prots);
        for(size_t id_L = 0; id_L < L.size(); id_L++){
            if(p->norm == 0) break;
            if(L[id_L]->norm == 0) continue;
            genRotations(L[id_L], lrots);
            for(int rot1 = 0; rot1 <= num_rots_; rot1++){
                if(prots[rot1]->norm < p->norm) CopyListPoint(p, prots[rot1]);
                for(int rot2 = 0; rot2 <= num_rots_; rot2++){
                    if(lrots[rot2]->norm < L[id_L]->norm) CopyListPoint(L[id_L], lrots[rot2]);
                    if(prots[rot1]->norm >= lrots[rot2]->norm){
                        if(reduceVector(prots[rot1], lrots[rot2])){
                            if(prots[rot1]->norm == 0){
                                if(p->norm >= L[id_L]->norm) CopyListPoint(p, prots[rot1]);
                                else CopyListPoint(L[id_L], prots[rot1]);
                                rot1 = num_rots_; rot2 = num_rots_;
                            }
                            else if(prots[rot1]->norm < p->norm){
                                CopyListPoint(p, prots[rot1]);
                                vec_change = true;
                            }
                        }
                    }
                    else{
                        if(reduceVector(lrots[rot2], prots[rot1])){
                            if(lrots[rot2]->norm == 0){
                                if(p->norm >= L[id_L]->norm) CopyListPoint(p, lrots[rot2]);
                                else CopyListPoint(L[id_L], lrots[rot2]);
                                rot1 = num_rots_; rot2 = num_rots_;
                            }
                            else if(lrots[rot2]->norm < L[id_L]->norm){
                                CopyListPoint(L[id_L], lrots[rot2]);
                                vec_change_L[id_L] = true;
                                vec_change = true;
                            }
                        }
                    }
                }
            }
        }
    }

    for(int i = 0; i < num_rots_+1; i++){
        DeleteListPoint(prots[i]);
        DeleteListPoint(lrots[i]);
    }

    // L -> S
    vector<ListPoint*> L_tmp;
    for(size_t i = 0; i < L.size(); i++){
        if(L[i]->norm == 0){
            DeleteListPoint(L[i]);
            collisions_++;
        }
        else{
            if(vec_change_L[i]) S.push(L[i]);
            else L_tmp.emplace_back(L[i]);
        }
    }
    L.swap(L_tmp);
    if(p->norm == 0){
        DeleteListPoint(p);
        return 0;
    }
    else return p->norm;
}

int64 IdealGSieve::IdealTripleReduce_rot2_2red_parallel(){
    vector<bool> vec_change_L(L.size(), false);
    vector<bool> vec_change_V(V.size(), false);
    bool loop = true;
    while(loop){
        loop = false;
        // V <- rot(L)
        #pragma omp parallel for
        for(size_t id_V = 0; id_V < V.size(); id_V++){
            if(vec_change_V[id_V]) continue;
            vector<ListPoint*> prots(num_rots_+1), lrots(num_rots_+1);
            for(int i = 0; i < num_rots_+1; i++){
                prots[i] = NewListPoint(m_);
                lrots[i] = NewListPoint(m_);
            }
            bool vec_change = true;
            while(vec_change){
                vec_change = false;
                genRotations(V[id_V], prots);
                for(size_t id_L = 0; id_L < L.size(); id_L++){
                    if(vec_change_L[id_L]) continue;
                    genRotations(L[id_L], lrots);
                    for(int rot1 = 0; rot1 <= num_rots_; rot1++){
                        if(prots[rot1]->norm < V[id_V]->norm) CopyListPoint(V[id_V], prots[rot1]);
                        for(int rot2 = 0; rot2 <= num_rots_; rot2++){
                            if(prots[rot1]->norm < lrots[rot2]->norm) continue;
                            if(reduceVector(prots[rot1], lrots[rot2])){
                                if(prots[rot1]->norm < V[id_V]->norm){
                                    CopyListPoint(V[id_V], prots[rot1]);
                                    vec_change = true;
                                }
                            }
                        }
                    }
                }
            }
            for(int i = 0; i < num_rots_+1; i++){
                DeleteListPoint(prots[i]);
                DeleteListPoint(lrots[i]);
            }
        }

        // V <- rot(V)
        vector<ListPoint*> prots1(num_rots_+1), prots2(num_rots_+1);
        for(int i = 0; i < num_rots_+1; i++){
            prots1[i] = NewListPoint(m_);
            prots2[i] = NewListPoint(m_);
        }
        bool loop2 = true;
        while(loop2){
            loop2 = false;
            for(size_t i = 0; i < V.size(); i++){
                if(vec_change_V[i]) continue;
                genRotations(V[i], prots1);
                for(size_t j = 0; j < V.size(); j++){
                    if(vec_change_V[j]) continue;
                    if(i == j) continue;
                    genRotations(V[j], prots2);
                    for(int rot1 = 0; rot1 <= num_rots_; rot1++){
                        if(prots1[rot1]->norm < V[i]->norm) CopyListPoint(V[i], prots1[rot1]);
                        for(int rot2 = 0; rot2 <= num_rots_; rot2++){
                            if(prots2[rot2]->norm < V[j]->norm) CopyListPoint(V[j], prots2[rot2]);
                            if(prots1[rot1]->norm >= prots2[rot2]->norm){
                                if(reduceVector(prots1[rot1], prots2[rot2])){
                                    if(prots1[rot1]->norm < V[i]->norm){
                                        CopyListPoint(V[i], prots1[rot1]);
                                        vec_change_V[i] = true;
                                        loop2 = true;
                                    }
                                }
                            }
                            else{
                                if(reduceVector(prots2[rot2], prots1[rot1])){
                                    if(prots2[rot2]->norm < V[j]->norm){
                                        CopyListPoint(V[j], prots2[rot2]);
                                        vec_change_V[j] = true;
                                        loop2 = true;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        for(int i = 0; i < num_rots_+1; i++){
            DeleteListPoint(prots1[i]);
            DeleteListPoint(prots2[i]);
        }

        // rot(V) -> L
        #pragma omp parallel for
        for(size_t id_L = 0; id_L < L.size(); id_L++){
            if(vec_change_L[id_L]) continue;
            vector<ListPoint*> prots(num_rots_+1), lrots(num_rots_+1);
            for(int i = 0; i < num_rots_+1; i++){
                prots[i] = NewListPoint(m_);
                lrots[i] = NewListPoint(m_);
            }
            bool vec_change = true;
            while(vec_change){
                vec_change = false;
                genRotations(L[id_L], lrots);
                for(size_t id_V = 0; id_V < V.size(); id_V++){
                    if(vec_change_V[id_V]) continue;
                    genRotations(V[id_V], prots);
                    for(int rot1 = 0; rot1 <= num_rots_; rot1++){
                        if(lrots[rot1]->norm < L[id_L]->norm) CopyListPoint(L[id_L], lrots[rot1]);
                        for(int rot2 = 0; rot2 <= num_rots_; rot2++){
                            if(lrots[rot1]->norm < prots[rot2]->norm) continue;
                            if(reduceVector(lrots[rot1], prots[rot2])){
                                if(lrots[rot1]->norm < L[id_L]->norm){
                                    CopyListPoint(L[id_L], lrots[rot1]);
                                    vec_change_L[id_L] = true;
                                    vec_change = true;
                                    loop = true;
                                }
                            }
                        }
                    }
                }
            }
            for(int i = 0; i < num_rots_+1; i++){
                DeleteListPoint(prots[i]);
                DeleteListPoint(lrots[i]);
            }
        }
    }

    // V -> S
    vector<ListPoint*> tmp_V;
    for(size_t i = 0; i < V.size(); i++){
        if(V[i]->norm == 0){
            DeleteListPoint(V[i]);
            collisions_++;
        }
        else{
            if(vec_change_V[i]) S.push(V[i]);
            else tmp_V.emplace_back(V[i]);
        }
    }
    V.swap(tmp_V);
    sort(V.begin(), V.end(), [](const ListPoint* i, const ListPoint* j){
        return i->norm < j->norm;
    });
    
    // L -> S
    vector<ListPoint*> L_tmp;
    for(size_t i = 0; i < L.size(); i++){
        if(L[i]->norm == 0){
            DeleteListPoint(L[i]);
            collisions_++;
        }
        else{
            if(vec_change_L[i]) S.push(L[i]);
            else L_tmp.emplace_back(L[i]);
        }
    }
    L.swap(L_tmp);

    if(!V.empty()) return V[0]->norm;
    else return 0;    
}
// ok

int64 IdealGSieve::IdealTripleReduce_anti_cyclic_2red(ListPoint* p, size_t &p_pos){
    vector<ListPoint*> rots(num_rots_+1);
    for(int i = 0; i < num_rots_+1; i++) rots[i] = NewListPoint(m_);
    
    // p <- rot(L)
    bool vec_change = true;
    size_t id_L;
    while(vec_change){
        vec_change = false;
        for(id_L = 0; id_L < L.size(); id_L++){
            if(p->norm < L[id_L]->norm) break;
            genRotations_anti_cyclic(L[id_L], rots);
            for(int rot = 0; rot <= num_rots_; rot++){
                if(reduceVector(p, rots[rot])){
                    vec_change = true;
                }
            }
        }
    }

    if(p->norm == 0){
        for(int i = 0; i < num_rots_+1; i++) DeleteListPoint(rots[i]);
        DeleteListPoint(p);
        return 0;
    }

    p_pos = id_L;

    // rot(p) -> L
    vector<bool> vec_change_L(L.size(), false);
    genRotations(p, rots);
    for(; id_L < L.size(); id_L++){
        for(int rot = 0; rot <= num_rots_; rot++){
            if(reduceVector(L[id_L], rots[rot])){
                vec_change_L[id_L] = true;
            }
        }
    }

    for(int i = 0; i < num_rots_+1; i++) DeleteListPoint(rots[i]);

    // L -> S
    vector<ListPoint*> tmp_L;
    for(size_t i = 0; i < L.size(); i++){
        if(L[i]->norm == 0){
            DeleteListPoint(L[i]);
            collisions_++;
        }
        else{
            if(vec_change_L[i]) S.push(L[i]);
            else tmp_L.emplace_back(L[i]);
        }
    }
    L.swap(tmp_L);

    return p->norm;
}

bool IdealGSieve::Ideal_check_red3(const ListPoint *p1, const ListPoint *p2, const ListPoint *p3, ListPoint *p_new){
    // input : p1, p2, p3 (p1 <= p2 <= p3)
    // output: 3-reduced p_new,
    //         if reduced 1 else 0.

    // l' = red(v, rot(l)) > l の場合reductionは行われないのでこのチェックは満たさない
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
    for(int i = 0; i < dims; i++){
        p_new->v[i] = p1->v[i] - sign12 * p2->v[i] - sign13 * p3->v[i];
        p_new->norm += p_new->v[i] * p_new->v[i];
    }

    if(p_new->norm < p3->norm) return true;
    else return false;
}

bool IdealGSieve::Ideal_check_red3_anti_cyclic(const ListPoint *p1, const ListPoint *p2, const ListPoint *p3, ListPoint *p_new){
    // input : p1, p2, p3 (p1 <= p2 <= p3)
    // output: 3-reduced p_new,
    //         if reduced 1 else 0.

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
    
    p_new->norm = p1->norm + p2->norm + p3->norm - sign12 * 2 * dot12 - sign13 * 2 * dot13 + sign12 * sign13 * 2 * dot23;
    if(p_new->norm < p3->norm){
        for(int i = 0; i < dims; i++) p_new->v[i] = p1->v[i] - sign12 * p2->v[i] - sign13 * p3->v[i];
        return true;
    }
    else return false;
}

int IdealGSieve::adjust_order_Ideal_check_red3(const ListPoint *p1, const ListPoint *p2, const ListPoint *p3, ListPoint *p_new){
    if(p1->norm <= p2->norm){
        if(p2->norm <= p3->norm) {if(Ideal_check_red3(p1, p2, p3, p_new)) return 3;}
        else{
            if(p1->norm <= p3->norm) {if(Ideal_check_red3(p1, p3, p2, p_new)) return 2;}
            else if(Ideal_check_red3(p3, p1, p2, p_new)) return 2;
        }
    }
    else{
        if(p1->norm <= p3->norm) {if(Ideal_check_red3(p2, p1, p3, p_new)) return 3;}
        else{
            if(p2->norm <= p3->norm) {if(Ideal_check_red3(p2, p3, p1, p_new)) return 1;}
            else if(Ideal_check_red3(p3, p2, p1, p_new)) return 1;
        }
    }
    return 0;
}

int64 IdealGSieve::IdealTripleReduce_rot2_rot1(ListPoint* p){
    int64 current_norm;
    size_t p_pos;
    ListPoint* lp = NewListPoint(m_);

    // p <- L
    // cout << "--- p <- L ---" << endl;
    bool ok = true;
    int repeat = 0;
    while(ok){
        ok = false;

        //  Ideal 2-reduce
        current_norm = IdealTripleReduce_rot2_2red(p);
        if(current_norm == 0){
            DeleteListPoint(lp);
            return 0;
        }
        
        //  Ideal minkowski reduce
        vector<ListPoint*> prots(num_rots_+1);
        for(int i = 0; i < num_rots_+1; i++) prots[i] = NewListPoint(m_);
        vector<bool> vec_change_L(L.size(), false);
        bool vec_change = true;
        while(vec_change){
            vec_change = false;
            if(p->norm == 0) break;
            genRotations(p, prots);
            for(size_t i = 0; i < L.size(); i++){
                if(p->norm == 0) break;
                for(size_t j = i+1; j < L.size(); j++){
                    if(p->norm == 0) break;
                    for(int rot = 0; rot < num_rots_; rot++){
                        int res = adjust_order_Ideal_check_red3(prots[rot], L[i], L[j], lp);
                        if(res == 1){
                            if(lp->norm == 0){
                                p->norm = 0;
                                rot = num_rots_;
                            }
                            else if(lp->norm < p->norm){
                                CopyListPoint(p, lp);
                                vec_change = true;
                                ok = true;
                            }
                        }
                        else if(res == 2){
                            if(lp->norm == 0){
                                CopyListPoint(L[i], p);
                                p->norm = 0;
                                rot = num_rots_;
                            }
                            else if(lp->norm < L[i]->norm){
                                CopyListPoint(L[i], lp);
                                vec_change_L[i] = true;
                                vec_change = true;
                            }
                        }
                        else if(res == 3){
                            if(lp->norm == 0){
                                CopyListPoint(L[j], p);
                                p->norm = 0;
                                rot = num_rots_;
                            }
                            else if(lp->norm < L[j]->norm){
                                CopyListPoint(L[j], lp);
                                vec_change_L[j] = true;
                                vec_change = true;
                            }
                        }
                    }
                }
            }
        }
        for(int i = 0; i < num_rots_+1; i++) DeleteListPoint(prots[i]);

        vector<ListPoint*> L_tmp;
        for(size_t i = 0; i < L.size(); i++){
            if(L[i]->norm == 0){
                DeleteListPoint(L[i]);
                collisions_++;
            }
            else{
                if(vec_change_L[i]) S.push(L[i]);
                else L_tmp.emplace_back(L[i]);
            }
        }
        L.swap(L_tmp);

        if(p->norm == 0) break;
        repeat++;
    }

    DeleteListPoint(lp);

    // insert p into L
    if(p->norm == 0){
        DeleteListPoint(p);
        return 0;
    }
    else{
        vector<ListPoint*>::iterator itr = lower_bound(L.begin(), L.end(), p, [](const ListPoint* i, const ListPoint* j){
            return i->norm < j->norm;
        });
        if(itr != L.end()) L.insert(itr, p);
        else L.emplace_back(p);
    }

    return p->norm;
}

int64 IdealGSieve::IdealTripleReduce_rot2_rot2(ListPoint* p){
    int64 current_norm;
    size_t p_pos;
    ListPoint* lp = NewListPoint(m_);

    // p <- L
    bool ok = true;
    int repeat = 0;
    while(ok){
        ok = false;

        //  Ideal 2-reduce
        current_norm = IdealTripleReduce_rot2_2red(p);
        if(current_norm == 0){
            DeleteListPoint(lp);
            return 0;
        }
        //  Ideal minkowski reduce
        vector<ListPoint*> prots(num_rots_+1), lrots(num_rots_+1);
        for(int i = 0; i < num_rots_+1; i++){
            prots[i] = NewListPoint(m_);
            lrots[i] = NewListPoint(m_);
        }
        vector<bool> vec_change_L(L.size(), false);
        bool vec_change = true;
        while(vec_change){
            vec_change = false;
            if(p->norm == 0) break;
            genRotations(p, prots);
            for(size_t i = 0; i < L.size(); i++){
                if(p->norm == 0) break;
                genRotations(L[i], lrots);
                for(size_t j = i+1; j < L.size(); j++){
                    if(p->norm == 0) break;
                    for(int rot1 = 0; rot1 < num_rots_; rot1++){
                        for(int rot2 = 0; rot2 < num_rots_; rot2++){
                            int res = adjust_order_Ideal_check_red3(prots[rot1], lrots[rot2], L[j], lp);
                            if(res == 1){
                                if(lp->norm == 0){
                                    p->norm = 0;
                                    rot1 = num_rots_; rot2 = num_rots_;
                                }
                                else if(lp->norm < p->norm){
                                    CopyListPoint(p, lp);
                                    vec_change = true;
                                    ok = true;
                                }
                            }
                            else if(res == 2){
                                if(lp->norm == 0){
                                    CopyListPoint(L[i], p);
                                    p->norm = 0;
                                    rot1 = num_rots_; rot2 = num_rots_;
                                }
                                else if(lp->norm < L[i]->norm){
                                    CopyListPoint(L[i], lp);
                                    vec_change_L[i] = true;
                                    vec_change = true;
                                }
                            }
                            else if(res == 3){
                                if(lp->norm == 0){
                                    CopyListPoint(L[j], p);
                                    p->norm = 0;
                                    rot1 = num_rots_; rot2 = num_rots_;
                                }
                                else if(lp->norm < L[j]->norm){
                                    CopyListPoint(L[j], lp);
                                    vec_change_L[j] = true;
                                    vec_change = true;
                                }
                            }
                        }
                    }
                }
            }
        }
        for(int i = 0; i < num_rots_+1; i++){
            DeleteListPoint(prots[i]);
            DeleteListPoint(lrots[i]);
        }

        vector<ListPoint*> L_tmp;
        for(size_t i = 0; i < L.size(); i++){
            if(L[i]->norm == 0){
                DeleteListPoint(L[i]);
                collisions_++;
            }
            else{
                if(vec_change_L[i]) S.push(L[i]);
                else L_tmp.emplace_back(L[i]);
            }
        }
        L.swap(L_tmp);

        if(p->norm == 0) break;
        repeat++;
    }

    DeleteListPoint(lp);

    // insert p into L
    if(p->norm == 0){
        DeleteListPoint(p);
        return 0;
    }
    else{
        vector<ListPoint*>::iterator itr = lower_bound(L.begin(), L.end(), p, [](const ListPoint* i, const ListPoint* j){
            return i->norm < j->norm;
        });
        if(itr != L.end()) L.insert(itr, p);
        else L.emplace_back(p);
    }

    return p->norm;
}

int64 IdealGSieve::IdealTripleReduce_rot2_rot3(ListPoint* p){
    int64 current_norm;
    size_t p_pos;
    ListPoint* lp = NewListPoint(m_);

    // p <- L
    bool ok = true;
    int repeat = 0;
    while(ok){
        ok = false;

        //  Ideal 2-reduce
        current_norm = IdealTripleReduce_rot2_2red(p);
        if(current_norm == 0){
            DeleteListPoint(lp);
            return 0;
        }

        //  Ideal minkowski reduce
        vector<ListPoint*> prots(num_rots_+1), lrots(num_rots_+1), lrots2(num_rots_+1);
        for(int i = 0; i < num_rots_+1; i++){
            prots[i] = NewListPoint(m_);
            lrots[i] = NewListPoint(m_);
            lrots2[i] = NewListPoint(m_);
        }
        vector<bool> vec_change_L(L.size(), false);
        bool vec_change = true;
        while(vec_change){
            vec_change = false;
            genRotations(p, prots);
            for(size_t i = 0; i < L.size(); i++){
                genRotations(L[i], lrots);
                for(size_t j = i+1; j < L.size(); j++){
                    genRotations(L[j], lrots2);
                    for(int rot1 = 0; rot1 < num_rots_; rot1++){
                        for(int rot2 = 0; rot2 < num_rots_; rot2++){
                            for(int rot3 = 0; rot3 < num_rots_; rot3++){
                                int res = adjust_order_Ideal_check_red3(prots[rot1], lrots[rot2], lrots2[rot3], lp);
                                if(res == 1){
                                    if(lp->norm < p->norm){
                                        CopyListPoint(p, lp);
                                        vec_change = true;
                                        ok = true;
                                    }
                                }
                                else if(res == 2){
                                    if(lp->norm < L[i]->norm){
                                        CopyListPoint(L[i], lp);
                                        vec_change_L[i] = true;
                                        vec_change = true;
                                    }
                                }
                                else if(res == 3){
                                    if(lp->norm < L[j]->norm){
                                        CopyListPoint(L[j], lp);
                                        vec_change_L[j] = true;
                                        vec_change = true;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        for(int i = 0; i < num_rots_+1; i++){
            DeleteListPoint(prots[i]);
            DeleteListPoint(lrots[i]);
            DeleteListPoint(lrots2[i]);
        }

        vector<ListPoint*> L_tmp;
        for(size_t i = 0; i < L.size(); i++){
            if(L[i]->norm == 0){
                DeleteListPoint(L[i]);
                collisions_++;
            }
            else{
                if(vec_change_L[i]) S.push(L[i]);
                else L_tmp.emplace_back(L[i]);
            }
        }
        L.swap(L_tmp);

        if(p->norm == 0) break;
        repeat++;
    }

    DeleteListPoint(lp);

    // insert p into L
    if(p->norm == 0){
        DeleteListPoint(p);
        return 0;
    }
    else{
        vector<ListPoint*>::iterator itr = lower_bound(L.begin(), L.end(), p, [](const ListPoint* i, const ListPoint* j){
            return i->norm < j->norm;
        });
        if(itr != L.end()) L.insert(itr, p);
        else L.emplace_back(p);
    }

    return p->norm;
}

int64 IdealGSieve::IdealTripleReduce_rot2_rot3_parallel(){
    int64 current_norm;

    // reduce v, by l1, l2
    bool loop1 = true;
    while(loop1){
        loop1 = false;

        //  Ideal 2-reduce
        current_norm = IdealTripleReduce_rot2_2red_parallel();
        if(current_norm == 0){
            return 0;
        }

        //  Ideal minkowski reduce
        // reduce V by l1, l2 in L
        #pragma omp parallel for
        for(size_t id_V = 0; id_V < V.size(); id_V++){
            vector<ListPoint*> prots(num_rots_+1), lrots1(num_rots_+1), lrots2(num_rots_+1);
            for(int i = 0; i < num_rots_+1; i++){
                prots[i] = NewListPoint(m_);
                lrots1[i] = NewListPoint(m_);
                lrots2[i] = NewListPoint(m_);
            }
            ListPoint* lp = NewListPoint(m_);
            bool vec_change = true;
            while(vec_change){
                vec_change = false;
                genRotations(V[id_V], prots);
                for(size_t i = 0; i < L.size(); i++){
                    genRotations(L[i], lrots1);
                    for(size_t j = i+1; j < L.size(); j++){
                        genRotations(L[j], lrots2);
                        for(int rot1 = 0; rot1 <= num_rots_; rot1++){
                            // if(prots[rot1]->norm < p->norm) CopyListPoint(p, prots[rot1]);
                            for(int rot2 = 0; rot2 <= num_rots_; rot2++){
                                // if(lrots[rot2]->norm < L[i]->norm) CopyListPoint(L[i], lrots[rot2]);
                                for(int rot3 = 0; rot3 <= num_rots_; rot3++){
                                    // if(lrots2[rot3]->norm < L[j]->norm) CopyListPoint(L[j], lrots2[rot3]);
                                    int res = adjust_order_Ideal_check_red3(prots[rot1], lrots1[rot2], lrots2[rot3], lp);
                                    if(res == 1){
                                        if(lp->norm < V[id_V]->norm){
                                            CopyListPoint(V[id_V], lp);
                                            vec_change = true;
                                            loop1 = true;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            for(int i = 0; i < num_rots_+1; i++){
                DeleteListPoint(prots[i]);
                DeleteListPoint(lrots1[i]);
                DeleteListPoint(lrots2[i]);
            }
            DeleteListPoint(lp);
        }
    }

    // reduce v1 by v2, v3
    vector<bool> vec_change_V(V.size(), false);
    vector<ListPoint*> copy_V;
    for(size_t i = 0; i < V.size(); i++){
        ListPoint* tp = NewListPoint(m_);
        CopyListPoint(tp, V[i]);
        copy_V.emplace_back(tp);
    }
    #pragma omp parallel for
    for(size_t id_V = 0; id_V < V.size(); id_V++){
        vector<ListPoint*> prots1(num_rots_+1), prots2(num_rots_+1), prots3(num_rots_+1);
        for(int i = 0; i < num_rots_+1; i++){
            prots1[i] = NewListPoint(m_);
            prots2[i] = NewListPoint(m_);
            prots3[i] = NewListPoint(m_);
        }
        ListPoint* lp = NewListPoint(m_);
        genRotations(V[id_V], prots1);
        for(size_t i = 0; i < copy_V.size(); i++){
            genRotations(copy_V[i], prots2);
            for(size_t j = i + 1; j < copy_V.size(); j++){
                genRotations(copy_V[j], prots3);
                for(int rot1 = 0; rot1 <= num_rots_; rot1++){
                    for(int rot2 = 0; rot2 <= num_rots_; rot2++){
                        for(int rot3 = 0; rot3 <= num_rots_; rot3++){
                            int res = adjust_order_Ideal_check_red3(prots1[rot1], prots2[rot2], prots3[rot3], lp);
                            if(res == 1){
                                if(lp->norm < V[id_V]->norm){
                                    CopyListPoint(V[id_V], lp);
                                    vec_change_V[id_V] = true;
                                }
                            }
                        }
                    }
                }
            }
        }
        for(int i = 0; i < num_rots_+1; i++){
            DeleteListPoint(prots1[i]);
            DeleteListPoint(prots2[i]);
            DeleteListPoint(prots3[i]);
        }
        DeleteListPoint(lp);
    }
    for(ListPoint* v: copy_V) DeleteListPoint(v);

    // V -> S
    vector<ListPoint*> tmp_V;
    for(size_t i = 0; i < V.size(); i++){
        if(V[i]->norm == 0){
            DeleteListPoint(V[i]);
            collisions_++;
        }
        else tmp_V.emplace_back(V[i]);
    }
    V.swap(tmp_V);
    sort(V.begin(), V.end(), [](const ListPoint* i, const ListPoint* j){
        return i->norm < j->norm;
    });

    // reduce l1 in L by V, l2 in L
    vector<bool> vec_change_L(L.size(), false);
    vector<ListPoint*> copy_L;
    for(size_t i = 0; i < L.size(); i++){
        ListPoint* tp = NewListPoint(m_);
        CopyListPoint(tp, L[i]);
        copy_L.emplace_back(tp);
    }
    #pragma omp parallel for
    for(size_t id_L = 0; id_L < L.size(); id_L++){
        vector<ListPoint*> prots(num_rots_+1), lrots1(num_rots_+1), lrots2(num_rots_+1);
        for(int i = 0; i < num_rots_+1; i++){
            prots[i] = NewListPoint(m_);
            lrots1[i] = NewListPoint(m_);
            lrots2[i] = NewListPoint(m_);
        }
        ListPoint* lp = NewListPoint(m_);
        bool vec_change = true;
        while(vec_change){
            vec_change = false;
            genRotations(L[id_L], lrots1);
            for(size_t id_V = 0; id_V < V.size(); id_V++){
                genRotations(V[id_V], prots);
                for(size_t j = id_L + 1; j < copy_L.size(); j++){
                    genRotations(copy_L[j], lrots2);
                    for(int rot1 = 0; rot1 <= num_rots_; rot1++){
                        for(int rot2 = 0; rot2 <= num_rots_; rot2++){
                            for(int rot3 = 0; rot3 <= num_rots_; rot3++){
                                int res = adjust_order_Ideal_check_red3(lrots1[rot1], prots[rot2], lrots2[rot3], lp);
                                if(res == 1){
                                    if(lp->norm < L[id_L]->norm){
                                        CopyListPoint(L[id_L], lp);
                                        vec_change_L[id_L] = true;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        for(int i = 0; i < num_rots_+1; i++){
            DeleteListPoint(prots[i]);
            DeleteListPoint(lrots1[i]);
            DeleteListPoint(lrots2[i]);
        }
        DeleteListPoint(lp);
    }
    for(ListPoint* v: copy_L) DeleteListPoint(v);

    // L -> S
    vector<ListPoint*> L_tmp;
    for(size_t i = 0; i < L.size(); i++){
        if(L[i]->norm == 0){
            DeleteListPoint(L[i]);
            collisions_++;
        }
        else{
            if(vec_change_L[i]) S.push(L[i]);
            else L_tmp.emplace_back(L[i]);
        }
    }
    L.swap(L_tmp);

    // L += V
    if(V.empty()) return min_norm_;
    vector<ListPoint*> tmp_L;
    size_t i = 0, j = 0;
    current_norm = V[0]->norm;
    while(!(i == L.size() && j == V.size())){
        while(i < L.size() && (j == V.size() || L[i]->norm <= V[j]->norm)){
            tmp_L.emplace_back(L[i]);
            i++;
        }
        while(j < V.size() && (i == L.size() || L[i]->norm > V[j]->norm)){
            tmp_L.emplace_back(V[j]);
            j++;
        }
    }
    V.clear();
    L.swap(tmp_L);

    return current_norm;
}

int64 IdealGSieve::IdealTripleReduce_anti_cyclic_rot1(ListPoint* p){
    int64 current_norm;
    size_t p_pos;
    ListPoint* lp = NewListPoint(m_);
    vector<ListPoint*> rots(num_rots_ + 1);
    for(int i = 0; i < num_rots_+1; i++) rots[i] = NewListPoint(m_);

    // p <- L
    // cout << "--- p <- L ---" << endl;
    bool ok = true;
    int repeat = 0;
    while(ok){
        // cout << "itr: " << repeat << ", L: " << L.size() << ", p: " << p->norm << endl;
        ok = false;

        //  Ideal 2-reduce
        current_norm = IdealTripleReduce_anti_cyclic_2red(p, p_pos);
        if(current_norm == 0){
            for(size_t i = 0; i < rots.size(); i++) DeleteListPoint(rots[i]);
            DeleteListPoint(lp);
            return 0;
        }

        //  Ideal minkowski reduce
        genRotations(p, rots);
        for(size_t i = 0; i < p_pos; i++){
            for(size_t j = i+1; j < p_pos; j++){
                for(int rot = 0; rot < num_rots_; rot++){
                    if(Ideal_check_red3_anti_cyclic(L[i], L[j], rots[rot], lp)){
                        CopyListPoint(p, lp);
                        ok = true;
                        break;
                    }
                }
                if(ok) break;
            }
            if(ok) break;
        }

        repeat++;
    }

    // p -> L
    vector<bool> vec_change_L(L.size(), false);
    genRotations_anti_cyclic(p, rots);
    for(size_t i = p_pos + 1; i < L.size(); i++){
        for(size_t j = 0; j < i; j++){
            if(vec_change_L[j]) continue;
            if(j == p_pos) continue;
            for(int rot = 0; rot < num_rots_; rot++){
                if(j < p_pos){
                    if(Ideal_check_red3_anti_cyclic(L[j], rots[rot], L[i], lp)){
                        CopyListPoint(L[i], lp);
                        vec_change_L[i] = true;
                    }
                }
                else{
                    if(Ideal_check_red3_anti_cyclic(rots[rot], L[j], L[i], lp)){
                        CopyListPoint(L[i], lp);
                        vec_change_L[i] = true;
                    }
                }
                if(vec_change_L[i]) break;
            }
            if(vec_change_L[i]) break;
        }
    }

    for(size_t i = 0; i < rots.size(); i++) DeleteListPoint(rots[i]);
    DeleteListPoint(lp);

    vector<ListPoint*> L_tmp;
    for(size_t i = 0; i < L.size(); i++){
        if(L[i]->norm == 0){
            DeleteListPoint(L[i]);
            collisions_++;
        }
        else{
            if(vec_change_L[i]) S.push(L[i]);
            else L_tmp.emplace_back(L[i]);
        }
    }
    L.swap(L_tmp);

    // insert p into L
    vector<ListPoint*>::iterator itr = lower_bound(L.begin(), L.end(), p, [](const ListPoint* i, const ListPoint* j){
        return i->norm < j->norm;
    });
    if(itr != L.end()) L.insert(itr, p);
    else L.emplace_back(p);

    return p->norm;
}

int64 IdealGSieve::IdealTripleReduce_anti_cyclic_rot2(ListPoint* p){
    int64 current_norm;
    size_t p_pos;
    ListPoint* lp = NewListPoint(m_);
    vector<ListPoint*> prots(num_rots_+1), lrots(num_rots_+1);
    for(size_t i = 0; i < num_rots_+1; i++){
        prots[i] = NewListPoint(m_);
        lrots[i] = NewListPoint(m_);
    }

    // p <- L
    bool ok = true;
    int repeat = 0;
    while(ok){
        ok = false;

        //  Ideal 2-reduce
        current_norm = IdealTripleReduce_anti_cyclic_2red(p, p_pos);
        if(current_norm == 0){
            for(int i = 0; i < num_rots_+1; i++){
                DeleteListPoint(prots[i]);
                DeleteListPoint(lrots[i]);
            }
            DeleteListPoint(lp);
            return 0;
        }

        //  Ideal minkowski reduce
        genRotations_anti_cyclic(p, prots);
        for(size_t i = 0; i < p_pos; i++){
            genRotations_anti_cyclic(L[i], lrots);
            for(size_t j = i+1; j < p_pos; j++){
                for(int rot1 = 0; rot1 < num_rots_; rot1++){
                    for(int rot2 = 0; rot2 < num_rots_; rot2++){
                        if(Ideal_check_red3_anti_cyclic(lrots[rot2], L[j], prots[rot1], lp)){
                            CopyListPoint(p, lp);
                            ok = true;
                            rot1 = num_rots_; rot2 = num_rots_;
                        }
                    }
                }
                if(ok) break;
            }
            if(ok) break;
        }

        repeat++;
    }

    // p -> L
    vector<bool> vec_change_L(L.size(), false);
    genRotations(p, prots);
    for(size_t i = p_pos + 1; i < L.size(); i++){
        genRotations(L[i], lrots);
        for(size_t j = 0; j < i; j++){
            if(vec_change_L[j]) continue;
            if(j == p_pos) continue;
            for(int rot1 = 0; rot1 < num_rots_; rot1++){
                for(int rot2 = 0; rot2 < num_rots_; rot2++){
                    if(j < p_pos){
                        if(Ideal_check_red3_anti_cyclic(L[j], prots[rot1], lrots[rot2], lp)){
                            CopyListPoint(L[i], lp);
                            vec_change_L[i] = true;
                            rot1 = num_rots_; rot2 = num_rots_;
                        }
                    }
                    else{
                        if(Ideal_check_red3_anti_cyclic(prots[rot1], L[j], lrots[rot2], lp)){
                            CopyListPoint(L[i], lp);
                            vec_change_L[i] = true;
                            rot1 = num_rots_; rot2 = num_rots_;
                        }
                    }
                }
            }
            if(vec_change_L[i]) break;
        }
    }

    for(int i = 0; i < num_rots_+1; i++){
        DeleteListPoint(prots[i]);
        DeleteListPoint(lrots[i]);
    }
    DeleteListPoint(lp);

    vector<ListPoint*> L_tmp;
    for(size_t i = 0; i < L.size(); i++){
        if(L[i]->norm == 0){
            DeleteListPoint(L[i]);
            collisions_++;
        }
        else{
            if(vec_change_L[i]) S.push(L[i]);
            else L_tmp.emplace_back(L[i]);
        }
    }
    L.swap(L_tmp);

    // insert p into L
    vector<ListPoint*>::iterator itr = lower_bound(L.begin(), L.end(), p, [](const ListPoint* i, const ListPoint* j){
        return i->norm < j->norm;
    });
    if(itr != L.end()) L.insert(itr, p);
    else L.emplace_back(p);

    return p->norm;
}

int64 IdealGSieve::IdealTripleReduce_anti_cyclic_rot3(ListPoint* p){
    int64 current_norm;
    size_t p_pos;
    ListPoint* lp = NewListPoint(m_);
    vector<ListPoint*> prots(num_rots_+1), lrots1(num_rots_+1), lrots2(num_rots_+1);
    for(size_t i = 0; i < num_rots_+1; i++){
        prots[i] = NewListPoint(m_);
        lrots1[i] = NewListPoint(m_);
        lrots2[i] = NewListPoint(m_);
    }

    // p <- L
    bool ok = true;
    int repeat = 0;
    while(ok){
        ok = false;

        //  Ideal 2-reduce
        current_norm = IdealTripleReduce_anti_cyclic_2red(p, p_pos);
        if(current_norm == 0){
            for(int i = 0; i < num_rots_+1; i++){
                DeleteListPoint(prots[i]);
                DeleteListPoint(lrots1[i]);
                DeleteListPoint(lrots2[i]);
            }
            DeleteListPoint(lp);
            return 0;
        }

        //  Ideal minkowski reduce
        genRotations_anti_cyclic(p, prots);
        for(size_t i = 0; i < p_pos; i++){
            genRotations_anti_cyclic(L[i], lrots1);
            for(size_t j = i+1; j < p_pos; j++){
                genRotations_anti_cyclic(L[j], lrots2);
                for(int rot1 = 0; rot1 < num_rots_; rot1++){
                    for(int rot2 = 0; rot2 < num_rots_; rot2++){
                        for(int rot3 = 0; rot3 < num_rots_; rot3++){
                            if(Ideal_check_red3_anti_cyclic(lrots1[rot2], lrots2[rot3], prots[rot1], lp)){
                                CopyListPoint(p, lp);
                            }
                        }
                    }
                }
            }
        }

        repeat++;
    }

    // p -> L
    vector<bool> vec_change_L(L.size(), false);
    genRotations(p, prots);
    for(size_t i = p_pos + 1; i < L.size(); i++){
        genRotations(L[i], lrots1);
        for(size_t j = 0; j < i; j++){
            if(vec_change_L[j]) continue;
            if(j == p_pos) continue;
            genRotations(L[j], lrots2);
            for(int rot1 = 0; rot1 < num_rots_; rot1++){
                for(int rot2 = 0; rot2 < num_rots_; rot2++){
                    for(int rot3 = 0; rot3 < num_rots_; rot3++){
                        if(j < p_pos){
                            if(Ideal_check_red3_anti_cyclic(lrots2[rot3], prots[rot1], lrots1[rot2], lp)){
                                CopyListPoint(L[i], lp);
                                vec_change_L[i] = true;
                            }
                        }
                        else{
                            if(Ideal_check_red3_anti_cyclic(prots[rot1], lrots2[rot3], lrots1[rot2], lp)){
                                CopyListPoint(L[i], lp);
                                vec_change_L[i] = true;
                            }
                        }
                    }
                }
            }
            if(vec_change_L[i]) break;
        }
    }

    for(int i = 0; i < num_rots_+1; i++){
        DeleteListPoint(prots[i]);
        DeleteListPoint(lrots1[i]);
        DeleteListPoint(lrots2[i]);
    }
    DeleteListPoint(lp);

    vector<ListPoint*> L_tmp;
    for(size_t i = 0; i < L.size(); i++){
        if(L[i]->norm == 0){
            DeleteListPoint(L[i]);
            collisions_++;
        }
        else{
            if(vec_change_L[i]) S.push(L[i]);
            else L_tmp.emplace_back(L[i]);
        }
    }
    L.swap(L_tmp);

    // insert p into L
    vector<ListPoint*>::iterator itr = lower_bound(L.begin(), L.end(), p, [](const ListPoint* i, const ListPoint* j){
        return i->norm < j->norm;
    });
    if(itr != L.end()) L.insert(itr, p);
    else L.emplace_back(p);

    return p->norm;
}

void IdealGSieve::IdealGaussSieve(){
    ListPoint* new_v;
    int64 current_norm;

    while(collisions_ < max_list_size_/10+200 && min_norm_ > goal_norm_){
        // cout << "iteration: " << iterations_ << endl;
        iterations_++;
        max_list_size_ = std::max(max_list_size_, (long)L.size());
        
        if(!S.empty()){
            new_v = S.top();
            S.pop();
        }else{
            new_v = sampler_->Sample();
            sample_vectors_++;
        }

        // current_norm = IdealGaussReduce_rot1(new_v);
        current_norm = IdealGaussReduce_rot2(new_v);
        // current_norm = IdealGaussReduce_rot2_modified(new_v);
        // current_norm = IdealGaussReduce_anti_cyclic(new_v);
        // current_norm = IdealTripleReduce_rot2_rot1(new_v);
        // current_norm = IdealTripleReduce_rot2_rot2(new_v);
        // current_norm = IdealTripleReduce_rot2_rot3(new_v);
        // current_norm = IdealTripleReduce_rot1_v2(new_v);
        // current_norm = IdealTripleReduce_rot1_v3(new_v);
        // current_norm = IdealTripleReduce_anti_cyclic_rot1(new_v);
        // current_norm = IdealTripleReduce_anti_cyclic_rot2(new_v);
        // current_norm = IdealTripleReduce_anti_cyclic_rot3(new_v);

        if(current_norm == 0){
            collisions_++;
        }else{
            if(current_norm < min_norm_){
                min_norm_ = current_norm;
                cout << min_norm_ << endl;
            }
            if(!L.empty() && L[0]->norm < min_norm_){
                min_norm_ = L[0]->norm;
                cout << min_norm_ << endl;
            }
        }

        // cout << L.size() << "." << S.size() << ":" << flush;
    }

    if(!(collisions_ < max_list_size_/10+200)) cout << "list size" << endl;
    if(!(min_norm_ > goal_norm_)) cout << "min norm" << endl;
}

void IdealGSieve::IdealGaussSieve_parallel(){
    ListPoint* new_v;
    int64 current_norm;
    vector<ListPoint*> lps(num_rots_+1);
    for(int i = 0; i < num_rots_+1; i++) lps[i] = NewListPoint(m_);

    while(collisions_ < max_list_size_/10+200 && min_norm_ > goal_norm_){
        // cout << "iteration: " << iterations_ << endl;
        iterations_++;
        max_list_size_ = std::max(max_list_size_, (long)L.size());

        while(V.size() < simu_samp_){
            if(!S.empty()){
                new_v = S.top();
                S.pop();
            }
            else{
                new_v = sampler_->Sample();
                sample_vectors_++;
            }
            V.emplace_back(new_v);
        }
        
        // current_norm = IdealGaussReduce_rot2_parallel();
        current_norm = IdealTripleReduce_rot2_rot3_parallel();
        // current_norm = IdealGaussReduce_anti_cyclic_parallel();

        if((current_norm > 0) && (current_norm < min_norm_)){
            min_norm_ = current_norm;
            cout << min_norm_ << endl;
        }
        if(!L.empty() && L[0]->norm < min_norm_){
            min_norm_ = L[0]->norm;
            cout << min_norm_ << endl;
        }

        #if DEBUG
        if(L.size() >= 2){
            for(size_t i = 0; i < L.size()-1; i++){
                if(L[i+1]->norm < L[i]->norm) cout << "error: L isn't asscending." << endl;
            }
        }
        if(V.size() >= 2){
            for(size_t i = 0; i < V.size()-1; i++){
                if(V[i+1]->norm < V[i]->norm) cout << "error: V isn't asscending." << endl;
            }
        }
        #endif
        // cout << L.size() << "." << S.size() << ":" << flush;
    }

    for(int i = 0; i < num_rots_+1; i++) DeleteListPoint(lps[i]);

    if(!(collisions_ < max_list_size_/10+200)) cout << "list size" << endl;
    if(!(min_norm_ > goal_norm_)) cout << "min norm" << endl;
}