#include "Idealgsieve.h"

void IdealGSieve::printL(){
    cout << "L:" << endl;
    for(auto v: L) cout << v->norm << " ";
    cout << endl;
}

void IdealGSieve::printV(){
    cout << "V:" << endl;
    for(auto v: V) cout << v->norm << " ";
    cout << endl;
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

void IdealGSieve::Init(const mat_ZZ &B, KleinSampler* sampler, long index, const vec_ZZ &modf){
    n_ = B.NumRows();
    m_ = B.NumCols();
    sampler_ = sampler;
    index_ = index;
    iterations_ = 0;
    collisions_ = 0;
    sample_vectors_ = 0;
    CleanUp();

    sampler_->Init(B);
    min_norm_ = to_long(B[0] * B[0]);

    modf_.SetLength(n_);
    for(int i = 0; i < n_; i++) modf_[i] = to_long(modf[i]);

    // mat_ZZ rot;
    // rot.SetDims(n_, n_);
    // clear(rot);
    // for(int i = 0; i < n_-1; i++){
    //     rot[i][i+1] = 1;
    // }
    // for(int i = 0; i < n_; i++) rot[n_-1][i] = -modf[n_-i];
    // MatInt64FromMatZZ(rot, rot_);
    
    ListPoint* p;
    int64 current_norm;
    for(int i = 0; i < n_; i++){
        p = NewListPoint(m_);
        VecZZToListPoint(B[i], p);
        current_norm = IdealGaussReduce2(p);
        if(current_norm < min_norm_){
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

int64 IdealGSieve::IdealGaussReduce(ListPoint* p){
    // p <- L
    // Lが昇順に並んでいるので、pより大きいLの要素はスルー
    // スルーした要素の中にrotationでpより短くなるものが存在する可能性はある

    bool vec_change = true;
    size_t id_L;
    while(vec_change){
        vec_change = false;
        for(id_L = 0; id_L < L.size(); id_L++){
            if(p->norm < L[id_L]->norm){
                break;
            }
            if(IdealreduceVector(p, L[id_L], modf_, index_)){
                vec_change = true;
            }
        }
    }

    if(p->norm == 0){
        DeleteListPoint(p);
        return 0;
    }

    // p -> L
    // pより大きいLの要素に対してreduceする
    // p以下のLの要素でもrotationでp以上になる可能性は十分ある
    // ここではすべてのLの要素に対してrotationしてreduceを行う  <- してない
    L.insert(L.begin()+id_L, p);
    id_L++;

    // p -> L
    for(; id_L < L.size(); id_L++){
        if(IdealreduceVector(L[id_L], p, modf_, index_)){
            S.push(L[id_L]);
            L.erase(L.begin()+id_L);
            id_L--;
        }
    }

    return p->norm;
}

// こっちのほうが時間はかかるがリストサイズは短い
int64 IdealGSieve::IdealGaussReduce2(ListPoint* p){
    // p <- L
    // Lの全要素に対してreduce
    // 次元ごとにrotationによるpのノルムの変化の上界・下界が分かれば最適化できる

    bool vec_change = true;
    while(vec_change){
        vec_change = false;
        for(size_t id_L = 0; id_L < L.size(); id_L++){
            if(IdealreduceVector(p, L[id_L], modf_, index_)){
                vec_change = true;
            }
        }
    }

    if(p->norm == 0){
        DeleteListPoint(p);
        return 0;
    }

    // insert p into L
    size_t p_pos;
    auto itr = lower_bound(L.begin(), L.end(), p, [](const ListPoint *i, const ListPoint *j){
        return i->norm < j->norm;
    });
    p_pos = itr - L.begin();
    L.insert(itr, p);

    // p -> L
    // ここではすべてのLの要素に対してrotationしてreduceを行う

    for(size_t id_L = 0; id_L < L.size(); id_L++){
        if(id_L == p_pos) continue;
        if(IdealreduceVector(L[id_L], p, modf_, index_)){
            S.push(L[id_L]);
            L.erase(L.begin()+id_L);
            id_L--;
            if(id_L < p_pos) p_pos--;
        }
    }

    return p->norm;
}

int64 IdealGSieve::GaussReduce_Parallel(){
    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    int64 current_norm = V[0]->norm;
    ////  L -> V  ////
    vector<int> vec_change_V(V.size(), -1);
    vector<thread> threads_LV;

    for(int i = 0; i < concurrency_; i++){
        threads_LV.emplace_back(
            [&vec_change_V, this](int i){
                for(size_t id_V = i; id_V < V.size(); id_V += concurrency_){
                    bool vec_change_parallel = true;
                    bool vec_change_flag = false;
                    size_t id_L;
                    while(vec_change_parallel){
                        vec_change_parallel = false;
                        for(id_L = 0; id_L < L.size(); id_L++){
                            if(V[id_V]->norm < L[id_L]->norm){
                                break;
                            }
                            if(reduceVector(V[id_V], L[id_L])){
                                vec_change_parallel = true;
                                vec_change_flag = true;
                            }
                        }
                    }

                    // vec_change_Vはvより最初に大きいlを指す
                    if(!vec_change_flag){
                        vec_change_V[id_V] = (int)id_L;
                    }
                }
            }, i
        );
    }
    for(thread &th : threads_LV){
        th.join();
    }

    // そもそもこの処理内でSに移す必要ある？
    for(size_t i = 0; i < V.size(); i++){
        if(vec_change_V[i] != -1) continue;
        if(V[i]->norm == 0){
            DeleteListPoint(V[i]);
            collisions_++;
        }
        else S.push(V[i]);
        // vec_change_V の値で区別してもいいかも <- この後のVのループで影響出る
        V.erase(V.begin()+i);
        vec_change_V.erase(vec_change_V.begin()+i);
        i--;
    }
    end = chrono::system_clock::now();
    timeL2V += (double)chrono::duration_cast<chrono::microseconds>(end-start).count()/1000;


    start = chrono::system_clock::now();
    ////  V <- V  ////
    bool vec_change_VV = true;
    while(vec_change_VV){
        vec_change_VV = false;
        for(size_t id_V1 = 0; id_V1 < V.size(); id_V1++){
            for(size_t id_V2 = 0; id_V2 < V.size(); id_V2++){
                if(id_V1 == id_V2) continue;
                if(V[id_V1]->norm < V[id_V2]->norm){
                    continue;
                }
                if(reduceVector(V[id_V1], V[id_V2])){
                    vec_change_VV = true;
                    vec_change_V[id_V1] = -1;
                }
            }
        }
    }
    
    for(size_t i = 0; i < V.size(); i++){
        if(vec_change_V[i] != -1) continue;
        if(V[i]->norm == 0){
            DeleteListPoint(V[i]);
            collisions_++;
        }
        else{
            if(current_norm > V[i]->norm) current_norm = V[i]->norm;
            S.push(V[i]);
        }
        V.erase(V.begin()+i);
        vec_change_V.erase(vec_change_V.begin()+i);
        i--;
    }
    end = chrono::system_clock::now();
    timeV2V += (double)chrono::duration_cast<chrono::microseconds>(end-start).count()/1000;


    start = chrono::system_clock::now();
    ////  V -> L  ////
    vector<bool> vec_change_L(L.size(), false);
    vector<thread> threads_VL;

    for(int i = 0; i < concurrency_; i++){
        threads_VL.emplace_back(
            [&vec_change_L, &vec_change_V, this](int i){
                for(size_t id_L = i; id_L < L.size(); id_L += concurrency_){
                    bool vec_change_parallel = true;
                    while(vec_change_parallel){
                        vec_change_parallel = false;
                        for(size_t id_V = 0; id_V < V.size(); id_V++){
                            if(vec_change_V[id_V] > (int)id_L) continue;
                            if(reduceVector(L[id_L], V[id_V])){
                                vec_change_parallel = true;
                                vec_change_L[id_L] = true;
                            }
                        }
                    }
                }
            }, i
        );
    }
    for(thread &th : threads_VL){
        th.join();
    }

    // l' -> S
    for(size_t i = 0; i < L.size(); i++){
        if(!vec_change_L[i]) continue;
        // 多分起きない
        if(L[i]->norm == 0){
            DeleteListPoint(L[i]);
            collisions_++;
        }
        else S.push(L[i]);
        L.erase(L.begin()+i);
        vec_change_L.erase(vec_change_L.begin()+i);
        i--;
    }
    end = chrono::system_clock::now();
    timeV2L += (double)chrono::duration_cast<chrono::microseconds>(end-start).count()/1000;

    return current_norm;
}

void IdealGSieve::IdealGaussSieve(){
    // chrono::system_clock::time_point start, end;
    ListPoint* new_v;
    int64 current_norm;

    // double time1 = 0, time2 = 0, time3 = 0;
    while(collisions_ < max_list_size_/10+200 && min_norm_ > goal_norm_){
        // cout << "iteration: " << iterations_ << endl;
        iterations_++;

        // time1
        // start = chrono::system_clock::now();
        max_list_size_ = std::max(max_list_size_, (long)L.size());

        if(!S.empty()){
            new_v = S.top();
            S.pop();
        }else{
            new_v = sampler_->Sample();
            sample_vectors_++;
        }
        // end = chrono::system_clock::now();
        // time1 += (double)chrono::duration_cast<chrono::microseconds>(end-start).count()/1000;


        // time2
        // start = chrono::system_clock::now();
        current_norm = IdealGaussReduce2(new_v);
        // end = chrono::system_clock::now();
        // time2 += (double)chrono::duration_cast<chrono::microseconds>(end-start).count()/1000;
        

        // time3
        // start = chrono::system_clock::now();
        if(current_norm == 0){
            collisions_++;
        }else{
            if(current_norm < min_norm_){
                min_norm_ = current_norm;
                cout << min_norm_ << endl;
            }
        }
        // end = chrono::system_clock::now();
        // time3 += (double)chrono::duration_cast<chrono::microseconds>(end-start).count()/1000;

        // cout << L.size() << "." << S.size() << ":" << flush;
    }

    // chk_time_.push_back(time1);
    // chk_time_.push_back(time2);
    // chk_time_.push_back(time3);
}

void IdealGSieve::GaussSieve_Parallel(){
    chrono::system_clock::time_point start, end;
    ListPoint* new_v;
    vector<ListPoint*>::iterator itr;
    int64 current_norm;

    double time1 = 0, time2 = 0, time3 = 0;
    while(collisions_ < max_list_size_/10+200 && min_norm_ > goal_norm_){
        // time1
        start = chrono::system_clock::now();
        max_list_size_ = std::max(max_list_size_, (long)(L.size()));
        
        // Vに昇順でベクトルを保存
        if(!S.empty()){
            while(!S.empty() && (int)V.size() < simu_samp_){
                new_v = S.top(); S.pop();
                if(V.empty()) V.emplace_back(new_v);
                else{
                    itr = lower_bound(V.begin(), V.end(), new_v, [](ListPoint* i, ListPoint* j){
                        return i->norm < j->norm;
                    });
                    if(itr != V.end()) V.insert(itr, new_v);
                    else V.emplace_back(new_v);
                }
            }
        }
        if(S.empty()){
            while((int)V.size() < simu_samp_){
                new_v = sampler_->Sample();
                if(V.empty()) V.emplace_back(new_v);
                else{
                    itr = lower_bound(V.begin(), V.end(), new_v, [](ListPoint* i, ListPoint* j){
                        return i->norm < j->norm;
                    });
                    if(itr != V.end()) V.insert(itr, new_v);
                    else V.emplace_back(new_v);
                }

                sample_vectors_++;
            }
        }
        end = chrono::system_clock::now();
        time1 += (double)chrono::duration_cast<chrono::microseconds>(end-start).count()/1000;

        // time2
        start = chrono::system_clock::now();
        current_norm = GaussReduce_Parallel();
        end = chrono::system_clock::now();
        time2 += (double)chrono::duration_cast<chrono::microseconds>(end-start).count()/1000;

        // time3
        start = chrono::system_clock::now();
        if(current_norm != -1 && current_norm < min_norm_){
            min_norm_ = current_norm;
        }

        for(int i = 0; i < (int)V.size(); i++){
            // lower_boundの区間の上限はvec_change_V[i]
            itr = lower_bound(L.begin(), L.end(), V[i], [](ListPoint* i, ListPoint* j){
                return i->norm < j->norm;
            });
            if(itr != L.end()) L.insert(itr, V[i]);
            else L.emplace_back(V[i]);
        }
        V.clear();
        end = chrono::system_clock::now();
        time3 += (double)chrono::duration_cast<chrono::microseconds>(end-start).count()/1000;

        iterations_++;
    }

    chk_time_.push_back(time1);
    chk_time_.push_back(time2);
    chk_time_.push_back(time3);
}