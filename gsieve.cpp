#include "gsieve.h"

void GSieve::printL(){
    cout << "L:" << endl;
    for(auto v: L) cout << v->norm << " ";
    cout << endl;
}

void GSieve::printV(){
    cout << "V:" << endl;
    for(auto v: V) cout << v->norm << " ";
    cout << endl;
}

void GSieve::CleanUp(){
    for(size_t i = 0; i < V.size(); i++) DeleteListPoint(V[i]);
    for(size_t i = 0; i < V_.size(); i++) DeleteListPoint(V_[i]);
    for(size_t i = 0; i < V__.size(); i++) DeleteListPoint(V__[i]);
    for(size_t i = 0; i < L.size(); i++) DeleteListPoint(L[i]);
    while(!S.empty()){
        DeleteListPoint(S.front());
        S.pop();
    }
}

void GSieve::Init(const mat_ZZ &B, KleinSampler* sampler){
    n_ = B.NumRows();
    m_ = B.NumCols();
    sampler_ = sampler;
    iterations_ = 0;
    collisions_ = 0;
    sample_vectors_ = 0;
    CleanUp();

    sampler_->Init(B);
    min_norm_ = to_long(B[0] * B[0]);

    ListPoint* p;
    int64 current_norm;
    for(int i = 0; i < n_; i++){
        p = NewListPoint(m_);
        VecZZToListPoint(B[i], p);
        current_norm = GaussReduce(p);
        if(current_norm < min_norm_){
            min_norm_ = current_norm;
        }
    }
    max_list_size_ = L.size();
    concurrency_ = 1;
    simu_samp_ = 1;
    goal_norm_ = 0;
}

// void GSieve::VectorReduce_Parallel(){
//     vector<LatticeVector*> Tmp;
//     for(int i = 0; i < V.size(); i++){
//         LatticeVector *lv = newLatticeVector(m_);
//         lv->vec = V[i]->vec;
//         lv->norm2 = V[i]->norm2;
//         Tmp.emplace_back(lv);
//     }
//     cout << 111 << endl;

//     vector<bool> vec_change(V.size(), false);
//     int num_parallel = V.size() / concurrency_;
//     int num_remain = V.size() % concurrency_;

//     // 各スレッドにnum_parallel個のベクトルを分配
//     vector<thread> threads;
//     for(int i = 0; i < concurrency_; i++){
//         threads.emplace_back(
//             [&vec_change, i, this, &Tmp](int num){
//                 for(int j = 0; j < num; j++){
//                     for(int k = 0; k < V.size(); k++){
//                         if(i*num+j == k) continue;
//                         if(Tmp[i*num+j]->norm2 < V[k]->norm2){
//                             continue;
//                         }
//                         if(reduceVector(Tmp[i*num+j], V[k])){
//                             vec_change[i*num+j] = true;
//                         }
//                     }
//                 }
//             }, num_parallel
//         );
//     }
//     for(thread &th : threads){
//         th.join();
//     }

//     // 余ったベクトルは本スレッドで
//     int tmp_sz = Tmp.size();
//     for(int tmp_id = concurrency_ * num_parallel; tmp_id < tmp_sz; tmp_id++){
//         for(int v_id = 0; v_id < V.size(); v_id++){
//             if(tmp_id == v_id) continue;
//             if(Tmp[tmp_id]->norm2 < V[v_id]->norm2){
//                 continue;
//             }
//             if(reduceVector(Tmp[tmp_id], V[v_id])){
//                 vec_change[tmp_id] = true;
//             }
//         }
//     }

//     // 更新したベクトルがあればSに移動
//     for(int i = 0; i < vec_change.size(); i++){
//         if(vec_change[i] == true){
//             S.push(Tmp[i]);
//             // delete V[i];
//             V.erase(V.begin()+i);
//             vec_change.erase(vec_change.begin()+i);
//             i--;
//         }
//         // else{
//         //     delete Tmp[i];
//         // }
//     }
// }

int64 GSieve::GaussReduce(ListPoint* p){
    // p <- L
    bool vec_change = true;
    size_t id_L;
    vector<ListPoint*> L_tmp;
    while(vec_change){
        vec_change = false;
        for(id_L = 0; id_L < L.size(); id_L++){
            if(p->norm < L[id_L]->norm){
                break;
            }
            if(reduceVector(p, L[id_L])){
                vec_change = true;
            }
        }
    }

    if(p->norm == 0){
        DeleteListPoint(p);
        return 0;
    }

    // vectorでinsert, deleteは時間かかるので改善が必要、listなら速い -> 元のGSよりちょっと早いんでこのまま
    L.insert(L.begin()+id_L, p);
    id_L++;

    // p -> L
    for(; id_L < L.size(); id_L++){
        if(reduceVector(L[id_L], p)){
            S.push(L[id_L]);
            L.erase(L.begin()+id_L);
            id_L--;
        }
    }

    return p->norm;
}

int64 GSieve::GaussReduce_Parallel(){
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
        else S.push(V[i]);
        V.erase(V.begin()+i);
        vec_change_V.erase(vec_change_V.begin()+i);
        i--;
    }


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

    if(!V.empty()) return V[0]->norm;
    else return -1;
}

void GSieve::GaussSieve(){
    // chrono::system_clock::time_point start, end;
    ListPoint* new_v;
    int64 current_norm;

    // double time1 = 0, time2 = 0, time3 = 0;
    while(collisions_ < max_list_size_/10+200 && min_norm_ > goal_norm_){
        iterations_++;

        // time1
        // start = chrono::system_clock::now();
        max_list_size_ = std::max(max_list_size_, (long)L.size());

        if(!S.empty()){
            new_v = S.front();
            S.pop();
        }else{
            new_v = sampler_->Sample();
            sample_vectors_++;
        }
        // end = chrono::system_clock::now();
        // time1 += (double)chrono::duration_cast<chrono::microseconds>(end-start).count()/1000;


        // time2
        // start = chrono::system_clock::now();
        current_norm = GaussReduce(new_v);
        // end = chrono::system_clock::now();
        // time2 += (double)chrono::duration_cast<chrono::microseconds>(end-start).count()/1000;
        

        // time3
        // start = chrono::system_clock::now();
        if(current_norm == 0){
            collisions_++;
        }else{
            if(current_norm < min_norm_){
                min_norm_ = current_norm;
            }
        }
        // end = chrono::system_clock::now();
        // time3 += (double)chrono::duration_cast<chrono::microseconds>(end-start).count()/1000;
    }

    // chk_time_.push_back(time1);
    // chk_time_.push_back(time2);
    // chk_time_.push_back(time3);
}

void GSieve::GaussSieve_Parallel(){
    // chrono::system_clock::time_point start, end;
    ListPoint* new_v;
    int64 current_norm;

    // double time1 = 0, time2 = 0, time3 = 0;
    while(collisions_ < max_list_size_/10+200 && min_norm_ > goal_norm_){
        // time1
        // start = chrono::system_clock::now();
        max_list_size_ = std::max(max_list_size_, (long)(L.size()));
        
        // Vに昇順でベクトルを保存
        if(!S.empty()){
            while(!S.empty() && (int)V.size() < simu_samp_){
                new_v = S.front(); S.pop();
                if(V.empty()) V.emplace_back(new_v);
                else{
                    size_t id = 0;
                    while(id < V.size() && V[id]->norm < new_v->norm) id++;
                    if(id < V.size()) V.insert(V.begin() + id, new_v);
                    else V.emplace_back(new_v);
                }
            }
        }
        if(S.empty()){
            while((int)V.size() < simu_samp_){
                new_v = sampler_->Sample();
                if(V.empty()) V.emplace_back(new_v);
                else{
                    size_t id = 0;
                    while(id < V.size() && V[id]->norm < new_v->norm) id++;
                    if(id < V.size()) V.insert(V.begin() + id, new_v);
                    else V.emplace_back(new_v);
                }

                sample_vectors_++;
            }
        }
        // end = chrono::system_clock::now();
        // time1 += (double)chrono::duration_cast<chrono::microseconds>(end-start).count()/1000;

        // time2
        // start = chrono::system_clock::now();
        current_norm = GaussReduce_Parallel();
        // end = chrono::system_clock::now();
        // time2 += (double)chrono::duration_cast<chrono::microseconds>(end-start).count()/1000;

        // time3
        // start = chrono::system_clock::now();
        if(current_norm != -1 && current_norm < min_norm_){
            min_norm_ = current_norm;
        }

        // vector<ListPoint*> tmp_L;
        // size_t i = 0, j = 0;
        // while(i < V.size() && j < L.size()){
        //     if(V[i]->norm < L[j]->norm){
        //         tmp_L.emplace_back(V[i]);
        //         i++;
        //     }
        //     else{
        //         tmp_L.emplace_back(L[j]);
        //         j++;
        //     }
        // }
        // while(i < V.size()){
        //     tmp_L.emplace_back(V[i]);
        //     i++;
        // }
        // while(j < L.size()){
        //     tmp_L.emplace_back(L[j]);
        //     j++;
        // }
        // L = tmp_L;
        for(int i = 0; i < (int)V.size(); i++){
            int id = 0;
            while(id < (int)L.size() && L[id]->norm < V[i]->norm) id++;
            if(id < (int)L.size()) L.insert(L.begin()+id, V[i]);
            else L.push_back(V[i]);
        }
        V.clear();
        // end = chrono::system_clock::now();
        // time3 += (double)chrono::duration_cast<chrono::microseconds>(end-start).count()/1000;

        iterations_++;
    }

    // chk_time.push_back(time1);
    // chk_time.push_back(time2);
    // chk_time.push_back(time3);
}