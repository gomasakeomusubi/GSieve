#include "gsieve.h"
#include "sample.h"
#include "tool.h"

NTL_CLIENT

void GSieve::Setup(){
    n_ = B.NumRows();
    m_ = B.NumCols();
    iterations_ = 0;
    collisions_ = 0;
    sample_vectors_ = 0;
    CleanUp();

    min_vector_ = newLatticeVector(m_);
    min_vector_->vec = B[0];
    min_vector_->norm2 = norm2(B[0]);
    LatticeVector* lv;
    for(int i = 0; i < n_; i++){
        lv = newLatticeVector(m_);
        lv->vec = B[i];
        lv->norm2 = norm2(B[i]);
        GaussReduce(lv);
        if(lv->norm2 < min_vector_->norm2){
            min_vector_->vec = lv->vec;
            min_vector_->norm2 = lv->norm2;
        }
    }
    cout << L.size() << "," << sample_vectors_ << "," << iterations_ << "," << collisions_ << endl;
    max_list_size_ = L.size();
    concurrency_ = 1;
    simu_samp_ = 1;
}

void GSieve::SampleReduce(LatticeVector *p){
    bool vec_change = true;
    while(vec_change){
        vec_change = false;
        for(int i = 0; i < L.size(); i++){
            if(p->norm2 < L[i]->norm2){
                continue;
            }
            if(reduceVector(p, L[i])){
                vec_change = true;
            }
        }
    }
}

void GSieve::SampleReduce_Parallel(){
    vector<bool> vec_change(V.size(), false);
    int num_parallel = ceil((double)V.size() / concurrency_);
    int num_vectors_V = V.size();
    int num_vectors = 0;

    // 各スレッドにnum_parallel個のベクトルを分配
    // 排他制御なし並列処理
    vector<thread> threads;
    for(int i = 0; i < concurrency_; i++){
        num_vectors = min(num_parallel, num_vectors_V);
        num_vectors_V -= num_vectors;

        threads.emplace_back(
            [&vec_change, i, num_parallel, this](int num){
                for(int j = 0; j < num; j++){
                    LatticeVector *lv = newLatticeVector(m_);
                    lv->vec = V[i*num_parallel+j]->vec;
                    lv->norm2 = V[i*num_parallel+j]->norm2;

                    SampleReduce(V[i*num_parallel+j]);

                    if(lv->vec != V[i*num_parallel+j]->vec){
                        vec_change[i*num_parallel+j] = true;
                    }

                    delete lv;
                }
            }, num_vectors
        );
    }
    for(thread &th : threads){
        th.join();
    }

    // 更新したベクトルがあればSに移動
    for(int i = 0; i < vec_change.size(); i++){
        if(vec_change[i] == true){
            S.push(V[i]);
            V.erase(V.begin()+i);
            vec_change.erase(vec_change.begin()+i);
            i--;
        }
    }
}

void GSieve::ListReduce(LatticeVector *p){
    for(int i = 0; i < L.size(); i++){
        if(p->norm2 >= L[i]->norm2){
            continue;
        }
        if(reduceVector(L[i], p)){
            S.push(L[i]);
            L.erase(L.begin()+i);
            i--;
        }
    }
}

void GSieve::ListReduce_Parallel(){
    vector<bool> vec_change(L.size(), false);
    int num_parallel = L.size() / concurrency_;
    int num_remain = L.size() % concurrency_;

    // 各スレッドにnum_parallel個のベクトルを分配
    vector<thread> threads;
    for(int i = 0; i < concurrency_; i++){
        threads.emplace_back(
            [&vec_change, i, this](int num){
                for(int j = 0; j < num; j++){
                    for(int k = 0; k < V.size(); k++){
                        if(L[i*num+j]->norm2 < V[k]->norm2){
                            continue;
                        }
                        if(reduceVector(L[i*num+j], V[k])){
                            vec_change[i*num+j] = true;
                        }
                    }
                }
            }, num_parallel
        );
    }
    for(thread &th : threads){
        th.join();
    }

    // 余ったベクトルは本スレッドで
    int l_sz = L.size();
    for(int l_id = concurrency_ * num_parallel; l_id < l_sz; l_id++){
        for(int v_id = 0; v_id < V.size(); v_id++){
            if(L[l_id]->norm2 < V[v_id]->norm2){
                continue;
            }
            if(reduceVector(L[l_id], V[v_id])){
                vec_change[l_id] = true;
            }
        }
    }

    // 更新したベクトルがあればSに移動
    for(int i = 0; i < vec_change.size(); i++){
        if(vec_change[i] == true){
            S.push(L[i]);
            L.erase(L.begin()+i);
            vec_change.erase(vec_change.begin()+i);
            i--;
        }
    }
}

void GSieve::VectorReduce_Parallel(){
    vector<LatticeVector*> Tmp;
    for(int i = 0; i < V.size(); i++){
        LatticeVector *lv = newLatticeVector(m_);
        lv->vec = V[i]->vec;
        lv->norm2 = V[i]->norm2;
        Tmp.emplace_back(lv);
    }
    cout << 111 << endl;

    vector<bool> vec_change(V.size(), false);
    int num_parallel = V.size() / concurrency_;
    int num_remain = V.size() % concurrency_;

    // 各スレッドにnum_parallel個のベクトルを分配
    vector<thread> threads;
    for(int i = 0; i < concurrency_; i++){
        threads.emplace_back(
            [&vec_change, i, this, &Tmp](int num){
                for(int j = 0; j < num; j++){
                    for(int k = 0; k < V.size(); k++){
                        if(i*num+j == k) continue;
                        if(Tmp[i*num+j]->norm2 < V[k]->norm2){
                            continue;
                        }
                        if(reduceVector(Tmp[i*num+j], V[k])){
                            vec_change[i*num+j] = true;
                        }
                    }
                }
            }, num_parallel
        );
    }
    for(thread &th : threads){
        th.join();
    }

    // 余ったベクトルは本スレッドで
    int tmp_sz = Tmp.size();
    for(int tmp_id = concurrency_ * num_parallel; tmp_id < tmp_sz; tmp_id++){
        for(int v_id = 0; v_id < V.size(); v_id++){
            if(tmp_id == v_id) continue;
            if(Tmp[tmp_id]->norm2 < V[v_id]->norm2){
                continue;
            }
            if(reduceVector(Tmp[tmp_id], V[v_id])){
                vec_change[tmp_id] = true;
            }
        }
    }

    // 更新したベクトルがあればSに移動
    for(int i = 0; i < vec_change.size(); i++){
        if(vec_change[i] == true){
            S.push(Tmp[i]);
            // delete V[i];
            V.erase(V.begin()+i);
            vec_change.erase(vec_change.begin()+i);
            i--;
        }
        // else{
        //     delete Tmp[i];
        // }
    }
}

void GSieve::GaussReduce(LatticeVector *p){
    // p <- L
    SampleReduce(p);

    if(p->norm2 == 0) return;

    // p -> L
    ListReduce(p);

    L.emplace_back(p);
}

void GSieve::GaussReduce_Parallel(){
    // V <- L
    SampleReduce_Parallel();

    // V <- V
    vector<bool> vec_change2(V.size(), false);
    bool vec_change = true;
    while(vec_change){
        vec_change = false;
        for(int i = 0; i < V.size(); i++){
            for(int j = 0; j < V.size(); j++){
                if(i != j){
                    if(V[i]->norm2 < V[j]->norm2){
                        continue;
                    }
                    if(reduceVector(V[i], V[j])){
                        vec_change = true;
                        vec_change2[i] = true;
                    }
                }
            }
        }
    }

    for(int i = 0; i < vec_change2.size(); i++){
        if(vec_change2[i] == true){
            S.push(V[i]);
            V.erase(V.begin()+i);
            vec_change2.erase(vec_change2.begin()+i);
            i--;
        }
    }
    // VectorReduce_Parallel();

    // V -> L
    ListReduce_Parallel();
}

LatticeVector *GSieve::GaussSieve(vector<double> &chk_time){
    chrono::system_clock::time_point start, end;

    double time1 = 0, time2 = 0, time3 = 0;
    while(collisions_ < max_list_size_/10+200 && min_vector_->norm2 > goal_norm2_){
        // time1
        start = chrono::system_clock::now();
        max_list_size_ = std::max(max_list_size_, (long)L.size());

        LatticeVector *new_v = newLatticeVector(m_);
        if(!S.empty()){
            new_v = S.front();
            S.pop();
        }else{
            sample::sample(B, new_v);
            sample_vectors_++;
        }
        end = chrono::system_clock::now();
        time1 += (double)chrono::duration_cast<chrono::microseconds>(end-start).count()/1000;


        // time2
        start = chrono::system_clock::now();
        GaussReduce(new_v);
        end = chrono::system_clock::now();
        time2 += (double)chrono::duration_cast<chrono::microseconds>(end-start).count()/1000;
        

        // time3
        start = chrono::system_clock::now();
        if(new_v->norm2 == 0){
            delete new_v;
            collisions_++;
        }else{
            if(new_v->norm2 < min_vector_->norm2){
                min_vector_->norm2 = new_v->norm2;
                min_vector_->vec = new_v->vec;
            }
        }
        end = chrono::system_clock::now();
        time3 += (double)chrono::duration_cast<chrono::microseconds>(end-start).count()/1000;

        iterations_++;
    }

    cout << "time1: " << time1 << endl;
    cout << "time2: " << time2 << endl;
    cout << "time3: " << time3 << endl;
    chk_time.push_back(time1);
    chk_time.push_back(time2);
    chk_time.push_back(time3);

    return min_vector_;
}

LatticeVector *GSieve::GaussSieve_Parallel(vector<double> &chk_time){
    chrono::system_clock::time_point start, end;

    ZZ K; K = 0;
    int max_list_size = 1;
    double time1 = 0, time2 = 0, time3 = 0;
    while(K < max_list_size/10+200 && min_vector_->norm2 > goal_norm2_){
        // time1
        start = chrono::system_clock::now();
        max_list_size = std::max(max_list_size, int(L.size()));
        LatticeVector *new_v = newLatticeVector(m_);
        if(!S.empty()){
            if(simu_samp_ > S.size()){
                sample::sample_set(B, V, simu_samp_-S.size());
                sample_vectors_ += simu_samp_-S.size();
                while(!S.empty()){
                    new_v = S.front();
                    V.push_back(new_v);
                    S.pop();
                }
            }
            else{
                for(int i = 0; i < simu_samp_; i++){
                    new_v = S.front();
                    V.push_back(new_v);
                    S.pop();
                }
            }
        }else{
            sample::sample_set(B, V, simu_samp_);
            sample_vectors_ += simu_samp_;
        }
        end = chrono::system_clock::now();
        time1 += (double)chrono::duration_cast<chrono::microseconds>(end-start).count()/1000;

        // cout << "[a]." << L.size() << "." << S.size() << "." << V.size() << "." << flush;
        

        // time2
        start = chrono::system_clock::now();
        GaussReduce_Parallel();
        end = chrono::system_clock::now();
        time2 += (double)chrono::duration_cast<chrono::microseconds>(end-start).count()/1000;
        
        // cout << "[b]." << L.size() << "." << S.size() << "." << V.size() << "." << flush;


        // time3
        start = chrono::system_clock::now();
        for(auto &itr_V: V){
            if(itr_V->norm2 <= 0){
                delete itr_V;
                K++;
            }else{
                L.emplace_back(itr_V);
                if(itr_V->norm2 < min_vector_->norm2){
                    min_vector_->norm2 = itr_V->norm2;
                    min_vector_->vec = itr_V->vec;
                }
            }
        }
        V.clear();
        end = chrono::system_clock::now();
        time3 += (double)chrono::duration_cast<chrono::microseconds>(end-start).count()/1000;

        // cout << "[c]." << L.size() << "." << S.size() << "." << V.size() << "." << K << "_" << endl;

        if(!(K < max_list_size/10+200)) cout << "cond 1:" << K << endl;
        if(!(min_vector_->norm2 > goal_norm2_)) cout << "cond 2:" << min_vector_->norm2 << endl;

        iterations_++;
    }

    chk_time.push_back(time1);
    chk_time.push_back(time2);
    chk_time.push_back(time3);

    return min_vector_;
}