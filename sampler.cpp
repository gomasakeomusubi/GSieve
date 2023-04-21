#include "sampler.h"

long KleinSampler::SampleZ(RR c_, RR s_square_){
    double c = to_double(c_);
    double s_square = to_double(s_square_);
    double s = sqrt(s_square);
    long minimun_c = floor(c - s * t_);
    long maximum_c = ceil(c + s * t_);
    long x;
    double rho;
    while(true){
        x = minimun_c + round((maximum_c - minimun_c) * double(rand()) / RAND_MAX);
        rho = exp(-M_PI * (x - c) * (x - c) / s_square);
        if((double(rand()) / RAND_MAX) <= rho){
            return x;
        }
    }
}

void KleinSampler::Init(const mat_ZZ &B){
    cout << "aaa" << endl;
    n_ = B.NumRows();
    cout << "bbb" << endl;
    m_ = B.NumCols();
    vec_RR Bstar_square;
    ComputeGS(B, mu_, Bstar_square);
    coef_.SetLength(n_);
    s_prime_square_.SetLength(n_);
    long max_star_sqr_norm = 0;
    for(int i = 0; i < n_; i++){
        max_star_sqr_norm = max(max_star_sqr_norm, to_long(Bstar_square[i]));
    }
    t_ = log(n_);
    double s_square = max_star_sqr_norm * log(n_);
    for(int i = 0; i < n_; i++){
        s_prime_square_[i] = s_square / to_double(Bstar_square[i]);
    }
}

LatticeVector* KleinSampler::Sample(){
    LatticeVector* lv = newLatticeVector(m_);
    for(int i = 0; i < n_; i++) coef_[i] = 0;
    for(int i = n_-1; i >= 0; i--){
        coef_[i] = to_RR(SampleZ(coef_[i], s_prime_square_[i]));
        for(int j = 0; j < i; j++){
            coef_[j] -= (coef_[i] * mu_[i][j]);
        }
    }
    lv->norm2 = 0;
    for(int i = 0; i < m_; i++){
        for(int j = 0; j < n_; j++){
            lv->vec[i] += to_ZZ(coef_[j] * to_RR(B_[j][i]));
        }
        lv->norm2 += lv->vec[i] * lv->vec[i];
    }
    return lv;
}