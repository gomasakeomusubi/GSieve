#include "sampler.h"

long KleinSampler::SampleZ(double c, double s_square){
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
    n_ = B.NumRows();
    m_ = B.NumCols();
    mat_RR mu;
    vec_RR Bstar_square;
    ComputeGS(B, mu, Bstar_square);
    MatInt64FromMatZZ(B, B_);
    MatDoubleFromMatRR(mu, mu_);
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

ListPoint* KleinSampler::Sample(){
    ListPoint* lp = NewListPoint(m_);
    for(int i = 0; i < n_; i++) coef_[i] = 0;
    for(int i = n_-1; i >= 0; i--){
        coef_[i] = SampleZ(coef_[i], s_prime_square_[i]);
        for(int j = 0; j < i; j++){
            coef_[j] -= (coef_[i] * mu_[i][j]);
        }
    }
    lp->norm = 0;
    for(int i = 0; i < m_; i++){
        for(int j = 0; j < n_; j++){
            lp->v[i] += coef_[j] * B_[j][i];
        }
        lp->norm += lp->v[i] * lp->v[i];
    }
    return lp;
}