#include "sample.h"
#include "tool.h"

void sample::gen_rref_mat(mat_ZZ& A, int seed){
    // generate a random upper-triangular matrix
    double alpha, beta;
    int dim;
    int i,j,k;
    ZZ sign;
  
    dim = A.NumCols();
  
    // generate rand number from 
    // beta distribution from gsl lib.

    const gsl_rng_type * T; 
    gsl_rng * r;
  
    gsl_rng_env_setup();

    T = gsl_rng_default; //default is mt19937
    r = gsl_rng_alloc (T);
  
    gsl_rng_set(r, seed);  

    default_random_engine generator(seed); //set another generator for uni_dist
    uniform_real_distribution<double> distribution(0.0,1.0);
  
    for(i=0;i<dim-1;i++){
        j = i+1;
        alpha = 6;
        beta = 4;
        for(k=0;k<i+1;k++){
            sign =  (to_ZZ)(2 * distribution(generator) -1);
            if(sign == 0) sign = 1;
            A[k][j] = sign * to_ZZ(gsl_ran_beta(r, alpha, beta) * (1 - (1.0*k)/dim) * 7 +0.5);
        }
    }
}

NTL::ZZ sample::max(mat_ZZ& A){
    //returns maximum absolute element
    int i,j,dim;
    ZZ max_ent;

    dim = A.NumCols();
    max_ent = abs(A[0][0]);
  
    for(i=0;i<dim;i++)
        for(j=0;j<dim;j++)
	        max_ent = max(max_ent,abs(A[i][j]));
    return max_ent;
}

NTL::ZZ sample::max(vec_ZZ& v){
    //returns maximum absolute element
    int i,dim;
    ZZ max_ent;

    dim = v.length();
    max_ent = abs(v[0]);
  
    for(i=1;i<dim;i++)
        max_ent = max(max_ent,abs(v[i]));
    return max_ent;
}

void sample::gen_random_unimodular2(mat_ZZ &L,int dim,int seed,int bits,int vl) {
    
    int i,j,k;
    int reCompute, tries, max_tries;
    max_tries = 5000;  // in sage the default max_tries=200
    tries = 0;
  
    mat_ZZ A, A_copy;
    ZZ max_ent,upper_bound;
    ident(A, dim);
    upper_bound = power2_ZZ(bits);

    int row_index;
    SetSeed(to_ZZ(seed));

    double time_start, time_start1;

    time_start1 = clock();
  
    if(bits == 0){
        // no upper bound for entries of unimodular matrix

        time_start = clock();
        /****** step 1. to generate a random upper triangular matrix ******/
        if (vl>=2) cout << "##### Step 1. generate a random upper triangular matrix. #####" << endl;
        
        gen_rref_mat(A, seed);

        if (vl>=2) cout << "The upper-triangular matrix is as:\n" << A << endl;
        /****** step 1 over. ******/
        
        if (vl>=1) cout << "elapsed time for step 1: " << (clock()- time_start) / CLOCKS_PER_SEC << "sec." << endl;

        /****** step 2. do rows transformations. ******/
        if (vl>=2) cout << "\n##### Step 2. do rows transformations. #####" << endl;

        for(k=dim-1;k>-1;k--){
            //cout << k << endl;
            row_index = 0;
            while(row_index < dim){
                if(row_index == k)
                    row_index++;
                if(k != row_index && row_index != dim){
                    A[row_index] = A[row_index] + A[k] * (RandomBnd(11)-5);
                    row_index++;
                }
            }
        }
        
        i=RandomBnd(dim-1)+1;
        j=RandomBnd(7)-3;
        if (vl>=2) cout << "Finally, to add the multiple of the " << i << "th row by " << j << " to the first row.\n"<< endl;
        A[0] += A[i] * j;

        if (vl>=1) cout << "elapsed time for step 2: " << (clock()- time_start) / CLOCKS_PER_SEC << "sec." << endl;
        if (vl>=2) cout << "##### The unimodular matrix is as: #####\n" << A << endl;
        if(determinant(A) == 1){
            if (vl>=2) cout << "check determinant: ok!" << endl;
        }
        else
            cout << "determinant ERROR!" << endl;
        
        }else{ // the entries of unimodular matrix are bounded in bits bits
        time_start = clock();
        /****** step 1. to generate a random upper triangular matrix ******/
        if (vl>=2) cout << "##### Step 1. generate a random upper triangular matrix. #####" << endl;
        do{
            reCompute = 0;
            gen_rref_mat(A, seed);
            max_ent = max(A);
            if(max_ent >= upper_bound){
                if (vl>=2) cout << "max entry in A: " << NumBits(max_ent) << " bits." << endl;;
                reCompute = 1;
                seed += 300;
            }
            tries++;
            if (vl>=2) cout << "tried " << tries << " times to generate the upper-tri-mat." << endl;
        }while(reCompute == 1 && tries < max_tries);

        if (vl>=2) cout << "The upper-triangular matrix is as:\n" << A << endl;
        /****** step 1 over. ******/
        
        if (vl>=1) cout << "elapsed time for step 1: " << (clock()- time_start) / CLOCKS_PER_SEC << "sec." << endl;
        
        /****** step 2. do rows transformations. ******/
        if (vl>=2) cout << "\n##### Step 2. do rows transformations. #####" << endl;

        A_copy = A;
        max_ent = max(A);
        vec_ZZ next_row;
        int nn = A.NumCols();
        int mm = A.NumRows();
        next_row.SetLength(nn);
        ZZ prev_max_ent;
        vec_ZZ rowmax;
        rowmax.SetLength(mm);
        for (k=0;k<mm;k++) rowmax[k] = max(A[k]);

        for(k=dim-1;k>-1;k--){
            row_index = tries = 0;
            while(row_index < dim){
            //cout << k << " ### " << row_index <<endl;
            prev_max_ent = max_ent;
            //cout << prev_max_ent << " " << max_ent << endl;
            if(k != row_index){
                next_row = A[row_index] + A[k] * (RandomBnd(11)-5);
                tries++;
            } else {
                next_row = A[row_index];
            }
            //cout << next_row << endl;
            max_ent = max(next_row);
            for (int j=0;j<dim;j++) {
                if (j!=row_index) max_ent = max(max_ent,rowmax[j]);
            }   //max_ent is the maximum of new matrix

            if(max_ent < upper_bound){
                //update the matrix
                A[row_index] = next_row;
                rowmax[row_index] = max(A[row_index]);
                row_index++;
                tries = 0;
            } else {
            }
            if(tries > max_tries){
                cout << "Error! set more tries." << endl;
                cout << A[k]<<endl;
                cout << next_row << endl;
                exit(1);
            }
            }
        }
        if (vl>=1) cout << "The output of step 2:\n" << A << endl;
        /****** setp 2 over ******/

        if (vl>=1) cout << "elapsed time for step 2: " << (clock()- time_start) / CLOCKS_PER_SEC << "sec." << endl;

        
        /****** step 3. to alter the first row vector.  ******/
        if (vl>=1) cout << "\n##### Step 3. alter the first row. #####" << endl;
        for (k=0;k<dim;k++) rowmax[k] = max(A[k]);
        max_ent = max(rowmax);
        k = 0;
        while(k<1){
            i=RandomBnd(dim-1)+1;
            j=(RandomBnd(7)-3);
            if (vl>=2) cout << "Finally, to add the multiple of the " << i << "th row by " << j << " to the first row.\n"<< endl;
            //A_copy[k] = A[k] + A[i] * j;
            next_row = A[k] + A[i] * j;
            max_ent = max(next_row);
            for (int j=0;j<dim;j++) {
                if (j!=k) max_ent = max(max_ent,rowmax[j]);
            }   //max_ent is the maximum of new matrix
            if(max_ent < upper_bound){
                A[k] = next_row;
                k++;
            }
        }

        /***** step 3 over *****/
        
        if (vl>=1) cout << "elapsed time for step 3: " << (clock()- time_start) / CLOCKS_PER_SEC << "sec." << endl;
        
        if (vl>=2) cout << "##### The unimodular matrix is as: #####\n" << A << endl;
        if(determinant(A) == 1) {
            if (vl>=2) cout << "check determinant: ok!" << endl;
        } else {
            cout << "determinant ERROR:" << determinant(A) << endl;
        }    
    }
    if (vl>=1) cout << "elapsed time: " << (clock()- time_start1) / CLOCKS_PER_SEC << "sec." << endl;

    L = A;
}

void sample::sample(mat_ZZ B, LatticeVector *ret) {
    mat_ZZ unim_mat;
    int seed = std::chrono::system_clock::now().time_since_epoch().count() + rand() % 1000;
    gen_random_unimodular2(unim_mat, B.NumCols(), seed , 16);
    ret->vec.SetLength(B.NumCols());
    mat_ZZ new_sample;
    mul(new_sample, unim_mat, B);
    auto temp = unim_mat * B;
    ret->vec = temp[rand() % B.NumCols()];
    ret->norm2 = norm2(ret->vec);
    //std::cout << "ret = " << ret << std::endl;
    //std::cout << "end_sample" << std::endl;
}

void sample::sample_set(mat_ZZ B, vector<LatticeVector*>& V, long num){
    mat_ZZ unim_mat;
    int seed = 0;
    LatticeVector* ret;
    for(int i = 1; i <= num; i++){
        seed = std::chrono::system_clock::now().time_since_epoch().count() + rand() % 1000;
        gen_random_unimodular2(unim_mat, B.NumCols(), seed , 16);
        ret = newLatticeVector(B.NumCols());
        mat_ZZ new_sample;
        mul(new_sample, unim_mat, B);
        auto temp = unim_mat * B;
        ret->vec = temp[rand() % B.NumCols()];
        ret->norm2 = norm2(ret->vec);
        V.push_back(ret);
    }
}
