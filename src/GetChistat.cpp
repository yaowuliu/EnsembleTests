// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
arma::mat Get_chistat (arma::mat X, arma::mat Sigma, int B, int m_select, bool is_Sigma_identity){
    const int p = X.n_rows;
    const int N_rep =  X.n_cols;
    arma::mat chistat(N_rep,B);
    
    //initialize
    Rcpp::IntegerVector ind_all = Rcpp::seq(0,p-1);
    Rcpp::IntegerVector ind0(m_select);
    arma::uvec ind(m_select);
    arma::vec v(p);
    arma::vec v_sub(m_select);
    arma::mat Sigma_sub(m_select,m_select);
    arma::mat Sigma_sub_inv(m_select,m_select);
    arma::mat val0(1,1);
    
    // loop
    for (int j = 0; j < N_rep; j++){
        v = vectorise(X.col(j));
        
        for (int k = 0; k < B; k++){
            ind0 = Rcpp::sample(ind_all,m_select,false);
            ind = Rcpp::as<arma::uvec>(ind0);
            
            v_sub = v.elem(ind);
            if (is_Sigma_identity){
                val0 = v_sub.t()*v_sub;
            }else{
                Sigma_sub = Sigma.submat(ind,ind);
                Sigma_sub_inv = inv(Sigma_sub);
                
                val0 = v_sub.t()*(Sigma_sub_inv*v_sub);
            }
            
            chistat(j,k) = val0(0,0);
        }
        
    }
    
    return(chistat);
}

