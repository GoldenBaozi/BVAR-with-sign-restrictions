#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List IRF_compute(mat beta, mat B, int hor, int nvar, int plag) {
    mat beta_1 = beta.t();
    mat beta_lag = beta_1.submat(0, 1, beta_1.n_rows-1, beta_1.n_cols-1);
    // Rcout << beta_lag.n_cols << "\n";
    mat large_eye = eye(nvar * (plag - 1), nvar * (plag - 1));
    mat large_zero = zeros((plag - 1) * nvar, nvar);
    // Rcout << large_eye.n_cols << "\n";
    // Rcout << large_zero.n_cols << "\n";
    mat beta_compact = join_cols(beta_lag, join_rows(large_eye, large_zero));
    mat irf_trans = join_rows(eye<mat>(nvar, nvar), zeros<mat>(nvar, nvar*(plag-1)));
    cube Psi(nvar, nvar, hor+1);
    cube IRF(nvar, nvar, hor+1);
    for (int h = 0; h < hor+1; h++)
    {
        Psi.slice(h) = irf_trans * powmat(beta_compact, h) * irf_trans.t();
        IRF.slice(h) = Psi.slice(h) * B;
    }
    List out = List::create(
        _["Psi"] = Psi,
        _["IRF"] = IRF
    );
    return out;
}

// [[Rcpp::export]]
cube FEVD_compute(mat Sigma, mat B, cube Psi, int nvar, int hor) {
    int nstep = hor + 1;
    cube MSE(nvar, nvar, nstep);
    cube MSE_shock(nvar, nvar, nstep);
    cube FEVD(nvar, nvar, nstep);
    for (int i = 0; i < nvar; i++)
    {
        MSE.slice(0) = Sigma;
        MSE_shock.slice(0) = B.col(i) * trans(B.col(i));
        for (int j = 1; j < nstep; j++)
        {
            MSE.slice(j) = MSE.slice(j-1) + Psi.slice(j) * Sigma * trans(Psi.slice(j));
            MSE_shock.slice(j) = MSE_shock.slice(j-1) + Psi.slice(j) * MSE_shock.slice(0) * trans(Psi.slice(j));
            FEVD.slice(j).col(i) = diagvec(MSE_shock.slice(j)) / diagvec(MSE.slice(j));
        }
        
    }
    return FEVD;
}

// [[Rcpp::export]]
double HDC(int start, int end, int shock, int res_var, cube IRF, mat eps) {
    int t = start - 1;
    int h = end - 1;
    int i = res_var - 1;
    int j = shock - 1;
    double HDC = 0.0;
    for (int k = t; k < h + 1; k++)
    {
        HDC += eps(k, j)*IRF(i, j, h-k);
    }
    return HDC;
}

// [[Rcpp::export]]
vec HDC_ts(int start, int end, int shock, int res_var, cube IRF, mat eps) {
    vec HDC_ts(end-start+1);
    for (int t = start; t < end+1; t++)
    {
        HDC_ts(t - start) = HDC(start, t, shock, res_var, IRF, eps);
    }
    return HDC_ts;
}