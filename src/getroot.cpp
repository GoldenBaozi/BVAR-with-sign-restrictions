#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
vec get_root(mat beta)
{
    int nvar = beta.n_cols;
    int plag = (beta.n_rows - 1) / nvar;
    mat beta_1 = beta.t();
    mat beta_lag = beta_1.submat(0, 1, beta_1.n_rows - 1, beta_1.n_cols - 1);
    mat large_eye = eye(nvar * (plag - 1), nvar * (plag - 1));
    mat large_zero = zeros((plag - 1) * nvar, nvar);
    mat beta_compact = join_cols(beta_lag, join_rows(large_eye, large_zero));
    cx_vec roots = eig_gen(beta_compact);
    vec root = abs(roots);
    return root;
}