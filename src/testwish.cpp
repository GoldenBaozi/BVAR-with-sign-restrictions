#include <RcppArmadillo.h>
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
mat test_invwish(mat mean, int nu) {
    mat out = iwishrnd(mean, nu);
    return out;
}