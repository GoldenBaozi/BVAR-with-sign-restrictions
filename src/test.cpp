#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
mat test_resize() {
    mat x(4, 1);
    x.col(0) = linspace(1,4,4);
    mat y = reshape(x, 2, 2);
    cube z(2, 2, 2);
    Rcout << z.n_rows << "\n";
    z.slice(0) = y;
    z.slice(1) = y - 1;
    cube haha = reshape(z, 2, 4, 1);
    mat out = haha.slice(0);
    return out;
}

// [[Rcpp::export]]
mat flatten_cube(cube obj) {
    int hor = obj.n_slices;
    int nvar = obj.n_rows;
    mat out(nvar*nvar, hor);
    for (int i = 0; i < nvar*nvar; i++)
    {
        int row = i / nvar;
        int col = i % nvar;
        cube tmp = reshape(obj.subcube(row, col, 0, row, col, hor-1), 1, hor, 1);
        out.row(i) = tmp.slice(0);
    }
    return out;
}

/*** R
input <- array(1:8, c(2,2,2))
input
flatten_cube(input)
*/