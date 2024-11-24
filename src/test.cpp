#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

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

// [[Rcpp::export]]
int add_progress(int range) {
    int report = range / 10;
    Progress p(10, true);
    Function f("Sys.sleep");
    for (int i = 0; i < range; i++)
    {
        f(1);
        if (i % report == 0)
        {
            p.increment();
            Rcout << i << " draws done.\n";
        }
    }
    return 0;
}

/*** R
# input <- array(1:8, c(2,2,2))
# input
# flatten_cube(input)
add_progress(10)
*/