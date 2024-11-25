#include <RcppArmadillo.h>
// #include <progress.hpp>
// #include <progress_bar.hpp>
#include <chrono>

// [[Rcpp::depends(RcppArmadillo)]]

/*
[[Rcpp::depends(RcppProgress)]]
*/

using namespace Rcpp;
using namespace arma;
using namespace std;
using namespace chrono;

// [[Rcpp::export]]
List IRF_compute(mat beta, mat B, int hor, int nvar, int plag)
{
    mat beta_1 = beta.t();
    mat beta_lag = beta_1.submat(0, 1, beta_1.n_rows - 1, beta_1.n_cols - 1);
    // Rcout << beta_lag.n_cols << "\n";
    mat large_eye = eye(nvar * (plag - 1), nvar * (plag - 1));
    mat large_zero = zeros((plag - 1) * nvar, nvar);
    // Rcout << large_eye.n_cols << "\n";
    // Rcout << large_zero.n_cols << "\n";
    mat beta_compact = join_cols(beta_lag, join_rows(large_eye, large_zero));
    mat irf_trans = join_rows(eye<mat>(nvar, nvar), zeros<mat>(nvar, nvar * (plag - 1)));
    cube Psi(nvar, nvar, hor + 1);
    cube IRF(nvar, nvar, hor + 1);
    for (int h = 0; h < hor + 1; h++)
    {
        Psi.slice(h) = irf_trans * powmat(beta_compact, h) * irf_trans.t();
        IRF.slice(h) = Psi.slice(h) * B;
    }
    List out = List::create(
        _["Psi"] = Psi,
        _["IRF"] = IRF);
    return out;
}

cube IRF_compute_1(mat beta, mat B, int hor, int nvar, int plag)
{
    mat beta_1 = beta.t();
    mat beta_lag = beta_1.submat(0, 1, beta_1.n_rows - 1, beta_1.n_cols - 1);
    // Rcout << beta_lag.n_cols << "\n";
    mat large_eye = eye(nvar * (plag - 1), nvar * (plag - 1));
    mat large_zero = zeros((plag - 1) * nvar, nvar);
    // Rcout << large_eye.n_cols << "\n";
    // Rcout << large_zero.n_cols << "\n";
    mat beta_compact = join_cols(beta_lag, join_rows(large_eye, large_zero));
    mat irf_trans = join_rows(eye<mat>(nvar, nvar), zeros<mat>(nvar, nvar * (plag - 1)));
    cube IRF(nvar, nvar, hor + 1);
    for (int h = 0; h < hor + 1; h++)
    {
        mat Psi = irf_trans * powmat(beta_compact, h) * irf_trans.t();
        IRF.slice(h) = Psi * B;
    }

    return IRF;
}

// [[Rcpp::export]]
double max_root(mat beta)
{
    int nvar = beta.n_cols;
    int plag = (beta.n_rows - 1) / nvar;
    mat beta_1 = beta.t();
    mat beta_lag = beta_1.submat(0, 1, beta_1.n_rows - 1, beta_1.n_cols - 1);
    mat large_eye = eye(nvar * (plag - 1), nvar * (plag - 1));
    mat large_zero = zeros((plag - 1) * nvar, nvar);
    mat beta_compact = join_cols(beta_lag, join_rows(large_eye, large_zero));
    cx_vec roots = eig_gen(beta_compact);
    double max_root = max(abs(roots));
    return max_root;
}

// [[Rcpp::export]]
cube FEVD_compute(mat Sigma, mat B, cube Psi, int nvar, int hor)
{
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
            MSE.slice(j) = MSE.slice(j - 1) + Psi.slice(j) * Sigma * trans(Psi.slice(j));
            MSE_shock.slice(j) = MSE_shock.slice(j - 1) + Psi.slice(j) * MSE_shock.slice(0) * trans(Psi.slice(j));
            FEVD.slice(j).col(i) = diagvec(MSE_shock.slice(j)) / diagvec(MSE.slice(j));
        }
    }
    return FEVD;
}

// [[Rcpp::export]]
double HDC(int start, int end, int shock, int res_var, cube IRF, mat eps)
{
    int t = start - 1;
    int h = end - 1;
    int i = res_var - 1;
    int j = shock - 1;
    double HDC = 0.0;
    for (int k = t; k < h + 1; k++)
    {
        HDC += eps(k, j) * IRF(i, j, h - k);
    }
    return HDC;
}

// [[Rcpp::export]]
vec HDC_ts(int start, int end, int shock, int res_var, cube IRF, mat eps)
{
    vec HDC_ts(end - start + 1);
    for (int t = start; t < end + 1; t++)
    {
        HDC_ts(t - start) = HDC(start, t, shock, res_var, IRF, eps);
    }
    return HDC_ts;
}

// [[Rcpp::export]]
mat flatten_cube(cube obj)
{
    int hor = obj.n_slices;
    int nvar = obj.n_rows;
    mat out(nvar * nvar, hor);
    for (int i = 0; i < nvar * nvar; i++)
    {
        int row = i / nvar;
        int col = i % nvar;
        cube tmp = reshape(obj.subcube(row, col, 0, row, col, hor - 1), 1, hor, 1);
        out.row(i) = tmp.slice(0);
    }
    return out;
}

// [[Rcpp::export]]
cube IRF_draws(cube beta_draw, cube B_draw, int hor)
{
    int nvar = B_draw.n_cols;
    int draw_num = B_draw.n_slices;
    int plag = (beta_draw.n_rows - 1) / nvar;
    cube IRF_draws(nvar * nvar, hor + 1, draw_num);
    for (int i = 0; i < draw_num; i++)
    {
        mat beta = beta_draw.slice(i);
        mat B = B_draw.slice(i);
        cube IRF_1draw = IRF_compute_1(beta, B, hor, nvar, plag);
        IRF_draws.slice(i) = flatten_cube(IRF_1draw);
    }
    return IRF_draws;
}

bool check_SR_and_EBR(List restrictions, cube IRF)
{
    int n = restrictions.length();
    int flag = 0;
    int num = 0;
    for (int i = 0; i < n; i++)
    {
        List res_1 = restrictions[i];
        string type = as<string>(res_1[0]);
        if (type == "SR")
        {
            num++;
            int shock = as<int>(res_1[1]) - 1;
            int var = as<int>(res_1[2]) - 1;
            int hor = as<int>(res_1[3]);
            double sign = as<double>(res_1[4]);
            mat IRF_test = IRF.slice(hor);
            double SR_test = IRF_test(var, shock) * sign;
            if (SR_test > 0)
            {
                flag++;
            }
        }
        if (type == "EBR")
        {
            num++;
            int shock = as<int>(res_1[1]) - 1;
            int var_1 = as<int>(res_1[2]) - 1;
            int var_2 = as<int>(res_1[3]) - 1;
            int hor = as<int>(res_1[4]);
            double mb = as<double>(res_1[5]);
            double lb = as<double>(res_1[6]);
            mat IRF_test = IRF.slice(hor);
            double EBR_test = IRF_test(var_1, shock) / IRF_test(var_2, shock);
            // Rcout << mb << " " << lb << " " << EBR_test << "\n";
            if (EBR_test <= mb && EBR_test >= lb)
            {
                flag++;
            }
        }
    }
    // Rcout << flag << "\n";
    bool out = (flag == num);
    return out;
}

bool check_NSR(List restrictions, cube IRF, mat eps)
{
    int n = restrictions.length();
    int flag = 0;
    int NSR_num = 0;
    for (int i = 0; i < n; i++)
    {
        List res_1 = restrictions[i];
        string type = as<string>(res_1[0]);
        if (type == "NSR")
        {
            NSR_num++;
            string NSR_type = as<string>(res_1[1]);
            if (NSR_type == "sign")
            {
                int shock = as<int>(res_1[2]) - 1;
                uvec period = as<uvec>(res_1[3]) - 1;
                double sign = as<double>(res_1[4]);
                mat eps_test = eps.submat(period(0), shock, period(period.n_elem - 1), shock) * sign;
                bool test = all(vectorise(eps_test) > 0);
                if (test)
                {
                    flag++;
                }
            }
            else
            {
                // here I don't set shock = ...-1, because the function HDC_ts will do this, to avoid messy
                int shock = as<int>(res_1[2]);
                int var = as<int>(res_1[3]);
                int nvar = eps.n_cols;
                uvec period = as<uvec>(res_1[4]);
                int start = period(0);
                int end = period(period.n_elem - 1);
                double sign = as<double>(res_1[5]);
                string inten = as<string>(res_1[6]);
                mat HDC_test(end - start + 1, nvar);
                for (int i = 0; i < nvar; i++)
                {
                    HDC_test.col(i) = HDC_ts(start, end, i + 1, var, IRF, eps);
                }
                mat HDC_target = repmat(HDC_test.col(shock - 1), 1, nvar);
                mat test_sign = HDC_test.col(shock - 1) * sign;
                mat test_inten = abs(HDC_test) - abs(HDC_target);
                bool cond0 = all(vectorise(test_sign) >= 0);
                bool cond1 = all(vectorise(test_inten) >= 0);
                bool cond2 = all(vectorise(test_inten) <= 0);
                if ((inten == "strong" && cond0 && cond1) || (inten == "weak" && cond0 && cond2))
                {
                    flag++;
                }
            }
        }
        else
        {
            continue;
        }
    }
    bool out = (flag == NSR_num);
    return out;
}

double compute_importance_weight(List restrictions, cube IRF, int M, int row, int col)
{
    int satisfy = 1; // to avoid zero weight
    for (int i = 0; i < M; i++)
    {
        mat eps_v = randn(row, col);
        bool cond = check_NSR(restrictions, IRF, eps_v);
        if (cond)
        {
            satisfy++;
        }
    }
    double weight = 1.0 / (satisfy * 1.0 / M);
    return weight;
}

// [[Rcpp::export]]
List sign_restrictions_main(List alpha_post, List Sigma_post, List restrictions, mat Y, mat X, int draw, int plag, int hor)
{
    int nvar = Y.n_cols;
    int T_est = Y.n_rows;
    // save draws
    cube B_draw(nvar, nvar, draw);
    cube beta_draw(nvar * plag + 1, nvar, draw);
    cube Sigma_draw(nvar, nvar, draw);
    // cube IRF_draw(nvar*nvar, hor+1, draw);
    vec weights(draw);
    // unpack priors
    vec alpha_post_mean = as<vec>(alpha_post[0]);
    mat alpha_post_cov = as<mat>(alpha_post[1]);
    mat Sigma_post_mean = as<mat>(Sigma_post[0]);
    int nu_post = as<int>(Sigma_post[1]);
    // flag and supervisor
    int flag = 0;
    int report = draw / 10;
    int try_per_beta_sigma = 100;
    auto start_time = high_resolution_clock::now();
    // Progress p(10, true);
    // loop on draw
    while (flag < draw)
    {
        mat Sigma_1draw = iwishrnd(Sigma_post_mean, nu_post);
        mat N0 = inv(X.t()*X);
        mat NT = kron(Sigma_1draw, N0);
        mat cov = chol(NT, "lower");
        vec alpha_1draw = alpha_post_mean + cov * randn(alpha_post_mean.n_elem, 1);
        mat beta_1draw = reshape(alpha_1draw, nvar * plag + 1, nvar);
        mat u_1draw = Y - X * beta_1draw;
        mat P = chol(Sigma_1draw, "lower");
        double maxroot = max_root(beta_1draw);
        if (maxroot <= 1)
        {
            for (int i = 0; i < try_per_beta_sigma; i++)
            {
                mat X = randn(nvar, nvar);
                mat Q, R;
                qr(Q, R, X);
                Q = Q * diagmat(sign(R.diag()));
                mat B_1draw = P * Q;
                mat B_inv_1draw = inv(trans(B_1draw));
                // mat eps_1draw = u_1draw * B_inv_1draw;
                cube IRF_1draw = IRF_compute_1(beta_1draw, B_1draw, hor, nvar, plag);
                // check sign restrictions, to decide save or not
                bool save_me = check_SR_and_EBR(restrictions, IRF_1draw);
                if (save_me)
                {
                    B_draw.slice(flag) = B_1draw;
                    beta_draw.slice(flag) = beta_1draw;
                    Sigma_draw.slice(flag) = Sigma_1draw;
                    // Rcout << IRF_1draw.n_slices << "\n";
                    // IRF_draw.slice(flag) = flatten_cube(IRF_1draw);
                    // weights(flag) = compute_importance_weight(restrictions, IRF_1draw, M, T_est, nvar);
                    // Rcout << weights(flag) <<"\n";
                    flag++;
                    if (flag % report == 0)
                    {
                        auto current_time = high_resolution_clock::now();
                        auto elapsed_time = duration_cast<seconds>(current_time - start_time).count();
                        Rcout << "-> " << flag << " draws saved, total time usage: " << elapsed_time << " seconds.\n";
                        // p.increment();
                    }
                }
            }
        }
        else
        {
            continue;
        }
        
    }
    Rcout << "* All draws are done\n";
    // re-weight using importance sampling
    // vec x = linspace(0, draw - 1, draw);
    // NumericVector xx = wrap(x);
    // NumericVector my_weights = wrap(weights);
    // NumericVector new_id = sample(xx, save, true, my_weights);
    // uvec save_id = as<uvec>(new_id);
    // cube B_saved = B_draw.slices(save_id);
    // cube beta_saved = beta_draw.slices(save_id);
    // cube Sigma_saved = Sigma_draw.slices(save_id);
    // cube IRF_saved = IRF_draw.slices(save_id);

    return List::create(
        _["B_saved"] = B_draw,
        _["beta_saved"] = beta_draw,
        _["Sigma_saved"] = Sigma_draw
        // _["IRF_saved"] = IRF_saved
    );
}
// [[Rcpp::export]]
List sign_restrictions_NSR(List restrictions, mat Y, mat X, cube B_draw, cube beta_draw, cube Sigma_draw, int plag, int hor, int M)
{
    int nvar = Y.n_cols;
    int T_est = Y.n_rows;
    int init_draw = B_draw.n_slices;
    // save draws
    cube B_NSR(nvar, nvar, init_draw);
    cube Sigma_NSR(nvar, nvar, init_draw);
    cube beta_NSR(nvar * plag + 1, nvar, init_draw);
    vec weights(init_draw);
    // flag and supervisor
    int flag = 0;
    // loop on draw
    for (int i = 0; i < init_draw; i++)
    {
        mat beta_1draw = beta_draw.slice(i);
        mat B_1draw = B_draw.slice(i);
        mat Sigma_1draw = Sigma_draw.slice(i);
        mat B_inv_1draw = inv(trans(B_1draw));
        mat eps_1draw = (Y - X * beta_1draw) * B_inv_1draw;
        cube IRF_1draw = IRF_compute_1(beta_1draw, B_1draw, hor, nvar, plag);
        bool save_me = check_NSR(restrictions, IRF_1draw, eps_1draw);
        if (save_me)
        {
            B_NSR.slice(flag) = B_1draw;
            beta_NSR.slice(flag) = beta_1draw;
            Sigma_NSR.slice(flag) = Sigma_1draw;
            weights(flag) = compute_importance_weight(restrictions, IRF_1draw, M, T_est, nvar);
            flag++;
        }
        else
        {
            continue;
        }
    }
    Rcout << "* " << flag << " draws in " << init_draw << " satisfy NSR, resampling using weights...\n";
    // re-weight using importance sampling
    vec x = linspace(0, flag - 1, flag);
    NumericVector xx = wrap(x);
    NumericVector my_weights = wrap(weights);
    NumericVector new_id = sample(xx, flag, true, my_weights);
    uvec save_id = as<uvec>(new_id);
    cube B_saved = B_NSR.slices(save_id);
    cube beta_saved = beta_NSR.slices(save_id);
    cube Sigma_saved = Sigma_NSR.slices(save_id);

    return List::create(
        _["B_NSR"] = B_saved,
        _["beta_NSR"] = beta_saved,
        _["Sigma_NSR"] = Sigma_saved);
}