# Estimation and Inference of VAR with Bayesian methods

- [x] basic VAR, OLS Estimation
- [x] IRF, FEVD, HDC and plots
- [ ] identification via recursive and IV approach (will be finished on about 11.20)
- [ ] inference via bootstrap (optional, not required with bayesian method)
- [x] bayesian VAR (`R` code is done, need to modify to `Rcpp` code)
- [x] identification via narrative sign restrictions (Antolín-Díaz and Rubio-Ramírez, 2018) with bayesian method, if time permits
- [ ] testing on examples (SW 2001, AR 2018, etc.)
- [ ] R package

The sign restriction procedure cannot be vectorized, must put it in Rcpp !!

For detailed illustration, see https://goldenbaozi.github.io/BayesVAR.html

# Progress

- most code finished, start testing using AR 2018 data
- current results not satisfied, IRF HPD set too wide.
- **time consuming**: use all 12 restrictions of example 1, 10 draws require about 2500 seconds.