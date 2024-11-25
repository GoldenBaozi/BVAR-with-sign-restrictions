# Estimation and Inference of VAR with Bayesian methods

- [x] basic VAR, OLS Estimation
- [x] IRF, FEVD, HDC and plots
- [ ] identification via recursive and IV approach (will be finished on about 11.20)
- [ ] inference via bootstrap (optional, not required with bayesian method)
- [x] bayesian VAR (`R` code is done, need to modify to `Rcpp` code)
- [x] identification via narrative sign restrictions (Antolín-Díaz and Rubio-Ramírez, 2018) with bayesian method, if time permits
- [x] testing on examples (SW 2001, AR 2018, etc.)
- [ ] R package

The sign restriction procedure cannot be vectorized, must put it in Rcpp !!

For detailed illustration, see https://goldenbaozi.github.io/BayesVAR.html

# Progress

- Main result get.
- **Note**: must check the max root of companion matrix, to ensure the stationary of the system !!!
- I re-read AR2018's data, and now strictly follow there procedure of checking restrictions. It greatly increases running efficiency
  - first, check NR and EBR, no NSR, and only compute IRF needed (don't save it)
  - for every drawing $\alpha$ and $\Sigma$, draw another 100 $Q$ trying to match them
  - second, impose NSR
  - third, compute IRF for `first` set and `second` set

# Rep result

![AR2018 figure 2 part 1](./out/IRF_oil_no_NSR.png)