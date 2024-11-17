# Estimation and Inference of VAR with Bayesian methods

- [x] basic VAR, OLS Estimation
- [x] IRF, FEVD, HDC and plots
- [ ] identification via recursive and IV approach (soon...)
- [ ] inference via bootstrap (optional, not required with bayesian method)
- [ ] bayesian VAR
- [ ] identification via narrative sign restrictions (Antolín-Díaz and Rubio-Ramírez, 2018) with bayesian method, if time permits
- [ ] testing on examples (SW 2001, A and R 2018, etc.)
- [ ] R package

## A reduced form Bayesian VAR

Let's consider a m-dimension VAR(p) model:

```math
Y_t = \alpha + \sum_{s=1}^{p}Y_{t-s} \beta_{s} + \epsilon_{t} \tag{1}
```

where we assume disturbance terms are iid along time series, $\epsilon_{t} \sim \mathcal{N}(0, \Sigma)$, but $\Sigma$ may not be diagonal. An illustrative example is a bi-variate VAR(2) model containing GDP $y_{t}$ and inflation $\pi_{t}$:

$$
\begin{align*}
    y_{t} &= \alpha_{0} + \alpha_{11} y_{t-1} + \alpha_{12} \pi_{t-1} + \alpha_{21} y_{t-2} + \alpha_{22} \pi_{t-2} + \epsilon_{t}^{y} \\
    \pi_{t} &= \beta_{0} + \beta_{11} y_{t-1} + \beta_{12} \pi_{t-1} + \beta_{21} y_{t-2} + \beta_{22} \pi_{t-2} + \epsilon_{t}^{y} \\
\end{align*}
$$

which can be written as

$$
\begin{bmatrix} 
    y_{t} \\
    \pi_{t} 
\end{bmatrix} = \begin{bmatrix} 
    \alpha_{0} \\
    \beta_{0} 
\end{bmatrix} + \begin{bmatrix} 
    \alpha_{11} & \alpha_{12} \\
    \beta_{21} & \beta_{22}
\end{bmatrix} \begin{bmatrix} 
    y_{t-1} \\
    \pi_{t-1}
\end{bmatrix} + \begin{bmatrix} 
    \alpha_{21} & \alpha_{22} \\
    \beta_{21} & \beta_{22}
\end{bmatrix} \begin{bmatrix} 
    y_{t-2} \\
    \pi_{t-2} 
\end{bmatrix} + \begin{bmatrix} 
    \epsilon_{t}^{y}\\
    \epsilon_{t}^{\pi} 
\end{bmatrix} \tag{2}
$$

The matrix form (2) can be used to understand the compact form (1).

## Estimating VAR via OLS

In model (1), we know innovation at time $t$ won't affect lagged terms of $Y$, so the OLS assumption $\mathbb{E}(\epsilon_{t} x_{t})=0$ and $\mathbb{E}(\epsilon_{t} \epsilon_{t-s})=0$ are naturally satisfied. Thus we can estimate model (1) via OLS method.

Still, we illustrate the problem using the bi-variate VAR(2) example. We can stack all observations vertically and write:

$$
\begin{bmatrix} 
    y_{3} & \pi_{3} \\
    y_{4} & \pi_{4} \\
    \vdots & \vdots \\
    y_{T} & \pi_{T}
\end{bmatrix} = \begin{bmatrix} 
    1 & y_{2} & \pi_{2} & y_{1} & \pi_{1} \\
    1 & y_{3} & \pi_{3} & y_{2} & \pi_{2} \\
    \vdots & \vdots & \vdots & \vdots & \vdots \\
    1 & y_{T-1} & \pi_{T-1} & y_{T-2} & \pi_{T-2} \\
\end{bmatrix} \begin{bmatrix} 
    \alpha_{0} & \beta_{0} \\ 
    \alpha_{11} & \beta_{11} \\ 
    \alpha_{12} & \beta_{12} \\ 
    \alpha_{21} & \beta_{21} \\ 
    \alpha_{22} & \beta_{22} \\ 
\end{bmatrix} + \begin{bmatrix} 
    \epsilon_{3}^{y} & \epsilon_{3}^{\pi} \\ 
    \epsilon_{4}^{y} & \epsilon_{4}^{\pi} \\ 
    \vdots & \vdots \\
    \epsilon_{T}^{y} & \epsilon_{T}^{\pi} \\ 
\end{bmatrix} 
$$

which can be written as

$$
Y = X \beta + \epsilon
$$

then we directly get

$$
\hat{\beta} = (X'X)^{-1}X'Y, \hat{\epsilon}=Y-X \hat{\beta}, \hat{\Sigma}=\hat{\epsilon}' \hat{\epsilon}
$$

this is done in the R code.

## Estimation VAR via Bayesian method

Why do we bother to use bayesian method where OLS would simply solve the problem? Because bayesian method would be of great help if we want to do inference, which I will talk about below. Here I first outline the basic idea of bayesian method (mainly from [this author's github repo](https://github.com/kthohr/BMR) and his book).

Let's write (1) in the following more compact form

$$
y = (\mathbb{I}_{m} \otimes X) \alpha + \varepsilon \tag{3}
$$

where $y=\text{vec}(Y)$, $\alpha=\text{vec}(\beta)$, $\varepsilon=\text{vec}(\epsilon) \sim \mathcal{N}(0, \Sigma \otimes \mathbb{I}_{m})$ are all of dimension $(m \times T) \times 1$ (in fact it should be $(m \times (T-lag)) \times 1$, for convenience I write $T$). Here $\text{vec}$ means to stack the matrix vertically by column and $\otimes$ means the Kronecker product. In R, it can be done by `matlib::vec` and `%x%`. Note $\mathbb{I}_{m}$ is the m-dimension identity matrix.

Again, we illustrate (3) using the VAR(2) example, it says

$$
\begin{align*}
    \begin{bmatrix} 
    y_{3} \\
    \vdots \\
    y_{T} \\
    \pi_{3} \\
    \vdots \\
    \pi_{T} 
\end{bmatrix}&= \left(\begin{bmatrix} 
    1 & 0 \\
    0 & 1 
\end{bmatrix} \otimes X\right) \alpha + \varepsilon \\
&= \begin{bmatrix} 
    X_{5 \times 5} & \mathbf{0} \\
    \mathbf{0} & X_{5 \times 5}
\end{bmatrix} \begin{bmatrix} 
    \alpha_{0} \\
    \vdots \\
    \alpha_{22} \\
    \beta_{0} \\
    \vdots \\
    \beta_{22}
\end{bmatrix} + \begin{bmatrix} 
    \epsilon_{3}^{y} \\
    \vdots \\
    \epsilon_{T}^{y} \\
    \epsilon_{3}^{\pi} \\
    \vdots \\
    \epsilon_{T}^{\pi} \\
\end{bmatrix} 
\end{align*}
$$

The parameters we want to estimate are $\alpha$ and $\Sigma$. Since we assume disturbance terms to be iid normal, with parameters given we know

$$
y \sim \mathcal{N} \left((\mathbb{I}_{m} \otimes X) \alpha, \Sigma \otimes \mathbb{I}_{m}\right)
$$

thus we can write the likelihood function as

$$
\begin{align*}
    \mathcal{L}(\alpha, \Sigma \mid y, X)&=(2 \pi)^{-mT / 2}\left|\Sigma \otimes \mathbb{I}_T\right|^{-1 / 2} \\
    \quad &\times \exp \left\{-\frac{1}{2}\left(y-\left(\mathbb{I}_m \otimes Z\right) \alpha\right)'\left(\Sigma^{-1} \otimes \mathbb{I}_T\right)\left(y-\left(\mathbb{I}_m \otimes Z\right) \alpha\right)\right\} \\
    & \equiv p(y \mid \alpha, \Sigma, X)
\end{align*}
$$

after some manipulation we yield

$$
\begin{align*}
\ln \mathcal{L} &= -\frac{T}{2} \ln |\Sigma|-\left\{\frac{1}{2}(\alpha-\widehat{\alpha})' \mathcal{X}'\left(\Sigma^{-1} \otimes \mathbb{I}_T\right) \mathcal{X}(\alpha-\widehat{\alpha})\right\} \\
& \quad -\left\{\frac{1}{2} \operatorname{tr}\left[\left[(y-\mathcal{X} \widehat{\alpha})(y-\mathcal{X} \widehat{\alpha})'\left(\Sigma^{-1} \otimes \mathbb{I}_T\right)\right]\right\}\right. \\
& \propto \ln \left[ \mathcal{N}(\alpha \mid \widehat{\alpha},\Sigma,\mathcal{X},y) \times \mathcal{IW}(\Sigma \mid \widehat{\alpha},\mathcal{X},y)\right]
\end{align*}
$$

where 

$$
\begin{align*}
    \widehat{\alpha} &\equiv (\Sigma^{-1} \otimes X'X)^{-1}(\Sigma^{-1} \otimes X)'y \\
    \mathcal{X} &\equiv \mathbb{I}_{m} \otimes X
\end{align*}
$$

which means the likelihood function is (conditionally) proportional to the product of **a normal distribution** for $\alpha$ and **a inverse-Wishart distribution** for $\Sigma$. It seems that we can sample the parameters independently from the posterior density using **Gibbs Sampler**, if we use appropriate priors. 

---
By the way, I find a code snippet to generate random draws from inverse-Wishart distribution from [this package](https://rdrr.io/cran/nicheROVER/src/R/rwish.R), which should be easily transformed to `Rcpp` code.

---

In fact it is the case if we use the Normal-inverse-Wishart Prior, which specifies the priors of $\alpha$ and $\Sigma$ to be **independently normal and inverse-Wishart**. The result is as follow:

First, use bayesian formula we know the joint posterior distribution satisfies

$$
p(\alpha, \Sigma \mid \mathcal{X}, y) \propto p(y \mid \alpha, \Sigma, X) p(\alpha,\Sigma) =  p(y \mid \alpha, \Sigma, X) p(\alpha) p(\Sigma)
$$

The priors are $p(\alpha) \sim \mathcal{N}(\bar{\alpha},\mathcal{E}_{\alpha})$ and $\Sigma \sim \mathcal{IW}(\mathcal{E}_{\Sigma},\nu)$. 

If given $\Sigma$, after some computation (still follow [this author's github repo](https://github.com/kthohr/BMR) or [this book](https://www.cambridge.org/core/books/structural-vector-autoregressive-analysis/DAF4217439EA585D10902D58A8849E06)'s chapter 5) we know the prior density of $\alpha$ satisfies

$$
p(\alpha \mid \Sigma, \mathcal{X}, y) \propto \exp \left(-\frac{1}{2}\left[(\alpha-\tilde{\alpha})' \tilde{\Sigma}_{\alpha}^{-1}(\alpha-\tilde{\alpha})+C \right]\right)=\mathcal{N}(\tilde{\alpha},\tilde{\Sigma}_{\alpha}) \tag{4}
$$

where

$$
\begin{align*}
    \tilde{\Sigma}_{\alpha}^{-1}&=\mathcal{E}_{\alpha}^{-1}+\mathcal{X}'(\Sigma^{-1} \otimes \mathbb{I}_{T}) \mathcal{X} \\
    \tilde{\alpha}&=\tilde{\Sigma}_{\alpha}(\mathcal{E}_{\alpha}^{-1}\bar{\alpha}+\mathcal{X}'(\Sigma^{-1}\otimes \mathbb{I}_{T})y) \\
    C&=y'(\Sigma^{-1}\otimes \mathbb{I}_{T})y + \bar{\alpha}'\mathcal{E}_{\alpha}^{-1} \bar{\alpha} - \tilde{\alpha}'\tilde{\Sigma}_{\alpha}^{-1}\tilde{\alpha}
\end{align*}
$$

If given $\alpha$, we have

$$
p(\Sigma \mid \beta, X, Y) = \mathcal{IW} \left(\mathcal{E}_{\Sigma}+(Y-X\beta)'(Y-X\beta),T+\nu \right) \tag{5}
$$

where $\alpha=\text{vec}(\beta)$.

On (4) and (5) we can apply the Gibbs sampler to estimate the parameters we are interested in. A possible selection of prior parameters is:

- $\bar{\alpha}$: set the prior mean of the first lag of the dependent variable to one and set the prior mean of all other slope coefficients to zero, i.e. mimic a multivariate random walk model
- $\mathcal{E}_{\alpha}$: set to $\eta \mathbb{I}_{m(mp+1)}$, where $\eta$ is a hyper parameter
- $\mathcal{E}_{\Sigma}$: set to $\mathbb{I}_{m}$
- $\nu$: set to $m+1$

## Identification and Inference

Note that in (1), the random shocks are serially uncorrelated, but different shocks at the same period may be correlated. In macroeconomic analysis, people are mainly interested in the response of economic variables on **structural shocks**, i.e. *exogenous and mutually independent shocks with explicit economic meaning*. To achieve this, we want to find a matrix $\mathcal{B}$ s.t. it transforms structural shocks $\xi_{t}$ to reduced-form shocks $\epsilon_{t}$

$$
\mathcal{B} \xi_{t} = \epsilon_{t}, \xi_{t} \sim \mathcal{N}(0,D)
$$

where $D$ is a diagonal matrix. Let $B=\mathcal{B}^{-1}$, we can write

$$
BY_t = B\alpha + \sum_{s=1}^{p}B\beta_{s} Y_{t-s} + \xi_{t}
$$

which is 

$$
Y_{t} = B\alpha + (I-B)Y_{t} + \sum_{s=1}^{p} B\beta_{s}Y_{t-s} + \xi_{t} \tag{6}
$$

Now it becomes inappropriate to estimate (6) using OLS since $\mathbb{E}(x_{t} \xi_{t}) \neq 0$, and OLS estimator is inconsistent.

To better understand the problem, let's go back to the reduced form VAR where structural shocks are specified

$$
Y_t = \alpha + \sum_{s=1}^{p} \beta_{s} Y_{t-s} + \mathcal{B} \xi_{t}
$$

We can directly using OLS to estimate this equation, get $\hat{\beta}_{s}$ and $\hat{\epsilon}_{t}$. To estimate $\mathcal{B}$, we utilize

$$
\text{Var}(\hat{\epsilon}_{t}) \equiv \Sigma = \mathcal{B} \mathcal{B}'
$$

where we normalize $D=\mathbb{I}_{m}$ to work with 1-unit-standard-error shock. Note we have $m\times m$ entries in $\mathcal{B}$ to estimate while the equality just give us $m(m-1)/2$ equations, this $\mathcal{B}$ is not identified, i.e. we lack enough assumptions to identify $\mathcal{B}$ matrix.

A possible solution is setting $\mathcal{B}$ to be lower triangle, which requires $b_{ij}=0$ if $i<j$. This assumption is widely used in old times because it's easy to implemented and has straightforward economic meaning: variable 1 won't affected by shock 2 to last in the short run, variable 2 won't affected by shock 3 to last in the shprt run, and so on.
But recently people find this assumption might be too strong and may not be realistic.

Another solution is using IV, which heavily relies on accessible data but easy to implement. It not the focus of this project.

We focus on identification by **sign restrictions**, which will be illustrated later.

After identification of $\mathcal{B}$, we can analyze the quantity we want, i.e. **the response of economic variables on structural shocks**, which is known as structural impulse-response functions (IRFs). The structural IRF of variable $i$ to shock $j$ after $h$ periods the shock took place is defined as

$$
\Phi_{ij}^{h} \equiv \frac{\partial y_{it+h}}{\partial \xi_{jt}} = \frac{\partial y_{it+h}}{\partial \epsilon_{t}} \frac{\partial \epsilon_{t}}{\partial \xi_{jt}}=\Psi_{i \cdot} \mathcal{B}_{\cdot j}
$$

where $\Psi$ is the reduced-form IRF which can be easily computed after $\hat{\beta}$ is estimated (and has been implemented in the R code).

To do inference, with frequentist methods we bootstrap on $\epsilon_{t}$ to reconstruct data, re-estimate model and re-compute IRF for, typically, 2000 times and use 5% and 95% quantiles of IRF computed as credible interval, with bayesian methods the sampling procedure naturally gives the credible interval of IRF.

## Sign Restrictions

TODO

### Basics

### Narrative Sign Restrictions
