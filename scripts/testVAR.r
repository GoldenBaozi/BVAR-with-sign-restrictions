library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)

plag <- 4
hor <- 24
sample_start <- "1960q1"
sample_end <- "2000q4"
var_name <- c("date", "infl", "unemp", "ffr")

# import data
data <- read.csv("./data/SW2001_Data.csv")
start_index <- which(data$date == sample_start)
end_index <- which(data$date == sample_end)
used <- data[start_index:end_index, var_name]

# test basic VAR
source("./scripts/VARclass.R")
VAR.ols <- VAR$new(data = used, p.lag = plag)
VAR.ols$fit() # fit model using OLS
VAR.ols$identify() # identify B matrix using recursive approach (cholesky decomposition)

# test VARtools in Cpp
## test IRF
sourceCpp("./src/VARtools.cpp")
VAR.ols$Sigma
IRF.out <- IRF_compute(VAR.ols$beta, VAR.ols$B, 24, VAR.ols$n.var, VAR.ols$p.lag)

VAR.ols$tools()
VAR.ols$IRF[,,25]
IRF.out$IRF[,,25] # note: the index will be automatically modified from Cpp to R

## test HDC
start <- which(VAR.ols$time == "1970q1") - VAR.ols$p.lag
end <- which(VAR.ols$time == "1980q4") - VAR.ols$p.lag
shock <- which(VAR.ols$var.names == "infl")
res <- which(VAR.ols$var.names == "ffr")
HDC.out <- HDC_ts(start, end, shock, res, VAR.ols$IRF, VAR.ols$eps)

## test FEVD
FEVD.out <- FEVD_compute(VAR.ols$Sigma, VAR.ols$B, VAR.ols$Psi, VAR.ols$n.var, 24)
FEVD.out
VAR.ols$FEVD[,,25]

t(VAR.ols$eps) %*% VAR.ols$eps # verify that epsilon is mutually independent shock

# test IRF plot
png("./out/IRF.png", width = 16, height = 9, units = "in", pointsize = 16, res = 300)
par(mfrow = c(3, 1))
for (i in var_name[2:4]) {
  VAR.ols$IRF.plot(shock = "ffr", response = i, hor = 24)
}
dev.off()

# test HD
HDC <- VAR.ols$HDC.compute("1970q1", "1980q4", "infl", "unemp")
plot(HDC, type = "l")

# test bayesian estimation

source("./scripts/BVAR.R")
y <- as.vector(VAR.ols$Y)
alpha.bar <- c(
  c(0,1,0,0,rep(0,9)),
  c(0,0,1,0,rep(0,9)),
  c(0,0,0,1,rep(0,9))
  )
Sigma.alpha <- 0.1 * diag(39)
Sigma.sigma <- diag(3)
nu <- 4
priors <- list(alpha.bar, Sigma.alpha, Sigma.sigma, nu)
bvar.out <- gsampler(VAR.ols$Y, VAR.ols$X, priors)
bvar.out$alpha.mean.post
bvar.out$Sigma.mean.post

alpha.fin <- rowMeans(bvar.out$alpha)
beta <- matrix(alpha.fin, ncol = 3)
beta - VAR.ols$beta # note that the beta result of OLS and bayesian are somehow different

# test Rcpp
# now, use OLS result as prior
sourceCpp('./src/BVARdraw.cpp')
# alpha.bar <- as.vector(VAR.ols$beta)
alpha.bar <- c(
  c(0,1,0,0,rep(0,9)),
  c(0,0,1,0,rep(0,9)),
  c(0,0,0,1,rep(0,9))
  )
Sigma.alpha <- 0.1 * diag(39)
# Sigma.sigma <- VAR.ols$Sigma
Sigma.sigma <- diag(3)
nu <- 4
priors <- list(alpha.bar, Sigma.alpha, Sigma.sigma, nu)
bvar.out.c <- gibbs_sampler(VAR.ols$Y, VAR.ols$X, priors, post_method = 1)

# test bayesian posterior beta
matrix(bvar.out.c$alpha_mean_post, 13,3) - VAR.ols$beta
VAR.ols$beta

# test bayesian posterior Sigma
bvar.out.c$Sigma_mean_post
VAR.ols$Sigma
rwish(1, bvar.out.c$Sigma_mean_post, bvar.out.c$nu_post, TRUE)

# test time
microbenchmark("gibbs R" = gsampler(VAR.ols$Y, VAR.ols$X, priors), "gibbs Cpp" = gibbs_sampler(VAR.ols$Y, VAR.ols$X, priors), times = 10)
