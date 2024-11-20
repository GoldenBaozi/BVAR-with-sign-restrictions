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
VAR.ols$tools()

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
