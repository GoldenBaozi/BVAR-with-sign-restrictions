library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
library(R.matlab)
library(vars)

sourceCpp("./src/VARtools.cpp")
sourceCpp("./src/BVARdraw.cpp")
source("./scripts/VARclass.R")

df <- readMat("./NS2018/data/data/Kilian_Data_Updated.mat")
varnames <- c(df$varNames[1, 1][[1]][[1]][1, 1], df$varNames[1, 2][[1]][[1]][1, 1], df$varNames[1, 3][[1]][[1]][1, 1])

dates <- as.Date(df$dates - 719529)

df.in <- data.frame(
  Date = dates,
  df$data
)
colnames(df.in)[2:4] <- varnames


plag <- 24

bvar <- bVAR$new(df.in, plag)
bvar$Sigma
# test conjugate method
t(chol(bvar$Sigma))
bvar$est(method = "conjugate")

# check S post is same as AR 2018
bvar$Sigma.post$mean
bvar$Sigma.post$nu
tail(bvar$alpha.post$mean) # alpha post is also the same
sourceCpp("./src/testwish.cpp") # check inv wish posterior

test_invwish(bvar$Sigma.post$mean, bvar$Sigma.post$nu)

SR <- list(
  list("Oil Production Growth", c("Oil Production Growth", "Economic Activity Index"), 0, -1),
  list("Oil Production Growth", "Real Oil Price", 0, 1),
  list("Economic Activity Index", c("Oil Production Growth", "Economic Activity Index", "Real Oil Price"), 0, 1),
  list("Real Oil Price", c("Oil Production Growth", "Real Oil Price"), 0, 1),
  list("Real Oil Price", "Economic Activity Index", 0, -1)
)

EBR <- list(
  list("Economic Activity Index", c("Oil Production Growth", "Real Oil Price"), 0, 0.05, NA),
  list("Real Oil Price", c("Oil Production Growth", "Real Oil Price"), 0, 0.05, NA)
)

NSR <- list(
  list("Economic Activity Index", "contribution", as.Date("1990-08-01"), as.Date("1990-08-01"), "Real Oil Price", -1, "strong")
)

which(bvar$time == NSR[[1]][[3]]) # correct time input
which(bvar$var.names == NSR[[1]][[1]]) # correct var name input
str(bvar$.__enclos_env__$private$get.restrictions(SR, EBR, NSR))
bvar$identify(SR = SR, EBR = EBR, draw = 100, save = 100, M = 10)

# plot(bvar$IRF.avg[7, ], type = "l", lty = 1, lwd = 2, col = 4, ylim = c(-5, 10))
# lines(bvar$IRF.lb[7, ], lty = 2, col = 2)
# lines(bvar$IRF.ub[7, ], lty = 2, col = 2)
# bvar$IRF.avg

# rm(list = c("bvar"))
# gc()

# test progress bar
sourceCpp("./src/test.cpp")
save.image("mydata.RData")