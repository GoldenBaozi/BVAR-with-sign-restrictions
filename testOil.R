library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
library(R.matlab)
library(MASS)
# library(vars)

load("./mydata.RData")
sourceCpp("./src/VARtools.cpp")
sourceCpp("./src/BVARdraw.cpp")
source("./scripts/VARclass.R")
sourceCpp("./src/getroot.cpp")

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
bvar$est(method = "conjugate")

# check S post is same as AR 2018
bvar$Sigma.post$mean
bvar$Sigma.post$nu
tail(bvar$alpha.post$mean) # alpha post is also the same
sourceCpp("./src/testwish.cpp") # check inv wish posterior

Sigma_draw <- test_invwish(bvar$Sigma.post$mean, bvar$Sigma.post$nu)
Sigma.1 <- Sigma_draw %x% solve(t(bvar$X) %*% bvar$X)
cov <- t(chol(Sigma.1))
alpha_draw <- bvar$alpha.post$mean + cov %*% rnorm(length(bvar$alpha.post$mean))
beta_draw <- matrix(alpha_draw, ncol = 3)
max_root(beta_draw)


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

bvar$identify(SR = SR, EBR = EBR, draw = 2000, M = 1000)
IRF.1 <- bvar$IRF.compute(18, 0.68)


var.names <- c("Oil production", "Economic Activity Index", "Real Oil Price")
shock.names <- c("Oil Supply", "Aggregate Demand", "Oil-specific Demand")
xaxis <- 0:18
png("./out/IRF_oil_no_NSR.png", width=12, height=9, units="in", res=300)
par(mfrow = c(3, 3))
par(family="serif", cex.main = 1.5, cex.lab = 1.2, cex.axis=1.2)
for (i in 0:8) {
  var <- i %% 3 + 1
  shk <- i %/% 3 + 1
  var.name <- var.names[var]
  shk.name  <- shock.names[shk]
  idx <- (i %% 3)*3 + (i %/% 3 + 1)
  median <- IRF.1$avg[idx,]
  ub <- IRF.1$ub[idx,]
  lb <- IRF.1$lb[idx,]

  plot(xaxis, median,
    type = "l", col = "blue", ylim = range(c(median, lb, ub,-2,1)),
    xlab = "Months", ylab = "Percent", main = paste(var.name,"to",shk.name,"shock")
  )

  polygon(c(xaxis, rev(xaxis)), c(ub, rev(lb)), col = "grey", border = NA)
  lines(xaxis, median, col = "blue", lty = 1, lwd = 3)
  abline(h = 0, lwd = 3, lty = 2)
}
dev.off()
save.image("mydata.RData")
