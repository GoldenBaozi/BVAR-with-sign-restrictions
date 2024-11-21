library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
library(R.matlab)

sourceCpp("./src/VARtools.cpp")
sourceCpp('./src/BVARdraw.cpp')
source("./scripts/VARclass.R")

df <- readMat("./NS2018/data/data/Kilian_Data_Updated.mat")
df$varNames[1,1][[1]][[1]][1,1]
varnames <- c(df$varNames[1,1][[1]][[1]][1,1], df$varNames[1,2][[1]][[1]][1,1], df$varNames[1,3][[1]][[1]][1,1])
dates <- as.Date(df$dates - 719529)
str(dates)
str(df$data)
df.in <- data.frame(
  Date = dates,
  df$data
)
colnames(df.in)[2:4] <- varnames
head(df.in)

plag <- 24
bvar <- bVAR$new(df.in, plag)
loglik <- bvar$est()
