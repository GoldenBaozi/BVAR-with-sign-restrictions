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

source("./scripts/VARclass.R")
VAR.ols <- VAR$new(data = used, p.lag = plag)
VAR.ols$fit() # fit model using OLS
VAR.ols$identify() # identify B matrix using recursive approach (cholesky decomposition)
t(VAR.ols$eps) %*% VAR.ols$eps # verify that epsilon is mutually independent shock

png("./out/IRF.png", width = 16, height = 9, units = "in", pointsize = 16, res = 300)
par(mfrow = c(3, 1))
for (i in var_name[2:4]) {
  VAR.ols$IRF.plot(shock = "ffr", response = i, hor = 24)
}
dev.off()

HDC <- VAR.ols$HDC.compute('1970q1','1980q4', 'infl', 'unemp')
plot(HDC, type='l')
