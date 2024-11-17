plag <- 4
hor <- 24
sample_start <- "1960q1"
sample_end <- "2000q4"
var_name <- c("date", "infl", "unemp", "ffr")

# import data
data <- read.csv("../data//SW2001_Data.csv")
start_index <- which(data$date == sample_start)
end_index <- which(data$date == sample_end)
used <- data[start_index:end_index, var_name]

source("./code/VARclass.R")
VAR.ols <- VAR$new(data = used, p.lag = plag, hor = hor)
VAR.ols$fit()
VAR.ols$identify()

png("IRF.png", width = 16, height = 9, units = "in", pointsize = 14, res = 300)
par(mfrow = c(3, 1))
for (i in var_name[2:4]) {
  VAR.ols$IRF.plot(shock = i, response = "ffr")
}
dev.off()
