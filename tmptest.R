data(Canada)
head(Canada)
var.1 <- vars::VAR(Canada, p = 2, type = "const")
lmtest()
var.1$varresult
resid(var.1)
cov(resid(var.1))

var.2 <- vars::VAR(df.in[, 2:4], p = 24, type = "const")
cov(resid(var.2))