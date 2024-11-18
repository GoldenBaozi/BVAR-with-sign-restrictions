library(R6)
library(expm)
library(MASS)
# TODO add bootstrap of IRF and historical decomposition
VAR <- R6Class(
  "VAR",
  public = list(
    T = NULL,
    T.est = NULL,
    n.var = NULL,
    p.lag = NULL,
    hor = NULL,
    data = NULL,
    var.names = NULL,
    time = NULL,
    Y = NULL,
    Y.start = NULL,
    X = NULL,
    beta = NULL,
    U = NULL,
    Sigma = NULL,
    B = NULL,
    eps = NULL,
    Psi = NULL,
    IRF = NULL,
    FEVD = NULL,
    initialize = function(data = NA, p.lag = NA) {
      # TODO check data is a matrix, p.lag is a positive integer
      self$T <- dim(data)[1]
      self$n.var <- dim(data)[2] - 1
      self$time <- data[, 1]
      self$data <- as.matrix(data[, 2:(self$n.var + 1)])
      self$p.lag <- p.lag
      self$T.est <- self$T - p.lag
      self$hor <- self$T
      self$var.names <- colnames(self$data)
      self$Y <- self$data[(1 + p.lag):self$T, ]
      Y.lag <- self$data[p.lag:(self$T - 1), ]
      for (i in 2:p.lag) {
        Y.lag <- cbind(Y.lag, self$data[(1 + p.lag - i):(self$T - i), ])
      }
      self$X <- cbind(rep(1, self$T.est), Y.lag)
      self$Y.start <- self$data[1:(self$T - p.lag), ]
      cat("* VAR class initialized.\n")
      cat("-> variables: ", self$var.names, "\n")
      cat("-> time period: ", self$time[1], " to ", self$time[self$T], "\n")
    },
    fit = function(method = "OLS") {
      if (method == "OLS") {
        # TODO check the OLS problem is appropriate
        self$beta <- private$est.OLS()
        self$U <- self$Y - self$X %*% self$beta
        self$Sigma <- t(self$U) %*% self$U / (self$T.est - self$n.var - 1)
      }
      cat("* parameters of VAR estimated using ", method, ", residual and Cov matrix yield.\n")
    },
    identify = function(method = "recursive", IV = NA, sign = NA) {
      # TODO add IV and sign identification methods
      if (method == "recursive") {
        self$B <- t(chol(self$Sigma))
        B.invT <- t(solve(self$B))
        self$eps <- self$U %*% B.invT # use B matrix and reduced-form shock to compute structural shock
      }
      cat("* model identified using ", method, " approach.\n")
      cat("==============================================\n")
      cat("* computing structural VAR tools ...\n")
      IRF.res <- private$IRF.compute()
      self$Psi <- IRF.res$Psi
      self$IRF <- IRF.res$IRF
      cat("-> IRF computed.\n")
      self$FEVD <- private$FEVD.compute()
      cat("-> FEVD computed.")
    },
    IRF.plot = function(shock = NA_character_, response = NA_character_, hor = self$hor) { # by default plot the whole horizon
      shock.id <- which(self$var.names == shock)
      res.id <- which(self$var.names == response)
      violation <- (identical(shock.id, integer(0)) | identical(res.id, integer(0)))
      if (violation) {
        stop("Please give appropriate variable names")
      }
      irf.to.plot <- self$IRF[res.id, shock.id, 1:(hor + 1)]
      x <- 0:hor
      if (hor <= self$hor) {
        plot(x, irf.to.plot,
          type = "l",
          lty = 1, lwd = 2, col = 4,
          xlab = "horizon", ylab = response
        )
        title(main = paste("IRF of", response, "to", shock, "shock"))
        abline(h = 0, lwd = 3, lty = 2)
      } else {
        cat("horizon > ", self$hor, " not allowed")
      }
    },
    HDC.compute = function(end1 = NA_character_, end2 = NA_character_, shock = NA_character_, res = NA_character_) {
      end1.id <- which(self$time[(self$p.lag + 1):self$T] == end1)
      end2.id <- which(self$time[(self$p.lag + 1):self$T] == end2)
      violation <- (identical(end1.id, integer(0)) | identical(end2.id, integer(0)) | end1.id > end2.id)
      if (violation) {
        stop("Please set appropriate time periods")
      } else if (end2.id > end1.id) {
        HDC.ij <- vector("numeric", (end2.id - end1.id))
        for (t in (end1.id + 1):end2.id) {
          HDC.ij[(t - end1.id)] <- private$HDC.core(t, shock, res) - private$HDC.core(end1.id, shock, res)
        }
        return(HDC.ij)
      } else {
        shock.id <- which(self$var.names == shock)
        res.id <- which(self$var.names == res)
        eps.used <- self$eps[end1.id, shock.id]
        Phi.used <- self$IRF[res.id, shock.id, 1]
        HDC.ij <- Phi.used * eps.used
        return(HDC.ij)
      }
    }
  ),
  private = list(
    # FIXME how to compute bootstrap of IRF ? every bootstrap need identification ?
    est.OLS = function(Y = self$Y, X = self$X) {
      beta <- solve(t(X) %*% X, t(X) %*% Y)
      return(beta)
    },
    IRF.compute = function(beta = self$beta, B = self$B, hor = self$hor, boot = FALSE) {
      nvar <- self$n.var
      plag <- self$p.lag
      beta.1 <- t(beta)
      beta.lag <- beta.1[, 2:(nvar * plag + 1)]
      beta.compact <- rbind(beta.lag, cbind(diag(nvar * (plag - 1)), matrix(0, (plag - 1) * nvar, nvar)))
      Psi <- array(0, c(nvar, nvar, hor + 1))
      irf_trans <- cbind(diag(nvar), matrix(0, nvar, nvar * (plag - 1)))
      if (boot == FALSE) {
        IRF <- array(0, c(nvar, nvar, hor + 1))
        for (h in 1:(hor + 1)) {
          Psi[, , h] <- irf_trans %*% (beta.compact %^% (h - 1)) %*% t(irf_trans)
          IRF[, , h] <- Psi[, , h] %*% B
        }
        return(
          list(
            "Psi" = Psi,
            "IRF" = IRF
          )
        )
      } else {
        for (h in 1:(hor + 1)) {
          Psi[, , h] <- irf_trans %*% (beta.compact %^% (h - 1)) %*% t(irf_trans)
        }
        return(Psi)
      }
    },
    FEVD.compute = function(Sigma = self$Sigma, B = self$B, Psi = self$Psi) {
      # initialize variables
      nvar <- self$n.var
      nstep <- self$hor + 1
      MSE <- array(0, c(nvar, nvar, nstep))
      MSE_shock <- array(0, c(nvar, nvar, nstep))
      FEVD <- array(0, c(nvar, nvar, nstep))
      # calculate FEVD
      for (mm in 1:nvar) {
        MSE[, , 1] <- Sigma
        MSE_shock[, , 1] <- B[, mm] %*% t(B[, mm])
        for (kk in 2:nstep) {
          MSE[, , kk] <- MSE[, , kk - 1] + Psi[, , kk] %*% Sigma %*% t(Psi[, , kk])
          MSE_shock[, , kk] <- MSE_shock[, , kk - 1] + Psi[, , kk] %*% MSE_shock[, , 1] %*% t(Psi[, , kk])
          FEVD[, mm, kk] <- diag(MSE_shock[, , kk]) / diag(MSE[, , kk])
        }
      }
      return(FEVD)
    },
    HDC.core = function(end = NA_integer_, shock = NA_character_, res = NA_character_) {
      # ensure 'end' is an id
      shock.id <- which(self$var.names == shock)
      res.id <- which(self$var.names == res)
      violation <- (identical(shock.id, integer(0)) | identical(res.id, integer(0)))
      if (violation) {
        stop("Please give appropriate variable names")
      } else {
        eps.used <- rev(self$eps[1:end, shock.id])
        Phi.used <- self$IRF[res.id, shock.id, 1:end]
        HDC.ijt <- Phi.used %*% eps.used
        return(HDC.ijt)
      }
    }
  )
)
# TODO write a gibbs sampler for VAR bayesian estimation
Gsampler <- R6Class(
  "Gibbs Sampler",
  public = list(

  ),
  private = list(

  )
)
