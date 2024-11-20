library(R6)
library(expm)
# TODO add bootstrap of IRF
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
    identify = function(method = "recursive", IV = NA) {
      # TODO add IV and sign identification methods
      if (method == "recursive") {
        self$B <- t(chol(self$Sigma))
        B.invT <- t(solve(self$B))
        self$eps <- self$U %*% B.invT # use B matrix and reduced-form shock to compute structural shock
      }
      cat("* model identified using ", method, " approach.\n")
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
    },
    tools = function(boot = FALSE) {
      cat("* computing structural VAR tools ...\n")
      IRF.res <- self$IRF.compute(boot = boot)
      self$Psi <- IRF.res$Psi
      self$IRF <- IRF.res$IRF
      cat("-> IRF computed.\n")
      self$FEVD <- self$FEVD.compute()
      cat("-> FEVD computed.")
    },
    IRF.plot = function(shock = NA_character_, response = NA_character_, hor = self$hor, CI = FALSE) { # by default plot the whole horizon
      shock.id <- which(self$var.names == shock)
      res.id <- which(self$var.names == response)
      violation <- (identical(shock.id, integer(0)) | identical(res.id, integer(0)))
      if (violation) {
        stop("Please give appropriate variable names")
      }
      irf.to.plot <- self$IRF[res.id, shock.id, 1:(hor + 1)]
      x <- 0:hor
      if (hor > self$hor) {
        msg <- paste("horizon >", self$hor, "not allowed")
        stop(msg)
      }
      plot(x, irf.to.plot,
        type = "l",
        lty = 1, lwd = 2, col = 4,
        xlab = "horizon", ylab = response
      )
      # TODO add CI option
      title(main = paste("IRF of", response, "to", shock, "shock"))
      abline(h = 0, lwd = 3, lty = 2)
    }
  ),
  private = list(
    est.OLS = function(Y = self$Y, X = self$X) {
      beta <- solve(t(X) %*% X, t(X) %*% Y)
      return(beta)
    }
    # TODO bootstraps
  )
)
# TODO write a gibbs sampler for VAR bayesian estimation

bVAR <- R6Class(
  "Bayes VAR",
  inherit = VAR,
  public = list(
    # prior and post parameters
    prior.type = NULL,
    alpha.prior = NULL,
    Sigma.prior = NULL,
    alpha.post = NULL,
    Sigma.post = NULL,
    # parameter draws
    alpha.draw = NULL, ## vector of beta
    beta.draw = NULL, ## matrix of alpha
    Sigma.draw = NULL, ## cov of Y-X \times beta, used for cholesky decomposition
    # related non-orthogonal/orthogonal shocks and IRFs
    u.draw = NULL, ## Y - X \times \beta
    eps.draw = NULL, ## B^{-1} u
    Psi.draw = NULL,
    IRF.draw = NULL,
    initialize = function(data = NA, p.lag = NA, prior.type = "independent", priors = NA) {
      # need to check data and p.lag input
      super$initialize(data, p.lag)
      self$prior.type <- prior.type
      if (prior.type == "independent") {
        if (is.na(priors)) {
          alpha.prior <- list(
            "mean" = as.vector(t(cbind(rep(0, self$n.var), diag(n.var), matrix(0, self$n.var, self$n.var * (self$p.lag - 1))))),
            "cov" = 0.1 * diag(self$n.var * (self$n.var * self$p.lag + 1))
          )
          Sigma.prior <- list(
            "mean" = diag(self$n.var),
            "nu" = self$n.var + 1
          )
        } else {
          alpha.prior <- list(
            "mean" = priors[[1]],
            "cov" = priors[[2]]
          )
          Sigma.prior <- list(
            "mean" = priors[[3]],
            "nu" = priors[[4]]
          )
        }
      }
      # TODO other priors, e.g. conjugate
    },
    est = function(Y = self$Y, X = self$X, prior.type = self$prior.type, burn.in = 100, draw = 1000, thin = 1, post.save = NA, post.method = "last") {
      if (prior.type == "independent") {
        y <- as.vector(Y)
        priors <- c(self$alpha.prior, self$Sigma.prior)
        cat("Gibbs sampler start...\n")
        res <- gsampler(y, Y, X, priors, burn.in = burn.in, draw = draw, thin, usage = "sign", post.save = post.save, post.method = post.method)
        ## save posterior parameters
        self$alpha.post <- list("mean" = out$alpha.mean.post, "cov" = out$alpha.cov.post)
        self$Sigma.post <- list("mean" = out$Sigma.mean.post, "nu" = out$nu.post)
        ## temporally save parameters estimated, if not using sign restrictions
        self$beta <- matrix(rowMeans(out$alpha), self$n.var * self$p.lag + 1, self$n.var)
        self$Sigma <- apply(out$Sigma, c(1, 2), mean)
        cat("posterior parameters estimated.") # TODO add a supervisor, and other priors
        ## return loglik for analysis
        return(out$loglik)
      }
    }
  ),
  private = list(
    get.restrictions = function(SR = NA, EBR = NA, NSR = NA) {
      restrictions <- list()
      if (!is.na(SR)) {
        for (i in 1:length(SR)) {
          len.var <- length(SR[[i]][[2]])
          shock <- SR[[i]][[1]]
          h <- SR[[i]][[3]]
          sign <- SR[[i]][[4]]
          for (j in 1:len.var) {
            var <- SR[[i]][[2]][[j]]
            restrictions <- c(restrictions, list(
              "type" = "SR",
              "shock" = which(self$var.names == shock),
              "var" = which(self$var.names == var),
              "h" = h,
              "sign" = sign
            ))
          }
        }
      }
      if (!is.na(EBR)) {
        for (i in 1:length(EBR)) {
          shock <- SR[[i]][[1]]
          h <- SR[[i]][[3]]
          mb <- SR[[i]][[4]]
          lb <- SR[[i]][[5]]
          if (is.na(mb)) mb <- Inf
          if (is.na(lb)) lb <- -Inf
          var.1 <- SR[[i]][[2]][[1]]
          var.2 <- SR[[i]][[2]][[2]]
          restrictions <- c(restrictions, list(
            "type" = "EBR",
            "shock" = which(self$var.names == shock),
            "var.1" = which(self$var.names == var.1),
            "var.2" = which(self$var.names == var.2),
            "h" = h,
            "max.bound" = mb,
            "low.bound" = lb
          ))
        }
      }
      if (!is.na(NSR)) {
        for (i in 1:length(NSR)) {
          shock <- NSR[[i]][[1]]
          type <- NSR[[i]][[2]]
          start <- which(self$time == NSR[[i]][[3]]) - self$p.lag
          end <- which(self$time == NSR[[i]][[4]]) - self$p.lag
          if (length(NSR[[i]]) == 5) {
            sign <- NSR[[i]][[5]]
            restrictions <- c(restrictions, list(
              "type" = "NSR",
              "NSR.type" = "sign",
              "shock" = which(self$var.names == shock),
              "period" = start:end,
              "sign" = sign,
            ))
          } else {
            var <- NSR[[i]][[5]]
            sign <- NSR[[i]][[6]]
            intensity <- NSR[[i]][[7]]
            restrictions <- c(restrictions, list(
              "type" = "NSR",
              "NSR.type" = "contribution",
              "shock" = which(self$var.names == shock),
              "var" = which(self$var.names == var),
              "period" = start:end,
              "sign" = sign,
              "intensity" = intensity
            ))
          }
        }
      }
    }
    # draw.1 = function() {

    # }
  )
)
