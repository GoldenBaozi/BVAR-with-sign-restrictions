library(R6)
library(expm)
library(Rcpp)
library(RcppArmadillo)
library(tictoc)
# TODO add bootstrap of IRF
VAR <- R6Class(
  "basic VAR",
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
      self$var.names <- colnames(self$data)
      self$Y <- self$data[(1 + p.lag):self$T, ]
      Y.lag <- self$data[p.lag:(self$T - 1), ]
      for (i in 2:p.lag) {
        Y.lag <- cbind(Y.lag, self$data[(1 + p.lag - i):(self$T - i), ])
      }
      self$X <- cbind(rep(1, self$T.est), Y.lag)
      self$Y.start <- self$data[1:(self$T - p.lag), ]
      cat("* VAR class initialized.\n")
      cat("-> variables: ", paste0(self$var.names, collapse = ", "), "\n")
      cat("-> time period: ", as.character(self$time[1]), " to ", as.character(self$time[self$T]), "\n")
      cat("-> model: ", as.character(self$p.lag), " lags of all variables\n")
      self$beta <- solve(t(self$X) %*% self$X, t(self$X) %*% self$Y)
      self$U <- self$Y - self$X %*% self$beta
      self$Sigma <- (t(self$U) %*% self$U) / (self$T.est)
      cat("* reduced form VAR parameters estimated using OLS.\n")
    },
    identify = function(method = "recursive", IV = NA) {
      # TODO add IV identification
      if (method == "recursive") {
        self$B <- t(chol(self$Sigma)) # to get lower triangular factor
        B.Tinv <- solve(t(self$B)) # Note here eps is actually eps' and U is U', and eps'=u'B'^{-1}
        self$eps <- self$U %*% B.Tinv # use B matrix and reduced-form shock to compute structural shock
      }
      cat("* model identified using ", method, " approach.\n")
    },
    IRF.compute = function(hor = NA_integer_) {
      #' @return IRF.list$Psi, IRF.list$IRF, 3-D array
      IRF.out <- IRF_compute(self$beta, self$B, hor, self$n.var, self$p.lag)
      self$hor <- hor
      self$Psi <- IRF.out$Psi
      self$IRF <- IRF.out$IRF
    },
    FEVD.compute = function(hor = self$hor) {
      #' @return FEVD 3-D array
      self$FEVD <- FEVD_compute(self$Sigma, self$B, self$Psi, self$n.var, self$hor)
    },
    HDC.compute = function(start = NA_character_, end = NA_character_, shock = NA_character_, res = NA_character_) {
      #' @return a matrix of (end-start) x 1
      shock.id <- which(self$var.names == shock)
      res.id <- which(self$var.names == res)
      start.id <- which(self$time == start) - self$p.lag
      end.id <- which(self$time == end) - self$p.lag
      HDC <- HDC_ts(start.id, end.id, shock.id, res.id, self$IRF, self$eps)
      return(HDC)
    }
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
    beta.draw = NULL, ## matrix of alpha
    Sigma.draw = NULL, ## cov of Y-X \times beta, used for cholesky decomposition
    B.draw = NULL,
    beta.NSR = NULL, ## matrix of alpha
    Sigma.NSR = NULL, ## cov of Y-X \times beta, used for cholesky decomposition
    B.NSR = NULL,
    # related non-orthogonal/orthogonal shocks and IRFs
    IRF.draw = NULL,
    IRF.NSR = NULL,
    initialize = function(data = NA, p.lag = NA, prior.type = "info", priors = NA) {
      # need to check data and p.lag input
      super$initialize(data, p.lag)
      self$prior.type <- prior.type
      if (!is.na(priors)) {
        self$alpha.prior <- list(
          "mean" = priors[[1]],
          "cov" = priors[[2]]
        )
        self$Sigma.prior <- list(
          "mean" = priors[[3]],
          "nu" = priors[[4]]
        )
        self$prior.type <- "user"
        cat("* Bayesian VAR priors set by user.\n")
      } else if (prior.type == "flat") {
        self$alpha.prior <- list(
          "mean" = as.vector(t(cbind(rep(0, self$n.var), diag(n.var), matrix(0, self$n.var, self$n.var * (self$p.lag - 1))))),
          "cov" = 0.1 * diag(self$n.var * (self$n.var * self$p.lag + 1))
        )
        self$Sigma.prior <- list(
          "mean" = diag(self$n.var),
          "nu" = 0
        )
        cat("* Bayesian VAR priors set as cheap guess.\n")
      } else if (prior.type == "info") {
        self$alpha.prior <- list(
          "mean" = as.vector(self$beta),
          "cov" = 0.1 * diag(self$n.var * (self$n.var * self$p.lag + 1))
        )
        self$Sigma.prior <- list(
          "mean" = self$Sigma,
          "nu" = 0
        )
        cat("* Bayesian VAR priors set as OLS estimates.\n")
      } else {
        stop("please follow the instructions to set priors")
      }
    },
    est = function(Y = self$Y, X = self$X, method = "gibbs", burn.in = 100, draw = 1000, thin = 1, post.save = NA, post.method = 0) {
      if (method == "gibbs") {
        if (is.na(post.save)) post.save <- draw / 2
        priors <- c(self$alpha.prior, self$Sigma.prior)
        cat("* Gibbs sampler start...\n")
        tic("Gibbs sampler time usage")
        out <- gibbs_sampler(Y, X, priors, burn.in, draw, thin, post.save, post.method)
        toc()
        ## save posterior parameters
        self$alpha.post <- list("mean" = out$alpha_mean_post, "cov" = out$alpha_cov_post)
        self$Sigma.post <- list("mean" = out$Sigma_mean_post, "nu" = out$nu_post)
        ## temporally save parameters estimated, if not using sign restrictions
        self$beta <- matrix(rowMeans(out$alpha), self$n.var * self$p.lag + 1, self$n.var)
        self$Sigma <- apply(out$Sigma, c(1, 2), mean)
        cat("* posterior parameters estimated using ", method, " method.\n")
        ## return loglik for potential analysis
        return(out$loglik)
      } else if (method == "conjugate") {
        # define quantities following Kilian (2017)
        V <- 0.1 * diag(self$n.var * self$p.lag + 1)
        A.star <- t(matrix(self$alpha.prior$mean, self$n.var * self$p.lag + 1, self$n.var))
        Sigma.mu <- self$Sigma
        S.star <- self$Sigma.prior$mean
        # conjugate estimation of posteriors, instead of Gibbs sampler, to save time
        conj.out <- private$conjugate.post(self$Y, self$X, V, A.star, Sigma.mu, S.star)
        self$alpha.post <- list("mean" = as.vector(t(conj.out$A.bar)), "cov" = conj.out$Sigma.alpha.bar)
        self$Sigma.post <- list("mean" = conj.out$S, "nu" = conj.out$tau)
        cat("* posterior parameters estimated using ", method, " method.\n")
      } else {
        stop("please set appropriate estimation method!")
      }
    },
    identify.seq = function(method = "sign", SR = NA, EBR = NA, NSR = NA, draw = 5000, M = 1000, IV = NA) {
      if (method == "recursive" | method == "IV") {
        super$identify(method, IV)
      } else if (method == "sign") {
        restrictions.out <- private$get.restrictions(SR, EBR, NSR)
        restrictions <- restrictions.out$restrictions
        max.irf <- restrictions.out$max.irf
        is.nsr <- restrictions.out$is.nsr
        res.num <- length(restrictions)
        cat("* ", res.num, " restrictions get.\n* First impose SR and EBR on drawn samples:\n")
        mystr <- paste("*", as.character(draw), "draws saved with sign and elasticity restrictions satisfied, total time usage")
        tic(msg = mystr)
        sign.out <- impose_SR_and_EBR(self$alpha.post, self$Sigma.post, restrictions, self$Y, self$X, draw, self$p.lag, max.irf)
        toc()
        # save drawn structural parameters and IRFs
        self$beta.draw <- sign.out$beta_saved
        self$B.draw <- sign.out$B_saved
        self$Sigma.draw <- sign.out$Sigma_saved
        if (is.nsr) {
          cat("* checking NSR...\n")
          mystr <- "* Narrative restrictions imposed, time usage"
          tic(msg = mystr)
          NSR.out <- impose_NSR(restrictions, self$Y, self$X, self$B.draw, self$beta.draw, self$Sigma.draw, self$p.lag, max.irf, M)
          toc()
          self$beta.NSR <- NSR.out$beta_NSR
          self$B.NSR <- NSR.out$B_NSR
          self$Sigma.NSR <- NSR.out$Sigma_NSR
        }
        cat("* identification done.\n")
      } else {
        stop("Please set appropriate identify method!")
      }
    },
    identify.all = function(SR = NA, EBR = NA, NSR = NA, draw = 5000, M = 1000) {
      restrictions.out <- private$get.restrictions(SR, EBR, NSR)
      restrictions <- restrictions.out$restrictions
      max.irf <- restrictions.out$max.irf
      res.num <- length(restrictions)
      cat("* ", res.num, " restrictions get.\n* impose all restrictions on drawn samples:\n")
      mystr <- paste("*", as.character(draw), "draws saved with all restrictions satisfied, total time usage")
      tic(msg = mystr)
      sign.out <- impose_all_restrictions(self$alpha.post, self$Sigma.post, restrictions, self$Y, self$X, draw, self$p.lag, max.irf, M)
      toc()
      # save drawn structural parameters and IRFs
      self$beta.draw <- sign.out$beta_saved
      self$B.draw <- sign.out$B_saved
      self$Sigma.draw <- sign.out$Sigma_saved
      cat("identification done.\n")
    },
    IRF.compute = function(hor, prob, set = "traditional") {
      msg <- paste("* IRF computed and", prob, "HDP get, time usage:")
      HDP <- c((1 - prob) / 2, 1 - (1 - prob) / 2)
      if (set == "traditional") {
        tic(msg)
        self$IRF.draw <- IRF_draws(self$beta.draw, self$B.draw, hor)
        IRF.avg <- apply(self$IRF.draw, c(1, 2), median)
        IRF.ub <- apply(self$IRF.draw, c(1, 2), function(x) quantile(x, probs = HDP[2]))
        IRF.lb <- apply(self$IRF.draw, c(1, 2), function(x) quantile(x, probs = HDP[1]))
        toc()
      } else if (set == "NSR") {
        tic(msg)
        self$IRF.NSR <- IRF_draws(self$beta.NSR, self$B.NSR, hor)
        IRF.avg <- apply(self$IRF.NSR, c(1, 2), median)
        IRF.ub <- apply(self$IRF.NSR, c(1, 2), function(x) quantile(x, probs = HDP[2]))
        IRF.lb <- apply(self$IRF.NSR, c(1, 2), function(x) quantile(x, probs = HDP[1]))
        toc()
      } else {
        cat("please specify correct IRF computation type")
      }
      return(
        list(
          "avg" = IRF.avg,
          "ub" = IRF.ub,
          "lb" = IRF.lb
        )
      )
    }
  ),
  private = list(
    conjugate.post = function(Y.data, X, V, A.star, Sigma.mu, S.star) {
      Y <- t(Y.data)
      Z <- t(X)
      T <- dim(Y)[2]
      n <- dim(Y)[1]
      V.inv <- solve(V)
      my.inv <- solve(V.inv + tcrossprod(Z))
      Sigma.alpha.bar <- my.inv %x% Sigma.mu
      A.bar <- (A.star %*% V.inv + Y %*% t(Z)) %*% my.inv
      A.hat <- tcrossprod(Y, Z) %*% solve(tcrossprod(Z))
      Sigma.mu.tilde <- tcrossprod(Y - A.hat %*% Z) / T
      S <- T * Sigma.mu.tilde + S.star + tcrossprod(A.hat %*% Z) + A.star %*% V.inv %*% t(A.star) - A.bar %*% (V.inv + tcrossprod(Z)) %*% t(A.bar)
      # tau <- T + n
      tau <- T
      return(
        list(
          "A.bar" = A.bar,
          "Sigma.alpha.bar" = Sigma.alpha.bar,
          "S" = S,
          "tau" = tau
        )
      )
    },
    get.restrictions = function(SR = NA, EBR = NA, NSR = NA) {
      restrictions <- list()
      max.irf <- 0
      flag <- 1
      is.NSR <- FALSE
      if (any(!is.na(SR))) {
        for (i in 1:length(SR)) {
          len.var <- length(SR[[i]][[2]])
          shock <- SR[[i]][[1]]
          h <- SR[[i]][[3]]
          if (h > max.irf) max.irf <- h
          sign <- SR[[i]][[4]]
          for (j in 1:len.var) {
            var <- SR[[i]][[2]][j]
            restrictions[[flag]] <- list(
              "type" = "SR",
              "shock" = which(self$var.names == shock),
              "var" = which(self$var.names == var),
              "h" = h,
              "sign" = sign
            )
            flag <- flag + 1
          }
        }
      }
      if (any(!is.na(EBR))) {
        for (i in 1:length(EBR)) {
          shock <- EBR[[i]][[1]]
          h <- EBR[[i]][[3]]
          if (h > max.irf) max.irf <- h
          mb <- EBR[[i]][[4]]
          lb <- EBR[[i]][[5]]
          if (is.na(mb)) mb <- Inf
          if (is.na(lb)) lb <- -Inf
          var.1 <- EBR[[i]][[2]][1]
          var.2 <- EBR[[i]][[2]][2]
          restrictions[[flag]] <- list(
            "type" = "EBR",
            "shock" = which(self$var.names == shock),
            "var.1" = which(self$var.names == var.1),
            "var.2" = which(self$var.names == var.2),
            "h" = h,
            "max.bound" = mb,
            "low.bound" = lb
          )
          flag <- flag + 1
        }
      }
      if (any(!is.na(NSR))) {
        is.NSR <- TRUE
        for (i in 1:length(NSR)) {
          shock <- NSR[[i]][[1]]
          type <- NSR[[i]][[2]]
          start <- which(self$time == NSR[[i]][[3]]) - self$p.lag
          end <- which(self$time == NSR[[i]][[4]]) - self$p.lag
          if ((end - start) > max.irf) max.irf <- end - start
          if (length(NSR[[i]]) == 5) {
            sign <- NSR[[i]][[5]]
            restrictions[[flag]] <- list(
              "type" = "NSR",
              "NSR.type" = "sign",
              "shock" = which(self$var.names == shock),
              "period" = start:end,
              "sign" = sign,
            )
            flag <- flag + 1
          } else {
            var <- NSR[[i]][[5]]
            sign <- NSR[[i]][[6]]
            intensity <- NSR[[i]][[7]]
            restrictions[[flag]] <- list(
              "type" = "NSR",
              "NSR.type" = "contribution",
              "shock" = which(self$var.names == shock),
              "var" = which(self$var.names == var),
              "period" = start:end,
              "sign" = sign,
              "intensity" = intensity
            )
            flag <- flag + 1
          }
        }
      }
      return(list(
        "restrictions" = restrictions,
        "max.irf" = max.irf,
        "is.nsr" = is.NSR
      ))
    }
  )
)
