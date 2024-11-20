# note: the values y, Y and X should come from a VARclass object.
# bayesian estimation of reduced-form VAR, with normal-inverse-Wishart prior and Gibbs sampler

# compute log-likelihood
loglik <- function(Sigma, alpha, y, X, m, T) {
  Sigma.inv <- solve(Sigma)
  alpha.hat <- solve(Sigma.inv %x% (t(X) %*% X), t(Sigma.inv %x% X)) %*% y
  X.cal <- diag(m) %x% X
  normal <- -T / 2 - (t(alpha - alpha.hat) %*% t(X.cal) %*% (Sigma.inv %x% diag(T)) %*% X.cal %*% (alpha - alpha.hat)) / 2
  inv.wishart <- -sum(diag((y - X.cal %*% alpha.hat) %*% t(y - X.cal %*% alpha.hat) %*% (Sigma.inv %x% diag(T)))) / 2
  loglik <- normal + inv.wishart
  return(loglik)
}

mvrnormR <- function(n, mu, sigma) {
  # generate multivariate normal variable
  # from https://gallery.rcpp.org/articles/simulate-multivariate-normal/
  # @return: draw = n (sample number) x m (dimension) matrix
  ncols <- ncol(sigma)
  mu <- rep(mu, each = n) ## not obliged to use a matrix (recycling)
  draw <- mu + matrix(rnorm(n * ncols), ncol = ncols) %*% chol(sigma)
  return(draw)
}

rwish <- function(n, Psi, nu, inv = FALSE) {
  # generate (inv-)wishart random variable
  # from https://rdrr.io/cran/nicheROVER/src/R/rwish.R
  # @return: a m x m x n dimension matrix, n is sample size
  if (inv) Psi <- solve(Psi)
  U <- chol(Psi)
  d <- nrow(Psi)
  ans <- array(0, dim = c(d, d, n))
  if (!is.null(dimnames(Psi))) dimnames(ans) <- c(dimnames(Psi), list(NULL))
  ans[rep(upper.tri(Psi), n)] <- rnorm(n * d * (d - 1) / 2)
  ans[rep(!lower.tri(Psi, diag = FALSE) &
    !upper.tri(Psi, diag = FALSE), n)] <- sqrt(rchisq(n * d, df = nu - 1:d + 1))
  for (ii in 1:n) {
    tmp <- ans[, , ii] %*% U
    if (inv) tmp <- backsolve(tmp, diag(d), transpose = TRUE)
    ans[, , ii] <- crossprod(tmp)
  }
  return(ans)
}

# sample on alpha
updt.alpha <- function(alpha, alpha.prior, sigma.prior.inv, Sigma, y, X, m, T) {
  Sigma.inv <- solve(Sigma)
  X.cal <- diag(m) %x% X
  Sigma.inv.cal <- Sigma.inv %x% diag(T)
  sigma.post <- solve(sigma.prior.inv + t(X.cal) %*% Sigma.inv.cal %*% X.cal)
  alpha.post <- sigma.post %*% (sigma.prior.inv %*% alpha.prior + t(X.cal) %*% Sigma.inv.cal %*% y)
  alpha.updt <- as.vector(mvrnormR(1, alpha.post, sigma.post))
  out <- list(
    "alpha" = alpha.updt,
    "mean.post" = alpha.post,
    "cov.post" = sigma.post
  )
  return(out)
}

# sample on sigma
updt.Sigma <- function(Sigma, Sigma.prior, nu.prior, alpha, Y, X, m, T) {
  beta <- matrix(alpha, ncol = m)
  residual <- Y - X %*% beta
  Sigma.post <- Sigma.prior + t(residual) %*% residual
  nu.post <- nu.prior + T
  Sigma.updt <- rwish(1, Sigma.post, nu.post, inv = TRUE)[, , 1]
  out <- list(
    "Sigma" = Sigma.updt,
    "Sigma.post" = Sigma.post,
    "nu.post" = nu.post
  )
  return(out)
}

# the main function
gsampler <- function(Y, X, priors, burn.in = 100, draw = 1000, thin = 1, post.save = NA, post.method = "last") {
  ## if used for estimating reduced form parameters, just return parameter draws
  ## if used for sign restriction, return parameter draws and posterior distribution parameters
  alpha.bar <- priors[[1]]
  Sigma.alpha <- priors[[2]]
  Sigma.alpha.inv <- solve(Sigma.alpha)
  Sigma.sigma <- priors[[3]]
  nu <- priors[[4]]
  m <- dim(Y)[2] # number of variables
  T.est <- dim(Y)[1] # estimate periods
  K <- length(alpha.bar) # number of parameters in alpha = m x (lag x m + 1)
  alpha.draw <- array(0, c(K, draw))
  Sigma.draw <- array(0, c(m, m, draw))
  loglik.draw <- rep(0, draw)

  if (is.na(post.save)) post.save <- draw / 2
  alpha.post <- array(0, c(K, post.save))
  Sigma.alpha.post <- array(0, c(K, K, post.save))
  Sigma.sigma.post <- array(0, c(m, m, post.save))
  nu.post <- rep(0, post.save)


  alpha <- mvrnormR(1, alpha.bar, Sigma.alpha)
  Sigma <- rwish(1, Sigma.sigma, nu, TRUE)[, , 1]
  y <- as.vector(Y)

  for (i in 1:burn.in) {
    alpha.out <- updt.alpha(alpha, alpha.bar, Sigma.alpha.inv, Sigma, y, X, m, T.est)
    alpha <- alpha.out$alpha
    Sigma.out <- updt.Sigma(Sigma, Sigma.sigma, nu, alpha, Y, X, m, T.est)
    Sigma <- Sigma.out$Sigma
  }

  for (i in 1:draw) {
    for (j in 1:thin) {
      alpha.out <- updt.alpha(alpha, alpha.bar, Sigma.alpha.inv, Sigma, y, X, m, T.est)
      alpha <- alpha.out$alpha
      Sigma.out <- updt.Sigma(Sigma, Sigma.sigma, nu, alpha, Y, X, m, T.est)
      Sigma <- Sigma.out$Sigma
    }
    alpha.draw[, i] <- alpha
    Sigma.draw[, , i] <- Sigma
    loglik.draw[i] <- loglik(Sigma, alpha, y, X, m, T.est)
    if (i >= (draw - post.save + 1)) {
      k <- i - (draw - post.save)
      alpha.post[, k] <- alpha.out$mean.post
      Sigma.alpha.post[, , k] <- alpha.out$cov.post
      Sigma.sigma.post[, , k] <- Sigma.out$Sigma.post
      nu.post[k] <- Sigma.out$nu.post
    }
  }
  if (post.method == "last") {
    alpha.mean.post <- alpha.post[, post.save]
    alpha.cov.post <- Sigma.alpha.post[, , post.save]
    Sigma.mean.post <- Sigma.sigma.post[, , post.save]
    nu.post <- nu.post[post.save]
  } else {
    alpha.mean.post <- colMeans(alpha.post)
    alpha.cov.post <- apply(Sigma.alpha.post, c(1, 2), mean)
    Sigma.mean.post <- apply(Sigma.sigma.post, c(1, 2), mean)
    nu.post <- mean(nu.post)
  }
  return(list(
    "alpha" = alpha.draw,
    "Sigma" = Sigma.draw,
    "loglik" = loglik.draw,
    "alpha.mean.post" = alpha.mean.post,
    "alpha.cov.post" = alpha.cov.post,
    "Sigma.mean.post" = Sigma.mean.post,
    "nu.post" = nu.post
  ))
}
