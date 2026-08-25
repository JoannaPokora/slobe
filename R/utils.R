
generate_data <- function(
    n = 100,
    p = 100,
    signals_num = 5,
    signals_strength = 1,
    rho = 0,
    p_missing = 0,
    response = "gaussian"
) {
  beta <- rep(0, p)
  beta[sample(1:p, signals_num)] <- signals_strength
  
  corr_matr <- toeplitz(rho^(0:(p - 1)))
  
  X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = corr_matr / n)
  
  y <- if (response == "gaussian") {
    X %*% beta + rnorm(n)
  } else {
    prob <- 1 / (1 + exp(-X %*% beta))
    rbinom(n, 1, prob)
  }
  
  if (p_missing != 0) {
    X_nomis <- X
    X[runif(n * p) < p_missing] <- NA 
  }
  
  list(
    X = X,
    X_nomis = {
      if (p_missing != 0)
        X_nomis
      else
        NULL
    },
    y = y,
    beta = beta
  )
}
