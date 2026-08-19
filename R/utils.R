generate_data <- function(
    n = 100,
    p = 100,
    signals_num = 5,
    signals_strength = 1,
    corr_matr = "identity",
    p_missing = 0,
    response = "gaussian"
) {
  beta <- rep(0, p)
  beta[sample(1:p, signals_num)] <- signals_strength
  
  if (corr_matr[1] == "identity") {
    corr_matr <- diag(nrow = p)
  }
  
  X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = corr_matr)
  
  y <- if (response == "gaussian") {
    X %*% beta + rnorm(n)
  } else {
    prob <- 1 / (1 + exp(-X %*% beta))
    rbinom(n, 1, prob)
  }
  
  if (p_missing != 0) {
    X_vals_n <- n * p
    X[sample(1:X_vals_n, X_vals_n * p_missing)] <- NA 
  }
  
  list(
    X = X,
    y = y,
    beta = beta
  )
}
