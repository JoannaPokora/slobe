#' rescale
#'
#' @param x A number
#' @param y A number
#'
#' @keywords internal

rescale <- function(y, x) {
  z <- data.frame(y = y, x = x)
  z <- stats::na.omit(z)

  stats::lm(x ~ y, data = z)$coef
}

#' rescale_all
#'
#' @param results A number
#' @param Xmis A number
#'
#' @keywords internal

rescale_all <- function(results, Xmis) {
  k <- ncol(results[["X"]])
  scales <- sapply(1:k, function(l) rescale(results[["X"]][, l], Xmis[, l]))
  results[["X"]] <- t(t(results[["X"]]) * scales[2, ] + scales[1, ])
  results[["Sigma"]] <- results[["Sigma"]] * (scales[2, ] %*% t(scales[2, ]))
  results[["mu"]] <- results[["mu"]] * scales[2, ] + scales[1, ]
  
  beta <- results[["coefficients"]] / scales[2,]
  intercept <- -sum(scales[1,] * results[["coefficients"]])
  coefs <- c(intercept, beta)
  coefs[is.na(coefs)] <- 0
  names(coefs) <- c("(Intercept)", paste0("V", 1:(length(beta))))

  results[["coefficients"]] <- coefs
  results
}

#' center_and_scale
#'
#' @param X A matrix
#'
#' @keywords internal

center_and_scale <- function(X) {
  X <- sweep(X, 2, colMeans(X), "-")
  sweep(X, 2, sqrt(colSums(X^2)), "/")
}

#' SLOBE
#'
#' @description Fit a gaussian model regularized with Adaptive Bayesian SLOPE
#'   and handle missing values by Stochastic Approximation of Expected
#'   Maximization (SAEM)
#'
#' @param start the initial vector of regression coefficients for the first
#' iteration. Default to the LASSO estimator obtained after
#'
#' @param Xmis design matrix
#' @param Y numeric. Response variable.
#' @param family model family - either "gaussian" (default) or "binomial"
#' @param start model coefficients used for initialization.
#' @param a_prior,b_prior non-negative parameters of the prior Beta distribution
#'   on theta.
#' @param Covmat numeric covariance matrix. Default to identity matrix.
#' @param sigma the variance of the noise. Default to 1.
#' @param FDR False Discovery Rate. Default to 0.05.
#' @param tol optimization tolerance.
#' @param max_iter the maximal number of iterations of the optimization
#'   algorithm. Default to 100.
#' @param verbose verbosity. Default to FALSE
#' @param BH logical. Indicates whether the Benjamini-Hochberg correction for
#'   multiple testing should be used.
#'
#' @details \code{ABSLOPE} is the combination of SLOPE and Spike-and-Slab
#'   LASSO (SSL). This approach relies on iterations of the weighted SLOPE
#'   algorithm and finds the solution by minimizing
#'   \deqn{
#'    G(Y,X,\beta) + \sum_{j = 1}^{p} w_j \lambda_{r(\beta,j)}|\beta_j|,
#'   }{
#'    G(Y,X,\beta) + \sum w_j |\beta_j|\lambda_r(\beta,j),
#'   }
#'   where \eqn{r(\beta,j) = {1,2,...,p}} is the rank of \eqn{\beta_j}  among
#'   elements in \eqn{\beta} in a descending order. The weight \eqn{w_j} depends
#'   on the posterior probability that a variable \eqn{X_j} is a true predictor
#'   and is calculated based on the prior knowledge and on the  estimator of
#'   \eqn{\beta_j}, the signal sparsity and its average strength from the
#'   previous iterations.
#'
#' @references
#'   Jiang, W., Bogdan, M., Josse, J., Majewski, S., Miasojedow, B.,
#'   Ročková, V., & TraumaBase® Group. (2021). Adaptive Bayesian SLOPE: Model
#'   Selection with Incomplete Data. Journal of Computational and Graphical
#'   Statistics, 1-25. \doi{10.1080/10618600.2021.1963263}
#'
#'   Bogdan, M., van den Berg, E., Sabatti, C., Su, W., & Candès, E. J. (2015).
#'   SLOPE -- adaptive variable selection via convex optimization. The Annals of
#'   Applied Statistics, 9(3), 1103–1140. \doi{10/gfgwzt}
#'
#'   Ročková, V., & George, E. I. (2018). The spike-and-slab lasso. Journal of
#'   the American Statistical Association, 113(521), 431-444.
#'   \doi{10.1080/01621459.2016.1260469}
#'
#' @examples
#' set.seed(17)
#' xy <- SLOPE:::randomProblem(1e2, 200, response = "gaussian")
#' X <- as.matrix(xy$x)
#' y <- xy$y
#' fit <- ABSLOPE(X, y)
#' @export slobe
#'

slobe <- function(X,
                  y,
                  family = "gaussian",
                  start = NULL,
                  a_prior = 0.01 * NROW(X),
                  b_prior = 0.01 * NROW(X),
                  Covmat = NULL,
                  sigma = NULL,
                  FDR = 0.05,
                  tol = 1e-04,
                  max_iter = 100L,
                  verbose = FALSE,
                  BH = TRUE) {

  checkmate::assert_number(a_prior)
  checkmate::assert_number(b_prior)
  checkmate::assert_number(FDR)
  checkmate::assert_number(tol)
  checkmate::assert_number(max_iter)
  checkmate::assert_logical(verbose)
  checkmate::assert_logical(BH)

  if(!(is.matrix(X) | is.data.frame(X))) {
    stop(paste0("X needs to be matrix or data.frame. You provided ",
                paste0(class(X), collapse = ", ")))
  }

  X <- as.matrix(X)

  ocall <- match.call()

  # if Covmat is null -> known_cov = FALSE
  known_cov <- !is.null(Covmat)
  # dummy value for a case with unknown covariance matrix
  if (is.null(Covmat)) {
    Covmat <- diag(NCOL(X))
  }

  # if sigma is null with gaussian family -> known_sigma = FALSE
  known_sigma <- (family == "binomial") || !is.null(sigma)
  if(!known_sigma || family == "binomial") {
    sigma <- 1
  }
  
  # if X has missing values -> column mean imputation
  if (any(is.na(suppressWarnings(as.numeric(X))))) {
    Xmis <- X
    Xinit <- apply(Xmis, 2, function(col) {
      missing_vals <- is.na(col)
      
      if (any(missing_vals)) {
        col[missing_vals] <- mean(col, na.rm = TRUE)
      }
      
      col
    })
  } else {
    Xmis <- NULL
    Xinit <- X
  }

  Xinit <- center_and_scale(Xinit)

  # if start is null -> LASSO gives starting coefficients
  if (is.null(start)) {
    lasso <- glmnet::cv.glmnet(
      Xinit, y, standardize = FALSE, intercept = FALSE, family = family)
    start <- stats::coefficients(lasso, s = "lambda.min")
    start <- start[2:(ncol(Xinit) + 1), 1]
  }

  out <- slobe_fit_cpp(
    start, Xinit, y, Xmis, family, a_prior, b_prior,
    Covmat, sigma, FDR, tol, known_sigma,
    max_iter, verbose, BH, known_cov
  )
  
  if(family == "binomial") {
    out$sigma <- NULL
  }

  out <- rescale_all(out, X)
  out[["call"]] <- ocall
  structure(
    out,
    class = "ABSLOPE"
  )
}
