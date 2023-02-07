RUV4_upgrade <- function(Y, X, ctl, k, Z = NULL, eta = NULL, include.intercept = TRUE, 
                 fullW0 = NULL, inputcheck = TRUE) 
{
  m = nrow(Y)
  n = ncol(Y)
  p = ncol(X)
  ctl = ctl2logi(ctl, n)
  Y = ruv::RUV1(Y, eta, ctl, include.intercept = include.intercept)
  Y0 = orthogonal_projection(Y, X)
  if (is.null(fullW0)) {
    full_U = svd(Y0 %*% t(Y0))$u[, 1:(m - p), drop = FALSE]
  }
  if (k > 0) {
    U = full_U[, 1:k, drop = FALSE]
    alpha = t(U) %*% Y0
    Y0c = Y0[, ctl, drop = FALSE]
    W = Y0c %*% t(Y0c) %*% U %*% solve(t(U) %*% Y0c %*% t(Y0c) %*% U)
    WA = W %*% alpha
  }
  newY = Y - W %*% alpha
  return(list(
    newY = newY,
    full_U = full_U,
    W = W,
    WA = W %*% alpha,
    alpha = alpha))
}

ctl2logi <- function(ctl, n)
{
  ctl2 = rep(FALSE, n)
  ctl2[ctl] = TRUE
  return(ctl2)
}

orthogonal_projection <- function (A, B) {
  return(A - B %*% solve(t(B) %*% B) %*% t(B) %*% A)
}