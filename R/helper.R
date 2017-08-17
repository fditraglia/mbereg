# ----------------------------------------------------------------------------
# Calculate the "MMM" test stat - "S1" from Andrews and Soares
# ----------------------------------------------------------------------------
# - m is the vector $\sqrt{n} * \bar{m}_n(theta_0)$
# - Sigma is the matrix $\widehat{\Sigma}_n(theta_0)
# - p is the number of inequality moment conditions, assumed to correpond to
#     the first p elements of m
# ----------------------------------------------------------------------------
get_test_stat <- function(m, Sigma, p){
  m[1:p] <- ifelse(m[1:p] < 0, m[1:p], 0)
  sum(m^2 / diag(Sigma))
}


#----------------------------------------------------------------------------
# Square root of symmetric, PSD matrix. Based on code from MASS::mvrnorm.
#----------------------------------------------------------------------------
#     Returns M_sqrt such that all.equal(M, M_sqrt %*% t(M_sqrt)) is TRUE.
#----------------------------------------------------------------------------
sqrtm <- function(M, tol = 1e-06){
  p <- nrow(M)
  eM <- eigen(M, symmetric = TRUE)
  ev <- eM$values
  stopifnot(all(ev >= -tol * abs(ev[1L])))
  M_sqrt <- eM$vectors %*% diag(sqrt(pmax(ev, 0)), p)
  return(M_sqrt)
}


get_coverage <- function(x, CIs, NAempty = TRUE){
  lower <- apply(CIs, 1, min)
  upper <- apply(CIs, 1, max)
  if(NAempty) {
    lower[is.na(lower)] <- Inf
    upper[is.na(upper)] <- -Inf
  } else {
    lower[is.na(lower)] <- -Inf
    upper[is.na(upper)] <- Inf
  }
  sapply(x, function(x)  mean((lower < x) & (upper > x)))
}

get_median_width <- function(CIs, NAempty = TRUE) {
  lower <- apply(CIs, 1, min)
  upper <- apply(CIs, 1, max)
  if(NAempty) {
    lower[is.na(lower)] <- 0
    upper[is.na(upper)] <- 0
  } else {
    lower[is.na(lower)] <- -Inf
    upper[is.na(upper)] <- Inf
  }
  median(upper - lower)
}
