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
  sapply(x, function(x)  100 * round(mean((lower < x) & (upper > x)), 2))
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
  round(median(upper - lower), 2)
}

get_twostep_CI <- function(CIs_GMM, CIs_bonf) {
 CIs_GMM[is.na(CIs_GMM[,1]),1] <- -Inf
 CIs_GMM[is.na(CIs_GMM[,2]),2] <- +Inf
 use_GMM <- (CIs_GMM[,1] > CIs_bonf[,1]) & (CIs_GMM[,2] < CIs_bonf[,2])
 lower <- ifelse(use_GMM, CIs_GMM[,1], CIs_bonf[,1])
 upper <- ifelse(use_GMM, CIs_GMM[,2], CIs_bonf[,2])
 out <- cbind(lower, upper)
 return(out)
}
