# This test uses the Frazis & Loewenstein / Mahajan moment equalities, which
# assume that Tstar and z are jointly exogenous. We use the *only* the weak
# bounds for a0 and a1.
GMS_test_alphas_FL <- function(a0, a1, dat, normal_sims){

  if((a0 > 0.99) || (a0 < 0)) return(0)
  if((a1 > 0.99) || (a1 < 0)) return(0)
  if(a0 + a1 > 0.99) return(0)

  Tobs <- dat$Tobs
  y <- dat$y
  z <- dat$z
  q <- mean(z) # treat this as fixed in repeated sampling
  n <- nrow(dat)
  kappa_n <- sqrt(log(n))

  # Preliminary estimates of theta1 and kappa1
  theta1 <- cov(z, y) / cov(z, Tobs)
  kappa1 <- mean(y) - theta1 * mean(Tobs)
  rho <- -theta1 * a0 * (1 - a1)
  eta <- theta1 * (1 + a0 - a1)

  # First moment inequalities
  mI1 <- cbind((1 - z) * (Tobs - a0),
               (1 - z) * (1 - Tobs - a1),
               z * (Tobs - a0),
               z * (1 - Tobs - a1))

  # Preliminary estimators and moment functions for equality conditions
  hE <- cbind(y - kappa1 - theta1 * Tobs,
              (y - kappa1 - theta1 * Tobs) * z)
  mE <- cbind((y - kappa1) * Tobs - rho - eta * Tobs,
              ((y - kappa1) * Tobs - rho - eta * Tobs) * z)

  # Stack all moment conditions
  h <- hE # no preliminary estimation for inequalities
  mI <- mI1 # only the weak bounds for (a0, a1)
  m <- cbind(mI, mE)

  # Derivative matrices to adjust variance for preliminary estimation
  p <- mean(Tobs)
  mu_Tz <- mean(Tobs * z)
  BE <- matrix(c(-p * mu_Tz, p^2,
                 -mu_Tz^2, mu_Tz * p), 2, 2, byrow = TRUE) / cov(Tobs, z)
  B <- rbind(matrix(0, 4, 2), BE)

  # Adjust the variance matrix
  V <- var(cbind(m, h))
  A <- cbind(diag(nrow(B)), B)
  Sigma_n <- A %*% V %*% t(A)
  s_n <- sqrt(diag(Sigma_n))

  # Calculate test statistic
  n_ineq <- ncol(mI)
  n_eq <- ncol(mE)
  m_bar <- colMeans(m)
  m_bar_ineq <- m_bar[1:n_ineq]
  m_bar_eq <- m_bar[(n_ineq + 1):(n_ineq + n_eq)]
  s_n_ineq <- s_n[1:n_ineq]
  s_n_eq <- s_n[(n_ineq + 1):(n_ineq + n_eq)]
  T_n_eq <- sum((sqrt(n) * m_bar_eq / s_n_eq)^2)
  T_n_ineq <- sum(ifelse(m_bar_ineq < 0, (sqrt(n) * m_bar_ineq / s_n_ineq)^2, 0))
  T_n <- T_n_ineq + T_n_eq

  # Calculate p-value of asymptotic test, only keeping inequalities that are
  # *far* from binding
  # Write the t-test as a *strict* inequality with s_n on the RHS
  # to correctly handle an inequality that does not bind but has zero variance:
  # it should be treated as *deterministic* and dropped
  keep_ineq <- (sqrt(n) * m_bar[1:n_ineq]) < (kappa_n * s_n[1:n_ineq])
  n_keep_ineq <- sum(keep_ineq)
  keep <- c(keep_ineq, rep(TRUE, n_eq))
  # If an inequality has zero variance and is satisfied, it was dropped in the
  # preceding step. If any of the inequalities that remain have zero variance
  # then they must bind. In this case, the p-value is zero: we have violated a
  # deterministic bound. Same holds for a moment equality
  if(isTRUE(all.equal(min(s_n[keep]), 0))) return(0)

  Omega_n_GMS <- cov2cor(Sigma_n[keep, keep])
  M_star_GMS <- sqrtm(Omega_n_GMS) %*% normal_sims[keep,]
  if(sum(n_keep_ineq) > 0){
    M_star_GMS_I <- M_star_GMS[1:n_keep_ineq,,drop=FALSE]
    M_star_GMS_I <- M_star_GMS_I * (M_star_GMS_I < 0)
    M_star_GMS_E <- M_star_GMS[-c(1:n_keep_ineq),,drop=FALSE]
    M_star_GMS <- rbind(M_star_GMS_I, M_star_GMS_E)
  }
  T_n_star <- colSums(M_star_GMS^2)
  return(mean(T_n_star > T_n))
}
