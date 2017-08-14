# This test uses the Frazis & Loewenstein / Mahajan moment equalities, which
# assume that Tstar and z are jointly exogenous. We use the *only* the weak
# bounds for a0 and a1.
GMS_test_alphas_noineq_FL <- function(a0, a1, dat, normal_sims){

  if((a0 > 0.99) || (a0 < 0)) return(0)
  if((a1 > 0.99) || (a1 < 0)) return(0)
  if(a0 + a1 > 0.99) return(0)

  Tobs <- dat$Tobs
  y <- dat$y
  z <- dat$z
  q <- mean(z) # treat this as fixed in repeated sampling
  n <- nrow(dat)

  # Preliminary estimates of theta1 and kappa1
  theta1 <- cov(z, y) / cov(z, Tobs)
  kappa1 <- mean(y) - theta1 * mean(Tobs)
  rho <- -theta1 * a0 * (1 - a1)
  eta <- theta1 * (1 + a0 - a1)

  # Preliminary estimators and moment functions for equality conditions
  hE <- cbind(y - kappa1 - theta1 * Tobs,
              (y - kappa1 - theta1 * Tobs) * z)
  mE <- cbind((y - kappa1) * Tobs - rho - eta * Tobs,
              ((y - kappa1) * Tobs - rho - eta * Tobs) * z)

  # Stack all moment conditions
  h <- hE # no preliminary estimation for inequalities
  m <- mE # no inequalities

  # Derivative matrices to adjust variance for preliminary estimation
  p <- mean(Tobs)
  mu_Tz <- mean(Tobs * z)
  BE <- matrix(c(-p * mu_Tz, p^2,
                 -mu_Tz^2, mu_Tz * p), 2, 2, byrow = TRUE) / cov(Tobs, z)
  B <- BE # no inequalities

  # Adjust the variance matrix
  V <- var(cbind(m, h))
  A <- cbind(diag(nrow(B)), B)
  Sigma_n <- A %*% V %*% t(A)
  s_n <- sqrt(diag(Sigma_n))

  # Calculate test statistic
  m_bar <- colMeans(m)
  T_n <- sum((sqrt(n) * m_bar / s_n)^2)

  # No inequalities, so no moment selection
  Omega_n <- cov2cor(Sigma_n)
  M_star <- sqrtm(Omega_n) %*% normal_sims
  T_n_star <- colSums(M_star^2)
  return(mean(T_n_star > T_n))
}
