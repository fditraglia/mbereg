# This test uses the Frazis & Loewenstein / Mahajan moment equalities, which
# assume that Tstar and z are jointly exogenous. We use the same moment
# inequalities as in the endogenous Tstar case.
GMS_test_alphas_nondiff_FL <- function(a0, a1, dat, normal_sims){

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

  # Calculate the quantiles, ensuring that they are always valid, and never
  # equal the sample max or min by "trimming" the r_tk but do so in a way that
  # is asymptotically negligible even when scaled up by sqrt(n)
  p0 <- mean(Tobs[z == 0])
  p1 <- mean(Tobs[z == 1])
  epsilon <- 2 / (sqrt(n) * log(n))
  r_min <- 0 + epsilon
  r_max <- 1 - epsilon
  r00 <- max(r_min, min(r_max, (a1 / (1 - p0)) * (p0 - a0) / (1 - a0 - a1)))
  r10 <- max(r_min, min(r_max, ((1 - a1) / p0) * (p0 - a0) / (1 - a0 - a1)))
  r01 <- max(r_min, min(r_max, (a1 / (1 - p1)) * (p1 - a0) / (1 - a0 - a1)))
  r11 <- max(r_min, min(r_max, ((1 - a1) / p1) * (p1 - a0) / (1 - a0 - a1)))
  q_under_00 <- quantile(y[Tobs == 0 & z == 0], r00)
  q_over_00 <- quantile(y[Tobs == 0 & z == 0], 1 - r00)
  q_under_10 <- quantile(y[Tobs == 1 & z == 0], r10)
  q_over_10 <- quantile(y[Tobs == 1 & z == 0], 1 - r10)
  q_under_01 <- quantile(y[Tobs == 0 & z == 1], r01)
  q_over_01 <- quantile(y[Tobs == 0 & z == 1], 1 - r01)
  q_under_11 <- quantile(y[Tobs == 1 & z == 1], r11)
  q_over_11 <- quantile(y[Tobs == 1 & z == 1], 1 - r11)
  quantiles <- c(q_under_00, q_over_00,
                 q_under_10, q_over_10,
                 q_under_01, q_over_01,
                 q_under_11, q_over_11)

  # Moment equalities for preliminary estimation of quantiles
  if(a1 != 0){
    hI_under_00 <- (y <= q_under_00) * (z == 0) * (1 - Tobs) -
      (a1 / (1 - a0 - a1)) * (z == 0) * (Tobs - a0)
    hI_over_00 <- (y <= q_over_00) * (z == 0) * (1 - Tobs) -
      ((1 - a0) / (1 - a0 - a1)) * (z == 0) * (1 - Tobs - a1)
  } else {
    hI_under_00 <- hI_over_00 <- rep(0, n)
  }

  if(a0 != 0){
    hI_under_10 <- (y <= q_under_10) * (z == 0) * Tobs -
      ((1 - a1) / (1 - a0 - a1)) * (z == 0) * (Tobs - a0)
    hI_over_10 <- (y <= q_over_10) * (z == 0) * Tobs -
      (a0 / (1 - a0 - a1)) * (z == 0) * (1 - Tobs - a1)
  } else {
    hI_under_10 <- hI_over_10 <- rep(0, n)
  }

  if(a1 != 0){
    hI_under_01 <- (y <= q_under_01) * (z == 1) * (1 - Tobs) -
      (a1 / (1 - a0 - a1)) * (z == 1) * (Tobs - a0)
    hI_over_01 <- (y <= q_over_01) * (z == 1) * (1 - Tobs) -
      ((1 - a0) / (1 - a0 - a1)) * (z == 1) * (1 - Tobs - a1)
  } else {
    hI_under_01 <- hI_over_01 <- rep(0, n)
  }

  if(a0 != 0){
    hI_under_11 <- (y <= q_under_11) * (z == 1) * Tobs -
      ((1 - a1) / (1 - a0 - a1)) * (z == 1) * (Tobs - a0)
    hI_over_11 <- (y <= q_over_11) * (z == 1) * Tobs -
      (a0 / (1 - a0 - a1)) * (z == 1) * (1 - Tobs - a1)
  } else {
    hI_under_11 <- hI_over_11 <- rep(0, n)
  }

  hI <- cbind(hI_under_00, hI_over_00,
              hI_under_10, hI_over_10,
              hI_under_01, hI_over_01,
              hI_under_11, hI_over_11)

  # Inequalities for non-differential measurement error
  if(a1 != 0){
    mI2_under_00 <- (y * (z == 0) * (Tobs - a0)) -
      ((1 - a0 - a1) / a1) * y * (y <= q_under_00) * (z == 0) * (1 - Tobs)
    mI2_over_00 <- (-y * (z == 0) * (Tobs - a0)) +
      ((1 - a0 - a1) / a1) * y * (y > q_over_00) * (z == 0) * (1 - Tobs)
  } else {
    mI2_under_00 <- mI2_over_00 <- rep(0, n)
  }

  if(a0 != 0){
    mI2_under_10 <- (y * (z == 0) * (Tobs - a0)) -
      ((1 - a0 - a1) / (1 - a1)) * y * (y <= q_under_10) * (z == 0) * Tobs
    mI2_over_10 <- (-y * (z == 0) * (Tobs - a0)) +
      ((1 - a0 - a1) / (1 - a1)) * y * (y > q_over_10) * (z == 0) * Tobs
  } else {
    mI2_under_10 <- mI2_over_10 <- rep(0, n)
  }

  if(a1 != 0){
    mI2_under_01 <- (y * (z == 1) * (Tobs - a0)) -
      ((1 - a0 - a1) / a1) * y * (y <= q_under_01) * (z == 1) * (1 - Tobs)
    mI2_over_01 <- (-y * (z == 1) * (Tobs - a0)) +
      ((1 - a0 - a1) / a1) * y * (y > q_over_01) * (z == 1) * (1 - Tobs)
  } else {
    mI2_under_01 <- mI2_over_01 <- rep(0, n)
  }

  if(a0 != 0){
    mI2_under_11 <- (y * (z == 1) * (Tobs - a0)) -
      ((1 - a0 - a1) / (1 - a1)) * y * (y <= q_under_11) * (z == 1) * Tobs
    mI2_over_11 <- (-y * (z == 1) * (Tobs - a0)) +
      ((1 - a0 - a1) / (1 - a1)) * y * (y > q_over_11) * (z == 1) * Tobs
  } else {
    mI2_under_11 <- mI2_over_11 <- rep(0, n)
  }

  mI2 <- cbind(mI2_under_00, mI2_over_00,
               mI2_under_10, mI2_over_10,
               mI2_under_01, mI2_over_01,
               mI2_under_11, mI2_over_11)

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
  h <- cbind(hI, hE)
  mI <- cbind(mI1, mI2)
  m <- cbind(mI, mE)

  # Derivative matrices to adjust variance for preliminary estimation
  p <- mean(Tobs)
  mu_Tz <- mean(Tobs * z)
  BE <- matrix(c(-p * mu_Tz, p^2,
                 -mu_Tz^2, mu_Tz * p), 2, 2, byrow = TRUE) / cov(Tobs, z)

  BI <- (1 - a0 - a1) * diag(quantiles / rep(c(rep(a1, 2), rep(1 - a1, 2)), 2))
  BI[is.infinite(BI)] <- 0
  B <- rbind(matrix(0, 4, 10),
             cbind(BI, matrix(0, 8, 2)),
             cbind(matrix(0, 2, 8), BE))

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
