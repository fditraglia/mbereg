# This test substitutes a preliminary estimator of theta1 and uses both the
# first moment and non-differential measurement error inequalities. This version
# should be faster since it avoids  unnecessary calculations by dropping
# inequalities that are far from binding rather than carrying them around.
GMS_test_alphas_nondiff <- function(a0, a1, dat, normal_sims){

  if((a0 > 0.99) || (a0 < 0)) return(0)
  if((a1 > 0.99) || (a1 < 0)) return(0)
  if(a0 + a1 > 0.99) return(0)

  Tobs <- dat$Tobs
  y <- dat$y
  z <- dat$z
  q <- mean(z) # treat this as fixed in repeated sampling
  n <- nrow(dat)
  kappa_n <- sqrt(log(n))

  # Preliminary estimate of theta1
  theta1 <- cov(z, y) / cov(z, Tobs)

  a2 <- (1 + (a0 - a1))
  a3 <- ((1 - a0 - a1)^2 + 6 * a0 * (1 - a1))
  theta2 <- theta1^2 * a2
  theta3 <- theta1^3 * a3

  # Functions of the data
  w <- cbind(Tobs, y, y * Tobs, y^2, y^2 * Tobs, y^3)
  w_bar <- colMeans(w)
  w_centered <- scale(w, scale = FALSE, center = TRUE)
  w_centered_z <- (w_centered * z)
  wz_bar <- colMeans(w * z)

  psi1 <- c(-theta1, 1, 0, 0, 0, 0)
  psi2 <- c(theta2, 0, -2 * theta1, 1, 0, 0)
  psi3 <- c(-theta3, 0, 3 * theta2, 0, -3 * theta1, 1)
  Psi <- cbind(psi1, psi2, psi3)

  # Derivatives of psi1, psi2, psi3 with respect to theta1
  D_psi1 <- c(-1, 0, 0, 0, 0, 0)
  D_psi2 <- c(2 * a2 * theta1, 0, -2, 0, 0, 0)
  D_psi3 <- c(-3 * a3 * theta1^2, 0, 6 * a2 * theta1, 0, -3, 0)
  D_Psi <- cbind(D_psi1, D_psi2, D_psi3)

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
  hE <- cbind(w_centered %*% Psi , w_centered_z %*% psi1)
  mE <- w_centered_z %*% Psi[,2:3]

  # Stack all moment conditions
  h <- cbind(hI, hE)
  mI <- cbind(mI1, mI2)
  m <- cbind(mI, mE)

  # Derivative matrices to adjust variance for preliminary estimation
  ME_kappa <- matrix(c(0, -q, 0,
                       0, 0, -q), nrow = 2, ncol = 3, byrow = TRUE)
  ME_theta1 <- crossprod(D_Psi[,2:3], wz_bar)
  ME <- cbind(ME_kappa, ME_theta1)

  p <- w_bar[1]
  p1_q <- wz_bar[1]
  d2_w <- drop(crossprod(D_psi2, w_bar))
  d3_w <- drop(crossprod(D_psi3, w_bar))
  HE_inv <- matrix(c(p1_q, 0, 0, -p,
                     -q * d2_w, -(p * q - p1_q), 0, d2_w,
                     -q * d3_w, 0, -(p * q - p1_q), d3_w,
                     -q, 0, 0, 1), byrow = TRUE, nrow = 4, ncol = 4) / (p * q - p1_q)
  BE <- -ME %*% HE_inv
  BI <- (1 - a0 - a1) * diag(quantiles / rep(c(rep(a1, 2), rep(1 - a1, 2)), 2))
  B <- rbind(matrix(0, 4, 12),
             cbind(BI, matrix(0, 8, 4)),
             cbind(matrix(0, 2, 8), BE))

  # Adjust the variance matrix
  V <- var(cbind(m, h))
  A <- cbind(diag(nrow(B)), B)
  Sigma_n <- A %*% V %*% t(A)

  # Calculate test statistic
  m_bar <- colMeans(m)
  n_ineq <- ncol(mI)
  n_eq <- ncol(mE)
  T_n <- get_test_stat(sqrt(n) * m_bar, Sigma_n, n_ineq)

  # Calculate p-value of asymptotic test, only keeping inequalities that are
  # *far* from binding
  s_n <- sqrt(diag(Sigma_n))
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
