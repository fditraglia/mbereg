# This test substitutes a preliminary estimator of theta1 and uses only the
# both the first and second moment inequalities.
GMS_test_alphas2 <- function(a0, a1, dat, normal_sims){

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

  # Preliminary estimation of parameters from moment inequalities
  # First moment inequalities
  h_p0 <- Tobs * (1 - z) / (1 - q)
  h_p1 <- Tobs * z / q
  p0 <- mean(h_p0)
  p1 <- mean(h_p1)
  h_p <- cbind(h_p0, h_p1)

  # Second moment inequalities
  h_nu_00 <- y * (1 - Tobs) * (1 - z) / sqrt(1 - q)
  h_nu_10 <- y * Tobs * (1 - z) / sqrt(1 - q)
  h_nu_01 <- y * (1 - Tobs) * z / sqrt(q)
  h_nu_11 <- y * Tobs * z / sqrt(q)
  nu_00 <- mean(h_nu_00)
  nu_10 <- mean(h_nu_10)
  nu_01 <- mean(h_nu_01)
  nu_11 <- mean(h_nu_11)
  h_nu <- cbind(h_nu_00, h_nu_10, h_nu_01, h_nu_11)

  # Full set of moment conditions for preliminary estimation
  hI <- scale(cbind(h_p, h_nu), scale = FALSE, center = TRUE)
  hE <- cbind(w_centered %*% Psi , w_centered_z %*% psi1)
  h <- cbind(hI, hE)

  # First moment inequality conditions
  mI_1 <- cbind(rep(p0 - a0, n),
                rep(1 - p0 - a1, n),
                rep(p1 - a0, n),
                rep(1 - p1 - a1, n))

  # Second moment inequality conditions
  y2_a0_z0 <- y^2 * (1 - z) * (Tobs - a0)
  d_a0_z0 <- (1 - a0) * nu_10 - a0 * nu_00

  y2_a1_z0 <- y^2 * (1 - z) * (1 - Tobs - a1)
  d_a1_z0 <- a1 * nu_10 - (1 - a1) * nu_00

  y2_a0_z1 <- y^2 * z * (Tobs - a0)
  d_a0_z1 <- (1 - a0) * nu_11 - a0 * nu_01

  y2_a1_z1 <- y^2 * z * (1 - Tobs - a1)
  d_a1_z1 <- a1 * nu_11 - (1 - a1) * nu_01

  mI_2 <- cbind((p0 - a0) * y2_a0_z0 - d_a0_z0^2,
                (1 - p0 - a1) * y2_a1_z0 - d_a1_z0^2,
                (p1 - a0) * y2_a0_z1 - d_a0_z1^2,
                (1 - p1 - a1) * y2_a1_z1 - d_a1_z1^2)

  # Full set of moment inequalities and equalities
  mI <- cbind(mI_1, mI_2)
  mE <- w_centered_z %*% Psi[,2:3]
  m <- cbind(mI, mE)

  # Derivative matrices to account for preliminary estimation
  R <- cbind(rbind(-1 * diag(3), c(-q, 0, 0)),
             rbind(crossprod(D_Psi, w_bar), crossprod(D_psi1, wz_bar)))
  n_h <- ncol(h)
  H_inv <- -diag(n_h)
  k <- ncol(cbind(h_p, h_nu)) + 1
  H_inv[k:n_h, k:n_h] <- solve(R)


  M_I1_p <- cbind(c(1, -1, 0, 0),
                  c(0, 0, 1, -1))
  M_I1 <- cbind(M_I1_p, matrix(0, 4, 8))
  M_I2_p <- matrix(c( mean(y2_a0_z0), 0,
                     -mean(y2_a1_z0), 0,
                      0, mean(y2_a0_z1),
                      0, -mean(y2_a1_z1)),
                   4, 2, byrow = TRUE)
  M_I2_nu <- matrix(c(2 * a0 * d_a0_z0, -2 * (1 - a0) * d_a0_z0, 0, 0,
                      2 * (1 - a1) * d_a1_z0, -2 * a1 * d_a1_z0, 0, 0,
                      0, 0, 2 * a0 * d_a0_z1, -2 * (1 - a0) * d_a0_z1,
                      0, 0, 2 * (1 - a1) * d_a1_z1, -2 * a1 * d_a1_z1),
                    4, 4, byrow = TRUE)
  M_I2 <- cbind(M_I2_p, M_I2_nu, matrix(0, 4, 4))
  M_I <- rbind(M_I1, M_I2)

  M_E_kappa <- matrix(c(0, -q, 0,
                       0, 0, -q), nrow = 2, ncol = 3, byrow = TRUE)
  M_E_theta1 <- crossprod(D_Psi[,2:3], wz_bar)
  M_E <- cbind(matrix(0, 2, 6), M_E_kappa, M_E_theta1)
  M <- rbind(M_I, M_E)

  # Adjust the variance matrix
  V <- var(cbind(m, h))
  B <- -1 * M %*% H_inv
  A <- cbind(diag(nrow(B)), B)
  Sigma_n <- A %*% V %*% t(A)

  # Calculate test statistic
  m_bar <- colMeans(m)
  n_ineq <- ncol(mI)
  n_eq <- ncol(mE)
  T_n <- get_test_stat(sqrt(n) * m_bar, Sigma_n, n_ineq)

  # Carry out moment selection on the inequality conditions
  s_n <- sqrt(diag(Sigma_n))
  phi <- c(ifelse((sqrt(n) * m_bar[1:n_ineq] / s_n[1:n_ineq]) > sqrt(log(n)),
                  Inf, 0), rep(0, n_eq))

  # Calculate the p-value of asymptotic test
  Omega_n <- cov2cor(Sigma_n)
  M_star <- sqrtm(Omega_n) %*% normal_sims
  M_star_GMS <- M_star + phi
  M_star_GMS_I <- M_star_GMS[1:n_ineq,]
  M_star_GMS_I[!is.finite(M_star_GMS_I)] <- 0
  M_star_GMS_I <- M_star_GMS_I * (M_star_GMS_I < 0)
  M_star_GMS_E <- M_star_GMS[(n_ineq + 1):(n_ineq + n_eq),]
  M_star_GMS <- rbind(M_star_GMS_I, M_star_GMS_E)
  T_n_star <- colSums(M_star_GMS^2)
  return(mean(T_n_star > T_n))

}
