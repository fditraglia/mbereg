# This test substitutes a preliminary estimator of theta1 and uses only the
# first moment inequalities.
GMS_test_alphas <- function(a0, a1, dat, normal_sims){

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

  # Moment functions
  h <- cbind(w_centered %*% Psi , w_centered_z %*% psi1)
  mI <- cbind(Tobs * (1 - z) / (1 - q) - a0,
              (1 - Tobs) * (1 - z) / (1 - q) - a1,
              Tobs * z / q - a0,
              (1 - Tobs) * z / q - a1)
  mE <- w_centered_z %*% Psi[,2:3]
  m <- cbind(mI, mE)

  # Derivative matrices to adjust variance for preliminary estimation
  MI <- matrix(0, nrow = 4, 4)
  ME_kappa <- matrix(c(0, -q, 0,
                       0, 0, -q), nrow = 2, ncol = 3, byrow = TRUE)
  ME_theta1 <- crossprod(D_Psi[,2:3], wz_bar)
  ME <- cbind(ME_kappa, ME_theta1)
  M <- rbind(MI, ME)

  H_kappa <- rbind(-1 * diag(3),
                   c(-q, 0, 0))
  H_theta1 <- rbind(crossprod(D_Psi, w_bar),
                    crossprod(D_psi1, wz_bar))
  H <- cbind(H_kappa, H_theta1)

  # Adjust the variance matrix
  V <- var(cbind(m, h))
  B <- -1 * M %*% solve(H)
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
