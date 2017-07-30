#----------------------------------------------------------------------------
# GMS test for our example but using only first moment inequalities
#----------------------------------------------------------------------------
#     No restrictions on a0 or a1, "weak" moment inequalities, along with
#     preliminary estimators of strongly identified parameters with appropriate
#     adjustment to the GMS variance matrix estimator. Based on the "asymptotic"
#     version of GMS, i.e. the version that makes standard normal draws to
#     simulate the limit experiment rather than bootstrapping the sample.
#----------------------------------------------------------------------------
# - a0 is the null hypothesis for $\alpha_0$
# - a1 is the null hypothesis for $\alpha_1$
# - beta is the null hypothesis for $
# - dat is a dataframe containing the data, with columns y, z, and Tobs
# - normal_sims is a (??? \times R) matrix of standard normal random draws
#----------------------------------------------------------------------------
GMS_test <- function(a0, a1, beta, dat, normal_sims){

  if((a0 > 0.99) || (a0 < 0)) return(0)
  if((a1 > 0.99) || (a1 < 0)) return(0)
  if(a0 + a1 > 0.99) return(0)

  Tobs <- dat$Tobs
  y <- dat$y
  z <- dat$z
  q <- mean(z) # treat this as fixed in repeated sampling
  n <- nrow(dat)
  kappa_n <- sqrt(log(n))

  theta1 <- beta / (1 - a0 - a1)
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

  # Moment functions
  h <- w_centered %*% Psi
  mI <- cbind(Tobs * (1 - z) / (1 - q) - a0,
              (1 - Tobs) * (1 - z) / (1 - q) - a1,
              Tobs * z / q - a0,
              (1 - Tobs) * z / q - a1)
  mE <- w_centered_z %*% Psi
  m <- cbind(mI, mE)

  # Derivative matrices to adjust variance for preliminary estimation
  MI <- matrix(0, nrow = 4, ncol = 3)
  ME <- -q * diag(3)
  M <- rbind(MI, ME)
  H <- -diag(3)

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


