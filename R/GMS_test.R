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


#----------------------------------------------------------------------------
# GMS test for our example:
#----------------------------------------------------------------------------
#     No restrictions on a0 or a1, "weak" and second moment bounds, and
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
  Tobs <- dat$Tobs
  y <- dat$y
  z <- dat$z
  q <- mean(z) # treat this as fixed in repeated sampling
  n <- nrow(dat)

  theta1 <- beta / (1 - a0 - a1)
  theta2 <- theta1^2 * (1 + (a0 - a1))
  theta3 <- theta1^3 * ((1 - a0 - a1)^2 + 6 * a0 * (1 - a1))

  # Moment functions for preliminary estimators
  h_p0 <- Tobs * (1 - z) / (1 - q)
  h_p1 <- Tobs * z / q
  h_nu_00 <- y * (1 - Tobs) * (1 - z) / sqrt(1 - q)
  h_nu_10 <- y * Tobs * (1 - z) / sqrt(1 - q)
  h_nu_01 <- y * (1 - Tobs) * z / sqrt(q)
  h_nu_11 <- y * Tobs * z / sqrt(q)
  h_k1 <- y - theta1 * Tobs
  h_k2 <- y^2 - theta1 * 2 * y * Tobs + theta2 * Tobs
  h_k3 <- y^3 - theta1 * 3 * y^2 * Tobs + theta2 * 3 * y * Tobs - theta3 * Tobs
  h <- cbind(h_p0, h_p1, h_nu_00, h_nu_10, h_nu_01, h_nu_11, h_k1, h_k2, h_k3)

  # Preliminary estimators
  p0 <- mean(h_p0)
  p1 <- mean(h_p1)
  nu_00 <- mean(h_nu_00)
  nu_10 <- mean(h_nu_10)
  nu_01 <- mean(h_nu_01)
  nu_11 <- mean(h_nu_11)
  k1 <- mean(h_k1)
  k2 <- mean(h_k2)
  k3 <- mean(h_k3)
  gamma <- c(p0, p1, nu_00, nu_10, nu_01, nu_11, k1, k2, k3)

  # Make sure I've ordered the parameters consistently
  stopifnot(all.equal(gamma, colMeans(h), check.attributes = FALSE))

  # Equality moment conditions
  u1 <- h_k1 - k1
  u2 <- h_k2 - k2
  u3 <- h_k3 - k3
  mE <- z * cbind(u1, u2, u3)

  # First-moment Inequalities
  mI_1 <- cbind(rep(p0 - a0, n),
                rep(1 - p0 - a1, n),
                rep(p1 - a0, n),
                rep(1 - p1 - a1, n))

  # Second-moment Inequalities
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

  # Full set of Inequalities
  mI <- cbind(mI_1, mI_2)

  # Full set of moment conditions
  m <- cbind(mI, mE)

  # Expected derivative matrix to account for preliminary estimation
  B_I1_p <- cbind(c(1, -1, 0, 0),
                  c(0, 0, 1, -1))
  B_I1 <- cbind(B_I1_p, matrix(0, 4, 4), matrix(0, 4, 3))
  B_I2_p <- matrix(c( mean(y2_a0_z0), 0,
                     -mean(y2_a1_z0), 0,
                      0, mean(y2_a0_z1),
                      0, -mean(y2_a1_z1)),
                   4, 2, byrow = TRUE)
  B_I2_nu <- matrix(c(2 * a0 * d_a0_z0, -2 * (1 - a0) * d_a0_z0, 0, 0,
                      2 * (1 - a1) * d_a1_z0, -2 * a1 * d_a1_z0, 0, 0,
                      0, 0, 2 * a0 * d_a0_z1, -2 * (1 - a0) * d_a0_z1,
                      0, 0, 2 * (1 - a1) * d_a1_z1, -2 * a1 * d_a1_z1),
                    4, 4, byrow = TRUE)
  B_I2 <- cbind(B_I2_p, B_I2_nu, matrix(0, 4, 3))
  B_E <- cbind(matrix(0, 3, 2), matrix(0, 3, 4), -q * diag(3))
  B <- rbind(B_I1, B_I2, B_E)

  # Estimate of asymptotic variance matrix of m
  A <- cbind(diag(nrow(B)), B)
  Sigma_n <- A %*% var(cbind(m, h)) %*% t(A)

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
  T_n_star <- apply(M_star, 2, function(x) get_test_stat(x + phi, Omega_n, n_ineq))
  return(mean(T_n_star > T_n))
}


