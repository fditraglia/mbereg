# ----------------------------------------------------------------------------
# Simple example of the GMS test from Andrews & Soares (2010)
# ----------------------------------------------------------------------------
# The differences between this implementation and the one in GMS_test.R are
# as follows:
# (1) We use preliminary estimates for the intercepts and have to adjust the
#     asymptotic variance matrix esimator accordingly.
# (2) We use the "asymptotic" version of the GMS test, which uses independent
#     standard normal draws in place of bootstrapping.
# ----------------------------------------------------------------------------

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
# Simple Example of GMS test for our example:
#----------------------------------------------------------------------------
#     Assumes that a0 = 0, uses only first and second moment information and
#     includes only the "weak" bounds for a1. Uses preliminary estimators for
#     the intercepts kappa1 and kappa2 and adjusts the variance matrix of the
#     GMS test accordingly. Based on the "asymptotic" version of GMS, i.e.
#     the version that makes standard normal draws to simulate the limit
#     experiment rather than bootstrapping the sample.
#----------------------------------------------------------------------------
# - a1 is the null hypothesis for $\alpha_1$
# - beta is the null hypothesis for $
# - dat is a dataframe containing the data, with columns y, z, and Tobs
# - normal_sims is a (4 \times R) matrix of standard normal random draws
#----------------------------------------------------------------------------
GMS_test2 <- function(a1, beta, dat, normal_sims){
  Tobs <- dat$Tobs
  y <- dat$y
  z <- dat$z
  q <- mean(z)
  n <- nrow(dat)
  theta1 <- beta / (1 - a1)
  theta2 <- beta^2 / (1 - a1)
  h1 <- y - theta1 * Tobs
  h2 <- y^2 - theta1 * 2 * y * Tobs + theta2 * Tobs
  h <- cbind(h1, h2)
  k1_hat <- mean(h1)
  k2_hat <- mean(h2)
  m1 <- (1 - a1) - Tobs * (1 - z) / (1 - q)
  m2 <- (1 - a1) - Tobs * z / q
  m3 <- (y - k1_hat - theta1 * Tobs) * z
  m4 <- (y^2 - k2_hat - theta1 * 2 * y * Tobs + theta2 * Tobs) * z
  m <- cbind(m1, m2, m3, m4)
  A <- cbind(diag(1, 4), -q * rbind(matrix(0, 2, 2), diag(1, 2)))
  Sigma_n <- A %*% var(cbind(m, h)) %*% t(A)
  nu <- sqrt(n) * colMeans(m)
  T_n <- get_test_stat(nu, Sigma_n, 2)
  s_n <- sqrt(diag(Sigma_n))
  phi <- c(ifelse((nu[1:2] / s_n[1:2]) > sqrt(log(n)), Inf, 0), rep(0, 2))
  Omega_n <- cov2cor(Sigma_n)
  M_star <- chol(Omega_n) %*% normal_sims
  T_n_star <- apply(M_star, 2, function(x) get_test_stat(x + phi, Omega_n, 2))
  return(mean(T_n_star > T_n))
}


set.seed(19295)
R <- 5000
sims <- matrix(rnorm(4 * R), 4, R)

check_GMS_test2 <- function(a1, beta){
  dat <- dgp(a0 = 0, a1 = a1, b = beta, n = 1000)
  GMS_test2(a1, beta, dat, sims)
}

p_values1 <- replicate(1000, check_GMS_test2(a1 = 0.0, beta = 1))
p_values2 <- replicate(1000, check_GMS_test2(a1 = 0.1, beta = 1))
p_values3 <- replicate(1000, check_GMS_test2(a1 = 0.2, beta = 1))
p_values4 <- replicate(1000, check_GMS_test2(a1 = 0.0, beta = 2))
p_values5 <- replicate(1000, check_GMS_test2(a1 = 0.1, beta = 2))
p_values6 <- replicate(1000, check_GMS_test2(a1 = 0.2, beta = 2))

