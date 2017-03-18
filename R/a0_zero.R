get_m_star <- function(d = 0.15, r = 0.5, c = 0){
  g0 <- qnorm(d)
  g1 <- qnorm(1 - d) - qnorm(d)

  f_U <- function(x, a) x * dnorm(x) * pnorm((a - r * x) / sqrt(1 - r^2))
  mstar00 <- c + integrate(function(x) f_U(x, -g0), -Inf, Inf)$value /
    pnorm(-g0)
  mstar01 <- c + integrate(function(x) f_U(x, -(g0 + g1)), -Inf, Inf)$value /
    pnorm(-(g0 + g1))

  f_L <- function(x, a) x * dnorm(x) * (1 - pnorm((a - r * x) / sqrt(1 - r^2)))
  mstar10 <- c + integrate(function(x) f_L(x, -g0), -Inf, Inf)$value /
    (1 - pnorm(-g0))
  mstar11 <- c + integrate(function(x) f_L(x, -(g0 + g1)), -Inf, Inf)$value /
    (1 - pnorm(-(g0 + g1)))
  out <- c(mstar00, mstar01, mstar10, mstar11)
  names(out) <- c('mstar00', 'mstar01', 'mstar10', 'mstar11')
  return(out)
}

# This is just a Monte Carlo version used to check the preceding function
get_m_star_sim <- function(d = 0.15, r = 0.5, c = 0, n_sims = 1000000){
  g0 <- qnorm(d)
  g1 <- qnorm(1 - d) - qnorm(d)
  S <- matrix(c(1, r, r, 1), 2, 2, byrow = TRUE)
  draws <- mvtnorm::rmvnorm(n_sims, sigma = S)
  epsilon <- draws[, 1]
  xi <- draws[, 2]
  mstar00 <- c + mean(epsilon[xi <= -g0])
  mstar01 <- c + mean(epsilon[xi <= -(g0 + g1)])
  mstar10 <- c + mean(epsilon[xi > -g0])
  mstar11 <- c + mean(epsilon[xi > -(g0 + g1)])
  out <- c(mstar00, mstar01, mstar10, mstar11)
  names(out) <- c('mstar00', 'mstar01', 'mstar10', 'mstar11')
  return(out)
}

# Check the theoretical calculation for the asymptotic variance of the reduced
# form estimator for the special case of our simulation DGP, in which the
# non-compliance is symmetric and q is fixed at 0.5
get_AVAR_rf <- function(b, d = 0.15, r = 0.5, c = 0){
  mstar <- get_m_star(d, r, c)
  m_diff <- mstar[4] - mstar[2]
  out <- 4 * (1 + d * (1 - d) * b * (b + 2 * m_diff))
  return(out)
}
get_AVAR_rf_sim <- function(b, d = 0.15, r = 0.5){
  sims <- dgp(a0 = 0, a1 = 0.1, b, 1000000, d, r)
  y_z0 <- sims$y[sims$z == 0]
  y_z1 <- sims$y[sims$z == 1]
  eta <- c(y_z0 - mean(y_z0), y_z1 - mean(y_z1))
  return(4 * var(eta))
}


get_V_true <- function(a1, b, d = 0.15, n_sims = 1000000){

  dat <- dgp(a0 = 0, a1, b, n_sims, d)
  n <- nrow(dat)
  q <- 0.5
  p1_star <- 1 - d
  p0_star <- d
  p_star <- p1_star * q + p0_star * (1 - q)
  p <- (1 - a1) * p_star
  mu <- b * p_star
  s <- with(dat, mean(y^2))
  r <- with(dat, mean(y * Tobs))

  s_T_z <- (1 - a1) * (p1_star - p0_star) * q * (1 - q)
  s_y_z <- b * (p1_star - p0_star) * q * (1 - q)
  s_y2_z <- with(dat, cov(y^2, z))
  s_yT_z <- with(dat, cov(y * Tobs, z))

  G_theta <- matrix(c(-s_T_z / (1 - a1),
                      -b * s_T_z / (1 - a1)^2,
                      2 * (b * s_T_z - s_yT_z) / (1 - a1),
                      (b^2 * s_T_z - 2 * b * s_yT_z) / (1 - a1)^2),
                    nrow = 2, ncol = 2, byrow = TRUE)

  G_gamma <- matrix(c(p * b / (1 - a1) - mu,
                      q * b / (1 - a1),
                      -q, 0, 0,
                      (b / (1 - a1)) * (2 * r - b * p) - s,
                      -q * b^2 / (1 - a1) ,
                      0, -q,
                      2 * b * q / (1 - a1)), nrow = 2, ncol = 5, byrow = TRUE)

  G_theta_inv <- solve(G_theta)
  FF_g <- cbind(G_theta_inv, G_theta_inv %*% G_gamma)
  FF_h <- cbind(matrix(0, nrow = 5, ncol = 2), diag(-1, 5))
  FF <- rbind(FF_g, FF_h)

  g1 <- with(dat, (z * y - q * mu) - (b / (1 - a1)) * (z * Tobs - q * p))
  g2 <- with(dat, (z * y^2 - q * s) -
               2 * (b / (1 - a1)) * (z * y * Tobs - q * r) +
               (b^2 / (1 - a1)) * (z * Tobs - q * p))
  h1 <- with(dat, z - q)
  h2 <- with(dat, Tobs - p)
  h3 <- with(dat, y - mu)
  h4 <- with(dat, y^2 - s)
  h5 <- with(dat, y * Tobs - r)
  f_hat <- cbind(g1, g2, h1, h2, h3, h4, h5)
  Omega <- crossprod(f_hat) / n

  V <- FF %*% Omega %*% t(FF)
  V_theta <- V[1:2, 1:2]
  return(V_theta)
}


est <- function(dat){
  n <- nrow(dat)
  q <- with(dat, mean(z))
  p <- with(dat, mean(Tobs))
  mu <- with(dat, mean(y))
  s <- with(dat, mean(y^2))
  r <- with(dat, mean(y * Tobs))

  # Inference for Reduced Form slope coefficient
  RF <- with(dat, mean(y[z == 1]) - mean(y[z == 0]))
  eta <- with(dat, y - z * mean(y[z == 1]) - (1 - z) * mean(y[z == 0]))
  s2_eta <- var(eta)
  E_z_eta2 <- mean(dat$z * eta^2)
  RF_AVAR <- s2_eta / (1 - q)^2 + ((1 - 2 * q) / (q^2 * (1 - q^2))) * E_z_eta2
  RF_SE <- sqrt(RF_AVAR / n)

  # Bounds for a1
  p0 <- with(dat, mean(Tobs[z == 0]))
  p1 <- with(dat, mean(Tobs[z == 1]))
  a1_upper <- min(1 - p0, 1 - p1)

  # Inference for IV slope

  s_T_z <- with(dat, cov(Tobs, z))
  s_y_z <- with(dat, cov(y, z))
  s_y2_z <- with(dat, cov(y^2, z))
  s_yT_z <- with(dat, cov(y * Tobs, z))

  b <- 2 * s_yT_z / s_T_z - s_y2_z / s_y_z
  a1 <- 1 - 2 * s_yT_z / s_y_z + s_y2_z * s_T_z / s_y_z^2

  G_theta <- matrix(c(-s_T_z / (1 - a1),
                      -b * s_T_z / (1 - a1)^2,
                      2 * (b * s_T_z - s_yT_z) / (1 - a1),
                      (b^2 * s_T_z - 2 * b * s_yT_z) / (1 - a1)^2),
                    nrow = 2, ncol = 2, byrow = TRUE)

  G_gamma <- matrix(c(p * b / (1 - a1) - mu,
                      q * b / (1 - a1),
                      -q, 0, 0,
                      (b / (1 - a1)) * (2 * r - b * p) - s,
                      -q * b^2 / (1 - a1) ,
                      0, -q,
                      2 * b * q / (1 - a1)), nrow = 2, ncol = 5, byrow = TRUE)

  G_theta_inv <- solve(G_theta)
  FF_g <- cbind(G_theta_inv, G_theta_inv %*% G_gamma)
  FF_h <- cbind(matrix(0, nrow = 5, ncol = 2), diag(-1, 5))
  FF <- rbind(FF_g, FF_h)

  g1 <- with(dat, (z * y - q * mu) - (b / (1 - a1)) * (z * Tobs - q * p))
  g2 <- with(dat, (z * y^2 - q * s) -
               2 * (b / (1 - a1)) * (z * y * Tobs - q * r) +
               (b^2 / (1 - a1)) * (z * Tobs - q * p))
  h1 <- with(dat, z - q)
  h2 <- with(dat, Tobs - p)
  h3 <- with(dat, y - mu)
  h4 <- with(dat, y^2 - s)
  h5 <- with(dat, y * Tobs - r)
  f_hat <- cbind(g1, g2, h1, h2, h3, h4, h5)
  Omega <- crossprod(f_hat) / n

  V <- FF %*% Omega %*% t(FF)
  V_theta <- V[1:2, 1:2]
  b_SE <- sqrt(V_theta[1, 1] / n)
  a1_SE <- sqrt(V_theta[2, 2] / n)
  return(c('b' = b, 'b_SE' = b_SE, 'a1' = a1, 'a1_SE' = a1_SE,
           'RF' = RF, 'RF_SE' = RF_SE, 'a1_upper' = a1_upper))
}

plot_dist <- function(a1, b, n, d = 0.15, rho = 0.5){
  V_true <- get_V_true(a1, b, d)
  b_SE_true <- sqrt(V_true[1, 1] / n)
  a1_SE_true <- sqrt(V_true[2, 2] / n)
  RF_SE_true <- sqrt(get_AVAR_rf(b, d, rho) / n)
  sim_draws <- as.data.frame(t(replicate(10000, est(dgp(a0 = 0, a1, b, n, d)))))
  par(mfrow = c(1, 3))
  par(oma = c(0, 0, 2, 0))

  MASS::truehist(sim_draws$b, col = 'lightskyblue', xlab = expression(beta))
  b_bias_text <- paste0('Bias = ', round(mean(sim_draws$b) - b, 3))
  b_SD_text <- paste0('SD = ', round(sd(sim_draws$b), 3))
  title(paste(b_bias_text, ',', b_SD_text), font.main = 1, cex.main = 1)
  b_seq <- seq(from = min(sim_draws$b), to = max(sim_draws$b), length.out = 1000)
  b_true_density <- dnorm(b_seq, mean = b, sd = b_SE_true)
  points(b_seq, b_true_density, type = 'l', lwd = 3, lty = 1)
  abline(v = b, lty = 1, lwd = 3, col = 'firebrick')

  MASS::truehist(sim_draws$a1, col = 'lightskyblue', xlab = expression(alpha[1]))
  a1_bias_text <- paste0('Bias = ', round(mean(sim_draws$a1) - a1, 3))
  a1_SD_text <- paste0('SD = ', round(sd(sim_draws$a1), 3))
  title(paste(a1_bias_text, ',', a1_SD_text), font.main = 1, cex.main = 1)
  a1_seq <- seq(from = min(sim_draws$a1), to = max(sim_draws$a1), length.out = 1000)
  a1_true_density <- dnorm(a1_seq, mean = a1, sd = a1_SE_true)
  points(a1_seq, a1_true_density, type = 'l', lwd = 3, lty = 1)
  abline(v = a1, lty = 1, lwd = 3, col = 'firebrick')

  MASS::truehist(sim_draws$RF, col = 'lightskyblue',
                 xlab = expression(beta * (1 - 2 * delta)))
  RF_bias_text <- paste0('Bias = ', round(mean(sim_draws$RF) - b * (1 - 2 * d), 3))
  RF_SD_text <- paste0('SD = ', round(sd(sim_draws$RF), 3))
  title(paste(RF_bias_text, ',', RF_SD_text), font.main = 1, cex.main = 1)
  RF_seq <- seq(from = min(sim_draws$RF), to = max(sim_draws$RF), length.out = 1000)
  RF_true_density <- dnorm(RF_seq, mean = b * (1 - 2 * d), sd = RF_SE_true)
  points(RF_seq, RF_true_density, type = 'l', lwd = 3, lty = 1)
  abline(v = b * (1 - 2 * d), lty = 1, lwd = 3, col = 'firebrick')

  mytitle <- substitute(paste(beta, ' = ', my_beta, ', ', alpha[1],
                               ' = ', my_alpha1, ', ',
                              delta, ' = ', my_delta,
                              ', ', n, ' = ', my_n),
                        list(my_beta = b, my_alpha1 = a1, my_n = n,
                             my_delta = d))


  title(main = mytitle, outer = T)
  par(mfrow = c(1, 1))
  par(oma = c(0, 0, 0, 0))
}

CI_sim <- function(a1, b, n, d, rho = 0.5, n_sims = 5000){
  sims <- replicate(n_sims, est(dgp(a0 = 0, a1, b, n, d, rho)))
  sims <- as.data.frame(t(sims))

  RF_lower <- with(sims, RF - qnorm(0.975) * RF_SE)
  RF_upper <- with(sims, RF + qnorm(0.975) * RF_SE)
  RF_true <- b * (1 - 2 * d)
  RF_bias <- with(sims, mean(RF - RF_true))
  RF_cover <- mean((RF_lower <= RF_true) & (RF_upper >= RF_true))
  RF_width <- median(RF_upper - RF_lower)

  # Check coverage of CIs
  #   GMM
  b_lower <- with(sims, b - qnorm(0.975) * b_SE)
  b_upper <- with(sims, b + qnorm(0.975) * b_SE)
  b_bias <- median(sims$b - b)
  GMM_cover <- mean((b_lower <= b) & (b_upper >= b))
  GMM_width <- median(b_upper - b)
  return(c('RF_cover' = RF_cover, 'GMM_cover' = GMM_cover,
           'RF_width' = RF_width, 'GMM_width' = GMM_width,
           'RF_bias' = RF_bias, 'b_bias' = b_bias))
}




