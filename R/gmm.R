gmm_a0_zero <- function(dat, params, gradient = FALSE){
  # Since we assume a0 = 0, denote a1 by a in this function

  # Unpack data
  y <- dat$y
  z <- dat$z
  Tobs <- dat$Tobs

  # Condition on P(z = 1)
  q <- mean(z)

  # Unpack parameters
  a <- params[1]
  b <- params[2]
  cc <- params[3]
  s_ee <- params[4]
  #lam1 <- params[5]
  #lam2 <- params[6]

  # Construct errors u(theta) and v(theta)
  u <- y - cc - (b / (1 - a)) * Tobs
  v <- y^2 - s_ee - cc^2 - (b / (1 - a)) * 2 * y * Tobs + (b^2 / (1 - a)) * Tobs

  # Construct sample analogues of moment conditions
  g1_n <- c(mean(u), mean(v))
  g2_n <- c(mean(u * z), mean(v * z))
  #h1_n <- (1 - a) - mean(Tobs * (1 - z)) / (1 - q)
  #h2_n <- (1 - a) - mean(Tobs * z) / q
  #psi_n <- c(g1_n, g2_n, h1_n - lam1, h2_n - lam2)
  psi_n <- c(g1_n, g2_n)

  if(gradient){
  # Calculate gradient
  k1_bar <- c(mean(Tobs),
              mean(2 * Tobs * y - b * Tobs),
              mean(Tobs * z),
              mean(2 * Tobs * y * z - b * Tobs * z))
  k2_bar <- c(mean(Tobs),
              mean(2 * Tobs * y - 2 * b * Tobs),
              mean(Tobs * z),
              mean(2 * Tobs * y * z - 2 * b * Tobs * z))
  R1_bar <- c(-1, -2 * cc, -mean(z), -2 * cc * mean(z))
  R2_bar <- c(0, -1, 0, -mean(z))

  G_n <- cbind((-b / ((1 - a)^2)) * k1_bar,
               (-1 / (1 - a)) * k2_bar,
               R1_bar, R2_bar)
  #H_n <- matrix(c(rep(-1, 2), rep(0, 6)), 2, 4)
  #F_n <- cbind(rbind(G_n, H_n))
  #Psi_n <- cbind(F_n, rbind(matrix(0, 4, 2), diag(nrow = 2, ncol = 2)))
  Psi_n <- G_n

  gradient <- drop(crossprod(Psi_n, psi_n)) # GMM criterion with 1/2 scaling
  names(gradient) <- NULL
  return(gradient)
  } else {
    return(0.5 * sum(psi_n^2)) # GMM criterion with 1/2 scaling
  }
}


solve_gmm_a0_zero <- function(dat){

  # Functions holding data *fixed* to pass to constrOptim
  #     gamma = (theta, lambda) is the full parameter vector
  criterion <- function(gamma) gmm_a0_zero(dat, gamma, gradient = FALSE)
  gradient <- function(gamma) gmm_a0_zero(dat, gamma, gradient = TRUE)

  # Unpack data to calculate starting values
  y <- dat$y
  Tobs <- dat$Tobs
  z <- dat$z

  # Starting value for a1 between zero and "weak" upper bound
  p0 <- mean(Tobs[z == 0])
  p1 <- mean(Tobs[z == 1])
  a1_upper <- 1 - max(p0, p1)
  a1_start <- 0.5 * a1_upper

  # Starting values for lambda1 and lambda2 based on a1_start
  #lam1_start <- 1 - p0 - a1_start
  #lam2_start <- 1 - p1 - a1_start

  # Starting value for beta between IV (wald) and Reduced Form (rf)
  wald <- cov(y, Tobs) / cov(z, Tobs)
  rf <- cov(y, z) / var(z)
  b_start <- 0.5 * (wald + rf)

  # Starting value for c and s_ee based on b_start
  c_start <- mean(y) - b_start * mean(Tobs)
  s_ee_start <- var(y - b_start * Tobs)

  # Collect all starting values
  #start_vals <- c(a1_start, b_start, c_start, s_ee_start, lam1_start, lam2_start)
  start_vals <- c(a1_start, b_start, c_start, s_ee_start)

  # Set up linear inequality constraints: M %*% params - v >= 0
  #M <- matrix(0, 4, 6)
  #M[1, 1] <- M[2, 4] <- M[3, 5] <- M[4, 6] <- 1
  #v <- rep(0, 4)
  M <- matrix(0, 2, 4)
  M[1, 1] <- M[2, 4] <- 1
  v <- rep(0, 2)

  out <- constrOptim(theta = start_vals, f = criterion, grad = NULL,#grad = gradient,
                     ui = M, ci = v)
  #names(out$par) <- c('a1', 'b', 'c', 's_ee', 'lam1', 'lam2')
  names(out$par) <- c('a1', 'b', 'c', 's_ee')

  # Also return the "naive" GMM estimates
  s_T_z <- cov(Tobs, z)
  s_y_z <- cov(y, z)
  s_y2_z <- cov(y^2, z)
  s_yT_z <- cov(y * Tobs, z)
  b_naive <- 2 * s_yT_z / s_T_z - s_y2_z / s_y_z
  a1_naive <- 1 - 2 * s_yT_z / s_y_z + s_y2_z * s_T_z / s_y_z^2
  naive <- c(a1 = a1_naive, b = b_naive)
  out$naive <- naive
  out$wald <- wald
  out$rf <- rf
  out$a1_upper <- a1_upper
  return(out)
}

sim_draw_a0_zero <- function(b, a1 = 0.1, n = 1000, d = 0.15, rho = 0.5){
  sim_dat <- dgp(a0 = 0, a1 = a1, b = b, n = n, d = d, rho = rho)
  results <- solve_gmm_a0_zero(sim_dat)
  out <- c(a1_naive = results$naive[1],
           a1_constr = results$par[1],
           b_naive = results$naive[2],
           b_constr = results$par[2],
           wald = results$wald,
           rf = results$rf,
           a1_upper = results$a1_upper)
  row.names(out) <- NULL
  return(out)
}

sim_study_a0_zero <- function(n_reps = 100, b, a1 = 0.1, n = 1000, d = 0.15,
                              rho = 0.5){
  out <- replicate(n_reps, sim_draw_a0_zero(b = b, a1 = a1, n = n, d = d, rho = rho))
  out <- as.data.frame(t(out))
  names(out) <- c('a1_naive', 'a1_constr', 'b_naive', 'b_constr', 'wald', 'rf',
                  'a1_upper')
  return(out)
}










