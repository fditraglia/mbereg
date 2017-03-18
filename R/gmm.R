gmm_a0_zero <- function(dat, params, gradient = FALSE){
  # Unpack data
  y <- dat$y
  z <- dat$z
  Tobs <- dat$Tobs

  # Unpack parameters
  a1 <- params[1]
  b <- params[2]
  q <- params[3]
  p <- params[4]
  mu <- params[5]
  s <- params[6]
  r <- params[7]

  # Calculate f_n
  W <- cbind(z, Tobs, y, y^2, y * Tobs)
  colnames(W) <- c('z', 'T', 'y', 'ysq', 'yT')
  h_n <- colMeans(W) - c(q, p, mu, s, r)
  cov_z_y <- mean(z * y) - q * mu
  cov_z_T <- mean(z * Tobs) - q * p
  cov_z_ysq <- mean(z * y^2) - q * s
  cov_z_yT <- mean(z * y * Tobs) - q * r
  g_n <- c(cov_z_y - (b / (1 - a1)) * cov_z_T,
          cov_z_ysq - 2 * (b / (1 - a1)) * cov_z_yT + (b^2 / (1 - a1)) * cov_z_T)
  f_n <- c(g_n, h_n)

  if(gradient){
  # Calculate F_n and GMM gradient
  D_g1 <- c(-b * cov_z_T / ((1 - a1)^2),
            -cov_z_T / (1 - a1),
            p * b / (1 - a1) - mu,
            q * b / (1 - a1),
            -q,
            0,
            0)
  D_g2 <- c((b^2 * cov_z_T - 2 * b * cov_z_yT) / ((1 - a1)^2),
            (2 * b * cov_z_T - 2 * cov_z_yT) / (1 - a1),
            2 * r * b / (1 - a1) - p * b^2 / (1 - a1) - s,
            -q * b^2 / (1 - a1),
            0,
            -q,
            2 * q * b / (1 - a1))
  F_n <- rbind(D_g1,
               D_g2,
               cbind(matrix(0, 5, 2), -diag(5)))

  gradient <- drop(2 * crossprod(F_n, f_n))
  names(gradient) <- NULL
  return(gradient)
  } else {
    return(sum(f_n^2))
  }
}


solve_gmm_a0_zero <- function(dat){

  # Functions holding data *fixed* to pass to constrOptim
  criterion <- function(theta) gmm_a0_zero(dat, theta, gradient = FALSE)
  gradient <- function(theta) gmm_a0_zero(dat, theta, gradient = TRUE)

  # Unpack data to calculate starting values
  y <- dat$y
  Tobs <- dat$Tobs
  z <- dat$z

  # Starting value for a1 between zero and "weak" upper bound
  p0 <- mean(Tobs[z == 0])
  p1 <- mean(Tobs[z == 1])
  a1_upper <- 1 - max(p0, p1)
  a1_start <- 0.5 * a1_upper

  # Starting value for beta between IV (wald) and Reduced Form (rf)
  wald <- cov(y, Tobs) / cov(z, Tobs)
  rf <- cov(y, z) / var(z)
  b_start <- 0.5 * (wald + rf)

  # Starting values for "nuisance" parameters
  q_start <- mean(z)
  p_start <- mean(Tobs)
  mu_start <- mean(y)
  s_start <- mean(y^2)
  r_start <- mean(y * Tobs)
  start_vals <- c(a1_start, b_start, q_start, p_start, mu_start, s_start, r_start)

  # Set up linear inequality constraints: M %*% params - v >= 0
  M <- matrix(0, 7, 7)
  M[1:2, 1] <- c(1, -1)
  M[3:4, 3] <- c(1, -1)
  M[5:6, 4] <- c(1, -1)
  M[7, 6] <- 1
  v <- c(0, -a1_upper, 0, -1, 0, -1, 0)

  out <- constrOptim(theta = start_vals, f = criterion, grad = gradient,
                     ui = M, ci = v)
  names(out$par) <- c('a1', 'b', 'q', 'p', 'mu', 's', 'r')

  # Also return the "naive" GMM estimates
  s_T_z <- cov(Tobs, z)
  s_y_z <- cov(y, z)
  s_y2_z <- cov(y^2, z)
  s_yT_z <- cov(y * Tobs, z)
  b_naive <- 2 * s_yT_z / s_T_z - s_y2_z / s_y_z
  a1_naive <- 1 - 2 * s_yT_z / s_y_z + s_y2_z * s_T_z / s_y_z^2
  naive <- c(a1 = a1_naive, b = b_naive, q = q_start, p = p_start, mu = mu_start,
             s = s_start, r = r_start)
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










