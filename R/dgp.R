dgp <- function(a0, a1, b = 1, n = 1000, d = 0.15, rho = 0.5, cc = 0){
  n_treat <- ceiling(n/2)
  n_control <- n - n_treat
  z <- c(rep(0, n_control), rep(1, n_treat)) # offer of treatment
  errors <- mvtnorm::rmvnorm(n, sigma = matrix(c(1, rho, rho, 1), 2, 2,
                                               byrow = TRUE))
  g0 <- qnorm(d)
  g1 <- qnorm(1 - d) - qnorm(d)
  Tstar <- as.numeric(g0 + g1 * z + errors[,2] > 0) #select into treatment
  y <- cc + b * Tstar + errors[,1]
  #mis-classification
  Tobs <- (1 - Tstar) * rbinom(n, 1, a0) + Tstar * rbinom(n, 1, 1 - a1)
  return(data.frame(Tstar, Tobs, y, z))
}
