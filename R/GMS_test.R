# Simple example of the GMS test from Andrews & Soares (2010) that imposes
# c = 0 and s_ee = 1, as is the case in our baseline simulation DGP.

# Calculate the "MMM" test stat - "S1" from Andrews and Soares
get_S1 <- function(m, Sigma, p){
  m[1:p] <- ifelse(m[1:p] < 0, m[1:p], 0)
  sum(m^2 / diag(Sigma))
}

# Calculate sample inequality and equality moment conditions under H0
get_m <- function(dat, a, b, q){
  y <- dat$y
  z <- dat$z
  Tobs <- dat$Tobs
  # Construct errors u(theta) and v(theta)
  u <- y - (b / (1 - a)) * Tobs
  v <- y^2 - 1 - (b / (1 - a)) * 2 * y * Tobs + (b^2 / (1 - a)) * Tobs

  # Construct sample analogues of moment conditions
  #     -- each row is an individual
  out <- cbind((1 - a) - Tobs * (1 - z) / (1 - q),
               (1 - a) - Tobs * z / q,
               u * z,
               v * z)
  return(out)
}

# Calculate sample analogue based on get_m
 get_m_bar <- function(dat, a, b, q){
   m <- get_m(dat, a, b, q)
   return(colMeans(m))
 }

# Estimate variance of moment conditions under H0
get_Sigma_hat <- function(dat, a, b, q){
  m <- get_m(dat, a, b, q)
  n <- nrow(m)
  return((n - 1) * var(m) / n) # var in R divides by (n - 1)
}


GMS_test <- function(dat, a_null, b_null, n_boot = 5000, get_S = get_S1){
  # Impose the null hypothesis and fix q
  q <- mean(dat$z)
  get_m_bar0 <- function(mydat) get_m_bar(mydat, a_null, b_null, q)
  get_Sigma_hat0 <- function(mydat) get_Sigma_hat(mydat, a_null, b_null, q)

  # Caclulate "true sample" test statistic (i.e. non-bootstrap)
  m_bar <- get_m_bar0(dat)
  Sigma_hat <- get_Sigma_hat0(dat)
  n <- nrow(dat)
  Tn <- get_S(sqrt(n) * m_bar, Sigma_hat, p = 2)

  # Determine which moment inequalities are far from binding (FB == TRUE)
  tn <- sqrt(n) * m_bar[1:2] / sqrt(diag(Sigma_hat)[1:2])
  FB <- tn > sqrt(log(n))

  keep_me <- which(c(!FB, rep(TRUE, 2))) # Indices of which moment conditions to keep
  p_GMS <- sum(!FB) # How many moment inequalities are *not* dropped?

  # Determine the bootstrap critical value
  Tn_boot <- rep(NA_real_, n_boot)
  for(r in 1:n_boot){

    dat_boot <- dat[sample(1:n, n, TRUE),]
    m_bar_star <- get_m_bar0(dat_boot)
    Sigma_hat_star <- get_Sigma_hat0(dat_boot)
    d_hat_star <- diag(Sigma_hat_star)
    M_star <- sqrt(n) * (m_bar_star - m_bar) / sqrt(diag(Sigma_hat_star))
    Omega_hat_star <- cov2cor(Sigma_hat_star)
    M_star_star <- M_star[keep_me]
    Omega_hat_star_star <- Omega_hat_star[keep_me, keep_me]
    Tn_boot[r] <- get_S(M_star_star, Omega_hat_star_star, p_GMS)
  }
  return(mean(Tn_boot > Tn))
}
