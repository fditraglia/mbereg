# This version assumes a0 = 0, uses only first moment inequality, and uses
# both first and second moment equalities.
BCS_test <- function(beta_null, dat, normal_draws){
  n <- nrow(dat)
  kappa_n <- sqrt(log(n))
  z <- with(dat, z)
  q <- mean(z)
  w1 <- with(dat, Tobs)
  w2 <- with(dat, y)
  w3 <- with(dat, 2 * y * Tobs)
  w4 <- with(dat, y^2)
  W <- cbind(w1, w2, w3, w4)
  W_tilde <- z * scale(W, center = TRUE, scale = FALSE)
  W_tilde_bar <- colMeans(W_tilde)
  w1_z0 <- w1 * (1 - z) / (1 - q)
  w1_z1 <- w1 * z / q
  w_z_bar <- c(mean(w1_z0), mean(w1_z1))
  Sigma_II <- var(cbind(w1_z0, w1_z1))
  s_II <- sqrt(diag(Sigma_II))
  M_EE <- var(W_tilde - q * W)
  M_IE <- rbind(cov(w1_z0, q * W) - cov(w1_z0, W_tilde),
                cov(w1_z1, q * W) - cov(w1_z1, W_tilde))
  
  get_Nu <- function(a1){
    theta1 <- beta_null / (1 - a1)
    theta2 <- beta_null * theta1
    nu1 <- c(-theta1, 1, 0, 0)
    nu2 <- c(theta2, 0, -theta1, 1)
    return(cbind(nu1, nu2))
  }
  
  get_Sigma_EE <- function(a1){
    Nu <- get_Nu(a1)
    crossprod(Nu, M_EE) %*% Nu
  }
  
  get_Sigma <- function(a1){
    Nu <- get_Nu(a1)
    Sigma_EE <- crossprod(Nu, M_EE) %*% Nu
    Sigma_IE <- M_IE %*% Nu
    rbind(cbind(Sigma_II, Sigma_IE), cbind(t(Sigma_IE), Sigma_EE))
  }
    
  Qn <- function(a1){
    m_bar_I <- (1 - a1) - w_z_bar
    Tn_I <- sum((m_bar_I < 0) * ((sqrt(n) * m_bar_I / s_II)^2))
    Nu <- get_Nu(a1)
    s_EE <- sqrt(diag(crossprod(Nu, M_EE) %*% Nu))
    m_bar_E <- drop(crossprod(Nu, W_tilde_bar))
    Tn_E <- sum((sqrt(n) * m_bar_E / s_EE)^2)
    return(Tn_I + Tn_E)
  }
  
  # We never optimize over this function, so I calculate the values for all
  # bootstrap draws simultaneously
  Qn_DR <- function(a1){
    # Moment selection
    m_bar_I <- (1 - a1) - w_z_bar
    phi <- ifelse(sqrt(n) * m_bar_I / s_II > kappa_n, Inf, 0)
    # Transform normal_draws to have appropriate covariance matrix
    Sigma <- get_Sigma(a1)
    Omega_sqrt <- sqrtm(cov2cor(Sigma))
    boot_draws <- tcrossprod(Omega_sqrt, normal_draws)
    # Row sums for Tn_DR_I since apply gives a transposed result
    Tn_DR_I <- rowSums(apply(boot_draws[1:2,] + phi, 1, function(x) pmin(x, 0)^2))
    Tn_DR_E <- colSums(boot_draws[3:4,]^2)
    return(Tn_DR_I + Tn_DR_E)
  }
  
  # Since we will optimize this for each bootstrap draw, I have written this 
  # differently from Qn_DR so that we evaluate for a fixed 
  Qn_PR <- function(a1, normal_draw){
    Sigma <- get_Sigma(a1)
    Omega_sqrt <- sqrtm(cov2cor(Sigma))
    boot_draw <- drop(Omega_sqrt %*% normal_draw)
    m_bar_I <- (1 - a1) - w_z_bar
    ell_I <- sqrt(n) * m_bar_I / (s_II * kappa_n)
    Tn_PR_I <- sum(pmin(boot_draw[1:2] + ell_I, 0)^2)
    Nu <- get_Nu(a1)
    m_bar_E <- drop(crossprod(Nu, W_tilde_bar))
    s_EE <- sqrt(c(Sigma[3,3], Sigma[4,4]))
    ell_E <- sqrt(n) * m_bar_E / (s_EE * kappa_n)
    Tn_PR_E <- sum((boot_draw[3:4] + ell_E)^2)
    return(Tn_PR_I + Tn_PR_E)
  }
  
  # If beta_null = 0, inference for beta does not depend on alpha: there is nothing
  # to profile, so we revert to a bootstrap version of the "standard" test.
  if(isTRUE(all.equal(0, beta_null))){
    Tn <- Qn(0)
    Sigma_EE <- get_Sigma_EE(0)
    Omega_EE_sqrt <- sqrtm(cov2cor(Sigma_EE))
    boot_draws_EE <- tcrossprod(Omega_EE_sqrt, normal_draws[,3:4])
    Tn_DR <- Tn_PR <- colSums(boot_draws_EE^2)
    
  } else {
    Qn_optimization <- optimize(Qn, lower = 0, upper = 0.99)
    Tn <- Qn_optimization$objective
    Tn_DR <- Qn_DR(Qn_optimization$minimum)
   
     
    Tn_PR <- apply(normal_draws, 1, function(x) optimize(function(a1) Qn_PR(a1, x),
                                                         lower = 0, upper = 0.99)$objective)
  }
  Tn_MR <- pmin(Tn_DR, Tn_PR)
  return(mean(Tn_MR >= Tn))
}


# # True parameter values
# a1_true <- 0.1
# beta_true <- 0.25
# d <- 0.15
# p0_star <- d
# p1_star <- 1 - d
# RF <- beta_true * (p1_star - p0_star)
# Wald <- beta_true / (1 - a1_true)
# b_lower <- RF * p1_star / (p1_star - p0_star)
# b_upper <- Wald
# 
# set.seed(72349)
# B <- 2001
# mynormals <- matrix(rnorm(4 * B), B, 4)
# 
# sim_test_null <- function(){
#   mydat <- dgp(a0 = 0, a1 = a1_true, b = beta_true)
#   p_lower <- BCS_test(b_lower, mydat, mynormals)
#   p_true <- BCS_test(beta_true, mydat, mynormals)
#   p_upper <- BCS_test(b_upper, mydat, mynormals)
#   out <- c(lower = p_lower, true = p_true, upper = p_upper)
#   return(out)
# }
# 
# # sim_test_alternative <- function(){
# #   mydat <- dgp(a0 = 0, a1 = a1_true, b = beta_true, n = n)
# #   p_0 <- BCS_test(0, mydat, myzeta)
# #   p_hundredth <- BCS_test(0.01, mydat, myzeta)
# #   p_tenth <- BCS_test(0.1, mydat, myzeta)
# #   out <- c(zero = p_0, hundredth = p_hundredth, tenth = p_tenth)
# #   return(out)
# # }
# 
# library(parallel)
# RNGkind("L'Ecuyer-CMRG")
# set.seed(12871)
# 
# n_reps <- 160
# #results <- replicate(n_reps, sim_test_alternative())
# #results <- replicate(n_reps, sim_test_null())
# system.time(results <- parallel::mclapply(1:n_reps, function(i) sim_test_null(), 
#                                           mc.cores = 8))
# results <- as.data.frame(do.call(rbind, results))
# # #results <- as.data.frame(t(results))
# # 
# apply(results, 2, function(pvalue) mean(pvalue <= 0.05))
