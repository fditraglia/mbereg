# Unconstrained GMM estimator of model with endogeneity
GMM_endog <- function(dat){
  y <- dat$y
  Tobs <- dat$Tobs
  z <- dat$z
  n <- nrow(dat)

  PI <- cov(Tobs, z)
  eta1 <- cov(y, z)
  eta2 <- cov(y^2, z)
  eta3 <- cov(y^3, z)
  tau1 <- cov(Tobs * y, z)
  tau2 <- cov(Tobs * y^2, z)
  theta1 <- eta1 / PI
  theta2 <- 2 * tau1 * eta1 / PI^2 - eta2 / PI
  theta3 <- eta3 / PI - 3 * (tau2 * eta1 + tau1 * eta2) / PI^2 +
    6 * tau1^2 * eta1 / PI^3

  b_squared <- 3 * (theta2 / theta1)^2 - 2 * (theta3 / theta1)

  est <- list(a0 = NA, a1 = NA, b = NA)
  SE <- list(a0 = NA, a1 = NA, b = NA)

  if(b_squared > 0){
    # Calculate GMM point estimates for beta, a0, a1
    b <- sign(theta1) * sqrt(b_squared)
    A <- theta2 / theta1^2
    B <- theta3 / theta1^3
    r <- Re(polyroot(c(A^2 - B, 2 * A, -2)))
    a0 <- min(r)
    a1 <- 1 - max(r)
    est$a0 <- a0
    est$a1 <- a1
    est$b <- b

    # Calculate GMM asymptotic standard errors
    w <- cbind(Tobs, y, y * Tobs, y^2, y^2 * Tobs, y^3)
    w_bar <- colMeans(w)
    w_centered <- scale(w, scale = FALSE, center = TRUE)
    w_centered_z <- (w_centered * z)
    wz_bar <- colMeans(w * z)

    psi1 <- c(-theta1, 1, 0, 0, 0, 0)
    psi2 <- c(theta2, 0, -2 * theta1, 1, 0, 0)
    psi3 <- c(-theta3, 0, 3 * theta2, 0, -3 * theta1, 1)
    Psi <- cbind(psi1, psi2, psi3)

    g <- w_centered_z %*% Psi
    h <- w_centered %*% Psi
    f <- cbind(g, h)

    D_theta1 <- c(theta1, theta1, 1) / (1 - a0 - a1)
    D_theta2 <- 2 * theta1 * (1 + a0 - a1) * D_theta1 + theta1^2 * c(1 , -1, 0)
    D_theta3 <- 3 * theta1^2 * ((1 - a0 - a1)^2 + 6 * a0 * (1 - a1)) * D_theta1 +
      theta1^3 * c(-2 * (1 - a0 - a1) + 6 * (1 - a1),
                   -2 * (1 - a0 - a1) - 6 * a0, 0)

    D_psi1 <- rbind(-D_theta1,
                    matrix(0, 5, 3))
    D_psi2 <- rbind(D_theta2,
                    matrix(0, 1, 3),
                    -2 * D_theta1,
                    matrix(0, 3, 3))
    D_psi3 <- rbind(-D_theta3,
                    matrix(0, 1, 3),
                    3 * D_theta2,
                    matrix(0, 1, 3),
                    -3 * D_theta1,
                    matrix(0, 1, 3))

    G <- rbind(t(wz_bar) %*% D_psi1,
               t(wz_bar) %*% D_psi2,
               t(wz_bar) %*% D_psi3)
    H <- rbind(t(w_bar) %*% D_psi1,
               t(w_bar) %*% D_psi2,
               t(w_bar) %*% D_psi3)
    q <- mean(z)

    Eff <- rbind(cbind(G, -q * diag(nrow(G))),
                 cbind(H, -diag(nrow(H))))
    if(kappa(Eff) * .Machine$double.eps < 1) {
      F_inv <- solve(Eff, tol = 1e-17)
      Omega <- var(f)
      Avar <- F_inv %*% Omega %*% t(F_inv)
      SE_full <- sqrt(diag(Avar) / n)
      SE$a0 <- SE_full[1]
      SE$a1 <- SE_full[2]
      SE$b <- SE_full[3]
    }
  }
  return(list(est = est, SE = SE))
}

# Unconstrained GMM estimator of model under the joint exogeneity assumption of
# Frazis & Loewenstein (2003), Mahajan (2006), etc.
GMM_exog <- function(dat){
  y <- dat$y
  Tobs <- dat$Tobs
  z <- dat$z
  n <- nrow(dat)

  PI <- cov(z, Tobs)
  theta1 <- cov(y, z) / PI
  kappa1 <- mean(y) - theta1 * mean(Tobs)
  eta <- cov((y - kappa1) * Tobs, z) / PI
  rho <- mean(y * Tobs) - (kappa1 + eta) * mean(Tobs)

  b_squared <- eta^2 + 4 * theta1 * rho

  est <- list(a0 = NA, a1 = NA, b = NA)
  SE <- list(a0 = NA, a1 = NA, b = NA)

  if(b_squared > 0){
    # Calculate GMM point estimates for beta, a0, a1
    b <- sign(theta1) * sqrt(b_squared)
    A <- eta / theta1
    B <- -rho / theta1
    r <- Re(polyroot(c(B, -A, 1)))
    a0 <- min(r)
    a1 <- 1 - max(r)
    est$a0 <- a0
    est$a1 <- a1
    est$b <- b

    # Calculate GMM asymptotic standard errors
    g <- cbind(y - kappa1 - theta1 * Tobs, (y - kappa1 - theta1 * Tobs) * z)
    h <- cbind((y - kappa1) * Tobs - rho -  eta * Tobs,
               ((y - kappa1) * Tobs - rho -  eta * Tobs) * z)
    f <- cbind(g, h)

    D_theta1 <- c(theta1, theta1, 1) / (1 - a0 - a1)
    D_rho <- -a0 * (1 - a1) * D_theta1 + c(-theta1 * (1 - a1), theta1 * a0, 0)
    D_eta <- (1 + a0 - a1) * D_theta1 + c(theta1, -theta1, 0)


    D_theta1 <- c(theta1, theta1, 1) / (1 - a0 - a1)
    D_theta2 <- 2 * theta1 * (1 + a0 - a1) * D_theta1 + theta1^2 * c(1 , -1, 0)
    D_theta3 <- 3 * theta1^2 * ((1 - a0 - a1)^2 + 6 * a0 * (1 - a1)) * D_theta1 +
      theta1^3 * c(-2 * (1 - a0 - a1) + 6 * (1 - a1),
                   -2 * (1 - a0 - a1) - 6 * a0, 0)

    q <- mean(z)
    p <- mean(Tobs)
    p1 <- mean(Tobs[z == 1])
    G <- rbind(cbind(-t(D_theta1) * p, -1),
               cbind(-t(D_theta1) * p1 * q, -q))
    H <- rbind(cbind(-t(D_rho + p * D_eta), -p),
               cbind(-t(q * D_rho + p1 * q * D_eta), -p1 * q))
    Eff <- rbind(G, H)


    if(kappa(Eff) * .Machine$double.eps < 1) {
      F_inv <- solve(Eff, tol = 1e-17)
      Omega <- var(f)
      Avar <- F_inv %*% Omega %*% t(F_inv)
      SE_full <- sqrt(diag(Avar) / n)
      SE$a0 <- SE_full[1]
      SE$a1 <- SE_full[2]
      SE$b <- SE_full[3]
    }
  }
  return(list(est = est, SE = SE))
}
