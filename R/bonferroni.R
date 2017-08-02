# Bonferroni Confidence Interval for beta constructed from the ordinary Wald
# interval for the reduced form parameter theta1 and a projection confidence
# set for s = (1 - a0 - a1) formed from a joint GMS confidence set for (a0, a1)
bonf_CI <- function(dat, normal_sims, test_alphas = GMS_test_alphas,
                    delta1, delta2, inc_coarse = 0.05,
                    inc_fine = 0.005){

  # Confidence interval for Wald estimator
  iv_slope <- with(dat, cov(y, z) / cov(Tobs, z))
  iv_intercept <- with(dat, mean(y - iv_slope * Tobs))
  epsilon_hat <- with(dat, y - iv_intercept - iv_slope * Tobs)
  s_epsilon <- sqrt(sum(epsilon_hat^2) / (length(epsilon_hat) - 2))
  se <- with(dat, (s_epsilon / sqrt(length(epsilon_hat)) * sd(z) / cov(Tobs, z)))
  UCL <- iv_slope + qnorm(1 - delta2 / 2) * se
  LCL <- iv_slope - qnorm(1 - delta2 / 2) * se
  theta1_CI <- c(LCL, UCL)

  coarse <- seq(0, 1, inc_coarse)
  coarse <- expand.grid(a0 = coarse, a1 = coarse)
  coarse <- subset(coarse, a0 + a1 < 1)
  p_coarse <- sapply(1:nrow(coarse),
                    function(i) test_alphas(coarse$a0[i], coarse$a1[i],
                                                dat, normal_sims))


  if(any(p_coarse > delta1)) {

    CI_coarse <- subset(coarse, p_coarse > delta1)
    a0_fine <- seq(max(min(CI_coarse$a0) - inc_coarse, 0),
                   max(CI_coarse$a0) + inc_coarse, inc_fine)
    a1_fine <- seq(max(min(CI_coarse$a1) - inc_coarse, 0),
                   max(CI_coarse$a1) + inc_coarse, inc_fine)
    fine <- expand.grid(a0 = a0_fine, a1 = a1_fine)
    fine <- subset(fine, a0 + a1 < 1)
    sum_coarse_max <- with(CI_coarse, max(a0 + a1))
    sum_coarse_min <- with(CI_coarse, min(a0 + a1))
    fine <- subset(fine, ((a0 + a1 ) < sum_coarse_min) |
                     ((a0 + a1) > sum_coarse_max))
    p_fine <- sapply(1:nrow(fine),
                      function(i) test_alphas(fine$a0[i], fine$a1[i],
                                                  dat, normal_sims))
    CI_fine <- subset(fine, p_fine >= delta1)
    CI_alphas <- rbind(CI_coarse, CI_fine)


    s_CI <- with(CI_alphas, range(1 - a0 - a1))
    bonf_LCL <- min(LCL * min(s_CI), LCL * max(s_CI))
    bonf_UCL <- max(UCL * min(s_CI), UCL * max(s_CI))
    bonf_CI <- c(bonf_LCL, bonf_UCL)

  } else {
    s_CI <- bonf_CI <- c(NA, NA)
  }
    out <- list(theta1 = theta1_CI,
                s = s_CI,
                b = bonf_CI)
    return(out)
}
