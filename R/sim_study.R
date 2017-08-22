sim_alphas_size <- function(true_params, normal_sims, ncores, nreps = 1000,
                              test_alphas = GMS_test_alphas){
  sim_rep <- function(){
    sim_dat <- dgp(a0 = true_params$a0,
                   a1 = true_params$a1,
                   b = true_params$b,
                   n = true_params$n,
                   d = true_params$d,
                   rho = true_params$rho,
                   cc = true_params$cc)
    pvalue <- test_alphas(a0 = true_params$a0,
                          a1 = true_params$a1,
                          dat = sim_dat,
                          normal_sims = normal_sims)
    return(pvalue)
  }
  return(parallel::mclapply(1:nreps, function(i) sim_rep(), mc.cores = ncores))
  #return(replicate(nreps, sim_rep()))
}

sim_bonf_CI <- function(true_params, normal_sims, delta1, delta2, ncores,
                        nreps = 1000, test_alphas = GMS_test_alphas){
  sim_rep <- function(){
    sim_dat <- dgp(a0 = true_params$a0,
                   a1 = true_params$a1,
                   b = true_params$b,
                   n = true_params$n,
                   d = true_params$d,
                   rho = true_params$rho,
                   cc = true_params$cc)
    CIs <- bonf_CI(dat = sim_dat,
                   normal_sims = normal_sims,
                   test_alphas = test_alphas,
                   delta1 = delta1,
                   delta2 = delta2)
    return(CIs)
  }
  return(parallel::mclapply(1:nreps, function(i) sim_rep(), mc.cores = ncores))
}

sim_GMM_endog <- function(true_params, delta, ncores, nreps = 1000) {
  sim_rep <- function(){
    sim_dat <- dgp(a0 = true_params$a0,
                   a1 = true_params$a1,
                   b = true_params$b,
                   n = true_params$n,
                   d = true_params$d,
                   rho = true_params$rho,
                   cc = true_params$cc)
    gmm <- GMM_endog(sim_dat)
    if((!is.na(gmm$est$b)) & (!is.na(gmm$SE$b))){
      b <- gmm$est$b
      SE <- gmm$SE$b
      ME <- SE * qnorm(1 - delta / 2)
      LCL <- b - ME
      UCL <- b + ME
      return(list(b = c(lower = LCL, upper = UCL)))
    } else {
      return(list(b = c(lower = NA, upper = NA)))
    }
  }
  return(parallel::mclapply(1:nreps, function(i) sim_rep(), mc.cores = ncores))
}


sim_GMM_exog <- function(true_params, delta, ncores, nreps = 1000) {
  sim_rep <- function(){
    sim_dat <- dgp(a0 = true_params$a0,
                   a1 = true_params$a1,
                   b = true_params$b,
                   n = true_params$n,
                   d = true_params$d,
                   rho = true_params$rho,
                   cc = true_params$cc)
    gmm <- GMM_exog(sim_dat)
    if((!is.na(gmm$est$b)) & (!is.na(gmm$SE$b))){
      b <- gmm$est$b
      SE <- gmm$SE$b
      ME <- SE * qnorm(1 - delta / 2)
      LCL <- b - ME
      UCL <- b + ME
      return(list(b = c(lower = LCL, upper = UCL)))
    } else {
      return(list(b = c(lower = NA, upper = NA)))
    }
  }
  return(parallel::mclapply(1:nreps, function(i) sim_rep(), mc.cores = ncores))
}


sim_bonf_vs_gmm <- function(true_params, normal_sims, delta1, delta2, ncores,
                            nreps = 1000, test_alphas = GMS_test_alphas,
                            get_gmm = GMM_endog){
  sim_rep <- function(){
    sim_dat <- dgp(a0 = true_params$a0,
                   a1 = true_params$a1,
                   b = true_params$b,
                   n = true_params$n,
                   d = true_params$d,
                   rho = true_params$rho,
                   cc = true_params$cc)
    bonf <- bonf_CI(dat = sim_dat,
                    normal_sims = normal_sims,
                    test_alphas = test_alphas,
                    delta1 = delta1,
                    delta2 = delta2)$b
    gmm <- get_gmm(sim_dat)
    if((!is.na(gmm$est$b)) & (!is.na(gmm$SE$b))){
      b <- gmm$est$b
      SE <- gmm$SE$b
      delta <- delta1 + delta2
      ME <- SE * qnorm(1 - delta / 2)
      LCL <- b - ME
      UCL <- b + ME
      gmm_CI <- c(LCL, UCL)
    } else {
      gmm_CI <- c(NA, NA)
    }
    CIs <- list(bonf = bonf, gmm = gmm_CI)
    return(CIs)
  }
  return(parallel::mclapply(1:nreps, function(i) sim_rep(), mc.cores = ncores))
}
