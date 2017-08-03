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
