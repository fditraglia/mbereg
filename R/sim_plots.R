CI_plot <- function(results, i, param = 'b', mult = 1, nx = 500, nominal = 0.9){
  stopifnot(param %in% c('b', 'theta1', 's'))

  true_params <- results$params[i,]
  b_true <- true_params$b
  a0_true <- true_params$a0
  a1_true <- true_params$a1
  s_true <- 1 - a0_true - a1_true
  n_true <- true_params$n
  theta1_true <- b_true / s_true
  d_true <- true_params$d
  RF_true <- b_true * (1 - 2 * d_true)

  CIs <- results$CIs_bonf[[i]]
  CIs <- lapply(CIs, function(x) getElement(x, param))
  CIs <- do.call(rbind, CIs)

  if(param == 's'){
    x <- seq(0, 0.999, length.out = nx / 4)
  } else {
    xmin <- min(RF_true, theta1_true)
    xmax <- max(RF_true, theta1_true)
    xwidth <- max(xmax - xmin, 1)
    x <- seq(xmin - mult * xwidth, xmax + mult * xwidth, length.out = nx)
  }
  coverage <- get_coverage(x, CIs)

  if(param == 'b') {
    myxlab <- '$\\beta$'
    param_true <- b_true
  } else if(param == 'theta1') {
    myxlab <- '$\\theta_1$'
    param_true <- theta1_true
  } else {
    myxlab <- '$s$'
    param_true <- s_true
  }
  mymain <- paste0('$\\beta = ', b_true,
                   ', \\alpha_0 = ', a0_true,
                   ', \\alpha_1 = ', a1_true,
                   ', n = ', n_true, '$')
  plot(x, coverage, type = 'l', xlab = myxlab, ylab = 'Coverage',
       main = mymain)
  abline(v = param_true, lty = 2)
  abline(h = nominal, lty = 2)

  if(param == 'b' & (RF_true != theta1_true)){
    abline(v = RF_true, col = 'red')
    abline(v = theta1_true, col = 'blue')
  } else if(param == 'theta1') {
    abline(v = b_true, col = 'blue')
  } else {
    #abline(v = 1 - (2 * d_true), col = 'red')
    #abline(v = 1, col = 'blue')
  }
}

CI_compare_plot <- function(results_weak, results_nondiff, i, xlims = NULL,
                            nx = 500, nominal = 0.9, alphaonly = FALSE){

  stopifnot(identical(results_weak$params[i,], results_nondiff$params[i,]))

  true_params <- results_weak$params[i,]
  b_true <- true_params$b
  a0_true <- true_params$a0
  a1_true <- true_params$a1
  s_true <- 1 - a0_true - a1_true
  n_true <- true_params$n
  theta1_true <- b_true / s_true
  d_true <- true_params$d
  RF_true <- b_true * (1 - 2 * d_true)

  CIs_weak <- results_weak$CIs_bonf[[i]]
  CIs_weak <- lapply(CIs_weak, function(x) x$b)
  CIs_weak <- do.call(rbind, CIs_weak)

  CIs_nondiff <- results_nondiff$CIs_bonf[[i]]
  CIs_nondiff <- lapply(CIs_nondiff, function(x) x$b)
  CIs_nondiff <- do.call(rbind, CIs_nondiff)

  if(is.null(xlims)){
    xmin <- min(quantile(CIs_nondiff[,1], 0.01, na.rm = TRUE),
                quantile(CIs_weak[,1], 0.01, na.rm = TRUE))
    xmax <- max(quantile(CIs_nondiff[,2], 0.99, na.rm = TRUE),
                quantile(CIs_weak[,2], 0.99, na.rm = TRUE))
  } else {
    xmin <- min(xlims)
    xmax <- max(xlims)
  }

  x <- seq(xmin, xmax, length.out = nx)
  coverage_weak <- get_coverage(x, CIs_weak)
  coverage_nondiff <- get_coverage(x, CIs_nondiff)

  myxlab <- '$\\beta$'
  if(alphaonly) {
    mymain <- paste0('$\\alpha_0 = ', a0_true,
                     ', \\alpha_1 = ', a1_true, '$')
  } else {
    mymain <- paste0('$\\beta = ', b_true,
                     ', \\alpha_0 = ', a0_true,
                     ', \\alpha_1 = ', a1_true,
                     ', n = ', n_true, '$')
  }
  plot(x, coverage_weak, type = 'l', xlab = '', ylab = '',
       main = mymain, lwd = 2, col = 'red', lty = 2)
  points(x, coverage_nondiff, type = 'l', lwd = 2, col = 'blue', lty = 1)
  abline(h = nominal, lty = 3)
  abline(v = RF_true, lty = 3)
  abline(v = theta1_true, lty = 3)
}
