summarize_CI_results <- function(results) {
  h <- function(i) {
    params <- results$params[i,]
    b_true <- params$b

    CIs_GMM <- do.call(rbind, lapply(results$CIs[[i]], function(x) x$gmm))
    cover_GMM <- get_coverage(b_true, CIs_GMM, NAempty = FALSE)
    width_GMM <- get_median_width(CIs_GMM, NAempty = FALSE)

    CIs_bonf <- do.call(rbind, lapply(results$CIs[[i]], function(x) x$bonf))
    cover_bonf <- get_coverage(b_true, CIs_bonf)
    width_bonf <- get_median_width(CIs_bonf)

    CIs_twostep <- get_twostep_CI(CIs_GMM, CIs_bonf)
    cover_twostep <- get_coverage(b_true, CIs_twostep)
    width_twostep <- get_median_width(CIs_twostep)

    GMM_na <- get_prop_CIs_na(CIs_GMM)

    out <- data.frame(cover_GMM = cover_GMM,
                      cover_bonf = cover_bonf,
                      cover_twostep = cover_twostep,
                      width_GMM = width_GMM,
                      width_bonf = width_bonf,
                      width_twostep = width_twostep,
                      GMM_na = GMM_na)
    return(out)
  }
  cbind(results$params, t(sapply(1:nrow(results$params), h)))
}



build_table <- function(results, TeX = TRUE, out_stat) {
  b_vals <- unique(results$b)
  tab <- subset(results, b == b_vals[1])[,c('a0', 'a1', out_stat)]
  names(tab)[names(tab) == out_stat] <- paste0('b=', b_vals[1])
  for(i in 2:length(b_vals)){
    temp <- subset(results, b == b_vals[i])[,c('a0', 'a1', out_stat)]
    names(temp)[names(temp) == out_stat] <- paste0('b=', b_vals[i])
    tab <- merge(tab, temp)
  }
  if(TeX){
    names(tab) <- c('\\alpha_0', '\\alpha_1', b_vals)
  }
  return(tab)
}

TeXtable <- function(tab) {
  n_a1 <- length(unique(tab[,2]))
  n_a0 <- length(unique(tab[,1]))
  keep_a0 <- n_a0 * 0:(n_a1 - 1) + 2
  nRow <- nrow(tab)
  nCol <- ncol(tab)
  body <- rbind(colnames(tab), format(tab))
  rownames(body) <- NULL
  colnames(body) <- NULL
  body <- apply(body, c(1,2), function(x) paste0('$', x, '$'))
  body[-c(1, keep_a0),1] <- ""
  body <- apply(body, 1, function(x) paste(x, collapse = ' & '))
  header <- body[1]
  body <- body[-1]
  body[keep_a0[-1] -1] <- paste('\\hline \n', body[keep_a0[-1] - 1])

  body <- paste(body, collapse = "\\\\ \n")

  header <- paste("\\hline\\hline\n", paste0("&& \\multicolumn{", nCol - 2,
                                             "}{c}{$\\beta$}\\\\\n"), header,
                  "\\\\ \n \\hline\n")
  header <- paste0("\\begin{tabular}{rr|",
                   paste(rep("r", nCol - 2), collapse = ""), "}\n",header)
  footer <- "\\\\ \n \\hline \n \\end{tabular}"
  return(paste0(header, body, footer))
}
