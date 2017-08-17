get_bonf_cover_raw <- function(results){
  g <- function(i) {
    params <- results$params[i,]
    b_true <- params$b
    CIs <- do.call(rbind, lapply(results$CIs_bonf[[i]], function(x) x$b))
    100 * round(get_coverage(b_true, CIs), 2)
  }
  as.data.frame(cbind(results$params,
                      stat = sapply(1:nrow(results$params), g)))
}
get_bonf_width_raw <- function(results){
  h <- function(i) {
    params <- results$params[i,]
    b_true <- params$b
    CIs <- do.call(rbind, lapply(results$CIs_bonf[[i]], function(x) x$b))
    round(get_median_width(CIs), 2)
  }
  as.data.frame(cbind(results$params,
                      stat = sapply(1:nrow(results$params), h)))
}

get_GMM_cover_raw <- function(results){
  g <- function(i) {
    params <- results$params[i,]
    b_true <- params$b
    CIs <- do.call(rbind, lapply(results$CIs_GMM[[i]], function(x) x$b))
    100 * round(get_coverage(b_true, CIs, NAempty = FALSE), 2)
  }
  as.data.frame(cbind(results$params,
                      stat = sapply(1:nrow(results$params), g)))
}
get_GMM_width_raw <- function(results){
  h <- function(i) {
    params <- results$params[i,]
    b_true <- params$b
    CIs <- do.call(rbind, lapply(results$CIs_GMM[[i]], function(x) x$b))
    round(get_median_width(CIs, NAempty = FALSE), 2)
  }
  as.data.frame(cbind(results$params,
                      stat = sapply(1:nrow(results$params), h)))
}

build_table <- function(results, TeX = TRUE) {
  b_vals <- unique(results$b)
  tab <- subset(results, b == b_vals[1])[,c('a0', 'a1', 'stat')]
  names(tab)[names(tab) == 'stat'] <- paste0('b=', b_vals[1])
  for(i in 2:length(b_vals)){
    temp <- subset(results, b == b_vals[i])[,c('a0', 'a1', 'stat')]
    names(temp)[names(temp) == 'stat'] <- paste0('b=', b_vals[i])
    tab <- merge(tab, temp)
  }
  if(TeX){
    names(tab) <- c('\\alpha_0', '\\alpha_1', paste0('\\beta=', b_vals))
  }
  return(tab)
}

TeXtable <- function(tab) {
  nRow <- nrow(tab)
  nCol <- ncol(tab)
  body <- rbind(colnames(tab), format(tab))
  body <- apply(body, c(1,2), function(x) paste0('$', x, '$'))
  rownames(body) <- NULL
  colnames(body) <- NULL
  body <- apply(body, 1, function(x) paste(x, collapse = ' & '))
  header <- body[1]
  body <- body[-1]
  body <- paste(body, collapse = "\\\\ \n")

  header <- paste("\\hline\\hline\n", header, "\\\\ \n \\hline\n")
  header <- paste0("\\begin{tabular}{rr|",
                   paste(rep("r", nCol - 2), collapse = ""), "}\n", header)
  footer <- "\\\\ \n \\hline \n \\end{tabular}"
  return(paste0(header, body, footer))
}
