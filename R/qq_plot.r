# This function samples k indices from 1:n starting with very dense sampling
# (every value) and then getting more and more sparse

#' Sampling for qq plots
#' 
#' Returns k indices between 1:n such that sampling is very dense at the start and much less dense at the end.
#' 
#' @export
thin <- function(n, k=2000) {
  if ( n < k ) {
    return (1:n)
  }
  
  # Generate k samples x such that exp(n-x) is uniformly distributed
  rx <- sort(log(runif(k)), decreasing=TRUE)
  x <- round(n-(rx-min(rx))/(max(rx)-min(rx))*(n-1))
  
  # Ensure each x is unique (Equivalently, diff(x) != 0)
  i <- 1
  j <- 1
  dx <- diff(x)
  
  while ( any(dx[i:(k-1)] == 0) ) {
    i <- i - 1 + which.max(dx[i:(k-1)] == 0)
    if ( dx[j] <= 1 )
      j <- j - 1 + which.max(dx[j:(k-1)] > 1)
    
    dx[i] <- dx[i] + 1
    dx[j] <- dx[j] - 1
  }
  
  cumsum(c(1, dx))
}

#' Fast QQ plots
#' 
#' QQ plots with lots of points can be slow. Since most of the points overlap
#' we don't actually have to plot all of them (and in many cases the output
#' will be identical.)
#' 
#' This function will either take a set of indices or a value k (default 2000)
#' and sample that many pvalues, then plot the log observed vs expected pvals
#' in a nice way.
#' 
#' @export
thin.qqplot <- function(pvals, thin.idx=NULL, k=2000, ...) {
  
  n.genotypes <- length(pvals)
  if ( missing(thin.idx) ) {
    if ( is.null(k) ) {
      thin.idx <- thin(n.genotypes)
    } else {
      thin.idx <- thin(n.genotypes, k)
    }
  }
  
  thin.logp.exp <- -log10(thin.idx/n.genotypes)
  
  thin.cint.95 <- sapply(thin.idx, function (x) { qbeta(0.95, x, n.genotypes - x + 1) })
  thin.cint.05 <- sapply(thin.idx, function (x) { qbeta(0.05, x, n.genotypes - x + 1) })
  
  
  thin.logp.obs <- -log10(pvals[order(pvals)[thin.idx]])
  
  plot(thin.logp.exp, thin.logp.obs, xlab=expression(-log[10](p[expected])), ylab=expression(-log[10](p[observed])), ...)
  abline(0, 1, col='gray', lty=2)
  lines(thin.logp.exp, -log10(thin.cint.95), lty=2, col='red')
  lines(thin.logp.exp, -log10(thin.cint.05), lty=2, col='red')
}