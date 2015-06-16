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
  
  # Figure out how many groups we need to divide n into
  # so that we end up with k values
  # We want each group to be twice as big as the previous group
  # We want to sample every index in the first group
  x <- 2
  while ( (2*x-2)/(2^x-2) > k/n ) { x <- x + 1 }
  
  # Each group is twice as big as the previous group
  end.idx <- c(round(2^(1:(x-1))*n/(2^x-2)), n)
  start.idx <- c(1, end.idx[-x]+1)
  
  # Sample the same number of values from each group
  sample.size <- ceiling(k/x)
  thin.idx <- NULL
  
  for ( i in 1:x ) {
    if ( end.idx[i] - start.idx[i] < sample.size) {
      # If there aren't enough values use all indices in group i
      thin.idx <- c(thin.idx, start.idx[i]:end.idx[i])
      # Adjust sample size accordingly
      sample.size <- ceiling((k-length(thin.idx))/(x-i))
    } else {
      # Sample from group i
      thin.idx <- c(thin.idx, sort(sample(start.idx[i]:end.idx[i], sample.size)))
    }
  }
  
  # Return the list of indices
  return ( thin.idx )
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