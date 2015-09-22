# This function samples k indices from 1:n starting with very dense sampling
# (every value) and then getting more and more sparse

#' Sampling for qq plots
#' 
#' Returns k indices between 1:max.idx (x) such that -log(x/n) is uniformly
#' distributed. This is useful to reduce overplotting when plotting
#' log-transformed uniform data (such as p-values in a qq plot).
#' 
#' @examples
#' n <- 1e5
#' x <- sort(runif(n))
#' thin.idx <- thin(n, 500)
#' par(mfcol=c(1, 2))
#' plot(-log10(1:n/n), -log10(x), ann=FALSE)
#' plot(-log10(thin.idx/n), -log10(x[thin.idx]), ann=FALSE)
#' 
#' @export
thin <- function(n, k=2000, max.idx=n) {
  
  if ( max.idx > n ) max.idx <- n
  if ( max.idx <= k ) return (1:max.idx)
  
  # Generate k samples between 1 and max.idx
  # such that -log10(x/n) is uniformly distributed
  x <- round(max.idx*exp(-(k:1)/k*log(n)))
  
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
#' we don't actually have to plot all of them (and in most cases the output
#' will be identical.)
#' 
#' This function will take a value k (default 2000) and sample that many
#' pvalues in an intelligent way such that we don't lose information, then
#' plot the log observed vs expected pvals in a nice way.
#' 
#' The optimal value of k depends on the resolution/size of the plot. Very high
#' resolution or very large plots may need larger k to look correct.
#' 
#' Additional arguments in ... are passed to plot().
#' 
#' @export
fastqq <- function(pvals, k=2000, ...) {
  
  np <- length(pvals)
  thin.idx <- thin(np, k)
  
  thin.cint.95 <- qbeta(0.95, thin.idx, np - thin.idx + 1)
  thin.cint.05 <- qbeta(0.05, thin.idx, np - thin.idx + 1)
  
  thin.logp.exp <- -log10(thin.idx/np)
  thin.logp.obs <- -log10(pvals[order(pvals)[thin.idx]])
  
  plot(thin.logp.exp, thin.logp.obs, xlab=expression(-log[10](p[expected])), ylab=expression(-log[10](p[observed])), ...)
  abline(0, 1, col='gray', lty=2)
  lines(thin.logp.exp, -log10(thin.cint.95), lty=2, col='red')
  lines(thin.logp.exp, -log10(thin.cint.05), lty=2, col='red')
}

#' @export
thin.qqplot <- function(pvals, thin.idx=NULL, k=2000, ...) {
  warning('This function is deprecated, use fastqq')
  fastqq(pvals, k, ...)
}