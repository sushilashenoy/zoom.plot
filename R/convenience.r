#' More useful defaults for plotting matrices
#' @export

image0 <- function(x, xlab, ylab, yaxt, ...) {

  if (missing(xlab)) 
    xlab <- ""
  if (missing(ylab)) 
    ylab <- ""
  
  image(1:ncol(x), 1:nrow(x), t(x[nrow(x):1, ]), xlab=xlab, ylab=ylab, yaxt='n', ...)
  if (missing(yaxt) || yaxt != 'n')
    axis(2, nrow(x)-axTicks(2)+1, axTicks(2))
}