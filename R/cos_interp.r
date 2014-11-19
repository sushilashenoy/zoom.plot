#' Piecewise cosine interpolation
#' 
#' Returns a list of points which interpolate given data points. Similar to
#' \link[stats]{approx}.
#' 
#' @param x,y numeric vectors giving the coordinates of the points to be 
#'   interpolated.
#' @param xout an optional set of numeric values specifying where interpolation 
#'   is to take place.
#' @param n If xout is not specified, interpolation takes place at n equally 
#'   spaced points spanning the interval [\code{min(x)}, \code{max(x)}].
#' @param yleft the value to be returned when input \code{x} values are less 
#'   than \code{min(x)}. The default is defined by the value of rule given 
#'   below.
#' @param yright the value to be returned when input \code{x} values are greater
#'   than \code{max(x)}. The default is defined by the value of rule given 
#'   below.
#' @param rule an integer (of length 1 or 2) describing how interpolation is to 
#'   take place outside the interval [\code{min(x)}, \code{max(x)}]. If rule is 
#'   \code{1} then \code{NA}s are returned for such points and if it is
#'   \code{2}, the value at the closest data extreme is used. Use, e.g.,
#'   \code{rule = 2:1}, if the left and right side extrapolation should differ.
#' @param ties Handling of tied \code{x} values. Either a function with a single
#'   vector argument returning a single number result or the string 
#'   \code{"ordered"}.
#'   
#' @details This function is based on and functions similarly to
#' \link[stats]{approx} and \link[stats]{spline} .
#' 
#' @examples
#' y <- rnorm(10)
#' x <- 1:10
#' plot(x, y)
#' lines(approx.cos(x, y, n=200))
#' @export
approx.cos <- function(x, y=NULL, xout, n=50, yleft, yright, rule=1, ties=mean) {
  require('stats')
  
  stopifnot(is.numeric(rule), (lenR <- length(rule)) >= 1L, 
            lenR <= 2L)
  if (lenR == 1) 
    rule <- rule[c(1, 1)]
  
  x <- stats:::regularize.values(x, y, ties)
  y <- x$y
  x <- x$x
  
  nx <- as.integer(length(x))
  if (is.na(nx)) 
    stop("invalid length(x)")
  if (nx <= 1) {
    stop("need at least two non-NA values to interpolate")
  }
  
  if (missing(yleft)) 
    yleft <- if (rule[1L] == 1) 
      NA
  else y[1L]
  if (missing(yright)) 
    yright <- if (rule[2L] == 1) 
      NA
  else y[length(y)]
  stopifnot(length(yleft) == 1L, length(yright) == 1L)
  if (missing(xout)) {
    if (n <= 0) 
      stop("'approx' requires n >= 1")
    xout <- seq.int(x[1L], x[nx], length.out = n)
  }
  x <- as.double(x)
  y <- as.double(y)
  
  # Select values of x we can actually interpolate
  xnew <- xout[xout >= min(x) & xout <= max(x)]
  
  # For each new x find which original x it is next to
  xi <- findInterval(xnew, x, TRUE)
  
  # calculate relative distance to original x
  dx <- (x[xi+1]-xnew)/(x[xi+1]-x[xi])
  
  # Interpolate with piecewise cosine function
  ynew <- (0.5 + cos(pi*dx)/2) * (y[xi+1] - y[xi]) + y[xi]
  
  # Add y values for any x values outside interval
  yout <- c(rep(yleft, sum(xout<min(x))), ynew, rep(yright, sum(xout>max(x))))
  
  list(x = xout, y = yout)
}