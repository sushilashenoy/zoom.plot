#' Make a stacked multi-figure plot
#' 
#' \code{fig.parts} takes a vector of percents and returns a convenient function
#' to help with plotting. 
#' 
#' @param percents vector of percents, totalling 1 (or less). Any NA values will
#' be substituted with equal values such that the total is 1
#' 
#' @examples
#' fig.part <- fig.parts(c(NA, 30)/100) # Make a figure with 2 stacked plots
#' fig.part(1, new=TRUE)
#' plot(-10:10, dnorm(-10:10), type='l', main='Part 1')
#' fig.part(2)
#' plot(-10:10, dnorm(-10:10), main='Part 2')
#' 
#' @export
fig.parts <- function(percents) {
  n <- length(percents)
  n.fill <- sum(is.na(percents))
  percents[is.na(percents)] <- (1-sum(percents, na.rm=TRUE))/n.fill
  fig.coords <- cbind(0, 1, 1-c(cumsum(percents)), 1-c(0, cumsum(percents)[-n]))
  fig.part <- function (i, new=FALSE) {
    par(fig=fig.coords[i, ], new=!new)
  }
  return ( fig.part )
}