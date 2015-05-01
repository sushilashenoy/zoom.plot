#' @export
extend.color.range <- function(colors, n, weight=rep(1, length(colors)-1)) {
  if ( n < length(colors) ) return ( colors )
  if ( length(weight) != length(colors)-1 ) stop('Must be one fewer weights than colors.')
  
  red.part <- strtoi(paste('0X', substring(colors, 2, 3), sep=''))/2^8
  grn.part <- strtoi(paste('0X', substring(colors, 4, 5), sep=''))/2^8
  blu.part <- strtoi(paste('0X', substring(colors, 6, 7), sep=''))/2^8
  
  new.red <- approx(c(0, cumsum(weight)/sum(weight)), red.part, n=n)$y
  new.grn <- approx(c(0, cumsum(weight)/sum(weight)), grn.part, n=n)$y
  new.blu <- approx(c(0, cumsum(weight)/sum(weight)), blu.part, n=n)$y
  
  return ( rgb(new.red, new.grn, new.blu) )
}

# # Example!
# library('RColorBrewer')
# wee.colors <- brewer.pal(9, 'Spectral')
# wee.colors <- c(rgb(0, 0:3/3, 0), rgb(1:3/3, 1, 1:3/3))
# wee.more.colors <- extend.color.range(wee.colors, 100)
# wee.nonlinear.colors <- extend.color.range(wee.colors, 100, c(1:8))
# 
# par(mfcol=c(1, 3))
# image(matrix(1:100, 1), col=wee.colors)
# image(matrix(1:100, 1), col=wee.more.colors)
# image(matrix(1:100, 1), col=wee.nonlinear.colors)

#' @export
lighten.colors <- function(colors, pct) {
  if ( length(pct) == 1 ) pct <- rep(pct, length(colors))
  
  red.part <- strtoi(paste('0X', substring(colors, 2, 3), sep=''))/2^8
  grn.part <- strtoi(paste('0X', substring(colors, 4, 5), sep=''))/2^8
  blu.part <- strtoi(paste('0X', substring(colors, 6, 7), sep=''))/2^8
  
  new.red <- pct+(1-pct)*red.part
  new.grn <- pct+(1-pct)*grn.part
  new.blu <- pct+(1-pct)*blu.part
  
  return ( rgb(new.red, new.grn, new.blu) )
}


#' @export
darken.colors <- function(colors, pct) {
  if ( length(pct) == 1 ) pct <- rep(pct, length(colors))
  
  red.part <- strtoi(paste('0X', substring(colors, 2, 3), sep=''))/2^8
  grn.part <- strtoi(paste('0X', substring(colors, 4, 5), sep=''))/2^8
  blu.part <- strtoi(paste('0X', substring(colors, 6, 7), sep=''))/2^8
  
  new.red <- (1-pct)*red.part
  new.grn <- (1-pct)*grn.part
  new.blu <- (1-pct)*blu.part
  
  return ( rgb(new.red, new.grn, new.blu) )
}