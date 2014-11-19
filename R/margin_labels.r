
shorten <- function (x0, y0, x1, y1, rad) {
  line.lengths <- sqrt((x1-x0)^2+(y1-y0)^2)
  pct.short <- rad/line.lengths
  
  new.pt <- list()
  new.pt$x1 <- x0+(x1-x0)*pct.short
  new.pt$x0 <- x1+(x0-x1)*pct.short
  new.pt$y1 <- y0+(y1-y0)*pct.short
  new.pt$y0 <- y1+(y0-y1)*pct.short
  
  return (new.pt)
}

#' Label points from the margin
#' 
#' @export
marginlabels <- function(x, y = NULL, labels=seq_along(x), margin=4,
                         col='black', lty=1, lwd=1, pch=1, pch.cex=1, las=2,
                         rad=0.15, ...) {
  
  len <- length(labels)
  
  if ( missing(y) || is.null(y) ) {
    y <- seq_along(labels)
  }
  
  if ( length(x) != len ) x <- rep(x, len)
  if ( length(y) != len ) y <- rep(y, len)
  
  if ( margin == 1 || margin == 3 ) {
    new.order <- order(x)
    x <- x[new.order]
    y <- y[new.order]
    if ( !missing(col) && length(col) == len ) col <- col[new.order]
    label.x <- seq(par('usr')[1], par('usr')[2], length.out=len+2)[-c(1, len+2)]
    label.y <- par('usr')[if ( margin ==  2) 3 else 4]
    tick.pos <- label.x
  } else {
    new.order <- order(y)
    x <- x[new.order]
    y <- y[new.order]
    if ( !missing(col) && length(col) == len ) col <- col[new.order]
    label.y <- seq(par('usr')[3], par('usr')[4], length.out=len+2)[-c(1, len+2)]
    label.x <- par('usr')[if ( margin ==  1) 1 else 2]
    tick.pos <- label.y
  }
  
  points(x, y, pch=pch, col=col, cex=pch.cex)
  
  connect.lines <- shorten(x, y, label.x, label.y, pch.cex*rad)
  connect.lines$lty <- lty
  connect.lines$lwd <- lwd
  connect.lines$x0 <- label.x
  connect.lines$y0 <- label.y
  
  do.call(segments, connect.lines)
  
  axis(margin, at=tick.pos, labels, las=las, lwd=0, line=-0.5)
}