#' Y axis scale bar for stacked plots
#' 
#' In stacked plots the standard y axis doesn't really work - it is vertical
#' lengths of line segments that correspond to data values, not the absolute y
#' positions. This function draws a scale bar in the space where the yaxis
#' would normally be.
#' 
#' @param length (Optional) if not supplied, the length will be the distance
#' between two ticks if the y axis had been plotted in the standard way (uses
#' \code{par('yaxp')}). 
#' @param pos Vertical location for scale bar where 0 is at the bottom, 1 is at
#' the top, and 0.5 (default) is in the middle.
#' @param las Controls orientation of axis label text (See \code{par('las')})
#' @param ... Other arguments to pass to axis (See \code{par()})
#' 
#' @seealso \code{\link{stack.plot}} for drawing stacked plots 
#' @export
y.scale.bar <- function(length, pos=0.5, abs.pos, las=par('las'), ...) {
  if ( missing(length) )
    length <- diff(par('yaxp')[1:2])/par('yaxp')[3]
  
  if ( missing(abs.pos) )
    start <- par('usr')[3]+(diff(par('usr')[3:4])-length)*pos
  else
    start <- abs.pos
  
  end <- start+length
  mid <- start+length/2
  
  axis(2, c(start, end), label=FALSE, lwd=0, lwd.ticks=1, mgp=c(0, 0, par('tcl')), ...)
  axis(2, c(start, end), label=FALSE, lwd=1, lwd.ticks=0, mgp=c(0, 0, 0), ...)
  axis(2, mid, label=length, lwd=0, las=las, mgp=c(0, -par('tcl'), 0), ...)
}


#' Stacked plots
#' 
#' Stacked plots allow plotting multiple vectors of data without overlap in
#' minimal space.
#' 
#' @param y a matrix of values. Each row of the matrix will be one line in the
#' plot.
#' @param x (Optional) x axis values - length should correspond to columns in y matrix
#' @param spacing Spacing between lines
#' @param style normal or zero
#' 
#' @export
stackplot <- function(y, x, spacing=0.5, style='normal', las=1, y.scale=TRUE,
                         lwd, lty, col, bty='n', xaxt=par('xaxt'),
                         yaxt='n', yaxs='i', xlab='', ylab='') {
  if ( missing(x) )
    x <- 1:ncol(y)
  
  recognized.styles <- c('normal', 'zero')
  if ( ! style %in% recognized.styles ) {
    warning('style ', style, ' not recognized. Using default (normal).')
  }
  
  nl <- nrow(y)
  y <- y[nl:1, ]
  
  if ( missing(col) ) {
    col <- hcl(1:nl/nl*300, 60, 40)
  } else if ( length(col) != nl ) {
    col <- rep(col, length.out=nl)
  } else {
    col <- col[nl:1]
  }
  
  if ( missing(lwd) ) {
    lwd <- rep(par('lwd'), nl)
  } else if ( length(lwd) != nl ) {
    lwd <- rep(lwd, length.out=nl)
  } else {
    lwd <- lwd[nl:1]
  }
  
  if ( missing(lty) ) {
    lty <- rep(par('lwd'), nl)
  } else if ( length(lty) != nl ) {
    lty <- rep(lty, length.out=nl)
  } else {
    lty <- lty[nl:1]
  }
  
  if ( style=='zero' ) {
    ymins <- pmin(0, apply(y, 1, min))
    ymaxs <- pmax(0, apply(y, 1, max))
  } else {
    ymins <- apply(y, 1, min)
    ymaxs <- apply(y, 1, max)
  }
  
  yheight <- ymaxs-ymins
  
  add.space <- sum(yheight)*spacing/nl
  offset <- c(0, cumsum(yheight) + cumsum(rep(add.space, nl))) - c(ymins, 0)
  
  plot(rep(x, each=nl), as.vector(y)+offset[1:nl], type='n', bty=bty, yaxt=yaxt, xlab=xlab, ylab=ylab, xaxt=xaxt, yaxs=yaxs)
  if ( style=='zero'  ) {
    sapply(1:nl, function (i) segments(x, offset[i], x, offset[i]+y[i, ], col=col[i], lty=lty[i], lwd=lwd[i]))
    segments(min(x), offset[1:nl], max(x), offset[1:nl], col=col, lty=lty, lwd=lwd)
  } else {
    sapply(1:nl, function (i) lines(x, y[i, ]+offset[i], col=col[i], lty=lty[i], lwd=lwd[i]))
  }
  
  if ( y.scale ) y.scale.bar(las=las)
  invisible(offset)
}


y <- matrix(rchisq(10*150, 1), 10)*rchisq(10, 1)

stackplot(y, style='zero', lwd=2, spacing=0.2, y.scale=FALSE)
y.scale.bar(pos=1.1, xpd=NA)
# 
# 
# simple.stack(y, style='zero', lwd=2, spacing=0.2)
