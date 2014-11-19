#' Print a bp length nicely
#' 
#' \code{prettybp} returns a string representation of a base pair length with optional rounding
#' 
#' @param n base pairs
#' @param round.digits (optional) number of digits to round to. Default is 2. Set
#' to NULL to disable rounding. This value is ignored if signif.digits is not \code{NULL}.
#' @param signif.digits (optional) set this if rounding to significant digits
#' is preferred.
#' @param space Set to false to eliminate the space between the number and unit.
#' @seealso \code{\link{draw.chrom.axis}} for drawing a bp scaled x-axis 
#' @export
prettybp <- function (n, round.digits=2, signif.digits=NULL, space=TRUE) {
  if ( space ) {
    use.sep <- ' '
  } else {
    use.sep <- ''
  }
  if ( log10(mean(n)) >= 6 ) {
    if ( !is.null(signif.digits) ) {
      return ( paste(signif(n/1e6, signif.digits), 'Mb', sep=use.sep))
    } else if ( !is.null(round.digits) ) {
      return ( paste(round(n/1e6, round.digits), 'Mb', sep=use.sep))
    }
    return ( paste(n/1e6, 'Mb', sep=use.sep))
  } else if ( log10(n) >= 3) {
    if ( !is.null(signif.digits) ) {
      return ( paste(signif(n/1e3, signif.digits), 'kb', sep=use.sep))
    } else if ( !is.null(round.digits) ) {
      return ( paste(round(n/1e3, round.digits), 'kb', sep=use.sep))
    }
    return ( paste(n/1e3, 'kb', sep=use.sep))
  }
  if ( !is.null(signif.digits) ) {
    return ( paste(signif(n, signif.digits), 'bp', sep=use.sep))
  } else if ( !is.null(round.digits) ) {
    return ( paste(round(n, round.digits), 'bp', sep=use.sep))
  }
  return ( paste(n, 'bp', sep=use.sep) )
}

#' Draw a horizontal chromosome axis
#' 
#' \code{draw.chrom.axis} draws a chromosome axis with appropriate scale (Mb, kb, or bp).
#' 
#' @param start.pos Starting position on the chromosome (in bp)
#' @param end.pos Ending position on the chromosome (in bp)
#' @param ... (optional) addition options to pass to \code{text} for drawing
#' labels
#' @seealso \code{\link{draw.scale}} for drawing a color scale bar
#' @export
draw.chrom.axis <- function(start.pos, end.pos, chrom=NULL, label.chrom=TRUE, label.scale=TRUE, ...) {
  plot(0, type='n', ylim=c(-1, 1), xlim=c(start.pos, end.pos),
       axes=FALSE, bty='n', xlab='', ylab='', yaxs='i')
  
  abline(h=0, xpd=FALSE)
  
  ticks.at <- axTicks(1)
  plot.height <- par('pin')[2]
  sapply(ticks.at, function (x) { lines(c(x, x), c(-0.1, 0.1)/plot.height)})
  text(ticks.at, rep(0, length(ticks.at)), ticks.at/1e6, pos=1, xpd=NA, ...)
  if ( label.scale ) {
    text(par('usr')[2], 0, 'Mb', pos=4, xpd=NA, ...)
  }
  if ( !is.null(chrom) && label.chrom ) {
    text(par('usr')[1], 0, chrom, pos=2, xpd=NA, ...)
  }
}


#' Assign colors to value according to a scale
#' 
#' \code{draw.scale} takes a vector of values, a color vector, and a range
#' vector and returns colors for plotting those values
#' 
#' 
#' @param x values to be plotted
#' @param scale.colors a vector of colors
#' @param scale.range range for scale (vector with 2 numeric elements)
#' @seealso \code{\link{draw.scale}} for drawing a color scale
#' @export
assign.scale.colors <- function(x, scale.colors, scale.range) {
  x.missing <- is.na(x)
  
  x[x.missing] <- mean(x, na.rm=TRUE)
  
  if ( missing(scale.range) ) scale.range <- range(x)
  x.scaled <- (x-scale.range[1])/diff(scale.range)
  color.idx <- 1+floor(x.scaled*(length(scale.colors)-1))
  color.idx <- pmax(pmin(color.idx, length(scale.colors)), 1)
  
  x.colors <- scale.colors[color.idx]
  x.colors[x.missing] <- NA
  return ( x.colors )
}


#' Draw a color scale bar
#' 
#' \code{draw.scale} adds a color scale to the corner of the current plot
#' with specified colors and range
#' 
#' 
#' Currently only a horizontal scale bar is supported. To adjust position in
#' plot coordinates use x.offset and y.offset. To adjust relative position use
#' x.shift and y.shift.
#' 
#' Keep in mind that if you shift the scale bar away from the plotting area
#' (i.e. into the margins), you may need to specify the additional parameter
#' xpd=NA which will be passed to the plotting commands and allows drawing
#' in the margins.
#' 
#' @param scale.colors a vector of colors
#' @param scale.range range for scale (vector with 2 numeric elements)
#' @param num.labs (optional) number of labels to draw
#' @param position (optional) either 'topleft' 'topright' 'bottomleft' or
#' 'bottomright'
#' @param size (optional) approximate length in inches
#' @param width.to.height (optional) ratio of scale bar width to height
#' @param x.offset (optional) adjust x position (uses \code{par('usr')} scale).
#' @param y.offset (optional) adjust y position (uses \code{par('usr')} scale).
#' @param x.shift (optional) adjust x position relative to scale size.
#' @param y.shift(optional) adjust y position relative to scale size.
#' @param ... (optional) addition options to pass to plotting commands for
#' drawing labels, lines and shapes
#' 
#' @seealso \code{\link{draw.chrom.axis}} for drawing x-axis 
#' @export
draw.scale <- function(scale.colors, scale.range, num.labs=6,
                       position='topleft', size=3, width.to.height=20,
                       x.offset, y.offset, x.shift, y.shift, ...) {
  par.usr <- par('usr')
  par.pin <- par('pin')
  
  if ( substr(position, 1, 3) == 'top' ) {
    y1 <- par.usr[4]
    y2 <- par.usr[4] - (par.usr[4]-par.usr[3])/(par.pin[2])*size/width.to.height
    top <- TRUE
  } else {
    y1 <- par.usr[3]
    y2 <- par.usr[3] + (par.usr[4]-par.usr[3])/(par.pin[2])*size/width.to.height
    top <- FALSE
  }
  
  if ( substr(position, 7-3*top, 11-3*top) == 'left' ) {
    x1 <- par.usr[1]
    x2 <- par.usr[1] + (par.usr[2]-par.usr[1])/(par.pin[1])*size
  } else {
    x2 <- par.usr[2]
    x1 <- par.usr[2] - (par.usr[2]-par.usr[1])/(par.pin[1])*size
  }
  
  if ( missing(x.offset) ) x.offset <- 0
  if ( missing(y.offset) ) y.offset <- 0
  
  if ( !missing(x.shift) ) {
    x.offset <- x.offset + x.shift*(y2-y1)
  }
  
  if ( !missing(y.shift) ) {
    y.offset <- y.offset + y.shift*(y2-y1)
  }
  
  y1 <- y1 + y.offset
  y2 <- y2 + y.offset
  x1 <- x1 + x.offset
  x2 <- x2 + x.offset
  
  x.points <- seq(x1+(x2-x1)*0.1, x1+(x2-x1)*0.9, length.out=length(scale.colors)+1)
  y.mid <- y1+(y2-y1)/2
  

  
  x.ticks <- seq(x1+(x2-x1)*0.1, x1+(x2-x1)*0.9, length.out=num.labs)
  arrows(x.ticks, rep(y1, num.labs), x.ticks, rep(y2, num.labs), length=0, ...)
  
  if ( substr(position, 1, 3) == 'top' ) {
    label.pos = 1
  } else {
    label.pos = 3
  }
  
  sapply(1:length(scale.colors), function (i) {
    x1 <- x.points[i]
    x2 <- x.points[i+1]
    polygon(c(x1, x1, x2, x2), c(y1, y.mid, y.mid, y1), border=scale.colors[i], col=scale.colors[i], ...)
  })
  
  text(x.ticks, rep(y2, num.labs), seq(scale.range[1], scale.range[2], length.out=num.labs),
       pos=label.pos, offset=0.25, ...)
  
}

#' @export
center.out.order <- function(n) {
  if ( n <= 2 )
    return (1:n)
  if ( n %% 2 == 1 ) {
    odd <- seq(1, n, 2)
    even <- seq((n-1), 2, -2)
  } else {
    odd <- seq(1, (n-1), 2)
    even <- seq(n, 2, -2)
  }
  return ( order(c(even, odd)) )
}

#' Calculate a range with margins (padding)
#' 
#' @export
mrange <- function(x, m=0.1) {
  if ( length(x) < 2 ) {
    cat('expected x with length of at least 2.\n')
    return ( x )
  }
  if ( length(x) > 2 ) {
    x <- range(x)
  }
  
  dx <- diff(x)
  x[1] <- x[1] - m*dx
  x[2] <- x[2] + m*dx
  
  return ( x )
}