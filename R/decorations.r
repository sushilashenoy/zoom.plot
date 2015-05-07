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
draw.chrom.axis <- function(start.pos, end.pos, chrom=NULL, label.chrom=TRUE,
                            label.scale=TRUE, tick.length=0.1, ...) {
  plot(0, type='n', ylim=c(-1, 1), xlim=c(start.pos, end.pos),
       axes=FALSE, bty='n', xlab='', ylab='', yaxs='i')
  
  abline(h=0, xpd=FALSE)
  
  ticks.at <- axTicks(1)
  plot.height <- par('pin')[2]
  sapply(ticks.at, function (x) { lines(c(x, x), c(-tick.length, tick.length)/plot.height)})
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
#' Adds a color scale legend to the current plot
#' with specified colors and range
#' 
#' 
#' @param scale.colors a vector of colors
#' @param scale.range range for scale (vector with 2 numeric elements)
#' @param num.labs (optional) number of labels to draw
#' @param pos position for scale, either a keyword such as "top", "topleft" or a numeric vector length 2
#' @param adj controls the anchoring of the scale in respect to \code{pos}
#' @param horiz  TRUE for horizontal (default) or FALSE for vertical bar
#' @param outside TRUE to display legend within the plotting area (default) or FALSE to place it outside
#' @param size approximate length in inches
#' @param ratio ratio of scale bar width to height
#' @param tick.length number between 0 and 1, controls how long the ticks are vs the color boxes. 
#' @param scale.offset inset (outset) amount for legend that is inside (outside)
#' @param label.offset spacing between ticks and text labels
#' @param box this controls whether a black border is drawn around the colors or not.
#' @param x.shift (optional) adjust x position in units of scale width.
#' @param y.shift (optional) adjust y position in units of scale height.
#' 
#' @seealso \code{\link{draw.chrom.axis}} for drawing x-axis 
#' @examples
#' plot(1:5, 1:5, col=gray(0:4/5), pch=15)
#' draw.scale(gray(0:4/5), c(0, 1), pos='top')
#' draw.scale(gray(0:4/5), c(0, 1), pos='right', horiz=FALSE, outside=TRUE)
#' @export
draw.scale <- function(scale.colors, scale.range, num.labs=6,
                       pos='topleft', adj=NULL, horiz=TRUE,
                       outside=FALSE, size=2, ratio=12, tick.length=0.25,
                       scale.offset=0.5, label.offset=0.1,
                       box=TRUE, x.shift=0, y.shift=0) {
  
  
  if ( outside ) {
    old.par <- par(xpd=NA)
  }
  
  # recognized keywords and corresponding positions
  pos.key <- list()
  pos.key$middle <- c(0.5, 0.5)
  pos.key$center <- c(0.5, 0.5)
  pos.key$topleft <- c(0, 1)
  pos.key$topright <- c(1, 1)
  pos.key$top <- c(0.5, 1)
  pos.key$bottomleft <- c(0, 0)
  pos.key$bottomright <- c(1, 0)
  pos.key$bottom <- c(0.5, 0)
  pos.key$left <- c(0, 0.5)
  pos.key$right <- c(1, 0.5)
  pos.key <- data.frame(pos.key, row.names=c('x', 'y'))
  
  if ( is.character(pos) ) pos <- tolower(pos)
  if ( is.character(pos) && length(pos)==1 && pos %in% names(pos.key) ) {
    pos <- pos.key[[pos]]
  }
  
  if ( !is.numeric(pos) || length(pos) != 2) {
    warning('Invalid pos Please specify a keyword or xy pair.')
    pos <- c(0, 1)
  }
  
  # If adj is not specified it is the same as pos for inside
  # opposite for outside
  if ( is.null(adj) ) {
    if ( outside )
      adj <- 1-pos
    else
      adj <- pos
  }
  
  par.usr <- par('usr')
  par.pin <- par('pin')
  
  plot.width <- par.usr[2]-par.usr[1]
  plot.height <- par.usr[4]-par.usr[3]
  
  # Calculate dimensions of scale legend
  if ( horiz ) {
    legend.width <- plot.width/par.pin[1]*size
    legend.height <- plot.height/par.pin[2]*size/ratio
  } else {
    legend.height <- plot.height/par.pin[2]*size
    legend.width <- plot.width/par.pin[1]*size/ratio
  }
  
  x1 <- par.usr[1]+pos[1]*plot.width-adj[1]*legend.width
  y1 <- par.usr[3]+pos[2]*plot.height-adj[2]*legend.height
  
  # Apply scale.offset
  if ( horiz ) {
    y1 <- y1-scale.offset*legend.height*2*(adj[2]-0.5)
    x1 <- x1-scale.offset*legend.height*plot.width/plot.height*2*(adj[1]-0.5)
  } else {
    x1 <- x1-scale.offset*legend.width*2*(adj[1]-0.5)
    y1 <- y1-scale.offset*legend.width*plot.height/plot.width*2*(adj[2]-0.5)
  }
  x2 <- x1+legend.width
  y2 <- y1+legend.height
  
  # Flip depending on orientation
  if ( horiz & adj[2] > 0.5 ) {
    tm <- y1
    y1 <- y2
    y2 <- tm
  } else if ( !horiz & adj[1] > 0.5 ) {
    tm <- x1
    x1 <- x2
    x2 <- tm
  }
  #     points(x1, y1, pch=8, col='red')
  #     rect(x1, y1, x2, y2, col='#ff6633')
  
  if ( !missing(x.shift) ) {
    x1 <- x1 + x.shift*legend.width
    x2 <- x2 + x.shift*legend.width
  }
  if ( !missing(y.shift) ) {
    y1 <- y1 + y.shift*legend.height
    y2 <- y2 + y.shift*legend.height
  }
  
  if ( horiz ) {
    x.points <- seq(x1, x2, length.out=length(scale.colors)+1)
    y.mid <- y1*tick.length+y2*(1-tick.length)
    y.text <- y2+label.offset*(y2-y1)
    # draw ticks
    x.ticks <- seq(x1, x2, length.out=num.labs)
    segments(x.ticks, rep(y1, num.labs), x.ticks, rep(y2, num.labs))
    # draw colored boxes
    sapply(1:length(scale.colors), function (i) {
      rect(x.points[i], y1, x.points[i+1], y.mid, border=scale.colors[i], col=scale.colors[i])
    })
    if ( box ) rect(x1, y1, x2, y.mid)
  } else {
    y.points <- seq(y1, y2, length.out=length(scale.colors)+1)
    x.mid <- x1*tick.length+x2*(1-tick.length)
    x.text <- x2+label.offset*(x2-x1)
    # draw ticks
    y.ticks <- seq(y1, y2, length.out=num.labs)
    segments(rep(x1, num.labs), y.ticks, rep(x2, num.labs), y.ticks)
    # draw colored boxes
    sapply(1:length(scale.colors), function (i) {
      rect(x1, y.points[i], x.mid, y.points[i+1], border=scale.colors[i], col=scale.colors[i])
    })
    if ( box ) rect(x1, y1, x.mid, y2, border='#000000')
  }
  
  
  tick.labels <- seq(scale.range[1], scale.range[2], length.out=num.labs)
  sigdig <- 0
  while ( length(unique(round(tick.labels, sigdig))) < length(tick.labels) ) sigdig <- sigdig + 1
  if ( horiz ) {
    text(x.ticks, rep(y.text, num.labs), format(round(tick.labels, sigdig)),
         adj=c(0.5, round(adj[2])))
  } else {
    text(rep(x.text, num.labs), y.ticks, format(round(tick.labels, sigdig)),
         adj=c(round(adj[1]), 0.5)) 
  }
  
  if ( outside ) {
    par(old.par)
  }
}



#' Draw a color scale bar (old version)
#' 
#' \code{draw.old.scale} adds a color scale to the corner of the current plot
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
#' @param x.shift (optional) adjust x position in units of scale width.
#' @param y.shift (optional) adjust y position in units of scale height.
#' @param ... (optional) addition options to pass to plotting commands for
#' drawing labels, lines and shapes
#' 
#' @seealso \code{\link{draw.chrom.axis}} for drawing x-axis 
#' @export
draw.old.scale <- function(scale.colors, scale.range, num.labs=6,
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
    x.offset <- x.offset + x.shift*(x2-x1)
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
  
  
  tick.labels <- seq(scale.range[1], scale.range[2], length.out=num.labs)
  sigdig <- 1
  while ( length(unique(signif(tick.labels, sigdig))) < length(tick.labels) ) sigdig <- sigdig + 1
  
  text(x.ticks, rep(y2, num.labs), signif(tick.labels, sigdig),
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