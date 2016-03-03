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
#' @param all.inside If TRUE, out-of-range points are colored with the ends of
#' the color scale. If FALSE, they will be assigned no color (NA). 
#' @seealso \code{\link{draw.scale}} for drawing a color scale
#' @export
assign.scale.colors <- function(x, scale.colors, scale.range, all.inside=TRUE) {
  x.missing <- is.na(x)
  
  x[x.missing] <- mean(x, na.rm=TRUE)
  
  if ( missing(scale.range) ) scale.range <- range(x)
  xbins <- seq(min(scale.range), max(scale.range), length.out=length(scale.colors)+1)
  color.idx <- findInterval(x, xbins, rightmost.closed=TRUE, all.inside=all.inside)
  if ( scale.range[1] > scale.range[2] )
    color.idx <- 1+length(scale.colors)-color.idx
  
  
  x.colors <- scale.colors[color.idx]
  x.colors[x.missing | color.idx == 0 | color.idx > length(scale.colors)] <- NA
  return ( x.colors )
}




#' Draw a color scale bar
#' 
#' Adds a color scale legend to the current plot with specified colors and range.
#' 
#' 
#' @param scale.colors a vector of colors
#' @param scale.range range for scale (vector with 2 numeric elements, e.g. as returned by \code{range()})
#' @param num.labs (default 6) number of text labels to draw
#' @param pos (default topright)position for scale, either a keyword such as "top", "topleft" or a numeric vector length 2
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
#' @details This function draws a color scale bar in the current plot. Position argument
#' \code{pos} can be specified by keyword: \code{'center'}, \code{'middle'},
#' \code{'topleft'}, \code{'topright'}, \code{'top'}, \code{'bottomleft'},
#' \code{'bottomright'}, \code{'bottom'}, \code{'left'}, or \code{'right'}.
#' Alternately \code{pos} can be set by specifying a numeric vector of length 2,
#' describing the x and y locations of legend, where each is a number between 0
#' (left/bottom) and 1 (right/top). The \code{adj} argument controls the anchor
#' point for the legend itself relative to the plot. If unspecified, it will
#' match \code{pos}, if \code{outside} is \code{FALSE} (default) or \code{1-pos}
#' if \code{outside} is \code{TRUE}. This alignment is offset by \code{scale.offset}
#' so that the legend does not overlap the plot border.
#' 
#' The position can be further tweaked using the \code{x.shift} and \code{y.shift}
#' arguments. 1 will shift the legend right/up by 1 scale width/height.
#' 
#' The size and shape of the color legend are controlled by \code{size}, which
#' is the approximate length of the longer dimension in inches, and \code{ratio},
#' which is the ratio between width and height (or vice versa for vertical scale bars)
#' and should generally be > 1.
#' 
#' There are several arguments that control how the legend is drawn. The
#' \code{tick.length} (between 0 and 1) argument controls how much of the shorter
#' dimension of the legend is taken up by the color scale vs the ticks. \code{label.offset}
#' controls the spacing between the ticks and the text labels. Parameters to modify
#' text size should be set via \code{par()} prior to calling this function. Setting
#' \code{box = FALSE} will remove the black outline around the colors.
#' 
#' 
#' @seealso \code{\link{draw.chrom.axis}} for drawing x-axis 
#' @examples
#' plot(1:5, 1:5, col=gray(0:4/5), pch=15)
#' # Place a horizontal scale bar a the top, inside the plot.
#' draw.scale(gray(0:4/5), c(0, 1), pos='top')
#' # Place a vertical scale bar to the right, outside the plot.
#' draw.scale(gray(0:4/5), c(76.1, 76.92), pos='right', horiz=FALSE, outside=TRUE)
#' # Place a horizontal  scale bar at the top left, outside the plot, but
#' aligned with the left edge of the plot
#' draw.scale(gray(0:4/5), c(1, 10), pos='topleft', adj=c(0, 0), outside=TRUE)
#' @export
draw.scale <- function(scale.colors, scale.range, num.labs=min(length(scale.colors)+1, 6),
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
  if ( length(x) < 2 || length(unique(x)) < 2 ) {
    warning('expected x with at least 2 unique values.\n')
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


#' Calculate offsets for plotting y ~ x scatter plots where y is continuous and x is categorical
#' 
#' 
#' @examples
#' rx <- ceiling(runif(250, 0, 5))
#' ry <- rnorm(250, 0, 100)+runif(5, 0, 2000)[rx]
#' offset <- splitter(ry, rx)
#' plot(offset+as.numeric(rx), ry, pch=20, col=rainbow(5)[rx])
#' @export
#' @export
splitter <- function (y, x=NULL, rad=.025, scale=TRUE) {
  zx <- rep(0, length(y))
  if ( length(y) < 2 ) return (zx)
  
  if ( !is.null(x) ) {
    if ( length(unique(x)) > length(x)/2 ) warning('x does not appear to be categorical')
    
    
    if ( scale )
      z <- (y-min(y))/diff(range(y))
    
    subs <- tapply(z, x, splitter, rad=rad, scale=FALSE, simplify=FALSE)
    subidx <- tapply(1:length(y), x)
    
    for ( s in 1:length(subs) ) {
      zx[subidx==s] <- subs[[s]]*length(subs)/2
    }
    return (zx)
  }
  
  y.order <- order(y)
  z <- y[y.order]
  if ( scale )
    z <- (z-min(z))/diff(range(z))
  
  for  ( i in 2:length(z) ) {
    
    dz <- z[i]-z[1:(i-1)]
    if ( any(dz < rad) ) {
      nbs <- (1:(i-1))[dz < rad]
      nbd <- sqrt((z[nbs]-z[i])^2+(zx[nbs]-zx[i])^2)
      if ( any(nbd < rad) ) {
        dz <- z[i]-z[nbs]
        ax <- sin(acos(dz/rad))*rad*1.01
        if ( mean(zx[nbs]) < 0 )
          nbx <- zx[nbs]+ax
        else
          nbx <-  zx[nbs]-ax
        
        for ( j in order(abs(nbx)) ) {
          zx[i] <- nbx[j]
          nbd <- sqrt((z[nbs]-z[i])^2+(zx[nbs]-zx[i])^2)
          if ( all(nbd >= rad) ) {
            break
          }
        }
        
        if ( any(nbd < rad) )
          cat('!')
        
      }
      
      
      zx[1:i] <- zx[1:i]-mean(zx[1:i])
    }
  }
  zx[y.order] <- zx
  return (zx)
}















#' Integer data binning
#' 
#' Sometimes with integer data which is not uniformly distributed, a linear
#' color scale is not ideal. The function \code{bin.ints} converts an integer
#' vector into a factor vector where each level corresponds to a range of
#' integers and the elements of x are evenly distributed (as much as possible)
#' across the factor levels.
#' 
#' This is a convenience function that makesuse of \code{find.bins} and
#' \code{label.bins}.
#' 
#' @param x data to be binned
#' @param num.bins number of bins to aim for
#'   
#'   
#' @return
#' \code{bin.ints} returns an ordered factor vector with appropriately
#' labeled levels.
#' @examples
#' x <- rpois(100, 10)
#' x.bin <- bin.ints(x)
#' 
#' 
#' @seealso \code{\link{find.bins}} for determining bin breakpoints
#' @seealso \code{\link{label.bins}} for bin labels suitable for a legend
#' @export
bin.ints <- function(x, num.bins=10) {
  bps <- find.bins(x, num.bins)
  xbin <- findInterval(x, bps)
  binlab <- label.bins(x, xbin)
  factor(xbin, labels=binlab, ordered=TRUE)
}



#' Find breakpoints for binning numeric values
#' 
#' @param x data to be binned
#' @param num.bins number of bins to aim for
#'   
#' @examples
#' x <- rpois(100, 10)
#' bins <- find.bins(x)
#' x.bin <- findInterval(x, bins)
#' 
#' 
#' @seealso \code{\link{label.bins}} for bin labels suitable for a legend
#' @seealso \code{\link{bin.ints}} for generating a factor with labels
#' @export
find.bins <- function(x, num.bins=10) {
  bins <- NULL
  x.table <- table(x)
  names(x.table) <- NULL
  x.levels <- sort(unique(x))
  x.levels <- (x.levels + c(x.levels[-1], 1+max(x)))/2
  
  i <- 1
  while ( length(bins) < (num.bins-1)  & i <= length(x.table)) {
    if ( sum(x.table[1:i])/sum(x.table) > 1/(num.bins-length(bins)) ) {
      i <- max(1, i - 1)
      bins <- c(bins, x.levels[i])
      x.table <- x.table[-(1:i)]
      x.levels <- x.levels[-(1:i)]
      i <- 1
    } else {
      i <- i + 1
    }
  }
  bins <- c(bins, max(x)+1)
  return (bins)
}


#' Label bins generated by \code{find.bins()}
#' 
#' @param x data to be binned
#' @param bin binned data
#' @param greedy (default TRUE) whether to include in labels values of x which
#' may not appear in the data. e.g. instead of \code{0, 1, 2-3, 5-6} we return
#' \code{0, 1, 2-4, 5-6}
#' 
#' @examples
#' x <- rpois(100, 10)
#' bins <- find.bins(x)
#' x.bin <- findInterval(x, bins)
#' bin.labels <- label.bins(x, x.bin)
#' 
#' @seealso \code{\link{find.bins}} for determining bin breakpoints
#' @seealso \code{\link{bin.ints}} for generating a factor with labels
#' @export
label.bins <- function(x, bins, greedy=TRUE) {
  if ( missing(bins)  ) {
    x.bins <- find.bins(x)
    bins <- findInterval(x, x.bins)
  } else if ( length(bins) == 1 && is.numeric(bins) ) {
    x.bins <- find.bins(x, bins)
    bins <- findInterval(x, x.bins)
  } else if  ( length(bins) != length(x) ) {
    stop('x and bins must have same length.')
  }
  
  x.ranges <- simplify2array(tapply(x, bins, range))
  
  if ( greedy && min(diff(sort(unique(x)))) >= 1 )
    x.ranges[2, -ncol(x.ranges)] <- x.ranges[1, -1]-1
  
  apply(x.ranges, 2, function (r) 
    if (diff(r)) paste0(r[1], if (r[2]==max(x)) '+' else paste0('-', r[2])) else r[1]
  )
}





#' @export
gw.snp.pos <- function(chromosome, position, spacing=0.1) {
  snp.order <- gtools::mixedorder(chromosome, position)
  
  co <- chromosome[snp.order]
  po <- position[snp.order]
  
  chr.ranges <- do.call(rbind, tapply(po, co, range))
  chr.names <- rownames(chr.ranges)
  chr.sizes <- apply(chr.ranges, 1, diff)
  space <- mean(chr.sizes)*spacing
  chr.bounds <- unname(cumsum(c(0, chr.sizes+spacing)))
  
  chr.offsets <- chr.bounds[1:nrow(chr.ranges)] - chr.ranges[, 1]
  
  gwpos <- chr.offsets[match(co, chr.names)] + po
  
  attr(gwpos, 'chr.names') <- chr.names
  attr(gwpos, 'chr.bounds') <- chr.bounds
  
  return (gwpos)
}

#' @export
gwaxis <- function(names, bounds) {
  midpts <- (bounds[-1] + bounds[-length(bounds)])/2
  axis(1, bounds, labels=FALSE)
  axis(1, midpts, labels=names, lwd=0)
}

#' @export
gwplot <- function (x, ...) {
  if ( is.null(attr(x, 'chr.names')) || is.null(attr(x, 'chr.bounds')) ) {
    warning('No gw.snp.pos attributes found.')
    return ( plot(x, ...) )
  }
  plot(x, xaxt='n', ... )
  gwaxis(attr(x, 'chr.names'), attr(x, 'chr.bounds'))
}

