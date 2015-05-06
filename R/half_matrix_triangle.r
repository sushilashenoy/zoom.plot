
#' @export
triangleplot <- function(mat, positions, colors, color.scale, drawscale=TRUE, scale.cex=0.6, ...) {
  
  if ( ncol(mat) != nrow(mat) ) stop('Matrix must be square.')
  
  n <- ncol(mat)
  
  if ( missing(positions) ) {
    positions <- 0:(n+1)
  } else if ( length(positions) != n && length(positions) != n+1 ) stop('positions must match matrix size')
  
  if ( is.unsorted(positions) ) warning('Positions are not in order!')
  
  if ( missing(colors) ) {
    if ( 'RColorBrewer' %in% rownames(installed.packages()) ) {
      require('RColorBrewer')
      colors <- extend.color.range(brewer.pal(11, 'RdYlGn'), 100)
    } else {
      colors <- grey(0:99/99)
    }
  }
  
  if ( missing(color.scale) ) color.scale <- range(mat, na.rm=TRUE)
  
  plot(0, type='n', ylim=c(min(positions)-max(positions), 0), 
       xlim=range(positions),
       axes=FALSE, xlab='', ylab='', yaxs='i', ...)
  
  # Calculate positions for SNPs
  if ( length(positions) == n) {
    x.pos <- positions
    x.pos.mids <- (x.pos[-1] + x.pos[-n])/2
    x.pos.left <- c(x.pos[1], x.pos.mids)
    x.pos.right <- c(x.pos.mids, x.pos[n])
  } else {
    x.pos.left <- positions[-(n+1)]
    x.pos.right <- positions[-1]
  }
  
  # Loop through each position in the matrix and draw polygon of appropriate shade
  for ( i in 1:n ) {
    for ( j in 1:i ) {
      fill.color <- assign.scale.colors(mat[i, j], colors, color.scale)
      # Draw a diamond shape 
      polygon(c(x.pos.left[i]+x.pos.left[j],
                x.pos.right[i]+x.pos.left[j],
                x.pos.right[i]+x.pos.right[j],
                x.pos.left[i]+x.pos.right[j])/2,
              -c(x.pos.left[i]-x.pos.left[j],
                 x.pos.right[i]-x.pos.left[j],
                 x.pos.right[i]-x.pos.right[j],
                 x.pos.left[i]-x.pos.right[j]),
              col=fill.color, border=fill.color, ljoin='bevel')
      # NOTE: not drawing polygon borders (i.e. border=NA) causes jagged,
      # pixelated edges. So we draw a border in the same color as the fill
    }
  }
  
  if ( drawscale ) {
    draw.scale(colors, color.scale, pos='bottomleft')
  }
}
