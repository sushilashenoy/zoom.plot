

# Function to print out maximum intervals where vec==1
find.intervals <- function(vec) {
  ints <- NULL
  j <- 1
  n <- length(vec)
  for ( k in 1:sum(vec) ) {
    i <- which.max(vec[j:n])+j-1
    j <- which.min(vec[i:n])+i-1
    
    if ( i == j )
      break
    
    ints <- rbind(ints, c(i, j-1))
  }
  if ( vec[n] == 1 ) {
    ints <- rbind(ints, c(i, n))
  }
  return ( ints )
}



# Interpolate missing values..
interpolate <- function (vec, pos=1:length(vec)) {
  
  missing.intervals <- find.intervals(is.na(vec))
  if ( missing.intervals[1, 1] == 1 ) {
    missing.intervals <- missing.intervals[-1, ]
  }
  if ( missing.intervals[nrow(missing.intervals), 2] == length(vec) ) {
    missing.intervals <- missing.intervals[-nrow(missing.intervals), ]
  }
  
  interpolate.interval <- function (interval) {
    idxs <- (interval[1]-1):(interval[2]+1)
    context.values <- vec[idxs]
    context.position <- pos[idxs]
    nrc <- length(idxs)
    
    if ( is.na(context.values[1]) || is.na(context.values[nrc]) ) {
      warning("adjacent value is missing!")
      return ( FALSE )
    }
    
    left.val <- context.values[1]
    right.val <- context.values[nrc]
    left.pos <- context.position[1]
    right.pos <- context.position[nrc]
    
    for ( i in 2:(nrc-1) ) {
      i.pos <- context.position[i]
      new.i <- (i.pos - left.pos)/(right.pos - left.pos) * (right.val - left.val) + left.val
      
      if ( is.na(context.values[i]) ) { vec[idxs[i]] <<- new.i }
    }
    
    return ( TRUE )
  }
  
  apply(missing.intervals, 1, interpolate.interval)
  return ( vec )
}

#' @export
extend.color.range <- function(colors, n) {
  if ( n < length(colors) ) return ( colors )
  
  red.part <- strtoi(paste('0X', substring(colors, 2, 3), sep=''))
  grn.part <- strtoi(paste('0X', substring(colors, 4, 5), sep=''))
  blu.part <- strtoi(paste('0X', substring(colors, 6, 7), sep=''))
  
  new.red <- new.grn <- new.blu <- rep(NA, n)
  
  old.idx <- round(seq(1, n, length.out=length(colors)))
  new.red[old.idx] <- red.part
  new.grn[old.idx] <- grn.part
  new.blu[old.idx] <- blu.part
  
  new.red <- interpolate(new.red)/255
  new.grn <- interpolate(new.grn)/255
  new.blu <- interpolate(new.blu)/255
  
  return ( rgb(new.red, new.grn, new.blu) )
}

# # Example!
# library('RColorBrewer')
# wee.colors <- brewer.pal(9, 'Spectral')
# wee.more.colors <- extend.color.range(wee.colors, 25)
# 
# par(mfcol=c(1, 2))
# image(matrix(1:100, 1), col=wee.colors)
# image(matrix(1:100, 1), col=wee.more.colors)
