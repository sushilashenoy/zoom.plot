

#' @export
get.encode.region <- function(chrom, start.pos, end.pos) {
  data('encode_tfbs')
  encode <- encode.tfbs
  
  if ( is.numeric(chrom) || nchar(chrom) == 1 ) {
    chrom <- paste('chr', chrom, sep='')
  }
  
  encode.rows <- encode$chrom == chrom & encode$end > start.pos & encode$start < end.pos
  
  cat('Found', sum(encode.rows), 'encode elements in region.\n')
  
  encode.result <- encode[encode.rows, ]
  
  attr(encode.result, 'start.pos') <- start.pos
  attr(encode.result, 'end.pos') <- end.pos
  
  return ( encode.result )
}

#' @export
arrange.encode.rows <- function(encode.region, start.pos, end.pos, labels=TRUE) {
  if ( is.null(encode.region) ) return ( NULL )
  
  if ( missing(start.pos) ) start.pos <- attr(encode.region, 'start.pos')
  if ( missing(end.pos) ) end.pos <- attr(encode.region, 'end.pos')
  
  char.width <- 0.5 * (end.pos-start.pos)/(par('pin')[1]/par('cin')[1])
  
  row.occupied <- rep(0, nrow(encode.region))
  encode.rows <- rep(0, nrow(encode.region))
  # Loop over elements in order based on start position
  for ( i in order(encode.region$start) ) {
    # Starting from row 1, search for a row where this gene will fit
    row <- 1
    label.width <- char.width * (nchar(encode.region$name[i])+1) * labels
    while ( row.occupied[row] + label.width > encode.region$start[i] ) {
      row <- row + 1
    }
    # This row is now occupied up to the end of this gene
    row.occupied[row] <- encode.region$end[i]
    encode.rows[i] <- row
    
#     print(round(100*(row.occupied[1:max(encode.rows)]-start.pos)/(end.pos-start.pos), 1))
    
#     if ( i > 9 ) { break }
  }
  
  return ( encode.rows )
}


#' @export
plotencode <- function(encode.region, encode.rows, start.pos, end.pos,
                        abs.scale=TRUE, labels=TRUE) {
  if ( is.null(encode.region) ) return ( invisible(NULL) )
  
  if ( missing(start.pos) ) start.pos <- attr(encode.region, 'start.pos')
  if ( missing(end.pos) ) end.pos <- attr(encode.region, 'end.pos')

  
  plot(0, type='n', ylim=0.5+c(0, max(encode.rows)), xlim=c(start.pos, end.pos),
       axes=FALSE, bty='n', xlab='', ylab='', yaxs='i')
  if ( abs.scale ) {
    encode.color <- gray(1-encode.region$score/1000)
  } else {
    max.score <- max(encode.region$score)
    min.score <- min(encode.region$score)
    encode.color <- gray(1-(encode.region$score - min.score)/(max.score - min.score))
  }
  
  sapply(1:nrow(encode.region), function (i) {
    lines(c(encode.region$start[i], encode.region$end[i]),
          c(encode.rows[i], encode.rows[i]),
          col=encode.color[i],
          lwd=8, lend=1)
    if ( labels ) {
      text(encode.region$start[i], encode.rows[i],
           encode.region$name[i],
           pos=2, offset=0.25, family='mono', cex=0.5)
    }
  })
  
  invisible(NULL)
}
