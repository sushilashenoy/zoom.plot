derez <- function(x, new.nrow, new.ncol, row.pos, col.pos) {
  # This function only decreases resolution, doesn't increase
  new.nrow <- min(new.nrow, nrow(x))
  new.ncol <- min(new.ncol, ncol(x))
  
  # Divide rows and cols into approximately equal bins
  row.end <- round(seq(0, nrow(x), length.out=new.nrow+1))[-1]
  row.beg <- c(1, row.end[-new.nrow]+1)
  col.end <- round(seq(0, ncol(x), length.out=new.ncol+1))[-1]
  col.beg <- c(1, col.end[-new.ncol]+1)
  
  # Vectorized function to calculate mean for each bin
  avg.block <- function(i, j) {
    sapply(1:length(i), function(k) {
      mean(x[row.beg[i[k]]:row.end[i[k]], col.beg[j[k]]:col.end[j[k]]])
    })
  }

  # Outer calls the vectorized function to generate new matrix
  new.x <- outer(1:new.nrow, 1:new.ncol, avg.block)

  if ( missing(row.pos) && missing(col.pos) ) {
    return ( new.x )
  }

  new.row.pos <- NULL
  new.col.pos <- NULL

  if ( !missing(row.pos) ) {
    new.row.pos <- sapply(1:new.nrow, function (i) mean(row.pos[row.beg[i]:row.end[i]]))
  }
  if ( !missing(col.pos) ) {
    new.col.pos <- sapply(1:new.ncol, function (i) mean(col.pos[col.beg[i]:col.end[i]]))
  }
  
  return ( list(x=new.x, row.pos=new.row.pos, col.pos=new.col.pos) )

}