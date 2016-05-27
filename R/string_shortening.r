#' Remove longest common prefix and suffix
#' @export
rlcps <- function(x) {
  require('Biostrings')
  xo <- order(x)
  s1 <- x[xo[1]]
  s2 <- x[xo[length(x)]]
  plen <- Biostrings::lcprefix(s1, s2)
  slen <- Biostrings::lcsuffix(s1, s2)
  result <- substring(x, 1+plen, nchar(x)-slen)
  attr(result, 'prefix') <- substring(s1, 1, plen)
  attr(result, 'suffix') <- substring(s1, nchar(s1)-slen+1)
  return (result)
}

#' Find a minimal unique representation
#' @export
mur <- function(x) {
  require('Biostrings')
  require('stringr')
  if ( any(duplicated(x))) stop('Duplicate values not allowed.')
  
  clean <-  stringr::str_replace_all(x, '[._\\s]', '')
  if ( !any(duplicated(clean)) ) x <- clean
  
  n <- length(x)
  xo <- order(x)
  sx <- x[xo]
  mx <- rep('', n)
  
  mxs <- NULL

  i <- 1
  j <- n
  while ( j < n || i < n ) {
    lcp <- Biostrings::lcprefix(sx[i], sx[j])
    if ( lcp > 0 ) {
      mx[i:j] <- paste0(mx[i:j], substring(sx[i:j], 1, 1))
      sx[i:j] <- substring(sx[i:j], lcp+1)
      mxs <- rbind(mxs, sx)
    }
    
    if ( i <= 1 ) {
      i <- 2+n-j
      j <- n
    } else {
      i <- i - 1
      j <- j - 1
    }
  }
  
  for ( j in 2:n ) {
    i <- j-1
    if ( mx[i] == mx[j] ) {
      if ( nchar(sx[i]) > 0 ) {
        mx[i] <- paste0(mx[i], substring(sx[i], 1, 1))
        sx[i] <- substring(sx[i], 2)
      }
      if ( nchar(sx[j]) > 0 ) {
        mx[j] <- paste0(mx[j], substring(sx[j], 1, 1))
        sx[j] <- substring(sx[j], 2)
      }
    }
  }
  
  result <- mx
  result[xo] <- mx
  return (mx)
}


