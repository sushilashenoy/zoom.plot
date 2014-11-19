
require(parallel)


# Reverse the coding of the major and minor alleles
complement <- function( marker ) {
  return ( 2 - marker )
}

# MLE method for calculating D'
LD.D.prime.matrix <- function(X, mc.cores) {
  
  if ( min(X, na.rm=TRUE) < 0 && max(X, na.rm=TRUE) < 2 ) {
    warn('Converting -1/0/1 coded genotypes to 0/1/2')
    X <- X+1
  }
  
  # X need to be coded as 0, 1, 2
  stopifnot(min(X, na.rm=TRUE) >= 0, max(X, na.rm=TRUE) <= 2)
  
  # Check for columns with NA values
  na.indices <- lapply(1:ncol(X), function (j) { which(is.na(X[, j])) })
  
  # Complement all neccessary columns 
  for (j in 1:ncol(X)) {
    if ( length(na.indices[[j]]) == 0) {
      p.A <- sum(X[, j])/ (2*nrow(X))
    } else {  
      p.A <- sum(X[-na.indices[[j]], j])/ (2*nrow(X))
    }			
    
    if ( p.A < 0.5 ) {
      X[, j] <- complement( X[, j] )
    }
  }
  
  n.geno <- ncol(X)
  # Make a list of all of all possible pairs of markers
  idx.pairs <- data.frame('a'=rep(1:n.geno, n.geno), 'b'=rep(1:n.geno, each=n.geno))
  # Keep only the pairs where a <= b
  idx.pairs <- idx.pairs[idx.pairs$a <= idx.pairs$b, ]
  
  # Remove NA values and calls linkage.disequilibrium.D.prime() for pair i
  do.pair <- function(i) {
    a <- idx.pairs$a[i]
    b <- idx.pairs$b[i]
    # exclude indices that are NA in either marker a or b
    exclude <- union(na.indices[[a]], na.indices[[b]])
    
    if ( length(exclude) == nrow(X))
      return ( 0 )
    
    if ( length(exclude) == 0) {			
      marker1 <- X[, a]
      marker2 <- X[, b]
    } else {
      marker1 <- X[-exclude, a]
      marker2 <- X[-exclude, b]
    }
    return ( linkage.disequilibrium.D.prime(marker1, marker2) )
  }
  
  # use MC.CORES to calculate D' for every pair
  D.prime.values <- mclapply(1:nrow(idx.pairs), do.pair, mc.cores=mc.cores)
  
  # Convert list into a matrix
  D.prime <- matrix(-1, n.geno, n.geno)
  for ( i in 1:nrow(idx.pairs) ) {
    D.prime[idx.pairs$a[i], idx.pairs$b[i]] <- D.prime.values[[i]]
    D.prime[idx.pairs$b[i], idx.pairs$a[i]] <- D.prime.values[[i]]
  }
  
  return ( D.prime ) 
}

# The log likelihood function for pAB for calculating D' 
# adaptived from LD function in library(genetics)
log.L.linkage.disequilibrium <- function(pAB, p.A, p.B, counts) {
  pAb <- p.A - pAB
  paB <- p.B - pAB
  pab <- 1 - p.A - p.B + pAB
  
  return (counts[1] * log(pab) +
            counts[2] * log(paB) +
            counts[3] * log(pAb) +
            counts[4] * log(pAB) +
            counts[5] * log(pAB*pab + pAb*paB))
}

# calculates D' much faster than the LD function in library(genetics)
# code is adapted from library(genetics)
linkage.disequilibrium.D.prime <- function(X1, X2) {
  
  if ( all(X1==X2) ) {
    return ( 1 )
  }
  
  p.A <- sum(X1) / (2*length(X1))
  p.B <- sum(X2) / (2*length(X2))
  
  if ( p.A < 0.5 ) {
    p.A <- 1 - p.A
    X1 <- complement(X1)
  }
  
  if ( p.B < 0.5 ) {
    p.B <- 1 - p.B
    X2 <- complement(X2)
  }
  
  p.a <- 1 - p.A
  p.b <- 1 - p.B
  
  # Make a contingency table for makers X1 and X2
  N <- table(factor(X1, levels=0:2), factor(X2, levels=0:2))
  
  counts <- c()
  counts[1] <- 2*N[1, 1] + N[1, 2] + N[2, 1] # p(ab)
  counts[2] <- 2*N[1, 3] + N[1, 2] + N[2, 3] # p(aB)
  counts[3] <- 2*N[3, 1] + N[2, 1] + N[3, 2] # p(Ab)
  counts[4] <- 2*N[3, 3] + N[3, 2] + N[2, 3] # p(AB)
  counts[5] <- N[2, 2] # p(AB)p(ab) + p(Ab)*p(aB)
  
  Dmin <- max(-p.A*p.B, -p.a*p.b)
  Dmax <- min(p.A*p.b, p.B*p.a);
  
  pmin <- p.A*p.B + Dmin;
  pmax <- p.A*p.B + Dmax;
  
  # Find p(AB) by maximizing log-likelihood function
  solution <- tryCatch(optimize(log.L.linkage.disequilibrium,
                                lower=pmin + .Machine$double.eps,
                                upper=pmax - .Machine$double.eps,
                                maximum=TRUE,
                                p.A=p.A, p.B=p.B, counts=counts),
                       error=function(e) { return ( NA ) })
  
  if ( is.na(solution) ) {
    return ( NA )
  }
  
  pAB <- solution$maximum
  
  D <- pAB - p.A*p.B
  
  if ( D > 0 ) {
    D.prime <- D / Dmax
  } else {
    D.prime <- D / Dmin
  }
  
  return ( D.prime )
}

#' @export
ld.matrix <- function(genotypes, mc.cores) {
  require('parallel')
  stopifnot(max(genotypes, na.rm=TRUE) <= 2,
            min(genotypes, na.rm=TRUE) >= 0)
  
  if ( missing(mc.cores) ) {
    mc.cores <- detectCores()
  }
  cat('Calculating LD for', ncol(genotypes), 'genotypes using', mc.cores, 'cores...\n')
  LDM <- LD.D.prime.matrix(genotypes, mc.cores)
  
  return ( LDM )
}

#' @export
ld.plot <- function(ld.matrix, snp.positions, start.pos, end.pos, ld.colors, drawscale=TRUE, scale.cex=0.6) {
  
  if ( missing(ld.colors) ) {
    ld.colors <- rgb(1, 99:0/99, 99:0/99)
  }
  stopifnot(ncol(ld.matrix)==length(snp.positions),
            nrow(ld.matrix)==length(snp.positions))
  
  if ( is.unsorted(snp.positions) ) warning('SNPs positions are not in order!')
  
  plot(0, type='n', ylim=c(min(snp.positions)-max(snp.positions), 0), 
       xlim=range(snp.positions),
       axes=FALSE, bty='n', xlab='', ylab='', yaxs='i')
  
  # Calculate positions for SNPs
  num.snps <- length(snp.positions)
  x.pos <- snp.positions
  x.pos.mids <- (x.pos[-1] + x.pos[-num.snps])/2
  x.pos.left <- c(x.pos[1], x.pos.mids)
  x.pos.right <- c(x.pos.mids, x.pos[num.snps])
  
  # Loop through each position in the LDM array and draw polygon of appropriate shade
  for ( i in 1:num.snps ) {
    for ( j in 1:i ) {
      fill.color <- ld.colors[round(abs(ld.matrix[i, j]) * (length(ld.colors) - 1)) + 1]
      # Draw a diamond shape 
      polygon(c(x.pos.left[i]+x.pos.left[j],
                x.pos.right[i]+x.pos.left[j],
                x.pos.right[i]+x.pos.right[j],
                x.pos.left[i]+x.pos.right[j])/2,
              -c(x.pos.left[i]-x.pos.left[j],
                 x.pos.right[i]-x.pos.left[j],
                 x.pos.right[i]-x.pos.right[j],
                 x.pos.left[i]-x.pos.right[j]),
              col=fill.color, border=fill.color)
      # NOTE: not drawing polygon borders (i.e. border=NA) causes jagged,
      # pixelated edges. So we draw a border in the same color as the fill
    }
  }
  
  if ( drawscale ) {
    draw.scale(ld.colors, c(0, 1), num.labs=6, position='bottomleft', cex=scale.cex)
  }
}
