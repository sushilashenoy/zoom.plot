
#' @export
load.tped <- function(prefix) {
  tped.file <- paste0(prefix, '.tped')
  tfam.file <- paste0(prefix, '.tfam')
  stopifnot(file.exists(tped.file), file.exists(tfam.file))
  
  geno.samples <- read.table(tfam.file)
  n.samples <- nrow(geno.samples)
  
  geno.data <- scan(tped.file, character())
  n.snps <- length(geno.data)/(2*n.samples+4)
  
  dim(geno.data) <- c(2*n.samples+4, n.snps)
  
  geno.raw <- geno.data[-(1:4), ]
  geno.info <- geno.data[1:4, ]
  rm('geno.data')
  gc(FALSE)
  
  colnames(geno.info) <- c('chromosome', 'id', 'distance', 'position')
  
  
  alleles <- sort(unique(as.vector(geno.raw)))
  # Put '0' (missing) at the end
  alleles <- c(alleles[alleles !='0'], '0')
  n.alleles <- length(alleles)-1
  rev.order <- c(n.alleles:1, n.alleles+1)

  # This is the most time consuming step
  has.what <- t(apply(geno.raw, 2, function (x) alleles %in% x))
  
  # A allele is first alphabetically
  # B allele is last alphabetically
  a.allele <- alleles[apply(has.what, 1, which.max)]
  b.allele <- (alleles[rev.order])[apply(has.what[, rev.order], 1, which.max)]
  b.allele[a.allele==b.allele] <- NA
  
  geno.both <- apply(geno.raw, 1, function (x) x != a.allele)
  geno.both[geno.raw=='0'] <- NA
  
  # Add even and odd columns to get diploid genotype code (# of b alleles)
  geno <- geno.both[, 2*(1:n.samples)-1] + geno.both[, 2*(1:n.samples)]
  
  geno.info$a <- a.allele
  geno.info$b <- b.allele
  
  colnames(geno) <- make.names(geno.samples$V2)
  rownames(geno) <- make.names(geno.info$id)

  return (list(geno=geno, geno.info=geno.info))
}