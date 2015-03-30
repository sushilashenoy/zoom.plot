
#' Zoom in manhattan plot with gene models and LD plotted below.
#' 
#' The zoom.plot function creates a png file with a zoom in manhattan plot
#' and LD plot.
#' 
#' @param pvals a vector where each pval corresponds to a single SNP, ordered
#' by chromosome and position (same order as \code{all.snps})
#' @param gene.id either a gene or phenotype name. If this matches a cis gene,
#' that gene will be drawn in a different color. The name is also used for the
#' plot title.
#' @param all.snps a data.frame containing one row for each snp (same order as
#' pvals). The data.frame should have three named columns (in any order): id,
#' chrom, and position. The snp order should match \code{pvals}.
#' @param snp.id The snp that the plot will center on. If not specified we find the
#' most significant one
#' @param genotypes a matrix of genotypes with 0/1/2 coding for all snps, with snps
#' in columns
#' @param window window size, in bp, of chromosomal region to show in zoom plot.
#' If there are more than a few hundred SNPs it can be very slow to calculate LD.
#' @param min.snps This argument is included for compatibility with previous
#' versions but DOES NOT DO ANYTHING.
#' @param color This can be \code{'maf'} if the manhattan plot should be colored by
#' minor allele frequency or \code{'het'} to color by percent heterozygotes.
#' Any other value will give black points.
#' @param img.height height of saved image
#' @param img.width width of saved image
#' @param img.prefix Any additional image prefix desired, in addition to the
#' gene.id which is always included.
#' @param save.image (default TRUE) set to false to plot with current device. Can
#' be used to plot within R or to save with a different image format.
#' 
#' This function was previously implemented in a standalone script
#' (\code{'zoom_ld_biomart.r'}) which was largely based on code
#' from G. Hoffman.) By default the most recent human gene ensembl database is used for
#' gene models. To use a different database please use the \code{load.gene.mart()}
#' function.
#' 
#' @export
zoom.plot <- function(pvals, gene.id, all.snps, snp.id=NULL, genotypes,
                      window=1e6, min.snps=10, color='maf',
                      img.height=700, img.width=1200, img.prefix="", save.image=TRUE) {
  # Colors for LD plot - ramp from white to red
  ld.colors <- c(rgb(1, seq(1, 0.7, length.out=120), seq(1, 0.7, length.out=120)),
                 rgb(1, seq(0.7, 0.4, length.out=41)[-1], seq(0.7, 0.4, length.out=41)[-1]),
                 rgb(1, seq(0.4, 0, length.out=6)[-1], seq(0.4, 0, length.out=6)[-1]))
  # maf.colors <- c(rgb(65:0/65, 0, 0:65/65), rgb(0, 0, 32:0/32))
  maf.colors <- c(rgb(20:0/20, 0, 0), rgb(0, 1:20/20, 0))
  
  het.colors <- c(rgb(0, 0:20/20, 0), rgb(1:20/20, 20:1/20, 0), rgb(1, 0, 1:20/20))

  # Some input checking
  if ( length(pvals) != nrow(all.snps) ) {
    stop('The length of the pvals vector does not match the rows of all.snps data.frame')
  }  else if ( ! all(c('id', 'chrom', 'position') %in% colnames(all.snps)) ) {
    stop('The all.snps data.frame needs to have named columns: id, chrom, and position')
  }
  
  # If snp.id is not provided, just pick the most significant SNP
  if ( is.null(snp.id) ) {
    snp.id <- all.snps$id[which.min(pvals)]
    cat('Centering around most significant SNP:', snp.id, '\n')
  }
  

  snp.index <- which(all.snps$id==snp.id)
  if ( snp.index == 0 ) {
    stop('Could not find specified SNP!')
  }
  
  chromosome <- all.snps$chrom[snp.index]
  
  position <- all.snps$position[snp.index]
  
  include.snps <- NULL

  include.snps <- intersect(which(all.snps$chrom == chromosome),
                            intersect(which(all.snps$position > position - window),
                                      which(all.snps$position < position + window)))

  
  cat('Selected ', length(include.snps), ' SNPs in window (', window/1e6, 'Mb)\n', sep='')
  
  # The positions of the SNPs determine the scale for the x-axis
  snp.positions <- all.snps$position[include.snps]
  
  snp.range <- range(snp.positions)
  x.positions <- snp.positions/1e6
  x.range <- range(x.positions)
  
  # Maximum for y-axis depends on most significant p-value in window
  max.y <- ceiling(max(-log10(pvals[include.snps]), na.rm=TRUE))
  
  # These are the SNPs we will look up to calculate MAF and later to make the LD matrix
  snp.ids <- all.snps$id[include.snps]
  
  if ( !all(snp.ids %in% colnames(genotypes)) ) {
    cat('Could not find genotypes!\n')
    stop('IDs in all.snps data frame must match columns of genotype matrix!')
  }
  
  # Get genotypes and calculate maf
  window.genotypes <- genotypes[, snp.ids] 
  if ( color == 'maf' ) {
    allele.freq <- apply(window.genotypes, 2, mean, na.rm=TRUE)/2
    maf <- pmin(allele.freq, 1-allele.freq)
    snp.colors <- maf.colors[1+floor(maf / 0.5 * (length(maf.colors)-1))]
    
  } else if ( color == 'het' ) {
    pct.het <- apply(window.genotypes, 2, function (x) mean(x==1))
    snp.colors <- het.colors[1+floor(pct.het * (length(het.colors)-1))]
  } else {
    snp.colors <- 'black'
  }
  
  if ( is.character(chromosome) ) {
    char.chroms <- paste('chr', c(1:22, 'X'), sep='')
    if ( chromosome %in% char.chroms ) {
      chromosome <- match(chromosome, char.chroms)
    }
  }
  
  zoom.genes <- get.regional.genes(chromosome, snp.range[1], snp.range[2])
  LDM <- ld.matrix(window.genotypes)
  
  
  if ( save.image ) {
    # Create the image file, name after the gene of interest
    if ( img.prefix == '' ) {
      img.filename <- sprintf('%s_zoom.png', gene.id)
    } else {
      img.filename <- sprintf('%s_%s_zoom.png', img.prefix, gene.id)
    }
    png(filename=img.filename, height=img.height, width=img.width, bg="white", pointsize=24)
  }
  
  layout(matrix(1:4, nrow=4), heights=c(4, 1, 1, 4))
  par(mar=c(0.5, 4, 2, 0.5))
  
  plot(1, xlab='', ylab=expression(-log[10](p)),
       ylim=c(0, max.y), xlim=x.range, bty='n', type='n',
       main = sprintf('%s hit region, chromosome %s', gene.id, chromosome),
       xaxt='n', las=1, xpd=NA, mgp=c(2.5, 1, 0))
  
  points(x.positions, -log10(pvals[include.snps]), pch=20, col=snp.colors)
  
  if ( color == 'maf' ) {
    draw.scale(maf.colors, c(0, 0.5), num.labs=6, position='topright', cex=0.6)
  } else if ( color == 'het' ) {
    draw.scale(het.colors, c(0, 1), num.labs=6, position='topright', cex=0.6)
  }
  
  par(mar=c(0.5, 4, 0, 0.5))
  plotgenes(zoom.genes, highlight.gene=gene.id, label.size=0.7)
  
  par(mar=c(0, 4, 0, 0.5))
  draw.chrom.axis(snp.range[1], snp.range[2], paste('chr', chromosome))
  
  par(mar=c(0.5, 4, 1.5, 0.5))
  ld.plot(LDM, snp.positions, snp.range[1], snp.range[2])
  
  if ( save.image ) {
    dev.off()
  }
}