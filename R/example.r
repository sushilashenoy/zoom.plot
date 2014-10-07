
if ( interactive() ) {
  
  library('zoom.plot')
  
  # Define colors for manhattan plot
  maf.colors <- c(rgb(20:0/20, 0, 0), rgb(0, 1:20/20, 0))
  
  # Set a location
  chrom <- 5
  loc <- 96252589
  window <- 5e5
  
  start.pos <- loc - window
  end.pos <- loc + window
  
  # Make up SNP data for a manhattan plot
  snp.coords <- floor(runif(100, start.pos, end.pos))
  snp.logp <- -log10(runif(100))
  snp.mafs <- runif(100, 0.1, 0.49)
  
  # Get colors for SNPs according to maf
  snp.colors <- assign.scale.colors(snp.mafs, maf.colors, c(0, 0.5))

  # Get genes in location
  test.genes <- get.regional.genes(chrom, start.pos, end.pos)
  
  # Make a stacked plot with middle portion size of 15% (for axis)
  fig <- fig.parts(c(NA, 0.15, 0.25))
  
  fig(2, new=TRUE) # Start a new figure, region 2 (middle)
  par(mar=c(0, 4, 0, 4))
  draw.chrom.axis(start.pos, end.pos, paste('chr', chrom))
  use.xlim <- par('usr')[1:2]
  
  fig(3) # Continue figure, region 3 (bottom)
  plotgenes(test.genes, highlight.gene='ERAP2')
  
  fig(1) # Continue figure, region 1 (top)
  par(mar=c(0, 4, 3, 4))
  plot(snp.coords, snp.logp, pch=20, col=snp.colors,
       xaxt='n', bty='n', xlim=use.xlim, xaxs='i',
       ylab=expression(-log[10](p)), xlab='')
  draw.scale(maf.colors, c(0, 0.5), pos='topright', y.offset=0.5, xpd=NA)

  
  
  # Instead of fig.parts you could just use mfcol for equally sized regions:
  par(mfcol=c(3, 1))
  par(mar=c(0, 4, 0, 4), oma=c(0, 0, 4, 0))
  plot(snp.coords, snp.logp, pch=20, col=snp.colors,
       xaxt='n', bty='n', xlim=c(start.pos, end.pos), xaxs='i',
       ylab=expression(-log[10](p)), xlab='')
  draw.scale(maf.colors, c(0, 0.5))
  draw.chrom.axis(start.pos, end.pos, paste('chr', chrom))
  use.xlim <- par('usr')[1:2]
  plotgenes(test.genes, highlight.gene='ERAP2')
  
  
}
