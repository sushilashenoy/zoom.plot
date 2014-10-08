
if ( interactive() ) {
  
  library('zoom.plot')
  # To use older chromosome coordinates!! Only have to do this once per package load.
  load.gene.mart('hg37')
  
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
  
  # Make up a fake LD matrix based on a subset of SNPs
  ld.mat <- matrix(runif(20*20), 20)
  ld.coords <- sort(sample(snp.coords, 20))
  
  # Get colors for SNPs according to maf
  snp.colors <- assign.scale.colors(snp.mafs, maf.colors, c(0, 0.5))

  # Get genes in region
  test.genes <- get.regional.genes(chrom, start.pos, end.pos)
  
  
  # To remake plot with visual adjustments you only need to rerun code below!
  # ( You do NOT have to look up genes or calculate LD every time! )
  
  # Set margins and layout for 4 plots
  par(mar=c(0, 4, 0, 4), oma=c(2, 1, 2, 0))
  layout(matrix(1:4, 4), height=c(30, 10, 10, 40))
#   layout.show(4)

  # Manhattan plot
  plot(snp.coords, snp.logp, pch=20, col=snp.colors,
       xaxt='n', bty='n', xlim=c(start.pos, end.pos), xaxs='i',
       ylab=expression(-log[10](p)), xlab='')
  draw.scale(maf.colors, c(0, 0.5), y.shift=-1, xpd=NA)

  # Genes plot
  plotgenes(test.genes, highlight.gene='ERAP2', label.size=1)

  # Chromosome axis plot
  draw.chrom.axis(start.pos, end.pos, paste('chr', chrom))
  
  # LD plot
  ld.plot(ld.mat, ld.coords, start.pos, end.pos)
  

}
