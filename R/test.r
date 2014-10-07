
if ( FALSE ) {
  library('zoom.plot')

  loc <- 96252589
  window <- 1e6
  
  par(mar=c(4, 4, 0, 4))
  
  test.genes <- get.regional.genes(5, loc - window, loc + window)
  test.gene.rows <- arrange.gene.rows(test.genes)
  plotgenes(test.genes, test.gene.rows, highlight.gene='ERAP2')
  
  test.gene.rows <- arrange.gene.rows(test.genes, label.size=0.75)
  plotgenes(test.genes, test.gene.rows, highlight.gene='ERAP2', label.size=0.75)
  
  test.enc <- get.encode.region(5, loc - window, loc + window)
  test.filtered.enc <- test.enc[test.enc$score > 900, ]
  test.enc.rows <- arrange.encode.rows(test.filtered.enc)
  plotencode(test.filtered.enc, test.enc.rows)
  
  
  stacked <- fig.parts(c(NA, NA, NA))
  
  par(mar=c(0, 4, 0, 4))
  stacked(1, new=TRUE)
  plotgenes(test.genes, test.gene.rows, loc - window, loc + window, highlight.gene='ERAP2')
  stacked(2)
  plotencode(test.filtered.enc, test.enc.rows, loc - window, loc + window)
  stacked(3)
  draw.chrom.axis(par('usr')[1], par('usr')[2], 'chr5')
  
  
  
  
  test.genes <- get.regional.genes(12, 82452, 1132743)
  test.gene.rows <- arrange.gene.rows(test.genes)
  plotgenes(test.genes, test.gene.rows)
  
}
