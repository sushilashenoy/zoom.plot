
# listMarts(host='www.ensembl.org')
# all.datasets <- listDatasets(useMart('ENSEMBL_MART_ENSEMBL', host='www.ensembl.org'))
# all.datasets <- listDatasets(useMart('ENSEMBL_MART_ENSEMBL', host='feb2014.archive.ensembl.org'))
# print(all.datasets[all.datasets$dataset=='hsapiens_gene_ensembl', ])

# test.mart <- useMart("ENSEMBL_MART_ENSEMBL", host='www.ensembl.org', dataset='hsapiens_gene_ensembl')
# test.mart <- useMart("ENSEMBL_MART_ENSEMBL", host='feb2014.archive.ensembl.org', dataset='hsapiens_gene_ensembl')
# head(listAttributes(test.mart))

BIOMART_OPTIONS <- new.env(parent=emptyenv())

#' @export
load.gene.mart <- function(version='hg38',
                           biomart.host,
                           biomart.fields,
                           biomart.filters,
                           biomart.mart,
                           biomart.dataset) {
  require('biomaRt')
  
  BIOMART_OPTIONS$host <- 'www.ensembl.org'
  BIOMART_OPTIONS$fields <- c('ensembl_transcript_id',
                              'external_gene_name',
                              'transcript_start',
                              'transcript_end',
                              'exon_chrom_start',
                              'exon_chrom_end')
  BIOMART_OPTIONS$filters <- c("chromosomal_region", 'with_ccds')
  BIOMART_OPTIONS$mart <- 'ENSEMBL_MART_ENSEMBL'
  BIOMART_OPTIONS$dataset <- 'hsapiens_gene_ensembl'
  
  if (!missing(biomart.host) ) {
    BIOMART_OPTIONS$host <- biomart.host
  }
  if ( !missing(biomart.fields) ) {
    BIOMART_OPTIONS$fields <- biomart.fields
  }
  if ( !missing(biomart.mart) ) {
    BIOMART_OPTIONS$mart <- biomart.mart
  }
  if ( !missing(biomart.dataset) ) {
    BIOMART_OPTIONS$dataset <- biomart.dataset
  }
  
  if ( version=='hg37' ) {
    BIOMART_OPTIONS$host <- 'feb2014.archive.ensembl.org'
    BIOMART_OPTIONS$fields <- c('ensembl_transcript_id',
                         'external_gene_id',
                         'transcript_start',
                         'transcript_end',
                         'exon_chrom_start',
                         'exon_chrom_end')
  } 
  
  BIOMART_OPTIONS$ens.mart <- useMart(BIOMART_OPTIONS$mart,
                                      host=BIOMART_OPTIONS$host,
                                      dataset=BIOMART_OPTIONS$dataset)
}

#' @export
get.regional.genes <- function(chrom, start.pos, end.pos) {
  if ( !exists('ens.mart', envir=BIOMART_OPTIONS) ) {
    load.gene.mart()
  }
  
  # Use biomaRt to get a list of genes in our region
  if ( chrom == 23 ) {
    chrom <- 'X'
  } else if ( chrom == 26 ) {
    chrom <- 'M' 
  }
  bm.range <- paste(chrom, start.pos, end.pos, sep=':')
  
  cat('Finding genes in region:', bm.range, '...\n')
  bm.exons <- getBM(attributes=BIOMART_OPTIONS$fields,
                    filters=BIOMART_OPTIONS$filters,
                    values=list(bm.range, TRUE), mart=BIOMART_OPTIONS$ens.mart)
  
  # Replace long names with shorter ones
  colnames(bm.exons) <- c('tx.id', 'name', 'start', 'end',
                          'exon.start', 'exon.end')
  
  if ( nrow(bm.exons) == 0 ) {
    cat('No genes found!\n')
    return ( NULL )
  }
  
  # Choose one transcript per gene (just use first??)
  use.tx <- aggregate(tx.id ~ name, bm.exons, head, 1)
  bm.exons <- bm.exons[bm.exons$tx.id %in% use.tx$tx.id, ]
  
  # Make a data.frame with info for each gene in region
  regional.genes <- unique(bm.exons[, c('tx.id', 'name', 'start', 'end')])
  
  cat('Found', nrow(regional.genes), 'genes in region.\n')
  
  # use pmax/pmin to determine left/right of each gene
  regional.genes$left <- pmin(regional.genes$start, regional.genes$end) 
  regional.genes$right <- pmax(regional.genes$start, regional.genes$end) 
  
  # condense exon information
  regional.genes$exon.starts <- sapply(regional.genes$tx.id, function (id) {
    paste(bm.exons$exon.start[bm.exons$tx.id==id], collapse=',')
  })
  regional.genes$exon.ends <- sapply(regional.genes$tx.id, function (id) {
    paste(bm.exons$exon.end[bm.exons$tx.id==id], collapse=',')
  })
  
  attr(regional.genes, 'start.pos') <- start.pos
  attr(regional.genes, 'end.pos') <- end.pos
  return ( regional.genes )
}

#' @export
arrange.gene.rows <- function(regional.genes, start.pos, end.pos, label.size=0.75) {
  if ( is.null(regional.genes) ) return ( NULL )
  
  if ( missing(start.pos) ) start.pos <- attr(regional.genes, 'start.pos')
  if ( missing(end.pos) ) end.pos <- attr(regional.genes, 'end.pos')
  
  char.width <- label.size * (end.pos-start.pos)/(par('pin')[1]/par('cin')[1])
  min.label.width <- (end.pos-start.pos)/(par('pin')[1]/par('cin')[1])
  row.occupied <- rep(-Inf, nrow(regional.genes))
  gene.rows <- rep(0, nrow(regional.genes))
  
  # Loop over genes in order based on left position
  for ( i in order(regional.genes$left) ) {
    # Starting from row 1, search for a row where this gene will fit
    row <- 1
    label.width <- max(min.label.width, char.width * (nchar(regional.genes$name[i])+1))
    while ( row.occupied[row] + label.width > regional.genes$left[i] ) {
      row <- row + 1
    }
    # This row is now occupied up to the end of this gene
    row.occupied[row] <- regional.genes$right[i]
    gene.rows[i] <- row
  }
  
  return ( gene.rows )
}



#' @export
plotgenes <- function(regional.genes, gene.rows, start.pos, end.pos, highlight.gene, label.size=0.75, max.rows=5) {
  if ( is.null(regional.genes) ) return ( invisible(NULL) )
  
  if ( missing(start.pos) ) start.pos <- attr(regional.genes, 'start.pos')
  else {
    regional.genes <- regional.genes[regional.genes$end > start.pos, ]
  }
  if ( missing(end.pos) ) end.pos <- attr(regional.genes, 'end.pos')
  else {
    regional.genes <- regional.genes[regional.genes$start < end.pos, ]
  }
  
  if ( missing(gene.rows) ) gene.rows <- arrange.gene.rows(regional.genes, start.pos, end.pos, label.size)
  
  if ( max(gene.rows) > max.rows ) {
    hide.labels <- TRUE
    gene.rows <- arrange.gene.rows(regional.genes, start.pos, end.pos, label.size=0)
  } else {
    hide.labels <- FALSE
  }
  
  plot(0, type='n', ylim=0.5+c(0, max(gene.rows)), xlim=c(start.pos, end.pos),
       axes=FALSE, bty='n', xlab='', ylab='', yaxs='i')
  
  gene.color <- rep('blue', nrow(regional.genes))
  if ( !missing(highlight.gene) ) {
    gene.color[regional.genes$name==highlight.gene] <- 'red'
  }
  
  sapply(1:nrow(regional.genes), function (i) {
    
    lines(c(regional.genes$start[i], regional.genes$end[i]),
          c(gene.rows[i], gene.rows[i]),
          col=gene.color[i],
          lwd=4, lend=1)
    
    exon.starts <- unlist(strsplit(regional.genes$exon.starts[i], ',', fixed=TRUE))
    exon.ends <- unlist(strsplit(regional.genes$exon.ends[i], ',', fixed=TRUE))

    sapply(1:length(exon.starts), function (e) {
      lines(c(exon.starts[e], exon.ends[e]), c(gene.rows[i], gene.rows[i]),
            lwd=8, lend=1, col=gene.color[i])
      })
    if ( ! hide.labels ) {
      text(max(par('usr')[1], regional.genes$left[i]), gene.rows[i],
           regional.genes$name[i],
           pos=2, offset=0.25, cex=label.size, xpd=NA) # family='mono', 
    }
  })
  
  invisible(NULL)
}
