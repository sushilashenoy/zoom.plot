% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zoom_plot.r
\name{zoom.plot}
\alias{zoom.plot}
\title{Zoom in manhattan plot with gene models and LD plotted below.}
\usage{
zoom.plot(pvals, gene.id, all.snps, snp.id = NULL, genotypes,
  window = 1e+06, min.snps = 10, color = "maf", ld.opts = NULL,
  ld.colors = NULL, img.height = 700, img.width = 1200, img.prefix = "",
  save.image = TRUE)
}
\arguments{
\item{pvals}{a vector where each pval corresponds to a single SNP, ordered
by chromosome and position (same order as \code{all.snps})}

\item{gene.id}{either a gene or phenotype name. If this matches a cis gene,
that gene will be drawn in a different color. The name is also used for the
plot title.}

\item{all.snps}{a data.frame containing one row for each snp (same order as
pvals). The data.frame should have three named columns (in any order): id,
chrom, and position. The snp order should match \code{pvals}.}

\item{snp.id}{The snp that the plot will center on. If not specified we find the
most significant one}

\item{genotypes}{a matrix of genotypes with 0/1/2 coding for all snps, with snps
in columns}

\item{window}{window size, in bp, of chromosomal region to show in zoom plot.
If there are more than a few hundred SNPs it can be very slow to calculate LD.}

\item{min.snps}{This argument is included for compatibility with previous
versions but DOES NOT DO ANYTHING.}

\item{color}{This can be \code{'maf'} if the manhattan plot should be colored by
minor allele frequency or \code{'het'} to color by percent heterozygotes.
Any other value will give black points.}

\item{ld.opts}{list of optional arguments for ld.matrix}

\item{ld.colors}{color scale for LD plot (default white->red)}

\item{img.height}{height of saved image}

\item{img.width}{width of saved image}

\item{img.prefix}{Any additional image prefix desired, in addition to the
gene.id which is always included.}

\item{save.image}{(default TRUE) set to false to plot with current device. Can
be used to plot within R or to save with a different image format.

This function was previously implemented in a standalone script
(\code{'zoom_ld_biomart.r'}) which was largely based on code
from G. Hoffman.) By default the most recent human gene ensembl database is used for
gene models. To use a different database please use the \code{load.gene.mart()}
function.}
}
\description{
The zoom.plot function creates a png file with a zoom in manhattan plot
and LD plot.
}

