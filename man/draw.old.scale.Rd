% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/decorations.r
\name{draw.old.scale}
\alias{draw.old.scale}
\title{Draw a color scale bar (old version)}
\usage{
draw.old.scale(scale.colors, scale.range, num.labs = 6,
  position = "topleft", size = 3, width.to.height = 20, x.offset,
  y.offset, x.shift, y.shift, ...)
}
\arguments{
\item{scale.colors}{a vector of colors}

\item{scale.range}{range for scale (vector with 2 numeric elements)}

\item{num.labs}{(optional) number of labels to draw}

\item{position}{(optional) either 'topleft' 'topright' 'bottomleft' or
'bottomright'}

\item{size}{(optional) approximate length in inches}

\item{width.to.height}{(optional) ratio of scale bar width to height}

\item{x.offset}{(optional) adjust x position (uses \code{par('usr')} scale).}

\item{y.offset}{(optional) adjust y position (uses \code{par('usr')} scale).}

\item{x.shift}{(optional) adjust x position in units of scale width.}

\item{y.shift}{(optional) adjust y position in units of scale height.}

\item{...}{(optional) addition options to pass to plotting commands for
drawing labels, lines and shapes}
}
\description{
\code{draw.old.scale} adds a color scale to the corner of the current plot
with specified colors and range
}
\details{
Currently only a horizontal scale bar is supported. To adjust position in
plot coordinates use x.offset and y.offset. To adjust relative position use
x.shift and y.shift.

Keep in mind that if you shift the scale bar away from the plotting area
(i.e. into the margins), you may need to specify the additional parameter
xpd=NA which will be passed to the plotting commands and allows drawing
in the margins.
}
\seealso{
\code{\link{draw.chrom.axis}} for drawing x-axis
}

