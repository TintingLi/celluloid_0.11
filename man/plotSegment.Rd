\name{plotSegment}
\alias{plotSegment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot segments 
}
\description{
Plot segments
}
\usage{
plotSegment(tumourrangedata, segments, ar = NULL, n.rc.seg = NULL, 
            columns = NULL, maskmap = NULL, file = "Rplot\%03d", 
            device="pdf", main = NULL, line=4.5, chr = NULL, perpage = 4, 
            layoutmat = NULL, width = 8.5, height = 11 , ylim = c(-0.5, 8), 
            tlwd = 5, tlty = 1, tcol = NULL, nlwd = 3, nlty = 3, 
            ncol = gray(0.5), alpha=50, cex.axis=1, cex.lab=1,cex.main=1, 
            annotation=NULL, cex.annotation=1,marCN=c(0.5,5,6,1), 
            marAR=c(3.5,5,0.5,1) , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tumourrangedata}{
Object of class IRanges containing tumour data, obtained from a call to gcCorrect
}
  \item{segments}{
segments data.frame obtained from a call to annotateSegments
}
  \item{ar}{
The ar data frame}
  \item{n.rc.seg}{
Segments in normal tissues. If provided, these segments will also be plotted.
}
  \item{columns}{
Vector of size 2 indicating the names of the columns to plot from tumourrangedata and segments. If NULL, defaults to c("icopy", "imean" ).
}
  \item{maskmap}{
Only plot a segment if its mean mappability is above maskmap. Plot all segments if set to NULL.
}
  \item{file}{
File prefix; suffix is appended based on device. Set to NULL to display on screen.
}
\item{device}{One of "pdf" (default),"bmp","jpeg", "png","tiff".}
  \item{main}{
main title for the graph (as in title(main=)). Defaults to the chromosome number.
}
\item{line}{Placement of main: this many lines outwards from the plot edge. See ?title.}
  \item{chr}{
Which chromosomes to plot. If NULL defaults to all chromosomes.
}
  \item{perpage}{
How many chromosomes to plot per page
}
  \item{layoutmat}{
See details.}
  \item{width}{
Argument passed to the graphic device to define the width of a page. The default value was set for a pdf device, in inch units (the default of pdf). 
For other devices, the user needs to define width and height in pixels, or else needs to define values for units and res (by adding arguments such as units="in", res=300; see e.g. ?png).  
}
  \item{height}{
Argument passed to the graphic device to define the height of a page. The default value was set for a pdf device, in inch units (the default of pdf). 
For other devices, the user needs to define width and height in pixels, or else needs to define values for units and res (by adding arguments such as units="in", res=300; see e.g. ?png).  
}
  \item{ylim}{
ylim in copy number graphs, same for all chromosomes.
}
  \item{tlwd}{
Line width of tumour segments
}
  \item{tlty}{
Line type of tumour segments
}
  \item{tcol}{
Color of tumour segments. See details. 
}
  \item{nlwd}{
Line width of segments in normal.
}
  \item{nlty}{
Line type of segments in normal.
}
  \item{ncol}{
Color of segments in normal. 
}
\item{alpha}{alpha argument passed to rgb to define transparency of the colors of the scaled read counts of the 1kb bins. Takes value between 0 and 255. Alpha value of 255 is opaque. Alpha value of 0 prevents points from being plotted. See Details. }
\item{cex.axis}{Magnification to be used for text annotation of the axis values. Passed to axis().}
\item{cex.lab}{Magnification to be used for text annotation of the axis labels. Passed to title().}
\item{cex.main}{Magnification to be used for text annotation of title (the main. Passed to title().}
\item{annotation}{Data frame containing: chr position label.  Labels are added to axis(1) of the segment plots
}
\item{cex.annotation}{ cex.axis for the annotation }
\item{marCN}{Margins for the copy-number component of the graph}
\item{marAR}{Margins for the allelic-ratio component of the graph}
\item{...}{Additional arguments passed to device. See ?pdf, ?png, etc.}
}
\details{
The default layout prints 4 chromosomes per page. For each chromosome, the graph displays 
two items: 1) the copy number segments on top and 2) allelic ratios on bottom.  
Underneath the copy number segments, the graph displays each individual 1kb bin.

The two items are plotted according to a 3:1 row height 
ratio by default. A full page of 4 uses the layout 

layout( matrix( c( 1,1,1,2,3,3,3,4,5,5,5,6,7,7,7,8 ), ncol=1 ) ) 

See ?layout for more details. 

tcol can take either a single color value or a vector of 10 colors. Colors are assigned based on segment height.
Segments with copy number in (-Inf,0.5] will be colored tcol[1], (0.5,1.5] will be colored tcol[2], ...,
up to segments with copy number in (8.5,Inf) that will be colored tcol[10].

When set to NULL, defaults to:
  rb<-rainbow(24)
  tcol<-c( rb[1], rb[4], "black", rb[6], rb[9], rb[13], rb[15], rb[17], rb[19] ,rb[20] )

The 1kb bins are plotted with semi-transparency by defaults (see ?rgb and the alpha argument). Not all
devices can plot transparencies, see help page for the device you want if you get a warning.  An alpha 
value of 255 is opaque.

}

\author{
Mathieu Lemire}

\examples{
\dontrun{
prepCN(12,1,NULL)
ePP<-ePeakPos( S=0.624, t=c(0.073, 0.923), cn=cn  )
tcs<- scaleReadCounts( tc , ePP )
segments<-scaleSegments(t.ar.seg ,  ePP )
segments<-annotateSegments(segments, ePP)
plotSegment( tcs,segments, ar , file="segments_page%1d",device="png",
             width=960,height=1320, cex.axis=2, cex.main=2, cex.lab=2, 
             type="cairo" , tlwd=8 ) 

# for sex chromosomes
segmentsXY<-scaleSegments(t.seg ,  ePP )
plotSegment( tcs,segmentsXY, ar=NULL , file="segments_pageXY",device="png",
             width=960,height=1320, cex.axis=2, cex.main=2, cex.lab=2, 
             type="cairo" , chr=c("chrX","chrY") , tlwd=8) 


}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
