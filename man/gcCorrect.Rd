\name{gcCorrect}
\alias{gcCorrect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Correction for GC content
}
\description{
This function performs GC correction of read counts in bins, by means of loess. 
}
\usage{
gcCorrect(rangedata, sampletype="normal",span = 0.3, mappability = 0.9, 
            samplesize = 50000, bprange = 6, maskmap=T ,... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rangedata}{
  Object of class IRanges created by HMMcopy::wigsToRangedData
}
\item{sampletype}{One of "tumor" or "normal"}
  \item{span}{
span parameter of loess().
}
  \item{mappability}{
Bins with mappability lower than this are ignored for GC correction
}
  \item{samplesize}{
That many random bins are selected to derive GC correction
}
  \item{bprange}{
  Bins with read counts that are outliers are ignored in GC correction, where an outlier is defined as in boxplot(..., range=bprange )
}
\item{maskmap}{Sets reads.gc and map columns to NA if bin's mappability is less than the value of mappability argument}
\item{...}{Arguments to be passed to copynumber::pcf through segmentSeqData}
}
\details{
If sampletype is "tumor", then a first segmentation round is done to adjust the read counts for the segments' height.}
\value{
Returns the rangedata, with an additional column named reads.gc that contains the (scaled) GC-corrected read counts. More precisely the column takes the value:

rangedata$reads/predict( loess(rangedata$reads~rangedata$gc, ...), rangedata$gc)
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
Mathieu Lemire
}

\examples{
\dontrun{
# normal data
n<- wigsToRangedData( normalWigFile, gcWigFile, mapWigFile )
nc<-gcCorrect(n) 
# tumour data
t <- wigsToRangedData( tumourWigFile, gcWigFile, mapWigFile )
tc<-gcCorrect( t, sampletype="tumor" )
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
