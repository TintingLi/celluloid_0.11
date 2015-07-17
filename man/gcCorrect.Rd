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
gcCorrect(rangedata, span = 0.3, mappability = 0.9, samplesize = 50000, bprange = 6, maskmap=T )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rangedata}{
  Object of class IRanges created by HMMcopy::wigsToRangedData
}
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
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Same object rangedata, with an additional column named reads.gc that contains the (scaled) GC-corrected read counts. More prexcisely the column takes the value rangedata$reads/predict( loess(rangedata$reads~rangedata$gc, ...), rangedata$gc)
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Mathieu Lemire
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (rangedata, span = 0.3, mappability = 0.9, samplesize = 50000, 
    bprange = 6) 
{
    reads <- rangedata$reads
    rangedata$gc[rangedata$gc < 0] <- NA
    reads[reads < 0] <- NA
    rangedata$ignore <- F
    rangedata$ignore[is.na(reads) | is.na(rangedata$gc)] <- T
    rangedata$ignore[!is.element(rangedata$space, paste("chr", 
        1:22, sep = ""))] <- T
    cat("Identifying and ignoring outlier read counts...")
    bp <- boxplot(reads, rangedata = bprange, plot = F)
    rangedata$ignore[reads < bp$stats[1, 1] | reads > bp$stats[5, 
        1]] <- T
    rangedata$ignore[rangedata$map < mappability] <- T
    keep <- which(!rangedata$ignore)
    keep <- sample(keep, samplesize)
    cat(" done.\nApplying loess correction...")
    lo <- loess(reads[keep] ~ rangedata$gc[keep], span = span)
    if (FALSE) {
        plot(rangedata$gc[keep], reads[keep], pch = ".")
        or <- order(rangedata$gc[keep])
        points(rangedata$gc[keep][or], predict(lo, rangedata$gc[keep][or]), 
            type = "l", col = "red")
    }
    rangedata$reads.gc <- rangedata$reads/predict(lo, rangedata$gc)
    rangedata$reads.gc[rangedata$reads.gc < 0] <- NA
    rangedata$reads.gc <- rangedata$reads.gc/mean(rangedata$reads.gc, 
        na.rm = T)
    rangedata$copy <- log(rangedata$reads.gc, 2)
    cat(" done.\nFind the result in $reads.gc\n")
    return(rangedata)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
