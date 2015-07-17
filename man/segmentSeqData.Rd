\name{segmentSeqData}
\alias{segmentSeqData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Copy number segmentation
}
\description{
Wrapped for function copynumber::pcf to do segmentation.
}
\usage{
segmentSeqData(rangedata, gamma = 500, kmin = 100, maskmap = 0.8, maskadj = FALSE, skipmeanmap=FALSE )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rangedata}{
Object of class IRange produced by function gcCorrect.
}
  \item{gamma}{
Parameter of copynumber::pcf
}
  \item{kmin}{
Parameter of copynumber::pcf
}
  \item{maskmap}{
Ignore bins with mappability less than this
}
  \item{maskadj}{
Ignore bins adjacent to bins with mappability less than maskmap. NOT USED.
}
\item{skipmeanmap}{Skip call to meanMapInSeg}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Data fram returned by copynumber::pcf, to which is added a column named meanmap, that contains the mean mappability calculated in each segment.
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
Mathieu Lemire}
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
function (rangedata, gamma = 500, kmin = 100, maskmap = 0.8, 
    maskadj = FALSE) 
{
    binlength <- end(rangedata[1, ]) - start(rangedata[1, ]) + 
        1
    tmp <- data.frame(chrom = rangedata$space, pos = start(rangedata), 
        rc = rangedata$reads.gc)
    tmp$chrom <- factor(as.character(tmp$chrom), levels = c(paste("chr", 
        1:22, sep = ""), "chrX", "chrY"))
    tmp$chrom <- as.numeric(tmp$chrom)
    sel <- rangedata$map < maskmap
    if (maskadj) {
        sell <- c(F, sel[-length(sel)])
        selr <- c(sel[-1], F)
        sel <- sel | sell | selr
    }
    tmp <- tmp[!sel, ]
    tmp <- tmp[!is.na(tmp$rc), ]
    require(copynumber)
    tmp.win <- winsorize(tmp)
    tmp.seg <- pcf(data = tmp.win, gamma = gamma, kmin = kmin, 
        digits = 4)
    tmp.seg$end.pos <- tmp.seg$end.pos + binlength - 1
    tmp.seg$chrom <- paste("chr", tmp.seg$chrom, sep = "")
    tmp.seg.meanmap <- meanMapInSeg(tmp.seg, rangedata)
    return(tmp.seg.meanmap)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
