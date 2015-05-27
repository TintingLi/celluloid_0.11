\name{scaleSegments}
\alias{scaleSegments}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scale segment data.frame}
\description{
Scale the mean gc- and mappability-corrected read counts in each segment to absolute copy number.
}
\usage{
scaleSegments(seg, epp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seg}{
segment data.frame obtained from a call to segmentSeqData or intersectSegments
}
  \item{epp}{
Expected peak positions. Data frame returned by ePeakPos. 
}
}
\details{
A column seg$mean is expected. The data.frame is returned with a column seg$imean.}
\value{
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
%%  ~~who you are~~
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
function (seg, epp) 
{
    seg$imean <- seg$mean
    mat <- epp[, seq(3, ncol(epp) - 2, 2)]
    pat <- epp[, seq(4, ncol(epp) - 2, 2)]
    matpat <- data.frame(mat + pat)
    sel <- apply(matpat == 0, 1, all)
    zero <- epp[sel, "x"][1]
    sel <- apply(matpat == 1, 1, all)
    one <- epp[sel, "x"][1]
    sel <- apply(matpat == 2, 1, all)
    two <- epp[sel, "x"][1]
    cntransform <- function(x) {
        (x - zero)/(two - zero)
    }
    seg$imean <- 2 * cntransform(seg$mean)
    return(seg)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
