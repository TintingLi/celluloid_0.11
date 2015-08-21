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
\dontrun{
prepCN(12)
ePP <- ePeakPos( S= 0.6260699, t=c(0.02222, 0.97778), cn=cn  )
segments<-scaleSegments(t.ar.seg ,  ePP )
segments<-annotateSegments(segments, ePP)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
