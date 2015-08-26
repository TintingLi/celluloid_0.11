\name{estimateLOHcurve}
\alias{estimateLOHcurve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to estimate LOH curve
%%  ~~function to do ... ~~
}
\description{Function to estimate LOH curve
}
\usage{
estimateLOHcurve(segments, cntr = NULL, plot = F, manual = F, minsize = 1e+07, minmap = 0.9, maxmean = 1.25)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segments}{
Segment data frame, such as t.ar.seg. Expecting to find columns "mean" and "p".
}
  \item{cntr}{
Object returned by showTumourProfile. Used if plot = T or manual = T.
}
  \item{plot}{
Display the tumour contour plot and the LOH curve.
}
  \item{manual}{
Manually select segments to enter curve estimation.
}
  \item{minsize}{
Minimum size of segments to be used to estimate the LOH curve.
}
  \item{minmap}{
Minumum mean mappability of segments to be used to estimate the LOH curve.
}
  \item{maxmean}{
Maximum mean value (x-axis) of segments to be used to estimate the LOH curve.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
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
\dontrun{
Sn<- estimateLOHcurve(t.ar.seg)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
