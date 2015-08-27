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

\author{
Mathieu Lemire
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
