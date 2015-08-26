\name{intersectSegments}
\alias{intersectSegments}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Intersect two segment data.frames
}
\description{
Intersect two segment data.frames}
\usage{
intersectSegments(seg1, seg2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seg1}{
  A segment data.frame obtained from copynumber::pcf
}
  \item{seg2}{
  A segment data.frame obtained from copynumber::pcf
}
}
\details{
NA is used to fill up data from seg1 not found in seg2 or vice versa.  Columns mean and meanmap are from seg1, columns mean.1 and meanmap.1 are from seg2.  These values are extracted from the larger segments in seg1 and seg2 and are not recomputed in the smaller, intersecting segments. 
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
\notrun{
t.n.seg<-intersectSegments( t.seg, n.seg )
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
