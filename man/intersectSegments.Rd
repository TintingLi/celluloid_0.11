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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (seg1, seg2) 
{
    interseg <- c()
    for (ch in union(seg1$chrom, seg2$chrom)) {
        for (arm in c("p", "q")) {
            sub1 <- seg1[seg1$chrom == ch & seg1$arm == arm, 
                ]
            sub2 <- seg2[seg2$chrom == ch & seg2$arm == arm, 
                ]
            if (nrow(sub1) + nrow(sub2) > 0) {
                bounds <- sort(unique(c(sub1$start.pos, sub1$end.pos + 
                  1, sub2$start.pos, sub2$end.pos + 1)))
                for (i in 1:(length(bounds) - 1)) {
                  s <- bounds[i]
                  e <- bounds[i + 1]
                  sel1 <- sub1$start.pos <= s & sub1$end.pos + 
                    1 >= e
                  sel2 <- sub2$start.pos <= s & sub2$end.pos + 
                    1 >= e
                  if (sum(sel1)) {
                    line1 <- as.data.frame(sub1[sel1, 7:ncol(sub1)])
                  }
                  else {
                    line1 <- as.data.frame(t(rep(NA, ncol(sub1) - 
                      6)))
                  }
                  if (sum(sel2)) {
                    line2 <- as.data.frame(sub2[sel2, 7:ncol(sub2)])
                  }
                  else {
                    line2 <- as.data.frame(t(rep(NA, ncol(sub2) - 
                      6)))
                  }
                  colnames(line1) <- names(sub1)[7:ncol(sub1)]
                  colnames(line2) <- names(sub2)[7:ncol(sub2)]
                  interseg <- rbind(interseg, data.frame(sampleID = NA, 
                    chrom = ch, arm = arm, start.pos = s, end.pos = e - 
                      1, size = e - s, line1, line2))
                }
            }
        }
    }
    return(interseg)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
