\name{arInSeg}
\alias{arInSeg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate the proportion of copies of chromosomes that have the same parental origin, in each segment.
}
\description{
Estimate the proportion of copies of chromosomes that have the same parental origin, in each segment.
}
\usage{
arInSeg(seg, ar, minhet = 50, maxhet = 5000, control = list(maxit = 1000), tumourrangedata = NULL, maskmap = 0.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seg}{
data.frame obtained from intersecSegments, between a read-count-derived seg1 (from segmentSeqData) and a allelic-ratio-derived seg2 (from segmentAR)
}
  \item{ar}{
The ar data.frame passed to segmentAR.
}
  \item{minhet}{
If the number of heterozygous positions in a segment is less than this, skip the segment.
}
  \item{maxhet}{
If the number of heterozygous positions in a segment is greater than this, randomly select that many.
}
  \item{control}{
Not used, placeholder for other optimization function
}
  \item{tumourrangedata}{
Object of class IRanges containing tumour data, obtained from a call to gcCorrect}
  \item{maskmap}{
Heterozygous positions falling in bins with mappability less than this are ignored.
}
}
\details{
Semantically speaking, if p is the proportion of copies that are of, say, paternal origin, the function returns an estimate of either p or 1-p.  In other words, the actual proportion of copies that are derived from the paternal chromosome is unknown.  Moreover, these proportions are not comparable between segments (the parental origin "attached" to a proportion can change across segments). 
}
\value{
The seg data.frame is returned, with an extra column named "p" containing the estimates.
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
\examples{\dontrun{
arInSeg( t.ar.seg, ar,  tumourrangedata=tc , minhet = 50 )
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
