\name{scaleReadCounts}
\alias{scaleReadCounts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scale read counts
}
\description{
Scale gc- and mappability-corrected read counts to absolute copy number.
}
\usage{
scaleReadCounts(tumourrangedata, epp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tumourrangedata}{
Object of class IRanges containing tumour data, obtained from a call to gcCorrect
}
  \item{epp}{
Expected peak positions. Data frame returned by ePeakPos. 
}
}
\details{
An additional column named icopy is added to tumourrangedata (tumourrangedata$icopy). Scaling is done by using the relationship given in epp: each integer copy number configurations in tumour (m1,p1,m2,p2,etc) is associated with the expected gc- and mappability-corrected read count (in x).  Values for icopy are in the copy-number scale (as opposed to tumourrangedata$copy which is in a log scale).
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
prepCN(12)
ePP<-ePeakPos( S=0.624, t=c(0.073, 0.923), cn=cn  )
tcs<- scaleReadCounts( tc , ePP )
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
