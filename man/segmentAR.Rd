\name{segmentAR}
\alias{segmentAR}
\title{
Function to segment allelic ratio data}
\description{
Function to segment allelic ratio data
}
\usage{
segmentAR(ar, tumourrangedata = NULL, gamma = 500, kmin = 100, 
          maskmap = 0.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ar}{
Allelic ratio data frame. See Details. 
}
  \item{tumourrangedata}{
Object of class IRange produced by function gcCorrect.
}
  \item{gamma}{
Parameter of copynumber::pcf
}
  \item{kmin}{
Parameter of copynumber::pcf
}
 \item{maskmap}{
Ignore SNPs falling in bins with mappability less than maskmap. Only used if tumourrangedata is provided.
}
}
\details{
ar must have the following columns:
\verb{

head(ar)
   CHR    POS REF_COUNT VAR_COUNT
1 chr1 754730        31        17
2 chr1 754813        10        23
3 chr1 754840        30        10
4 chr1 754873        16        20
5 chr1 755955        26        14
6 chr1 758555        11        11
}
}
\author{
Mathieu Lemire
}

\examples{
\dontrun{
ar.seg<- segmentAR( ar, tc ) 
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
