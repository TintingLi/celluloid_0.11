\name{prepCopyAr}
\alias{prepCopyAr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create copyAr object.
}
\description{
This function creates copy-number-related objects. The main component contains, 
at each germline heterozygous positions, the allelic ratio in the tumour (proportion of 
reads supporting the allele found on the reference sequence) and the average read counts 
(gc- and mappability-corrected) of the segments they fall into.}
\usage{
prepCopyAr(seg, ar, tumourrangedata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seg}{
A segment data.frame obtained from functions such as arInSeg}
  \item{ar}{
The ar data.frame passed to segmentAR
}
  \item{tumourrangedata}{
Object of class IRanges containing tumour data, obtained from a call to gcCorrect
}
}

\author{
Mathieu Lemire}

\examples{
\dontrun{
copyAr<-  prepCopyAr( t.ar.seg[!t.ar.seg$mask,] , ar,  tc  )
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
