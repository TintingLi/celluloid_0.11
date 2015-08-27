\name{segmentSeqData}
\alias{segmentSeqData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Copy number segmentation
}
\description{
Wrapped for function copynumber::pcf to do segmentation.
}
\usage{
segmentSeqData(rangedata, gamma = 500, kmin = 100, maskmap = 0.8, maskadj = FALSE, skipmeanmap=FALSE )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rangedata}{
Object of class IRange produced by function gcCorrect.
}
  \item{gamma}{
Parameter of copynumber::pcf
}
  \item{kmin}{
Parameter of copynumber::pcf
}
  \item{maskmap}{
Ignore bins with mappability less than this
}
  \item{maskadj}{
Ignore bins adjacent to bins with mappability less than maskmap. NOT USED.
}
\item{skipmeanmap}{Skip call to meanMapInSeg}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Data fram returned by copynumber::pcf, to which is added a column named meanmap, that contains the mean mappability calculated in each segment.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
Mathieu Lemire}

\examples{
\dontrun{
n<- wigsToRangedData( normalWigFile, gcWigFile, mapWigFile )
nc<-gcCorrect(n) 
n.seg <- segmentSeqData( nc, k=50 , maskmap = 0.8 )
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
