\name{addLabels}
\alias{addLabels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add labels to a plot.
}
\description{
Add labels that describe the copy number configuration of peaks from a showTumourProfile graph. 
}
\usage{
addLabels(epp, selectedPoints = NULL, manual = F, cex = 0.5, bg = "white", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{epp}{
Expected peak location obtained from a call to ePeakPos
}
  \item{selectedPoints}{
Selected peaks obtained from a call to selectPeaks
}
  \item{manual}{
Flag to indicate if labels are to be manually added. 
}
  \item{cex}{
Character expansion.
}
  \item{bg}{
Backgound color for the labels.  Passed to legend.
}
  \item{\dots}{
Additional arguments passed to legend.
}
}
\details{
Labels are added with a call to legend.  If manual is TRUE, the user is prompted to left-click near expected peaks (epp) to add their associated copy number configurations (right click when done). The format for the copy number configuration is "ab/cd/..", where a and b are the number of maternal and paternal chromosomes in sub-clone 1, c and d the number of maternal and paternal chromosomes in sub-clone 2, etc. 
}

\author{
Mathieu Lemire
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{\dontrun{
showTumourProfile(copyAr, maxPoints=50000 , flatten=.25 , nlev=20, 
       seed=12345  , xlim=c(0,2) , nx=200, ny=50  )
prepCN(12)
ePP<-plotModelPeaks( S= 0.6260699, t=c(0.02222, 0.97778),
       selectedPoints=sp,cn=cn, epcol="red",epcex=1,eplwd=3 )
addLabels( ePP, manual=T)
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
