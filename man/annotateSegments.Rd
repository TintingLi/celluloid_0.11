\name{annotateSegments}
\alias{annotateSegments}
\title{
Adds annotations to a segment data.frame
}
\description{
For each segment in a segment data.frame, the function finds the closest point in the given ePeakPos solution, and annotates the segment with the corresponding number of (say) paternal/maternal copies, in each subclones.  
}
\usage{
annotateSegments(seg, epp, weights="quadratic")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seg}{
segment data.frame obtained from a call to scaleSegments.
}
  \item{epp}{
data.frame returned from a call to ePeakPos. 
}
\item{weights}{Function to assign weights given to segments based on distance to nearest epp points. One of linear or quadratic.}
}

\author{
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{\dontrun{
prepCN(12)
ePP <- ePeakPos( S= 0.6260699, t=c(0.02222, 0.97778), cn=cn  )
segments<-scaleSegments(t.ar.seg ,  ePP )
segments<-annotateSegments(segments, ePP)
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
