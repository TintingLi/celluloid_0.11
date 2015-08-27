\name{removePeaks}
\alias{removePeaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Remove peaks from data.frame
}
\description{
Remove peaks from set of points obtained from selectPeaks
}
\usage{
removePeaks(selectedPoints)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{selectedPoints}{
data.frame obtained from a call of selectPeaks.
}
}
\details{
Left click to remove a point, right click to end. 
}

\author{
Mathieu Lemire
}

\examples{
\dontrun{
cntr<-showTumourProfile(copyAr, maxPoints=50000 , flatten=.25 , nlev=20, 
        xlim=c(0,2) , nx=200, ny=50 )
sp <-  selectPeaks( cntr, copyAr , getLocal=T , 
                    percentMax=.33, manual=F, nrand=200 , filtersymm=T   ) 
sp <- removePeaks(sp)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
