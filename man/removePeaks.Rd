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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (selectedPoints) 
{
    cat("left click on peaks, right click to finish\n")
    id <- identify(selectedPoints[, 1], selectedPoints[, 2], 
        n = 1, labels = "")
    while (length(id) > 0) {
        points(selectedPoints[id, ], pch = "X", col = "black")
        remove <- c(remove, id)
        id <- identify(selectedPoints[, 1], selectedPoints[, 
            2], n = 1, labels = "")
    }
    return(selectedPoints[!is.element(1:nrow(selectedPoints), 
        remove), ])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
