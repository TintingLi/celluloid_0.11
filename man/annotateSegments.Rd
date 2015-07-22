\name{annotateSegments}
\alias{annotateSegments}
\title{
Adds annotations to a segment data.frame
}
\description{
For each segment in a segment data.frame, the function finds the closest point in the given ePeakPos solution, and annotates the segment with the corresponding number of (say) paternal/maternal copies, in each subclones.  
}
\usage{
annotateSegments(seg, epp, weigths="quadratic")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seg}{
segment data.frame obtained from a call to scaleSegments.
}
  \item{epp}{
data.frame returned from a call to ePeakPos. 
}
\item(weights){Function to assign weights given to segments based on distance to nearest epp points. One of linear or quadratic.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
%%  ~~who you are~~
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
function (seg, epp) 
{
    tmpseg <- seg
    tmpseg$labels <- NA
    for (s in 1:nrow(tmpseg)) {
        if (s == 1 | s\%\%10 == 0) {
            cat("annotating segment ", s, " of ", nrow(tmpseg), 
                "...\n", sep = "")
        }
        if (!is.na(tmpseg[s, "mean"]) & !is.na(tmpseg[s, "p"])) {
            tmp <- rbind(data.frame(x = tmpseg[s, "mean"], ar = tmpseg[s, 
                "p"]), epp[, c("x", "ar")])
            d <- as.matrix(dist(tmp))
            sel <- which(d[1, -1] == min(d[1, -1]))[1]
            m <- epp[sel, ][seq(3, ncol(epp) - 2, 2)]
            p <- epp[sel, ][seq(4, ncol(epp) - 2, 2)]
            lab <- paste(paste(m, p, sep = ""), collapse = "/")
            tmpseg$labels[s] <- lab
        }
    }
    return(tmpseg)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
