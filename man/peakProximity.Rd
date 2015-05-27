\name{peakProximity}
\alias{peakProximity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Objective function
}
\description{
This is the objective function to be minimized. It is not meant to be called directly. 
}
\usage{
peakProximity(v, selectedPeaks, wd, verbose = T, npeaks=NULL , Sn=NULL, xonly, notSeenPenalty, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{
The parameter vector c(S,t[1],t[2],...,t[nsubclones]) where t[1] is the percentage of normal cells and t[j] j>1 is the percentage of subclone j. The last subclone is ignored in this vector ( sum( c(t[1],t[2],...,t[nsubclones] ) ) < 1 ).
}
  \item{selectedPeaks}{
Object created by selectPeaks(), or subset of lines from it (if usesubset is TRUE in coverParamSpace)}
\item{wd}{
Weights given to each point in a cluster of points, see eoDist.
}
  \item{verbose}{
Shut up or not.
}
  \item{npeaks}{
If selectedPeaks is a subset of the original data frame created by selectPeaks(), then npeaks specifies the original number of peaks (if usesubset is TRUE in coverParamSpace)
}
\item{Sn}{{Fixed value of S*n where n is the proportion of normal cells.}
\item{xonly}{Distance to peak to be calculated on the x-axis only}
\item{notSeenPenalty}{IGNORE, dropped from development}
\item{\dots}{
Additional arguments passed to eoDist.
}
}
\details{
Not meant to be called by the user. 
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

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
