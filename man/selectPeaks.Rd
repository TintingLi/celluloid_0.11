\name{selectPeaks}
\alias{selectPeaks}
\title{
Select peaks from a showTumourProfile graph.
%%  ~~function to do ... ~~
}
\description{
Select peaks from a showTumourProfile graph, either manually (recommended) or by random calls to optim, or both.
}
\usage{
selectPeaks(cntr, copyAr, manual = T, getLocalMax = F, percentMax = 0.05, 
  nrand = 100, filtersymm = T,  autocol="red", autopch=19, autocex=1, 
  manucol="black", manupch=19, manucex=1 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cntr}{
Values returned from a call to showTumourProfile
}
  \item{copyAr}{
The copyAr object passed to showTumourProfile
}
  \item{manual}{
Flag to indicate manual selection of peaks is required. 
}
  \item{getLocalMax}{
Get peaks from nrand calls of the optim function.
}
  \item{percentMax}{
Only return a peak if it is at least percentMax times the value of the max peak.
}
  \item{nrand}{
Use nrand calls of optim to find local maxima, with random starting values each time.
}
 \item{filtersymm}{Flag to indicate that points that are mirrored around y=0.5 should be merged into a single point.}
 \item{autocol}{color of points added to graph when getLocalMax=T}
 \item{autopch}{pch of points added to graph when getLocalMax=T}
 \item{autocex}{cex of points added to graph when getLocalMax=T}
  \item{manucol}{color of points added to graph when manual=T}
 \item{manupch}{pch of points added to graph when manual=T}
 \item{manucex}{cex of points added to graph when manual=T}
}


\details{
When manually selecting peaks (manual = T), left-click to select peaks and right click to exit the function.  When getLocalMax = F, the function find local maxima by calling the function optim a number of times (nrand times) with random starting values each time. 

The contour plot obtained by showTumourProfile should be symmetrical around y=0.5.  The filtersymm=T flag indicates that pairs of points that are mirror of each other should be merged into a single point. Some peaks may be located on the line y=0.5; if filtering of mirrored points is not performed, these peaks will only contribute half the weight of two mirrored points in the objective function to be minimized (each selected peak contribute to the objective function), even though mirrored points do not provide additional information on the tumour copy number profile. 

}
\value{
A data.frame of x and y coordinates. 
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
