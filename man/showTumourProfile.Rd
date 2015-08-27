\name{showTumourProfile}
\alias{showTumourProfile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Display the copy-number profile of the tumour.
%%  ~~function to do ... ~~
}
\description{
This function displays the copy-number profile of the tumour summarized in copyAr$ars[,1:2].
}
\usage{
showTumourProfile(copyAr, nx = 200, ny = 50, maxPoints = 50000, 
 selected = NULL, flatten = 1, xlim = c(0, 2), nlev = 50, 
 xaxt = "n", seed = NULL, nopoints = F, hist2d = F, chr = NULL, 
 noise = NULL, col = terrain.colors(50), onlypoints=F , plot=T )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{copyAr}{
 Object created from calling loadCopyAr.
%%     ~~Describe \code{copyAr} here~~
}
  \item{nx}{
Number of grid points in the x direction. Passed to kde2d or hist2d. 
}
  \item{ny}{
Number of grid points in the y direction. Passed to kde2d or hist2d. 
}
  \item{maxPoints}{
Only select and display that many lines from copyAr$ars, randomly selected.
}
  \item{selected}{
Not used.
}
  \item{flatten}{
Elevate the estimated density or height (estimated from a call to kde2d or hist2d) to the power given by flatten. This is to "elevate" smaller peaks.
}
  \item{xlim}{
Region of the x-axis to be plotted.
}
  \item{nlev}{
Number of levels used in the contour plot
}
  \item{xaxt}{
Argument passed to plot
}
  \item{seed}{
Seed for random number generator.
}
  \item{nopoints}{
Underneat (hidden) the contour plot are the actual points copyAr$ars[,1:2], so that a call to identify() can be made. This flag suppresses the plotting of these points. 
}
  \item{hist2d}{
Use gplots::hist2d instead of MASS::kde2d.
}
  \item{chr}{
Restrict plot to chromosome chr (character vector, as in "chr22")
}
  \item{noise}{
Add some noise to the copyAr$ars$copy value: copyAr$ars$copy + rnorm( length(copyAr$ars$copy), 0, noise).  May be needed if a single chromosome that only has one segment is to be plotted.
}
  \item{col}{
Color palette to use. Passed to contour().
}
 \item{onlypoints}{
Only plot pairs of points in copyAr. Do not plot contour plot. 
}
\item{plot}{Set to FALSE if no plot is required. The function still returns the output of kde2d.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}

\value{
(Invisible) list of values that are returned from kde2d or hist2d ($x, $y and $z)
}


\author{
Mathieu Lemire
}

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
