\name{eoDist}
\alias{eoDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Total distance between observed and (closest) expected peaks. 
}
\description{
Calculates the total distance between observed peaks and (closest) expected peaks, with added penalty to force a specific copy number configuration to one of the peak.  Not meant to be called directly. 
}
\usage{
eoDist(selectedPeaks, S, t, wd, xonly = F, weightbydiff = F, weightbyarp = F, 
 penaltymultiplier = 1000, penaltymultipliery = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{selectedPeaks}{
Passed from peakProximity
}
  \item{S}{
The S parameter
}
  \item{t}{
The t parameter (sums to 1).
}
 \item{wd}{
Weights given to each point in a cluster of points
}
  \item{xonly}{
Ignore the y-axis, only calculate distance based on the x-axis
}
  \item{weightbydiff}{
Ignore, not extensively tested.
}
  \item{weightbyarp}{
Ignore, not extensively tested.
}
  \item{penaltymultiplier}{
If a peak (a line of selectedPeak) is to be forced to represent a particular copy number configuration, then the fit is forced by adding a penalty multiplier to the distance between that peak and the expected location of the copy number configuration.
}
  \item{penaltymultipliery}{
If the penalty multiplier is to differ in the x- and y- directions, then use this as an additional multiplier in the y- direction (the penalty multiplier in the y-direction is then penaltymultiplier*penaltymultipliery).
}
  \item{\dots}{
Arguments passed to prepAR (ie, preserveMatPatDiff)
}
}
\details{
To force a peak, a third column ($w) in selectedPeak is needed. That column contains a single non-zero integer value taking the format xy (for a 2-subclone model), where x is the copy number in the first sub-clone, y in the second subclone. In a 3-subclone system, the format is xyz, etc.  All integers x,y,z,etc must take non-negative values no larger than 9. 
}

\author{
Mathieu Lemire
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
