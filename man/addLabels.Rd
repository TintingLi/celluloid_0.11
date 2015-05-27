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
