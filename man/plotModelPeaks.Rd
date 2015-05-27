\name{plotModelPeaks}
\alias{plotModelPeaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add expected peak locations 
}
\description{
Add expected peak locations to a plot obtained from showTumourProfile.
}
\usage{
plotModelPeaks(par = NULL, S = NULL, t = NULL, selectedPoints = NULL, cn, 
               epcol = "red", epcex = 3, eplwd = 3, eppch = 21, spcol = "black",
               sppch = 19, spcex = 1, splwd = 1, addlabels = F, preserveMatPatDiff=T, 
               preserveMaxSubClDiff=T , plot=TRUE , ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{
$par values returned from a call of coverParamSpace (as in li[[1]]$par ).  One of par or (S and t) must be given. 
%%     ~~Describe \code{par} here~~
}
  \item{S}{
The S parameter.
%%     ~~Describe \code{S} here~~
}
  \item{t}{
The t parameter vector: % of normal cells and % of each subclones. Sums to 1.  
%%     ~~Describe \code{t} here~~
}
  \item{selectedPoints}{
Set of points selected with selectPeaks.
}
  \item{cn}{
Data frame of allowed  copy-number configurations, typically created with prepCN from within coverParamSpace.
}
  \item{epcol}{
Color for expected peaks in epp.
}
  \item{epcex}{
Symbol expansion for expected peaks in epp.
}
  \item{eplwd}{
Line width for expected peaks in epp.
}
  \item{eppch}{
Plotting character for expected peaks in epp.
}
  \item{spcol}{
Color for selected peaks in selectedPoints. 
}
  \item{sppch}{
Plotting character for selected peaks in selectedPoints. 
}
  \item{spcex}{
Symbol expansion for selected peaks in selectedPoints. 
}
  \item{splwd}{
Line width for selected peaks in selectedPoints. 
}
  \item{addlabels}{
Flag to indicate that a call to addLabels is requested. 
}
 \item{preserveMatPatDiff}{
Flag to indicate if the sign of the difference between (say) the paternal and maternal chromosome copy-number is to be consistent in each clones. 
}
\item{preserveMaxSubClDiff}{
maxsubcldiff is a parameter of prepCN that indicates the upper bound for the difference in copy number between any two subclones (when >=1).  However, each parental chromosomes can have copy number differences greater than maxsubcldiff.   This flag is used to indicate that maxsubcldiff is to be applied to each parental chromosomes.
}
\item{plot}{ If FALSE, no points are plotted }

  \item{\dots}{
Additional parameters passed to addLabels.
}
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
