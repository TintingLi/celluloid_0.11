\name{coverParamSpace}
\alias{coverParamSpace}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Search for a fit between observed peaks and expected peaks.
}
\description{
Function to set up the optimization.
}
\usage{
coverParamSpace( selectedPeaks=NULL, segments=NULL, verbose=T , addToParamSpace=F , 
      control=NULL ,  Sfrom=NULL, Sto=NULL,  Sn=NULL, maxc=NULL, 
      maxsubcldiff=NULL, optimFct=2, lowerF, upperF, nrep=NULL, usesubsets=NULL , 
      xonly=FALSE,modeat=NULL, weight=NULL, notSeenPenalty=TRUE, method=NULL, ...  )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{selectedPeaks}{
data.frame of selected peaks, obtained from calling selectPeaks. The peak-based objective function is then used.
}
\item{segments}{Tumour segment data frame, such as t.ar.seg. Must contain columns "mean","p" and "size". The segment-based objective function is then used.}
  \item{verbose}{
More or less details in the output
}
  \item{addToParamSpace}{
A call to coverParamSpace crates a global dataframe with the name paramSpace. This flag indicates if we want to append to it.  Useless for grid search (when length(optimFct)> 1) as it is overwritten.  
}
  \item{control}{
Arguments passed to the optimization function used. 
}
  \item{Sfrom}{
Lower bound for S. 
}
  \item{Sto}{
Upper bound for S.
}
\item{Sn}{Fixed value of S*n where n is the proportion of normal cells.}
  \item{maxc}{
Argument passed to prepCN.
}
  \item{maxsubcldiff}{
Argument passed to prepCN.
}
  \item{optimFct}{
The optim function used: optim (optimFct=1); GenSA (optimFct=2). If optimFct is of length > 1, then grids are used as starting points to optim. See Details.
}
  \item{lowerF}{
Lower bound for c(t[1],...,t[nsubclones]) (see Details). where t[1] represents the percentage of normal cells, t[j] (j>1) the percentage of tumour cells in subclone j-1 (nsubcl is the number of subclones). Note that the last subclone is
}
  \item{upperF}{
Upper bound for  c(t[1],...,t[nsubclones]), where t is the usual parameter vector t. 
}
  \item{nrep}{
That many calls to the optimization function will be made.  The list that is returned have length equal to nrep. 
}
  \item{usesubsets}{
NULL or integer value. Randomly select at least "usesubsets" rows of selectedPeaks. The objective function ignore the other lines. 
}
 \item{xonly}{
 Ignore the y-axis values of the peaks (the allelic ratio axis).
}
\item{modeat}{ 
Force the most "frequent" or "highest" peak to represent "modeat" number of copies. Not extensively tested and incomplete, ignore.
}
\item{weight}{
Related to modeat, ignore.
}
\item{notSeenPenalty}{
Only for the case of no subclones. Selected peaks and expected peaks are projected onto the x-axis (the read count/copy number axis). This imposes an additional penalty (distance) for each expected peak that are far from any selected peak, only within the range of selected peaks 
}
\item{method}{ Method argument for optim function. Defaults to "L-BFGS-B" if NULL. See ?optim}
  \item{\dots}{
Arguments passed to peakProximity (the objective function to be minimized).
}
}
\details{
In c(t[1],...,t[nsubclones]), t[1] represents the percentage of normal cells t[j] (j>1) the percentage of tumour cells in subclone j-1 (nsubclones is the number of subclones). Note that the last subclone is absent from that vector.  In other words, lowerF and upperF are the lower and upper bounds for the usual parameter vector t, ignoring its last element.  

If optimFct is of length > 1, then grid search is performed with optim() using a number of starting values equal to prod(optimFct). Starting values for parameters S and t[1],...,t[nsubclones] are: seq(Sfrom, Sto, len=optimFct[1]) (for S),  seq( lowerF[i], upperF[i], len=optimFct[i+1] ) (for t[i]).
}
\value{
A list of length nrep if length(optimFct)==1 or prod(optimFct) if length(optimFct)>1. 

Each component of that list is a list that includes at the minimum the values:

  \item{value }{The value of the objective function corresponding to par.}
  \item{par }{Vector: the optimal parameters found.}

Other values returned by the optimization algorithms are also returned. See the man pages of these functions. 
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
\examples{\dontrun{
cntr<-showTumourProfile(copyAr, maxPoints=50000 , flatten=.25 , nlev=20, 
       seed=12345  , xlim=c(0,2) , nx=200, ny=50 )
sp <-  selectPeaks( cntr, copyAr , getLocal=F ) 
li1<-coverParamSpace( selectedPeaks=sp , optimFct=2, lowerF=c(0), upperF=c(1),  
       Sfrom=.25, Sto=2 , maxc=12 , control=list( maxit=1000  ) )
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
