\name{getLocalSolutions}
\alias{getLocalSolutions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prune output of coverParamSpace 
}
\description{
Takes results from coverParamSpace and prune out solutions that are close one another, with pretty print. 
}
\usage{
getLocalSolutions(li, max=TRUE, pruneS = 0.05, pruneT = 0.025, ploidy = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{li}{
List returned from coverParamSpace
}
\item{max}{Set to TRUE if coverParamSpace was set to maximize, or FALSE if was set to minimize}
  \item{pruneS}{
Skip the solution if its S is within +/- pruneS of the S of a better solution 
}
  \item{pruneT}{
Skip the solution if its t parameter vector is within +/- pruneT of the t of a better solution 
}
  \item{ploidy}{
Vector of length 2. Only output solutions consistent with a tumour ploidy above or equal to ploidy[1] and below or equal to ploidy[2].
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
Mathieu Lemire}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
lm<- coverParamSpace( segments=segmentsubset, control=list( maxit=1000 ), Sn=Sn ,
                      maxc=12,optimFct=1 , nrep=50, 
                      lowerF=c(0), upperF=c( Sn/.25 ), addToParamSpace=T  )

localSolutions <-getLocalSolutions(lm, max=TRUE) 
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
