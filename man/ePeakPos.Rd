\name{ePeakPos}
\alias{ePeakPos}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Expected peak positions
%%  ~~function to do ... ~~
}
\description{
Given values of S and t parameters and a data frame of allowed copy-number configurations in each clone, the function returns the expected peak locations.
}
\usage{
ePeakPos(par = NULL, S = NULL, t = NULL, cn, preserveMatPatDiff = T, preserveMaxSubClDiff=T )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{
$par values returned from a call of coverParamSpace (as in li[[1]]$par ).  One of par or (S and t) must be given. 
}
  \item{S}{
The S parameter.
}
  \item{t}{
The t parameter vector: % of normal cells and % of each subclones. Sums to 1.  

}
  \item{cn}{
Data frame of allowed  copy-number configurations, typically created with prepCN.
}
  \item{preserveMatPatDiff}{
Flag to indicate if the sign of the difference between (say) the paternal and maternal chromosome copy-number is to be consistent in each clones. 
}
\item{preserveMaxSubClDiff}{
maxsubcldiff is a parameter of prepCN that indicates the upper bound for the difference in copy number between any two subclones (when >=1).  However, each parental chromosomes can have copy number differences greater than maxsubcldiff.   This flag is used to indicate that maxsubcldiff is to be applied to each parental chromosomes.
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
prepCN(6,2,.5)
epp<-ePeakPos( S=.5, t=c(.1,.7,.2) , cn=cn , preserveMatPatDiff=FALSE )
head(epp)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
