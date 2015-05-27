\name{prepCN}
\alias{prepCN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Allowed copy number values in each clones
}
\description{
This prepares a global object called cn, that contains, for normal cells and all subclones, all allowed copy number combinations. Normal cells are assumed to have 2 copies of all autosomal chromosomes. 
}
\usage{
prepCN(maxc, nsubcl=1, maxsubcldiff=NULL, maxlines=10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{maxc}{
Maximum allowed copy number in any subclone. 
}
  \item{nsubcl}{
Number of subclones.  ncol(cn)=nsubcl+1.
}
  \item{maxsubcldiff}{
If >= 1, the upper bound for the difference in copy number between any two subclones.  If < 1, the upper bound for the ratio between
the maximum and minimum copy number seen across all subclones.  See Details. 
}
  \item{maxlines}{
Stop with an error message if the number of configurations is greater than maxlines. 
}
}
\details{

In details, assume two subclones (nsubcl<-2).  First cn is formed by listing in a two-column data.frame all (maxc+1)^2 combinations of numbers from 0:maxc:

 cn<- expand.grid( rep( list( 0:maxc ), nsubcl ) )


Then lines are selected out based on the value of maxsubcldiff.  If maxsubcldiff>=1:

 sel<- apply( cn,1, function(x){ max( dist( x ) ) } ) <= maxsubcldiff 
 cn<-cn[sel,] 

If maxsubcldiff<1:

 sel<-apply( cn,1, function(x){ mn<-max(1,min(x)) ; mx<-max(x); return( mn/mx>=maxsubcldiff ) } ) 
 cn<-cn[sel,]

Note that 0 copies is treated as if it was 1. 

}
\value{
An global object with the name cn is created. This object is expected in some functions. 
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
\examples{
 prepCN(6,2,.5)
 head(cn)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
