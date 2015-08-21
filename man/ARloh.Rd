\name{ARloh}
\alias{ARloh}
\title{Expected allelic ratio in a LOH region  
%%  ~~function to do ... ~~
}
\description{Expected allelic ratio in a LOH region with scaled read count x

}
\usage{
ARloh(x, S, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Value on the x axis of a tumor contour plot
}
  \item{S}{
Parameter S
}
  \item{n}{
Proportion of normal cells in the sample.
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
\examples{\dontrun{
# select a few peaks or segments that are LOH
sp <-  selectPeaks( cntr, copyAr , getLocal=F ) 
x<-sp$x; y<-sp$y
nls.fit<- nls( y ~ ARloh( x, 1 ,n ) , start=list(n=0.01  ) , upper=list( n=1 ), 
               lower=list(n=0 ) , algo="port" )
Sn<-  summary(nls.fit)$coefficients[1,1]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
