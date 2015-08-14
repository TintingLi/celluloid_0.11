\name{ARloh}
\alias{ARloh}
Expected allelic ratio in a LOH region  
\title{
Expected allelic ratio in a LOH region with scaled read count x
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
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
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, S, n) 
{
    k <- 2 * (x - S * n)/(S - S * n)
    return(n/(2 * n + (1 - n) * k))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
