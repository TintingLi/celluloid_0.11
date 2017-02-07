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
ePeakPos(par = NULL, S = NULL, t = NULL, cn, preserveMatPatDiff=T, 
          preserveMaxSubClDiff=T , xonly=FALSE )
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
\item{xonly}{Ignore allelic ratios, in which case the ar column is filled with dummy 0.5 values and the number of "paternal" copies arbitrarily set to 0.}
}

\author{
Mathieu Lemire}

\examples{
prepCN(12)
ePP <- ePeakPos( S= 0.6260699, t=c(0.02222, 0.97778), cn=cn  )
head(ePP)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
