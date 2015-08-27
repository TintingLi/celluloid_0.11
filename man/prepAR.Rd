\name{prepAR}
\alias{prepAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate expected allelic ratios
}
\description{
Calculates expected allelic ratios given percentage of normal cells and subclones. Not meant to be called directly.
}
\usage{
prepAR(t, preserveMatPatDiff = T, preserveMaxSubClDiff=T )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
The t parameter, vector containing percentages of normal and subclones. Sums to 1. }
 \item{preserveMatPatDiff}{
Flag to indicate if the sign of the difference between (say) the paternal and maternal chromosome copy-number is to be consistent in each clones. 
}
\item{preserveMaxSubClDiff}{
maxsubcldiff is a parameter of prepCN that indicates the upper bound for the difference in copy number between any two subclones (when >=1).  However, each parental chromosomes can have copy number differences greater than maxsubcldiff.   This flag is used to indicate that maxsubcldiff is to be applied to each parental chromosomes.
}
}
\details{
Not meant to be called by the user.
}

\author{
Mathieu Lemire}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
