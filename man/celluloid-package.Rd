\name{celluloid-package}
\alias{celluloid-package}
\alias{celluloid}
\docType{package}
\title{
CELLULOID: cellularity and ploidy tool
}
\description{
This package contains functions to calculates the cellularity, the ploidy and the purity of a whole-genome sequenced tumour sample (high coverage); functions to scale gc- and mappability-corrected read counts into absolute copy number; functions to segment data and plot copy number graphs.
}
\details{
\tabular{ll}{
Package: \tab celluloid\cr
Type: \tab Package\cr
Version: \tab 0.11\cr
Date: \tab 2015-02-27\cr
License: Apache License 2.0
}

In the following, read counts is typically short for gc- and mappability-corrected read counts, scaled to the average read counts in the autosomal genome.

The sequenced sample consists of all cells that were sequenced, normal (if any) and tumour cells. 

S is a parameter representing the expected read count in segments where all cells display two copies.  Ploidy of the sequenced sample (the number of haploid autosomal genomes averaged over all sequenced cells) is 2/S. 

t is the parameter vector c(t[1],...,t[nsubclones+1]), sum(t)=1,  where components represent the percentage of normal cells in the sequenced sample (t[1]) and percentage of cells belonging to each subclone (t[j] for j>1; nsubclones is the number of subclones that needs to be pre-specified).

}
\author{
Mathieu Lemire <mathieu.lemire@oicr.on.ca>
}

\keyword{ package }
