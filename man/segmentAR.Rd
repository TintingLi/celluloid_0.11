\name{segmentAR}
\alias{segmentAR}
\title{
Function to segment allelic ratio data}
\description{
Function to segment allelic ratio data
}
\usage{
segmentAR(ar, tumourrangedata = NULL, gamma = 500, kmin = 100, maskmap = 0.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ar}{
Allelic ratio data frame. See Details. 
}
  \item{tumourrangedata}{
Object of class IRange produced by function gcCorrect.
}
  \item{gamma}{
Parameter of copynumber::pcf
}
  \item{kmin}{
Parameter of copynumber::pcf
}
 \item{maskmap}{
Ignore SNPs falling in bins with mappability less than maskmap. Only used if tumourrangedata is provided.
}
}
\details{
ar must have the following columns:

head(ar)
   CHR    POS REF_COUNT VAR_COUNT
1 chr1 754730        31        17
2 chr1 754813        10        23
3 chr1 754840        30        10
4 chr1 754873        16        20
5 chr1 755955        26        14
6 chr1 758555        11        11
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
function (ar, tumourrangedata = NULL, gamma = 500, kmin = 100, 
    maskmap = 0.8) 
{
    ar.tmp <- ar
    if (maskmap > 0 & !is.null(tumourrangedata)) {
        binlength <- end(tumourrangedata[1, ]) - start(tumourrangedata[1, 
            ]) + 1
        bin <- paste(ar[, 1], binlength * floor(ar[, 2]/binlength) + 
            1, sep = "-")
        sel <- tumourrangedata$map < maskmap
        mask <- paste(tumourrangedata$space, start(tumourrangedata), 
            sep = "-")[sel]
        ar.tmp <- ar[!is.element(bin, mask), ]
    }
    tmp <- data.frame(chrom = ar.tmp[, 1], pos = ar.tmp[, 2], 
        ar = ar.tmp[, 3]/(ar.tmp[, 3] + ar.tmp[, 4]))
    tmp$chrom <- factor(as.character(tmp$chrom), levels = c(paste("chr", 
        1:22, sep = ""), "chrX", "chrY", "chrM"))
    tmp$chrom <- as.numeric(tmp$chrom)
    tmp <- tmp[tmp$chrom < 23, ]
    tmp <- tmp[!is.na(tmp$chrom), ]
    sel <- tmp$ar > 0.5
    tmp$ar[sel] <- 1 - tmp$ar[sel]
    require(copynumber)
    tmp.win <- winsorize(tmp)
    ar.seg <- pcf(data = tmp.win, gamma = 5000, kmin = 100, digits = 4)
    ar.seg$chrom <- paste("chr", ar.seg$chrom, sep = "")
    colnames(ar.seg)[7] <- "meanar"
    return(ar.seg)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
