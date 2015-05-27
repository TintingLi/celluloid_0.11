\name{makeCopyAr}
\alias{makeCopyAr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create copyAr object.
}
\description{This function creates copy-number-related objects. The main component contains, 
at each germline heterozygous positions, the allelic ratio in the tumour (proportion of 
reads supporting the allele found on the reference sequence) and the average read counts 
(gc- and mappability-corrected) of the segments they fall into.
}
\usage{
makeCopyAr(seg, ar, tumourrangedata = NULL, maskmap = 0.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seg}{
A segment data.frame obtained from functions such as arInSeg}
  \item{ar}{
The ar data.frame passed to segmentAR
}
  \item{tumourrangedata}{
Object of class IRanges containing tumour data, obtained from a call to gcCorrect
}
  \item{maskmap}{
Heterozygous positions falling in bins with mappability less than this are ignored.
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
function (seg, ar, tumourrangedata = NULL, maskmap = 0.8) 
{
    binlength <- NA
    ar.tmp <- ar
    ar.tmp$ar <- ar[, 3]/(ar[, 3] + ar[, 4])
    if (maskmap > 0 & !is.null(tumourrangedata)) {
        binlength <- end(tumourrangedata[1, ]) - start(tumourrangedata[1, 
            ]) + 1
        bin <- paste(ar.tmp[, 1], binlength * floor(ar.tmp[, 
            2]/binlength) + 1, sep = "-")
        sel <- tumourrangedata$map < maskmap
        mask <- paste(tumourrangedata$space, start(tumourrangedata), 
            sep = "-")[sel]
        ar.tmp <- data.frame(bin = bin, ar.tmp)
        ar.tmp <- ar.tmp[!is.element(bin, mask), ]
    }
    else {
        binlength <- 1000
        bin <- paste(ar.tmp[, 1], binlength * floor(ar.tmp[, 
            2]/binlength) + 1, sep = "-")
        ar.tmp <- data.frame(bin = bin, ar.tmp)
    }
    meaninbins <- apply(data.frame(seg[, c(2, 4, 5, 7)]), 1, 
        function(x) {
            tmp <- paste(x[1], seq(as.numeric(x[2]), as.numeric(x[3]) - 
                1, binlength), sep = "-")
            tmp <- data.frame(tmp, as.numeric(x[4]))
            names(tmp) <- c("bin", "mean")
            return(tmp)
        })
    require(data.table)
    meaninbins <- rbindlist(meaninbins)
    setnames(meaninbins, c("bin", "copy"))
    rcar <- merge(ar.tmp, meaninbins, by = "bin")
    rcar <- rcar[rcar[, 5] + rcar[, 6] > 0, ]
    copyAr <- list(copyAr = rcar[, c("copy", "ar", "CHR", "POS")])
    return(copyAr)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
