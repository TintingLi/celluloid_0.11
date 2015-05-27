\name{prepCopyAr}
\alias{prepCopyAr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create copyAr object.
}
\description{
This function creates copy-number-related objects. The main component contains, 
at each germline heterozygous positions, the allelic ratio in the tumour (proportion of 
reads supporting the allele found on the reference sequence) and the average read counts 
(gc- and mappability-corrected) of the segments they fall into.}
\usage{
prepCopyAr(seg, ar, tumourrangedata)
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (seg, ar, tumourrangedata) 
{
    require(data.table)
    tmp <- tumourrangedata[1, ]
    binlength <- end(tmp) - start(tmp) + 1
    cat("...extracting segment mean and assigning to each bin\n")
    meaninbins <- apply(seg[, c("chrom", "start.pos", "end.pos", 
        "mean")], 1, function(x) {
        cat(".")
        fr <- floor(as.numeric(x[2])/binlength) * binlength + 
            1
        to <- floor(as.numeric(x[3])/binlength) * binlength
        if (to < fr) {
            to <- fr
        }
        tmp <- paste(x[1], seq(fr, to, binlength), sep = ":")
        tmp <- data.frame(tmp, as.numeric(x[4]))
        names(tmp) <- c("bin", "mean")
        return(tmp)
    })
    cat(".")
    meaninbins <- rbindlist(meaninbins)
    cat(".")
    binlabels <- data.frame(line = 1:nrow(tumourrangedata), bin = paste(tumourrangedata$space, 
        ":", start(tumourrangedata), sep = ""))
    cat(".")
    m <- merge(binlabels, meaninbins, by = "bin", all.x = T)
    tumourrangedata$smcopy <- NA
    cat(".")
    tumourrangedata$smcopy <- m$mean[order(m$line)]
    cat(".")
    tumourrangedata$smcopy[tumourrangedata$ignore] <- NA
    cat(".")
    ars <- data.frame(chpo = paste(tumourrangedata$space, start(tumourrangedata), 
        sep = "-"), copy = tumourrangedata$smcopy)
    ars <- ars[!tumourrangedata$ignore, ]
    cat(".")
    alratio <- ar
    alratio <- alratio[alratio[, 3] + alratio[, 4] > 0, ]
    cat(".")
    alratio$chpo <- paste(alratio[, 1], floor(alratio[, 2]/binlength) * 
        binlength + 1, sep = "-")
    cat("\n...merging with ref/alt read counts")
    copyAr <- merge(ars, alratio, by = "chpo")
    copyAr <- copyAr[, c(3, 4, 2, 5, 6)]
    names(copyAr) <- c("chr", "pos", "copy", "ref", "alt")
    copyAr <- copyAr[order(copyAr$chr, copyAr$pos), ]
    copyAr <- copyAr[!is.na(copyAr$copy), ]
    rownames(copyAr) <- 1:nrow(copyAr)
    cat("...done.\n")
    return(copyAr)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
