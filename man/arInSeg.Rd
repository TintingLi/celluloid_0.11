\name{arInSeg}
\alias{arInSeg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate the proportion of copies of chromosomes that have the same parental origin, in each segment.
}
\description{
Estimate the proportion of copies of chromosomes that have the same parental origin, in each segment.
}
\usage{
arInSeg(seg, ar, minhet = 50, maxhet = 5000, control = list(maxit = 1000), tumourrangedata = NULL, maskmap = 0.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seg}{
data.frame obtained from intersecSegments, between a read-count-derived seg1 (from segmentSeqData) and a allelic-ratio-derived seg2 (from segmentAR)
}
  \item{ar}{
The ar data.frame passed to segmentAR.
}
  \item{minhet}{
If the number of heterozygous positions in a segment is less than this, skip the segment.
}
  \item{maxhet}{
If the number of heterozygous positions in a segment is greater than this, randomly select that many.
}
  \item{control}{
Not used, placeholder for other optimization function
}
  \item{tumourrangedata}{
Object of class IRanges containing tumour data, obtained from a call to gcCorrect}
  \item{maskmap}{
Heterozygous positions falling in bins with mappability less than this are ignored.
}
}
\details{
Semantically speaking, if p is the proportion of copies that are of, say, paternal origin, the function returns an estimate of either p or 1-p.  In other words, the actual proportion of copies that are derived from the paternal chromosome is unknown.  Moreover, these proportions are not comparable between segments (the parental origin "attached" to a proportion can change across segments). 
}
\value{
The seg data.frame is returned, with an extra column named "p" containing the estimates.
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
function (seg, ar, minhet = 50, maxhet = 5000, control = list(maxit = 1000), 
    tumourrangedata = NULL, maskmap = 0.8) 
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
    dz <- function(z, C, p) {
        return(dpois(z, C * p) + dpois(z, C * (1 - p)) - ppois(z, 
            C * p) * ppois(z, C * (1 - p)) + ppois(z - 1, C * 
            p) * ppois(z - 1, C * (1 - p)))
    }
    zlik <- function(pr) {
        p <- pr
        return(-sum(log(dz(x, n, p))))
    }
    dz <- function(z, C, p) {
        return(dpois(z, C * p) + dpois(z, C * (1 - p)) - ppois(z, 
            C * p) * ppois(z, C * (1 - p)) + ppois(z - 1, C * 
            p) * ppois(z - 1, C * (1 - p)))
    }
    zlik <- function(pr, z, C) {
        p <- pr
        return(-sum(log(dz(z, C, p))))
    }
    tmpseg <- seg
    tmpseg$p <- NA
    for (i in 1:nrow(tmpseg)) {
        cat("#########", i, "\n")
        chr <- as.character(tmpseg$chrom[i])
        from <- tmpseg$start.pos[i]
        to <- tmpseg$end.pos[i]
        sel <- ar.tmp[, 1] == chr & ar.tmp[, 2] >= from & ar.tmp[, 
            2] <= to
        if (sum(sel) >= minhet) {
            x <- ar.tmp[sel, 3]
            y <- ar.tmp[sel, 4]
            bpx <- boxplot(x, plot = F)
            bpy <- boxplot(y, plot = F)
            rem <- is.element(x, bpx$out) | is.element(y, bpy$out) | 
                x + y == 0
            x <- x[!rem]
            y <- y[!rem]
            n <- x + y
            sel <- sample(1:length(x), min(length(x), maxhet), 
                replace = F)
            x <- x[sel]
            n <- n[sel]
            x <- apply(data.frame(x, n - x), 1, min)
            pr <- max(0.005, min(tmpseg$meanar[i], 1 - tmpseg$meanar[i]))
            cat(i, " ", pr, "\n")
            bestmo <- tryCatch(optimize(f = zlik, interval = c(0, 
                0.5), z = x, C = n), error = function(e) {
            })
            if (!is.null(bestmo)) {
                print((bestmo))
                est <- bestmo$minimum
                tmpseg$p[i] <- est
            }
        }
        print(tmpseg[i, ])
        if (0) {
            plot(tmpseg[1:i, c("meanar", "p")], xlim = c(0, 1), 
                ylim = c(0, 1))
            abline(0, 1)
            abline(1, -1)
            abline(h = 0.5)
            abline(v = 0.5)
        }
    }
    return(tmpseg)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
