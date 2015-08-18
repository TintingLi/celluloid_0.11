\name{plotSegment}
\alias{plotSegment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot segments 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotSegment(tumourrangedata, segments, ar = NULL, n.rc.seg = NULL, columns = NULL, maskmap = NULL, 
            file = "Rplot\%03d.pdf", title = NULL, chr = NULL, perpage = 4, layoutmat = NULL, 
            width = 8.5, height = 11 , ylim = c(-0.5, 8), normal = F, tlwd = 5, tlty = 1, 
            tcol = NULL, nlwd = 3, nlty = 3, ncol = gray(0.5), annotation=NULL, cex.annotation=1,... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tumourrangedata}{
Object of class IRanges containing tumour data, obtained from a call to gcCorrect
}
  \item{segments}{
segments data.frame obtained from a call to annotateSegments
}
  \item{ar}{
%%     ~~Describe \code{ar} here~~
}
  \item{n.rc.seg}{
%%     ~~Describe \code{n.rc.seg} here~~
}
  \item{columns}{
%%     ~~Describe \code{columns} here~~
}
  \item{maskmap}{
%%     ~~Describe \code{maskmap} here~~
}
  \item{file}{
%%     ~~Describe \code{file} here~~
}
  \item{title}{
%%     ~~Describe \code{title} here~~
}
  \item{chr}{
%%     ~~Describe \code{chr} here~~
}
  \item{perpage}{
%%     ~~Describe \code{perpage} here~~
}
  \item{layoutmat}{
%%     ~~Describe \code{layoutmat} here~~
}
  \item{width}{
%%     ~~Describe \code{width} here~~
}
  \item{height}{
%%     ~~Describe \code{height} here~~
}
  \item{ylim}{
%%     ~~Describe \code{ylim} here~~
}
  \item{normal}{
%%     ~~Describe \code{normal} here~~
}
  \item{tlwd}{
%%     ~~Describe \code{tlwd} here~~
}
  \item{tlty}{
%%     ~~Describe \code{tlty} here~~
}
  \item{tcol}{
%%     ~~Describe \code{tcol} here~~
}
  \item{nlwd}{
%%     ~~Describe \code{nlwd} here~~
}
  \item{nlty}{
%%     ~~Describe \code{nlty} here~~
}
  \item{ncol}{
%%     ~~Describe \code{ncol} here~~
}
\item{annotation}{Data frame containing: chr position label.  Labels are added to axis(1) of the segment plots
}
\item{cex.annotation}{ cex.axis for the annotation }
\item{...}{Additional arguments passed to pdf}
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
%%  ~~who you are~~
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
function (tc, segments, ar = NULL, n.rc.seg = NULL, columns = NULL, 
    maskmap = 0.8, file = "Rplot\%03d.png", title = NULL, chr = NULL, 
    perpage = 4, layoutmat = NULL, width = 2 * 480, height = 11 * 
        480/4, ylim = c(-0.5, 8), normal = F, tlwd = 5, tlty = 1, 
    tcol = NULL, nlwd = 3, nlty = 3, ncol = gray(0.5)) 
{
    if (!is.null(file)) 
        png(filename = file, width = width, height = height)
    if (is.null(columns)) {
        columns = c("icopy", "imean")
    }
    if (length(columns == 2)) {
        tccolumn <- colnames(tc) == columns[1]
        segcolumn <- colnames(segments) == columns[2]
    }
    else {
        stop("columns argument must be of length 2\n")
    }
    if (sum(tccolumn) != 1) {
        stop("could not find column named", columns[1], "in argument tc\n")
    }
    if (sum(segcolumn) != 1) {
        stop("could not find column named", columns[2], "in argument segments\n")
    }
    chrlist <- c(paste("chr", 1:22, sep = ""), "chrX", "chrY")
    if (!is.null(chr)) {
        chrlist <- chr
    }
    if (is.null(layoutmat)) {
        tmp <- c()
        for (i in 1:perpage) {
            tmp <- c(tmp, c(2 * i - 1, 2 * i - 1, 2 * i - 1, 
                2 * i))
        }
        layout(matrix(tmp, ncol = 1))
    }
    else {
        layout(layoutmat)
    }
    if (is.null(tcol)) {
        rb <- rainbow(24)
        ra <- c(rb[1], rb[4], "black", rb[6], rb[9], rb[13], 
            rb[15], rb[17], rb[19], rb[20])
    }
    else if (length(tcol) == 12) {
        ra <- tcol
    }
    else if (length(tcol) == 1) {
        ra <- rep(tcol, 12)
    }
    else {
        stop("tcol has to be NULL or have length of 1 or 12")
    }
    segcol <- ra[as.integer(cut(segments[, segcolumn], br = c(-Inf, 
        seq(0.5, 8.5, 1), Inf)))]
    chrlist <- intersect(chrlist, unique(segments$chrom))
    for (chr in chrlist) {
        cat("plotting", chr, "\n")
        sel <- tc$space == chr
        x <- (start(tc[sel, ]) + end(tc[sel, ]))/2
        XLIM <- c(0, max(x))
        y <- tc[sel, tccolumn][[1]]
        if (!is.null(maskmap)) {
            mask <- tc$map[sel] < maskmap
            x <- x[!mask]
            y <- y[!mask]
        }
        par(mar = c(0.5, 5, 6, 1))
        plot.new()
        plot.window(ylim = ylim, xlim = XLIM)
        title(ylab = "copy number", xlab = "")
        box()
        axis(2)
        main <- chr
        if (!is.null(title)) {
            main <- paste(title, chr, sep = "/")
        }
        title(main = main, line = 4.5)
        abline(h = 0:ceiling(max(ylim)), lty = 1, col = "black")
        sel <- segments$chrom == chr
        if (sum(sel) > 0) {
            subseg <- segments[sel, ]
            for (s in 1:nrow(subseg)) {
                selx <- x >= subseg[s, "start.pos"] & x <= subseg[s, 
                  "end.pos"]
                points(x[selx], y[selx], pch = ".", col = makeTransparent(segcol[sel][s], 
                  alpha = 50))
                shadow <- "black"
                if (segcol[sel][s] == "black") {
                  shadow <- "white"
                }
                if (!is.null(maskmap)) 
                  if (subseg[s, "meanmap"] > maskmap) {
                    lines(c(subseg[s, "start.pos"], subseg[s, 
                      "end.pos"]), c(subseg[s, segcolumn], subseg[s, 
                      segcolumn]), lwd = tlwd + 1, lty = tlty, 
                      col = shadow)
                    lines(c(subseg[s, "start.pos"], subseg[s, 
                      "end.pos"]), c(subseg[s, segcolumn], subseg[s, 
                      segcolumn]), lwd = tlwd, lty = tlty, col = segcol[sel][s])
                  }
                abline(v = c(subseg[s, "start.pos"], subseg[s, 
                  "end.pos"]), lty = 3)
                mid <- (subseg[s, "start.pos"] + subseg[s, "end.pos"])/2
                le <- (subseg[s, "end.pos"] - subseg[s, "start.pos"])
                if (any(colnames(segments) == "labels") & le > 
                  1e+06) {
                  axis(3, at = mid, label = subseg[s, "labels"], 
                    las = 2, cex.axis = 1)
                }
            }
        }
        if (!is.null(n.rc.seg)) {
            sel <- n.rc.segments$chrom == chr
            if (sum(sel) > 0) {
                subsegn <- n.rc.seg[sel, ]
                for (s in 1:nrow(subsegn)) {
                  lines(c(subsegn[s, "start.pos"], subsegn[s, 
                    "end.pos"]), 2 * c(subsegn[s, "mean"], subsegn[s, 
                    "mean"]), lwd = nlwd, col = ncol, lty = nlty)
                }
            }
        }
        if (!is.null(ar)) {
            sel <- ar$CHR == chr
            if (sum(sel) > 0) {
                subar <- ar[sel, ]
                subar$ar <- subar[, 3]/(subar[, 3] + subar[, 
                  4])
                hh <- hist2d(subar$POS, subar$ar, nbins = c(floor(nrow(subar)/50), 
                  20), show = F)
                par(mar = c(3, 5, 0.5, 1), xaxs = "i", yaxs = "i")
                plot.new()
                plot.window(xlim = XLIM, ylim = c(0, 1))
                title(ylab = "AR")
                for (s in 1:nrow(subseg)) {
                  from <- subseg[s, 4]
                  to <- subseg[s, 5]
                  sel <- hh$x >= from & hh$x <= to
                  if (sum(sel) > 1) {
                    if (sum(hh$counts[sel, ] != 0) > 0) 
                      image(hh$x[sel], hh$y, log(hh$counts[sel, 
                        ]), col = gray(1 - seq(0, 1, 0.05)), 
                        add = T)
                  }
                }
                axis(1, at = seq(0, 3e+08, 1e+07), labels = as.character(seq(0, 
                  3e+08, 1e+07)/1000))
                axis(2, at = c(0, 0.5, 1))
                box()
            }
            else {
                plot.new()
                plot.window(xlim = XLIM, ylim = c(0, 1))
                axis(1, at = seq(0, 3e+08, 1e+07), labels = as.character(seq(0, 
                  3e+08, 1e+07)/1000))
                axis(2, at = c(0, 0.5, 1))
                box()
            }
        }
        else {
            plot.new()
            plot.window(xlim = XLIM, ylim = c(0, 1))
            axis(1, at = seq(0, 3e+08, 1e+07), labels = as.character(seq(0, 
                3e+08, 1e+07)/1000))
            axis(2, at = c(0, 0.5, 1))
            box()
        }
    }
    if (!is.null(file)) 
        dev.off()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
