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
            file = "Rplot\%03d", device="pdf", main = NULL, chr = NULL, perpage = 4, layoutmat = NULL, 
            width = 8.5, height = 11 , ylim = c(-0.5, 8), tlwd = 5, tlty = 1, 
            tcol = NULL, nlwd = 3, nlty = 3, ncol = gray(0.5), alpha=50, cex.axis=1, cex.lab=1, 
            cex.main=1, annotation=NULL, cex.annotation=1,marCN=c(0.5,5,6,1), marAR=c(3,5,0.5,1) , ... )
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
The ar data frame}
  \item{n.rc.seg}{
Segments in normal tissues. If provided, these segments will also be plotted.
}
  \item{columns}{
Vector of size 2 indicating the names of the columns to plot from tumourrangedata and segments. If NULL, defaults to c("icopy", "imean" ).
}
  \item{maskmap}{
Only plot a segment if its mean mappability is above maskmap. Plot all segments if set to NULL.
}
  \item{file}{
File prefix; suffix is appended based on device. Set to NULL to display on screen.
}
\item{device}{One of "pdf" (default),"bmp","jpeg", "png","tiff".}
  \item{title}{
main title for the graph (as in title(main=)). Defaults to the chromosome number.
}
  \item{chr}{
Which chromosomes to plot. If NULL defaults to all chromosomes.
}
  \item{perpage}{
How many chromosomes to plot per page
}
  \item{layoutmat}{
See details.}
  \item{width}{
Argument passed to the graphic device to define the width of a page. The default value was set for a pdf device, in inch units (the default of pdf). 
For other devices, the user needs to define width and height in pixels, or else needs to define values for units and res (by adding arguments such as units="in", res=300; see e.g. ?png).  
}
  \item{height}{
Argument passed to the graphic device to define the height of a page. The default value was set for a pdf device, in inch units (the default of pdf). 
For other devices, the user needs to define width and height in pixels, or else needs to define values for units and res (by adding arguments such as units="in", res=300; see e.g. ?png).  
}
  \item{ylim}{
ylim in copy number graphs, same for all chromosomes.
}
  \item{tlwd}{
Line width of tumour segments
}
  \item{tlty}{
Line type of tumour segments
}
  \item{tcol}{
Color of tumour segments. See details. 
}
  \item{nlwd}{
Line width of segments in normal.
}
  \item{nlty}{
Line type of segments in normal.
}
  \item{ncol}{
Color of segments in normal. 
}
\item{alpha}{alpha argument passed to rgb to define transparency of the colors of the scaled read counts of the 1kb bins. Takes value between 0 and 255. Alpha value of 255 is opaque. Alpha value of 0 prevents points from being plotted. See Details. }
\item{cex.axis}{Magnification to be used for text annotation of the axis values. Passed to axis().}
\item{cex.lab}{Magnification to be used for text annotation of the axis labels. Passed to title().}
\item{cex.main}{Magnification to be used for text annotation of title (the main. Passed to title().}
\item{annotation}{Data frame containing: chr position label.  Labels are added to axis(1) of the segment plots
}
\item{cex.annotation}{ cex.axis for the annotation }
\item{marCN}{Margins for the copy-number component of the graph}
\item{marAR}{Margins for the allelic-ratio component of the graph}
\item{...}{Additional arguments passed to device. See ?pdf, ?png, etc.}
}
\details{
The default layout prints 4 chromosomes per page. For each chromosome, the graph displays 
two items: 1) the copy number segments on top and 2) allelic ratios on bottom.  
Underneath the copy number segments, the graph displays each individual 1kb bin.

The two items are plotted according to a 3:1 row height 
ratio by default. A full page of 4 uses the layout 

layout( matrix( c( 1,1,1,2,3,3,3,4,5,5,5,6,7,7,7,8 ), ncol=1 ) ) 

See ?layout for more details. 

tcol can take either a single color value or a vector of 10 colors. Colors are assigned based on segment height.
Segments with copy number in (-Inf,0.5] will be colored tcol[1], (0.5,1.5] will be colored tcol[2], ...,
up to segments with copy number in (8.5,Inf) that will be colored tcol[10].

When set to NULL, defaults to:
  rb<-rainbow(24)
  tcol<-c( rb[1], rb[4], "black", rb[6], rb[9], rb[13], rb[15], rb[17], rb[19] ,rb[20] )

The 1kb bins are plotted with semi-transparency by defaults (see ?rgb and the alpha argument). Not all
devices can plot transparencies, see help page for the device you want if you get a warning.  An alpha 
value of 255 is opaque.

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
