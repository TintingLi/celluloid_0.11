\name{estimateLOHcurve}
\alias{estimateLOHcurve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
estimateLOHcurve(segments, cntr = NULL, plot = F, manual = F, minsize = 1e+07, minmap = 0.9, maxmean = 1.25)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segments}{
Segment data frame, such as t.ar.seg. Expecting to find columns "mean" and "p".
}
  \item{cntr}{
Object returned by showTumourProfile. Used if plot = T or manual = T.
}
  \item{plot}{
Display the tumour contour plot and the LOH curve.
}
  \item{manual}{
Manually select segments to enter curve estimation.
}
  \item{minsize}{
Minimum size of segments to be used to estimate the LOH curve.
}
  \item{minmap}{
Minumum mean mappability of segments to be used to estimate the LOH curve.
}
  \item{maxmean}{
Maximum mean value (x-axis) of segments to be used to estimate the LOH curve.
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
function (segments, cntr = NULL, plot = F, manual = F, minsize = 1e+07, 
    minmap = 0.9, maxmean = 1.25) 
{
    ARloh <- function(x, S, n) {
        k <- 2 * (x - S * n)/(S - S * n)
        return(n/(2 * n + (1 - n) * k))
    }
    sel <- segments$size > minsize & !segments$mask & segments$meanmap > 
        minmap & segments$mean < maxmean & !is.na(segments$p)
    subseg <- segments[sel, ]
    if (!manual) {
        f <- function(n) {
            x <- subseg$mean
            ear <- ARloh(x, 1, n)
            diff <- subseg$p - ear
            diff[diff > 0.025] <- 0
            if (all(diff == 0) | sum(diff != 0) < 4) {
                return(sum(subseg$p - ear)^2)
            }
            else {
                return(sum(diff^2))
            }
        }
        Sn <- optimize(f, c(0, 1))$minimum
        if (plot) {
            if (is.null(cntr)) {
                stop("estimateLOHcurve: missing cntr argument")
            }
            image(cntr, col = terrain.colors(50))
            contour(cntr, nlev = 50, add = T)
            le <- subseg$end.pos - subseg$start.pos
            cxcut <- as.integer(cut(le, c(1e+05, 1e+06, 5e+06, 
                1e+07, 2e+07, 5e+07, Inf)))/3
            points(x <- subseg$mean, y <- subseg$p, pch = 21, 
                col = "blue", lwd = 3, cex = cxcut)
            points(subseg$mean, subseg$p, pch = 19, col = "white", 
                cex = cxcut - 0.5)
            points(subseg$mean, 1 - subseg$p, pch = 21, col = "blue", 
                lwd = 3, cex = cxcut)
            points(subseg$mean, 1 - subseg$p, pch = 19, col = "white", 
                cex = cxcut - 0.5)
            x <- seq(Sn, 2, 0.01)
            points(x, ARloh(x, 1, Sn), type = "l", col = "black", 
                lwd = 3)
            points(x, 1 - ARloh(x, 1, Sn), type = "l", col = "black", 
                lwd = 3)
        }
        return(Sn)
    }
    else {
        if (is.null(cntr)) {
            stop("estimateLOHcurve: missing cntr argument")
        }
        image(cntr, col = terrain.colors(50))
        contour(cntr, nlev = 50, add = T)
        le <- subseg$end.pos - subseg$start.pos
        cxcut <- as.integer(cut(le, c(1e+05, 1e+06, 5e+06, 1e+07, 
            2e+07, 5e+07, Inf)))/3
        points(x <- subseg$mean, y <- subseg$p, pch = 21, col = "blue", 
            lwd = 3, cex = cxcut)
        points(subseg$mean, subseg$p, pch = 19, col = "white", 
            cex = cxcut - 0.5)
        points(subseg$mean, 1 - subseg$p, pch = 21, col = "blue", 
            lwd = 3, cex = cxcut)
        points(subseg$mean, 1 - subseg$p, pch = 19, col = "white", 
            cex = cxcut - 0.5)
        cat("Select segments on the LOH curve; right click to finish\n")
        selectedPoints <- c()
        zz <- subseg[, c("mean", "p")]
        id <- identify(zz[, 1], zz[, 2], n = 1, labels = "")
        while (length(id) > 0) {
            points(zz[id, ], pch = 19, col = "black")
            selectedPoints <- rbind(selectedPoints, data.frame(x = zz[id, 
                1], y = zz[id, 2]))
            id <- identify(zz[, 1], zz[, 2], n = 1, labels = "")
        }
        x <- selectedPoints[, 1]
        y <- selectedPoints[, 2]
        nnllss <- nls(y ~ ARloh(x, 1, n), start = list(n = 0.01), 
            upper = list(n = 1), lower = list(n = 0), algo = "port")
        Sn <- summary(nnllss)$coefficients[1, 1]
        x <- seq(Sn, 2, 0.01)
        points(x, ARloh(x, 1, Sn), type = "l", col = "black", 
            lwd = 3)
        points(x, 1 - ARloh(x, 1, Sn), type = "l", col = "black", 
            lwd = 3)
        return(Sn)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
