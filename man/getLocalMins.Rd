\name{getLocalMins}
\alias{getLocalMins}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prune output of coverParamSpace 
}
\description{
Takes results from coverParamSpace and prune out solutions that are close one another, with pretty print. 
}
\usage{
getLocalMins(li, pruneS = 0.05, pruneT = 0.025, ploidy = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{li}{
List returned from coverParamSpace
}
  \item{pruneS}{
Skip the solution if its S is within +/- pruneS of the S of a better solution 
}
  \item{pruneT}{
Skip the solution if its t parameter vector is within +/- pruneT of the t of a better solution 
}
  \item{ploidy}{
Vector of length 2. Only output solutions consistent with a tumour ploidy above or equal to ploidy[1] and below or equal to ploidy[2].
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
function (li, pruneS = 0.05, pruneT = 0.025, ploidy = NULL) 
{
    allParamSpace <<- c()
    if (DEBUG) {
        print(li)
    }
    out <- data.frame()
    for (i in 1:length(li)) {
        tmp <- data.frame(t(c(li[[i]]$value, li[[i]]$par, li[[i]]$subset)))
        names(tmp) <- paste("X", 1:ncol(tmp), sep = "")
        out <- rbind(out, tmp)
        allParamSpace <<- rbind(allParamSpace, li[[i]]$paramSpace)
    }
    if (ncol(out) == 4) {
        out <- cbind(out[, 1:3], 1 - as.numeric(as.character(out[, 
            3])), out[, 4])
        names(out) <- c("value", "S", "N", "T1", "subset")
    }
    else {
        out <- cbind(out[, 1:(ncol(out) - 1)], 1 - apply(out[, 
            3:(ncol(out) - 1)], 1, function(x) {
            sum(as.numeric(x))
        }), out[, ncol(out)])
        names(out) <- c("value", "S", "N", paste("T", 1:(ncol(out) - 
            4), sep = ""), "subset")
    }
    for (i in 1:(ncol(out) - 1)) {
        out[, i] <- as.numeric(as.character(out[, i]))
    }
    if (!is.null(pruneS)) {
        dS <- as.matrix(dist(out[, 2]))
        dT <- as.matrix(dist(out[, 3:(ncol(out) - 2)]))
        rem <- rep(FALSE, nrow(dS))
        for (i in 1:(nrow(dS) - 1)) {
            rem <- rem | c(rep(FALSE, i), dS[i, (i + 1):ncol(dS)] < 
                pruneS & dT[i, (i + 1):ncol(dT)] < pruneT)
        }
        out <- out[!rem, ]
    }
    if (!is.null(ploidy)) {
        rem <- rep(FALSE, nrow(out))
        for (i in 1:nrow(out)) {
            tploi <- (2/out[i, "S"] - 2 * out[i, "N"])/(1 - out[i, 
                "N"])
            cat(tploi, ploidy, "\n")
            rem[i] <- tploi < ploidy[1] | tploi > ploidy[2]
        }
        out <- out[!rem, ]
    }
    return(out[order(out[, 1]), 1:(ncol(out) - 1)])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
