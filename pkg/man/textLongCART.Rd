\name{textLongCART}
\alias{textLongCART}
\title{Place text on LongCART tree}
\description{
Labels the current plot of the tree generated from LongCART object with text.
}
\usage{
textLongCART(x, splits = TRUE, all = FALSE, 
             use.n = FALSE, minlength = 1L, ...)
}
\arguments{
\item{x}{a fitted object of class \code{"LongCART"}, containing a linear mixed effects tree.}
\item{splits}{similar to plot.rpart; logical flag. If \code{TRUE} (default), 
then the splits in the tree are labeled with the criterion for the split.}
\item{all}{Logical. If TRUE, all nodes are labeled, otherwise just terminal nodes.}
\item{use.n}{Logical. If TRUE, adds \code{n} to label.}
\item{minlength}{the length to use for factor labels. A value of 1 causes them to be printed as 
'a', 'b', ..... Larger values use abbreviations of the label names. See the labels.rpart function 
for details.}
\item{...}{arguments to be passed to or from other methods.}
}

\author{ Madan Gopal Kundu \email{madan_g.kundu@yahoo.com}}

\references{
Kundu, M. G., and Harezlak, J. (2019). Regression trees for longitudinal data with baseline covariates. Biostatistics & Epidemiology, 3(1):1-22.
}
\seealso{
\code{\link{plotLongCART}}, \code{\link{LongCART}}, \code{\link{StabCat}}, \code{\link{StabCont}} 
}

\examples{
#--- Get the data
data(ACTG175)
                
#--- Run LongCART()
gvars=c("age", "gender", "wtkg", "hemo", "homo", "drugs", 
        "karnof", "oprior", "z30", "zprior", "race",
        "str2", "symptom", "treat", "offtrt")
tgvars=c(1, 0, 1, 0, 0, 0,
         1, 0, 0, 0, 0,
         0, 0, 0, 0)

\donttest{
out<- LongCART(data=ACTG175, patid="pidnum", fixed=cd4~time, 
               gvars=gvars, tgvars=tgvars, alpha=0.05, 
               minsplit=100, minbucket=50, coef.digits=2)}
               
#--- Plot tree
\donttest{
par(xpd = T)
plot(out, compress = T)
text(out, use.n = T)}
}