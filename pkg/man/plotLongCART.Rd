\name{plotLongCART}
\alias{plotLongCART}
\title{Plot an LongCART Object}
\description{
Plots an LongCART object on the current graphics device.
}
\usage{
plotLongCART(x, uniform = FALSE, branch = 1, nspace,
     margin = 0, minbranch = 0.3, ...)
}
\arguments{
\item{x}{a fitted object of class \code{"LongCART"}, containing a linear mixed effects tree.}
\item{uniform}{similar to \code{plot.rpart}; if TRUE, uniform vertical spacing
 of the nodes is used; this may be less cluttered when fitting a large plot 
onto a page. The default is to use a non-uniform spacing proportional to the 
error in the fit.}
\item{branch}{similar to \code{plot.rpart}; 	
controls the shape of the branches from parent to child node. Any number from 
0 to 1 is allowed. A value of 1 gives square shouldered branches, a value of 
0 give V shaped branches, with other values being intermediate.}
\item{nspace}{similar to \code{plot.rpart}; the amount of extra space between a node with children 
and a leaf, as compared to the minimal space between leaves. Applies to compressed trees only. 
The default is the value of \code{branch}.}
\item{margin}{similar to \code{plot.rpart}; an extra fraction of white space to 
leave around the borders of the tree. (Long labels sometimes get cut off by the default computation). }
\item{minbranch}{similar to \code{plot.rpart}; set the minimum length for a branch to \code{minbranch} 
times the average branch length. This parameter is ignored if \code{uniform=TRUE}. Sometimes a split will 
give very little improvement, or even (in the classification case) no improvement at all. A tree with 
branch lengths strictly proportional to improvement leaves no room to squeeze in node labels.}
\item{...}{arguments to be passed to or from other methods.}
}
\details{
This function is a method for the generic function plot, for objects of class 
\code{LongCART}. The y-coordinate of the top node of the tree will always be 1.}
\value{
The coordinates of the nodes are returned as a list, with components \code{x}
 and \code{y}.
}
\author{ Madan Gopal Kundu \email{madan_g.kundu@yahoo.com}}

\references{
Kundu, M. G., and Harezlak, J. (2019). Regression trees for longitudinal data with baseline covariates. Biostatistics & Epidemiology, 3(1):1-22.
}
\seealso{
\code{\link{textLongCART}}, \code{\link{LongCART}}, \code{\link{StabCat}}, \code{\link{StabCont}} 
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