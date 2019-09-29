\name{StabCont }
\alias{StabCont }
\title{parameter stability test for continuous partitioning variable}
\description{Performs parameter stability test (Kundu and Harezlak, 2019) 
with continuous partitioning variable to determine whether the parameters 
of linear mixed effects model remains same across all distinct values of 
given continuous partitioning variable.
}
\usage{
StabCont(data, patid, fixed, splitvar)
}
\arguments{
\item{data}{name of the dataset. It must contain variable specified for \code{patid} (indicating subject id) and all the variables specified in the formula and the StabCont(data, fixed, splitvar)partitioning variable of interest specified in \code{splitvar}.}
\item{patid}{name of the subject id variable.}
\item{fixed}{a two-sided linear formula object describing the fixed-effects part of the model, with the response on the left of a \code{~ }operator and the terms, separated by \code{+} operators, on the right. Model with \code{-1} to the end of right side indicates no intercept. For model with no fixed effect beyond intercept, please specify only \code{1} right to the \code{~} operator.}
\item{splitvar}{the continuous partitioning variable of interest. It's value should not change over time.}
}
\details{
The continuous partitioning variable of interest. It's value should not change over time.

\code{Y_i(t)= W_i(t)theta + b_i + epsilon_{it}}

where \code{W_i(t)} is the design matrix, \code{theta} is the parameter associated with 
\code{W_i(t)} and \code{b_i} is the random intercept. Also, \code{epsilon_{it} ~  N(0,sigma ^2)} 
and \code{b_i ~ N(0, sigma_u^2)}. Let \code{X} be the baseline continuous partitioning 
variable of interest. \code{StabCont()} performs the following omnibus test

\code{H_0:theta_{(g)}=theta_0} vs. \code{H_1: theta_{(g)} ^= theta_0},     for all g

where, \code{theta_{(g)}} is the true value of \code{theta} for subjects with \code{X=C_g} 
where \code{C_g} is the any value realized by \code{X}.
}
\value{
\item{p}{It returns the p-value for parameter instability test}
}
\author{ Madan Gopal Kundu \email{madan_g.kundu@yahoo.com}}

\references{
Kundu, M. G., and Harezlak, J. (2019). Regression trees for longitudinal data with baseline covariates. Biostatistics & Epidemiology, 3(1):1-22.
}
\seealso{
\code{\link{StabCont}}, \code{\link{LongCART}}, \code{\link{plotLongCART}}, \code{\link{textLongCART}}
}

\examples{
#--- Get the data
data(ACTG175)
               
#--- Run StabCont()                
out<- StabCont(data=ACTG175, patid="pidnum", fixed=cd4~time, splitvar="age")
out$pval
}