\name{plot.PDlifeHist}
\alias{plot.PDlifeHist}
\alias{summary.PDlifeHist}
\alias{print.PDlifeHist}

\title{
Summarizing and plotting life history variables from parametric demographic models.
}
\description{
These functions are all generic methods for class \code{PDlifeHist}.
}
\usage{
\method{plot}{PDlifeHist}(x, type = "rates", ...)

\method{summary}{PDlifeHist}(object, ...)

\method{print}{PDlifeHist}(x, ...)
}

\arguments{
  \item{x }{Object of class \dQuote{\code{PDlifeHist}} obtained from function \code{\link{CalcLifeHist}}}
  
  \item{object }{Object of class \dQuote{\code{PDlifeHist}} obtained from function \code{\link{CalcLifeHist}}}

  \item{type }{Character string specifying the type of plot to be drawn, options are \dQuote{\code{rates}}, \dQuote{\code{sensitivities}}, \dQuote{\code{all}}}
  
  \item{... }{Additional arguments passed to functions \code{plot}, \code{print}, and \code{summary}}
}

\details{
For objects of class \dQuote{\code{PDlifeHist}}, function \code{print} displays on the console the data frame of life history variables produced by function \code{\link{CalcLifeHist}}. 

Function \code{summary} prints to the console, first the settings as specified by the user (i.e., \code{theta} and \code{beta} parameters, and the types of models used), and then, as function \code{print}, the table of life history variables. 

Function \code{plot} with argument \code{type} = \dQuote{\code{rates}} produces plots of the age-specific demographic rates (i.e., \eqn{p_x} and \eqn{m_x}) and the stable age structure and reproductive values (i.e., \eqn{\omega_x} and \eqn{\nu_x}). If argument \code{type} = \dQuote{\code{sensitivity}}, then it produces plots of the sensitivities and elasticities of \eqn{\lambda} to \eqn{p_x} and \eqn{m_x}. 
}

\value{No return value, called for plotting objects of class \dQuote{\code{PDlifeHist}}}

\author{Fernando Colchero \email{fernando_colchero@eva.mpg.de}}

\seealso{
\code{\link{CalcLifeHist}} to calculate life history variables from parametric demographic models.
}
\examples{
# Calculate life histories on the defalult models:
test <- CalcLifeHist(theta = c(b0 = -3, b1 = 0.1), 
                     beta = c(b0 = 1, b1 = 0.01, b2 = 5),
                     ageMatur = 2)
                     
# Print results:
test

# Summary:
summary(test)

# Plot results:
plot(test, type = "all")

}
\keyword{hplot}
