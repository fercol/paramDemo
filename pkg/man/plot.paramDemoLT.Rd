\name{plot.paramDemoLT}
\alias{plot.paramDemoLT}
\alias{plot.paramDemoPLE}

\title{
 Functions to plot life tables and product limit estimators
}

\description{
 Draw a plot of demographic rates from a life table produced with function \code{CalcLifeTable} or of product limit estimator produced with function \code{CalcProductLimitEst}
}

\usage{
\method{plot}{paramDemoLT}(x, demorate = "lx", inclCIs = FALSE, ...)

\method{plot}{paramDemoPLE}(x, inclCIs = FALSE, ...)
}

\arguments{
	\item{x }{Object of class \code{paramDemoLT} or \code{paramDemoPLE} produced with functions \code{CalcLifeTable} or \code{CalcProductLimitEst}}

	\item{demorate }{Demographic rate to be plotted, choices are \dQuote{\code{lx}} for the survival, \dQuote{\code{px}} for the age-specific survival probability, \dQuote{\code{qx}} for the age-specific mortality probability, and \dQuote{\code{ex}} for the remaining life expectancy}

  \item{inclCIs }{Logical indicating whether confidence intervals should be ploted in case they were calculated with function \code{CalcLifeTable} or \code{CalcProductLimitEst}}
  
	\item{... }{Additional arguments passed to function \code{plot}}

}

\value{No return value, called for plotting objects of class \dQuote{\code{paramDemoLT}}}

\author{Fernando Colchero \email{fernando_colchero@eva.mpg.de}}

\seealso{
\code{\link{CalcLifeTable}} to calculate life tables. \code{\link{CalcProductLimitEst}} to calculate product limit estimators.
}

\examples{
# Simulate age at death data from Gompertz model:
ages <- SampleRandAge(n = 100, theta = c(b0 = -5, b1 = 0.1))

# Calculate life table:
lt <- CalcLifeTable(ageLast = ages, departType = rep("D", 100))

# Plot life table:
plot(lt)
}

\keyword{models}

