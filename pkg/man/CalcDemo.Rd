\name{CalcDemo}
\alias{CalcDemo}

\title{
 FILL UP
}

\description{
 FILL UP
}

\usage{CalcDemo(theta = NULL, beta = NULL, x = NULL, dx = NULL, 
          model = "GO",  shape = "simple", modelFert = "quadratic", 
          type = "both", minSx = 0.01,  summarStats = TRUE, 
          ageMatur = 0, maxAge = NULL, agesAR = NULL,  
          SxValsAR = NULL)}

\arguments{
	\item{theta }{Numerical vector of mortality parameters.}

	\item{beta }{Numerical vector of fertility parameters.}

	\item{x }{Numerical vector of ages. If \code{NULL} (default), ages are constructed from age 0 to the age when the survival function reaches 0.0001.}

	\item{dx }{Numerical value for age increments. If \code{NULL}, the default value is 0.01.}

	\item{model }{Mortality model, options are \dQuote{\code{EX}} = exponential (i.e., constant mortality with age), \dQuote{\code{GO}} = Gompertz (default), \dQuote{\code{WE}} = Weibull, or \dQuote{\code{LO}} = logistic.}

	\item{shape }{The overall shape of the mortality model. Values are: \dQuote{\code{simple}} = no extra parameters added; \dQuote{\code{Makeham}} = a constant parameter is added to the mortality; and \dQuote{\code{bathtub}} = a Gompertz declining mortality for early ages and a constant parameter are added to the mortality model (see \code{details}).}

	\item{modelFert }{Age-specific fertility model. Options are \dQuote{\code{quadratic}} (default), \dQuote{\code{PeristeraKostaki}}, \dQuote{\code{ColcheroMuller}}, \dQuote{\code{Hadwiger}}, \dQuote{\code{gamma}}, \dQuote{\code{beta}}, \dQuote{\code{skewNormal}}, \dQuote{\code{gammaMixture}}, \dQuote{\code{HadwigerMixture}}, \dQuote{\code{skewSymmetric}}, \dQuote{\code{skewLogistic}} (see details)}

	\item{type }{Character string specifying the demographic rates to be calculated. Options are \dQuote{\code{both}} (default), \dQuote{\code{survival}}, \dQuote{\code{fertility}}.}

	\item{minSx }{Numerical value specifying the lower bound of the survival function to find the maximum age. The default is 0.01.}

	\item{summarStats }{Logical specifying whether summary statistics such as life expectancy, lifespan equality and ageing rates at the specified ages should be calculated.}

	\item{ageMatur }{Numerical value specifying the age at sexual maturity, used as the lower bound for the age-specific fertility. Default is 0.}

	\item{maxAge }{Numerical values for the maximum age for survival and fertility demographic rates. If \code{NULL} (default), the maximum age is calculated when the survival function reaches \code{minSx}. }

	\item{agesAR }{Numerical vector of ages at which ageing rates should be calculated. If \code{NULL} (default) they are calculated at the ages when the survival function reaches 0.5, 0.2, and 0.05. }

	\item{ SxValsAR }{Numerical value for the levels of survival at which ageing rates should be calculated. If \code{NULL} (default) they are calculated at the ages when the survival function reaches 0.5, 0.2, and 0.05.}

}

\details{
FILL UP 
}

\value{
	\item{fill up }{FILL UP}
}

\author{Fernando Colchero \email{fernando_colchero@eva.mpg.de}}

\seealso{FILL UP}

\examples{FILL UP}

\keyword{FILL UP}

