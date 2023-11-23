\name{CalcProductLimitEst}
\alias{CalcProductLimitEst}

\title{
 Calculating the product limit estimator from Data.
}

\description{
 \code{CalcProductLimitEst} uses non-parametric methods to calculate the product limit estimator and confidence intervals.
}

\usage{CalcProductLimitEst(ageLast, ageFirst = NULL, departType,
                     calcCIs = FALSE, nboot = 1000, alpha = 0.05)}

\arguments{
	\item{ageLast }{Numerical vector with the ages at last detection (i.e., death and censoring) (see \code{details})}

	\item{ageFirst }{Numerical vector of ages at first detection (i.e., truncation). If \code{NULL} then all values are set to 0 (see \code{details}).}

	\item{departType }{Character string vector for the type of departure (i.e., last detection), with values \dQuote{\code{D}} for death and \dQuote{\code{C}} for censoring (see \code{details}).}

  \item{calcCIs }{Logical indicating whether confidence intervals should be calculated}
  
	\item{nboot }{Number of bootstrap iterations}

	\item{alpha }{Alpha level. Default is 0.05 for 95\% CIs}

}

\details{
\bold{1) Data structure:}

The function allows to calculate product limit estimator (Wann \emph{et al.} 1987) for data that includes the following types of records:
\itemize{ 
\item\bold{Uncensored}: individuals with known ages at death;

\item\bold{right-censored}: individuals last seen alive;

\item\bold{left-truncated}: individuals born before the start of the study and are truncated at the age of entry. 
}

The data required are the ages at last detection (i.e., \emph{uncensored} or \emph{right-censored}) passed through argument \dQuote{\code{ageLast}}, the type of departure via argument \dQuote{\code{departType}}, which takes two values, namely \dQuote{\code{D}} for death, and \dQuote{\code{C}} for censored (i.e., right-censored). 

In addition, if there is \emph{left-truncation}, it takes the ages at entry to the study by means of argument \dQuote{\code{ageFirst}}. If all the individuals were born during the study, the value of \dQuote{\code{ageFirst}} can be left as \code{NULL}, which will make them all equal to 0. 
 }

\value{
  \code{CalcProductLimitEst} returns an object of class \dQuote{\code{paramDemoPLE}} with consists of a data frame and, if indicated by the user, the confidence interval of the PLE. The data frame includes the following numerical columns:
  
	\item{Ages }{ages when individuals died}
	
	\item{ple }{product limit estimator}
	
	\item{Lower }{If indicated with argument \code{calcCIs}, the lower \eqn{\alpha / 2} confidence bound}
	
		\item{Lower }{If indicated with argument \code{calcCIs}, the upper \eqn{1 - \alpha / 2} confidence bound}

}

\author{Fernando Colchero \email{fernando_colchero@eva.mpg.de}}

\seealso{
\code{\link{CalcLifeTable}} to calculate life tables.
}

\examples{
# Simulate age at death data from Gompertz model:
ages <- SampleRandAge(n = 100, theta = c(b0 = -5, b1 = 0.1))

# Calculate life table:
lt <- CalcLifeTable(ageLast = ages, departType = rep("D", 100))

}

\keyword{survival}

