\name{CalcLifeTableCIs}
\alias{CalcLifeTableCIs}

\title{
 Non-Parametric Bootstrap to Calculate Life Table Confidence Intervals.
}

\description{
 \code{CalcLifeTableCIs} uses non-parametric bootstrap to calculate confidence intervals for a life table.
}

\usage{CalcLifeTableCIs (ageLast, ageFirst = NULL, departType, 
                  nboot = 2000,  alpha = 0.05, dx = 1)}

\arguments{
	\item{ageLast }{Numerical vector with the ages at last detection (i.e., death and censoring)}

	\item{ageFirst }{Numerical vector of ages at first detection (i.e., truncation). If \code{NULL} then all values are set to 0.}

	\item{departType }{Character string vector for the type of departure (i.e., last detection), with values \dQuote{\code{D}} for death and \dQuote{\code{C}} for censoring.}

	\item{nboot }{Number of bootstrap iterations}

	\item{alpha }{Alpha level. Default is 0.05 for 95\% CIs}

	\item{dx }{Age interval size, default set at 1.}

}

\details{
FILL UP 
}

\value{
	\item{fill up }{FILL UP}
}

\author{Fernando Colchero \email{colchero@imada.sdu.dk}}

\seealso{FILL UP}

\examples{}

\keyword{models}

