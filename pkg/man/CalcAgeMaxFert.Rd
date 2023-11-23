\name{CalcAgeMaxFert}
\alias{CalcAgeMaxFert}

\title{
 Calculating Age at Maximum Fertility and Maximum Fertility.
}

\description{
 \code{CalcAgeMaxFert} is used to calculate the age at maximum fertility, and the maximum fertility level from parametric models of fertility.
}

\usage{CalcAgeMaxFert (beta, modelFert = "quadratic", ageMatur = 0, 
                           maxAge = 100)}

\arguments{
	\item{beta }{Numerical vector of age-specific fertility parameters (see details).}

	\item{modelFert }{Age-specific fertility model. Options are \dQuote{\code{quadratic}} (default), \dQuote{\code{PeristeraKostaki}}, \dQuote{\code{ColcheroMuller}}, \dQuote{\code{Hadwiger}}, \dQuote{\code{gamma}}, \dQuote{\code{beta}}, \dQuote{\code{skewNormal}}, \dQuote{\code{gammaMixture}}, \dQuote{\code{HadwigerMixture}}, \dQuote{\code{skewSymmetric}}, \dQuote{\code{skewLogistic}} (see details)}

	\item{ageMatur }{Numerical value for the age at sexual maturity.}

	\item{maxAge }{Numerical value for the maximum possible age of reproduction.}

}

\details{
For a given function \eqn{g} of age-specific fertility, \code{CalcAgeMaxFert} uses either an analytical solution to 
\deqn{
\frac{dg}{dx} = 0
}

or a numerical approximation. For details of the parametric models available see \code{\link{CalcFert}}.
}

\value{
 The function outputs a vector with the calculated or estimated age at maximum fertility, the maximum level of fertility, the error in the estimation and the number of iterations (for numerical approximations),  and the age at maturity as specified by the user. 
}

\author{Fernando Colchero \email{fernando_colchero@eva.mpg.de}}

\seealso{
\code{\link{CalcSurv}} to calculate age-specific survival, \code{\link{CalcMort}} to calculate age-specific mortality, \code{\link{CalcFert}} to calculate age-specific fertility.

\code{\link{CalcAgeingRateMort}} to calculate ageing rates from parametric models of age-specific mortality. \code{\link{CalcRemainLifeExp}} to calculate remaining life expectancy from parametric models of age-specific mortality.
}

\examples{
# Calculate age at maximum fertility from quadratic model:
maxg <- CalcAgeMaxFert(beta = c(b0 = 0.5, b1 = 0.01, b2 = 10))

# Calculate age at maximum fertility from gamma model:
maxg <- CalcAgeMaxFert(beta = c(b0 = 13, b1 = 2, b2 = 0.15),
                     modelFert = "gamma")

}

\keyword{models}
