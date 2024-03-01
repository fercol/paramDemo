\name{CalcAgeingRateMort}
\alias{CalcAgeingRateMort}
\alias{CalcAgeingRateFert}

\title{
 Calculating actuarial and reproduction aging rates.
}

\description{
 \code{CalcAgeingRateMort} and \code{CalcAgeingRateFert} are used to calculate actuarial and reproductive ageing rates from parametric models of mortality and fertility.
}

\usage{
CalcAgeingRateMort(theta, x, model = "GO", shape = "simple", 
                   checkTheta = TRUE)

CalcAgeingRateFert(beta, x, modelFert = "quadratic", ageMatur = 0,
                   checkBeta = TRUE)
}

\arguments{
	\item{theta }{Numerical vector of age-specific mortality parameters (see details).}
	
	\item{beta }{Numerical vector of age-specific fertility parameters (see details.)}

	\item{x }{Numerical vector of ages at which to calculate mortality.}

	\item{model }{The underlying mortality model to be used. \code{"EX"} = exponential,\code{"GO"} = Gompertz, \code{"WE"} = Weibull and \code{"LO"} = logistic (see \code{details}).}

	\item{shape }{The overall shape of the model. Values are: \code{simple} = no extra parameters added; \code{Makeham} = a constant parameter is added to the mortality; and \code{bathtub} = a Gompertz declining mortality for early ages and a constant parameter are added to the mortality model (see \code{details}).}
	
		\item{modelFert }{Age-specific fertility model. Options are \dQuote{\code{quadratic}} (default), \dQuote{\code{PeristeraKostaki}}, \dQuote{\code{ColcheroMuller}}, \dQuote{\code{Hadwiger}}, \dQuote{\code{gamma}}, \dQuote{\code{beta}}, \dQuote{\code{skewNormal}}, \dQuote{\code{gammaMixture}}, \dQuote{\code{HadwigerMixture}}, \dQuote{\code{skewSymmetric}}, \dQuote{\code{skewLogistic}} (see details)}

	\item{ageMatur }{Numerical value for the age at sexual maturity.}

	\item{checkTheta }{Logical to verify that the \code{theta} parameters conform with the mortality model's specification (see details).}
	
	\item{checkBeta }{Logical to verify that the \code{beta} parameters conform with the fertility model's specification (see details).}

}

\details{
The function \code{CalcAgeingRateMort} uses parametric functions to calculate the actuarial (i.e., mortality) rate of ageing. The function follows the conventions from package BaSTA (Colchero and Clark 2012, Colchero \emph{et al.} 2012, Colchero \emph{et al.} 2021) to select the parametric model of mortality. The mortality function describes how the risk of mortality changes with age, and is defined as 
\deqn{
\mu(x | \theta) = \lim_{\Delta x \rightarrow 0} \frac{\Pr[x < X < x + \Delta x | X > x]}{\Delta x},
}
where \eqn{X} is a random variable for ages at death, \eqn{x \geq 0} are ages and \eqn{\theta} is the vector of mortality parameters. (For further details on the mortality and survival models see \code{\link{CalcMort}}). 

Given a vector of ages \eqn{x_1, x_2, \dots, x_n} specified by the user with argument \code{x}, the function calculates ageing rates at age \eqn{x_i} as
\deqn{
\frac{d}{dx}\ln [\mu(x)] |_{x = x_i}, 
}
for \eqn{i = 1, 2, \dots, n}.

Similarly, function \code{CalcAgeingRatesFert} calculate reproductive ageing rates from parametric models of age-specific fertility, \eqn{g(x)}. It uses a numerical approximation to 
\deqn{
\frac{d}{dx}\ln [g(x)] |_{x = x_i}, 
}
for \eqn{i = 1, 2, \dots, n}.

}

\value{
 The functions output a matrix with the ages at which ageing rates were calculated, the estimated actuarial rate of ageing, and the level of survival at that age. 
}

\author{Fernando Colchero \email{fernando_colchero@eva.mpg.de}}

\seealso{
\code{\link{CalcSurv}} to calculate age-specific survival, \code{\link{CalcMort}} to calculate age-specific mortality, \code{\link{CalcFert}} to calculate age-specific fertility.

\code{\link{CalcAgeMaxFert}} to calculate the age at maximum fertility from parametric models of age-specific fertility. \code{\link{CalcRemainLifeExp}} to calculate remaining life expectancy from parametric models of age-specific mortality.
}

\examples{
# Calculate actuarial ageing rate from Gompertz model:
arm <- CalcAgeingRateMort(theta = c(b0 = -5, b1 = 0.1), x = 10)

# Calculate reproductive ageing rate from quadratic model:
arf <- CalcAgeingRateFert(beta = c(b0 = 2, b1 = 0.0025, b2 = 2), x = 10,
                          modelFert = "quadratic")

}

\keyword{models}

