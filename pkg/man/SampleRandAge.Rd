\name{SampleRandAge}
\alias{SampleRandAge}

\title{
 Sampling Ages at Death from a Parametric Mortality model
}

\description{
 \code{SampleRandAge} is used to randomly sample ages at death from a parametric mortality model 
}

\usage{SampleRandAge(n, theta, dx = 0.001, model = "GO", 
               shape = "simple",  minSx = 1e-04)}

\arguments{
	\item{n }{Number of individual ages at death to sample}

	\item{theta }{Numerical vector of age-specific mortality parameters (see details).}

	\item{dx }{Numeric value for the age increments (see details).}

	\item{model }{The underlying mortality model to be used. \code{"EX"} = exponential,\code{"GO"} = Gompertz, \code{"WE"} = Weibull and \code{"LO"} = logistic (see \code{details}).}

	\item{shape }{The overall shape of the model. Values are: \code{simple} = no extra parameters added; \code{Makeham} = a constant parameter is added to the mortality; and \code{bathtub} = a Gompertz declining mortality for early ages and a constant parameter are added to the mortality model (see \code{details}).}

	\item{minSx }{Minimum value for the survival function to set the maximum age (see details)}

}

\details{
\code{SampleRandAge} arguments \dQuote{\code{model}} and \dQuote{\code{shape}}, to calculate the survival function, \eqn{S(x)} for \eqn{x \geq 0}, at discrete age intervals \eqn{[x, x + \Delta x)} where \eqn{\Delta x} is specfied with argument \dQuote{\code{dx}}. It then calculates the CDF of ages at death, \eqn{F(x)}, and uses inverse sampling to sample \eqn{n} ages at death by random sampling from a uniform distribution bound as 

\deqn{
u_i \sim Unif[F(x_m), F(x_M)],
}
where \eqn{x_m = min(x)} and \eqn{x_M = max(x)}. It then finds the corresponding ages at death as \deqn{
x_i = F^{-1}(u_i).
}

In most cases, this is acheived numerically and thus the lower the value for argument \dQuote{\code{dx}} the higher the precision.
}

\value{
\code{SampleRandAge} returns a vector of class \dQuote{\code{numeric}} with randmly generated ages at death.
}

\author{Fernando Colchero \email{fernando_colchero@eva.mpg.de}}

\examples{
# Simulate age at death data from Gompertz model:
ages <- SampleRandAge(n = 100, theta = c(b0 = -5, b1 = 0.1))

# Simulate age at death data from Siler model:
ages <- SampleRandAge(n = 100, theta = c(a0 = 0, a1 = 1, c = 0.001, 
                      b0 = -5, b1 = 0.1), model = "GO", 
                      shape = "bathtub")

}

\keyword{datagen}