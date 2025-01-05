\name{CalcDiscrDemo}
\alias{CalcDiscrDemo}

\title{
 Calculating Demographic Rates on Discrete Age Intervals from Parametric Age-Specific Mortality
}

\description{
   \code{CalcDiscrDemo} is used to calculate survival and fertility functions (i.e., \code{lx}, \code{px}, and \code{qx}) on discrete age intervals from an object of class \dQuote{\code{paramDemo}} created with function \code{CalcDemo}.
}

\usage{CalcDiscrDemo(demo, dx = 1)}

\arguments{
	\item{demo }{Object of class \dQuote{\code{paramDemo}} created with function \code{CalcDemo}}

	\item{dx }{Numeric value for the age increments. Default is \code{dx = 1}}

}

\details{
\code{CalcDiscrDemo} takes the continuous survival functions from an object of class \dQuote{\code{paramDemo}} and calculates the basic survival probabilities and, if available, fertility rates on discrete age intervals. 

Let \eqn{x \geq 0} and \eqn{\omega \in \mathbb{N}} be the next integer of the age when the survival function \eqn{S(x)} reaches the minimum. Given the age increments, \eqn{\Delta x}, specified with argument \code{dx}, the function creates a partition of \eqn{[0, \omega]} such that \eqn{[x_i, x_i + \Delta x) \subset [0, \omega]} for \eqn{i = 1, \dots, n}. At each age interval it calculates the discrete survival, \code{lx}, as 
\deqn{
l_{x_i} = S(x_i).
} 

The age- or stage-specific survival probability, \code{px}, is calculated as 
\deqn{
p_{x_i} = \frac{S(x_i + \Delta x)}{S(x_i)},
}
while the age- or stage-specific mortality probability, \code{qx}, is \eqn{q_{x_i} = 1 - p_{x_i}}.

If fertility was also calculated with function \code{CalcDemo}, then the function includes the discrete age- or stage-specific fertility, \code{bx}, as
\deqn{
b_{x_i} = b\left(x_i + \frac{\Delta x}{2}\right) \Delta x,
}
where \eqn{b(x)} is the continuous fertility function.

}

\value{
\code{CalcDiscrDemo} returns an object of class \dQuote{\code{discrDemo}} that consist of a matrix with the following columns:

	\item{age }{Partition of the full age subset into discrete intervals}
	
	\item{lx }{Survival or cumulative survival}
		
	\item{px }{The age- or stage-specific survival probability}
	
	\item{qx }{The age- or stage-specific mortality probability}
	
	\item{bx }{If available, the age- or stage-specific fertility}

}

\author{Fernando Colchero \email{fernando_colchero@eva.mpg.de}}

\seealso{
\code{\link{CalcDemo}} to create an object of class \dQuote{\code{paramDemo}}.
}

\examples{
# Create paramDemo object from Gompertz mortality and
# quadratic fertility:
dem <- CalcDemo(theta = c(b0 = -5, b1 = 0.1), 
                 beta = c(b0 = 0.5, b1 = 0.01, b2 = 10),
                 summarStats = TRUE, agesAR = c(5, 10))

# Create discrete demographic object:
demDisc <- CalcDiscrDemo(dem, dx = 1)

}

\keyword{models}

