\name{plot.paramDemo}
\alias{plot.paramDemo}

\title{
 Functions to plot parametric demographic functions
}

\description{
 Draw a plot of demographic functions produced with function \code{CalcDemo} 
}

\usage{
\method{plot}{paramDemo}(x, demofun = "all", ...)}

\arguments{
	\item{x }{Object of class \dQuote{\code{paramDemo}}}

	\item{demofun }{Character string for the demographic function to be ploted}

	\item{... }{Additional arguments passed to function \code{plot}}

}

\value{No return value, called for plotting objects of class \dQuote{\code{paramDemo}}}

\author{Fernando Colchero \email{fernando_colchero@eva.mpg.de}}

\seealso{
\code{\link{CalcDemo}} to calculate parametric demographic functions. }

\examples{
# Create paramDemo object from Gompertz mortality and
# quadratic fertility:
dem <- CalcDemo(theta = c(b0 = -5, b1 = 0.1), 
                 beta = c(b0 = 0.5, b1 = 0.01, b2 = 10),
                 summarStats = TRUE, agesAR = c(5, 10))

# Plot demographic object:
plot(dem)

}

\keyword{models}

