\name{CalcMort}
\alias{CalcMort}

\title{
 Calculating Parametric Age-Specific Mortality.
}

\description{
 \code{CalcMort} is used to calculate the age-specific mortality function or hazard rate from different parametric models.
}

\usage{CalcMort (theta, x, model = "GO", shape = "simple", 
          checkTheta = TRUE)}

\arguments{
	\item{theta }{Numerical vector of age-specific mortality parameters (see details).}

	\item{x }{Numerical vector of ages at which to calculate mortality.}

	\item{model }{The underlying mortality model to be used. \code{"EX"} = exponential,\code{"GO"} = Gompertz, \code{"WE"} = Weibull and \code{"LO"} = logistic (see \code{details}).}

	\item{shape }{The overall shape of the model. Values are: \code{simple} = no extra parameters added; \code{Makeham} = a constant parameter is added to the mortality; and \code{bathtub} = a Gompertz declining mortality for early ages and a constant parameter are added to the mortality model (see \code{details}).}

	\item{checkTheta }{Logical to verify that the \code{beta} parameters conform with the model's specification (see details).}

}

\details{
\code{1) AGE-SPECIFIC MORTALITY MODELS:}

The function \code{CalcMort} uses parametric functions to calculate age-specific mortality, following the conventions from package BaSTA (Colchero and Clark 2012, Colchero \emph{et al.} 2012, Colchero \emph{et al.} 2021). The mortality function describes how the risk of mortality changes with age, and is defined as \eqn{\mu(x | \theta)}, where \eqn{x} corresponds to age and \eqn{\theta} is the vector of mortality parameters. 

\code{1.1)}\code{Argument }\dQuote{\code{model}}\code{:}

The \code{model} argument allows the user to choose between four basic mortality functions, namely 

(a) \code{model = }\dQuote{\code{EX}}: The exponential model (Cox and Oakes 1974), with constant mortality with age, specified as

\deqn{\mu_b(x | b) = b,}

where \eqn{b > 0}. 

(b) \code{model = }\dQuote{\code{GO}}: The Gompertz mortality model (Gompertz 1925, Pletcher 1999), calculated as 

\deqn{\mu_b(x | b) = \exp(b_0 + b_1 x),}

where \eqn{-\infty < b_0, b_1 < \infty}. 

(c) \code{model = }\dQuote{\code{WE}}: The Weibull mortality model (Pinder III \emph{et al.} 1978) calculated as

\deqn{\mu_b(x | b) = b_0 b_1^{b_0} x^{b_0 -1},}

where \eqn{b_0, b_1 > 0}. 

(d) \code{model = }\dQuote{\code{LO}}: The logistic mortality model (Pletcher 1999), calculated as 

\deqn{\mu_b(x | b) = \frac{\exp(b_0 + b_1 x)}{1 + b_2 \frac{e^{b_0}}{b_1} \left(e^{b_1 x}-1\right)},}
where \eqn{b_0, b_1, b_2 > 0}.

\code{1.2)}\code{Argument }\dQuote{\code{shape}}\code{:}

The \code{shape} argument allows the user to extend these models in order to explore more complex mortality shapes. 

The default value is \dQuote{\code{simple}} which leaves the model as defined above, in other words
\deqn{
\mu(x | \theta) = \mu_b(x | \theta)
}

With value \dQuote{\code{Makeham}}, a constant is added to the mortality, making the model equal to

\deqn{
\mu(x | \theta) = c + \mu_b(x | \theta_1),
}
where \eqn{\theta = [c, \theta_1]}. 

With value \dQuote{\code{bathtub}}, concave shapes in mortality can be explored. This is achieved by adding a declining Gompertz term and a constant parameter to the basic mortality model, namely

\deqn{
\mu(x | \theta) = \exp(a_0 - a_1 x) + c + \mu_b(x | \theta_1),
}
where \eqn{a_0 \in \mathbb{R}}, \eqn{a_1, c \geq 0} and \eqn{\theta_1 \in \theta} are specified based on argument \code{model}.

\code{2) SPECIFYING THETA PARAMETERS:}

Argument \code{theta} requires a numerical vector of parameters. For instance, for a Gompertz model with simple shape and with parameter vector \eqn{\theta^{\top} = [-5, 0.1]}, the argument should be specified as \code{theta = c(b0 = -5, b1 = 0.1)}. Note that in this example the parameter names are specified directly, this is required when \code{checkTheta = FALSE}. Although assigning the names to each parameter is not necessary when \code{checkTheta = TRUE}, it is advisable to ensure that the right values are assigned to the right parameter.  

If argument \code{checkTheta} is set to \code{TRUE}, then the vector of \code{theta} parameters is verified for consistency with the requirements of the model and shape selected. 
}

\references{
Colchero, F. and J.S. Clark (2012) Bayesian inference on age-specific survival from capture-recapture data for censored and truncated data. \emph{Journal of Animal Ecology}. 81, 139-149.

Colchero, F., O.R. Jones and M. Rebke. (2012) BaSTA: an R package for Bayesian estimation of age-specific survival from incomplete mark-recapture/recovery data with covariates. \emph{Method in Ecology and Evolution}. 3, 466-470. 

Colchero, F., \emph{et al.}  (2021) The long lives of primates and the "invariant rate of aging" hypothesis. \emph{Nature Communications} 12:3666

Cox, D. R., and Oakes D. (1984) \emph{Analysis of Survival Data}. Chapman and Hall, London.

Gompertz, B. (1825) On the nature of the function expressive of the law of human mortality, and on a new mode of determining the value of life contingencies. \emph{Philosophical Transactions of the Royal Society of London}, 115, 513-583.

Pinder III, J.E., Wiener, J.G. and Smith, M.H. (1978) The Weibull distribution: a new method of summarizing survivorship data. \emph{Ecology}, 59, 175-179.
}

\author{Fernando Colchero \email{fernando_colchero@eva.mpg.de}}

\seealso{\code{\link{CalcSurv}} to calculate age-specific survival, \code{\link{CalcFert}} to calculate age-specific fertility.}

\examples{
# Gompertz age specific mortality (default):
mort <- CalcMort(theta = c(b0 = -5, b1 = 0.1), x = 10)

# Siler age specific mortality:
mort <- CalcMort(beta = c(a0 = -1, a1 = 1, c = 0.0001, 
                          b0 = -6, b1 = 0.15), x = 10,
                  model = "GO", shape = "bathtub")

}

\keyword{models}

