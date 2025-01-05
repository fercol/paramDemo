\name{CalcMort}
\alias{CalcMort}
\alias{CalcSurv}

\title{
 Calculating Parametric Age-Specific Mortality and survival.
}

\description{
 \code{CalcMort} and \code{CalcSurv} are used to calculate the age-specific mortality and survival functions from different parametric models.
}

\usage{
CalcMort(theta, x, model = "GO", shape = "simple", 
          checkTheta = TRUE)

CalcSurv(theta, x, model = "GO", shape = "simple", 
          checkTheta = TRUE)          
}

\arguments{
	\item{theta }{Numerical vector of age-specific mortality parameters (see details).}

	\item{x }{Numerical vector of ages at which to calculate mortality.}

	\item{model }{The underlying mortality model to be used. \code{"EX"} = exponential,\code{"GO"} = Gompertz, \code{"WE"} = Weibull and \code{"LO"} = logistic (see \code{details}).}

	\item{shape }{The overall shape of the model. Values are: \code{simple} = no extra parameters added; \code{Makeham} = a constant parameter is added to the mortality; and \code{bathtub} = a Gompertz declining mortality for early ages and a constant parameter are added to the mortality model (see \code{details}).}

	\item{checkTheta }{Logical to verify that the \code{beta} parameters conform with the model's specification (see details).}

}

\details{

\bold{1) Age-specific mortality and survival models:}

The function \code{CalcMort} uses parametric functions to calculate age-specific survival, defined as the probability of surviving to a given age. The function follows the conventions from package BaSTA (Colchero and Clark 2012, Colchero \emph{et al.} 2012, Colchero \emph{et al.} 2021). The mortality function or hazard rate describes how the risk of mortality changes with age, and is defined as 
\deqn{
\mu(x | \theta) = \lim_{\Delta x \rightarrow 0} \frac{\Pr[x < X < x + \Delta x | X > x]}{\Delta x},
}
where \eqn{X} is a random variable for ages at death, \eqn{x \geq 0} are ages and \eqn{\theta} is the vector of mortality parameters. From the mortality function, the survival function is then given by
\deqn{
S(x | \bm{\theta}) = \Pr[X > x] = \exp[-\int_0^x \mu(t | \bm{\theta}) dt],
}
with the cumulative distribution function \eqn{F(x | \bm{\theta}) = 1 - S(x | \bm{\theta})}. The probability density function of ages at death is 
\deqn{
f(x | \bm{\theta}) = \mu(x | \bm{\theta}) S(x | \bm{\theta}), 
}
for \eqn{x \geq 0}. 

\code{1.1)}\code{Argument }\dQuote{\code{model}}\code{:}

The \code{model} argument allows the user to choose between four basic mortality functions, namely 

(a) \code{model = }\dQuote{\code{EX}}: The exponential model (Cox and Oakes 1974), with constant mortality with age, specified as

\deqn{\mu_b(x | \bm{\theta}) = b,}

where \eqn{b > 0}, with survival 
\deqn{
S_b(x | \bm{\theta}) = \exp[-b x].
}

(b) \code{model = }\dQuote{\code{GO}}: The Gompertz mortality model (Gompertz 1925, Pletcher 1999), calculated as 

\deqn{\mu_b(x | \bm{\theta}) = \exp(b_0 + b_1 x),}

where \eqn{-\infty < b_0, b_1 < \infty}, with survival
\deqn{
S_b(x | \bm{\theta}) = \exp\left[\frac{e^{b_0}}{b_1}\left(1 - e^{b_1 x}\right)\right].
}

(c) \code{model = }\dQuote{\code{WE}}: The Weibull mortality model (Pinder III \emph{et al.} 1978) calculated as

\deqn{\mu_b(x | \bm{\theta}) = b_0 b_1^{b_0} x^{b_0 -1},}

where \eqn{b_0, b_1 > 0}, with survival
\deqn{
S_b(x | \bm{\theta}) = \exp\left[-(b_1 x)^{b_0}\right].
}

(d) \code{model = }\dQuote{\code{LO}}: The logistic mortality model (Pletcher 1999), calculated as 

\deqn{\mu_b(x | \bm{\theta}) = \frac{\exp(b_0 + b_1 x)}{1 + b_2 \frac{e^{b_0}}{b_1} \left(e^{b_1 x}-1\right)},}
where \eqn{b_0, b_1, b_2 > 0}, with survival
\deqn{
S_b(x | \bm{\theta}) = \left[1 + b_2 \frac{e^{b_0}}{b_1} \left(e^{b_1 x} - 1\right)\right]^{-1 / b_2}.
}

\code{1.2)}\code{Argument }\dQuote{\code{shape}}\code{:}

The \code{shape} argument allows the user to extend these models in order to explore more complex mortality shapes. 

(a) \code{shape = }\dQuote{\code{simple}}: (default) Leaves the model as defined above, with mortality given by
\deqn{
\mu(x | \bm{\theta}) = \mu_b(x | \bm{\theta})
}
and survival
\deqn{
S(x | \bm{\theta}) = S_b(x | \bm{\theta}.
}

(b) \code{shape = }\dQuote{\code{Makeham}}: A constant is added to the mortality, such that the mortality is given by
\deqn{
\mu(x | \bm{\theta}) = c + \mu_b(x | \bm{\theta}_1),
}
where \eqn{\bm{\theta} = [c, \bm{\theta}_1]}, and with survival
\deqn{
S(x | \bm{\theta}) = e^{-cx} S_b(x | \bm{\theta}_1)
}

The most common models with this shape is the Gompertz-Makeham model (Gompertz 1825, Makeham 1866).

(c) \code{shape = }\dQuote{\code{bathtub}}: produces a concave shapes in mortality by adding a declining Gompertz term and a constant parameter to the basic mortality model, where the mortality function is
\deqn{
\mu(x | \bm{\theta}) = \exp(a_0 - a_1 x) + c + \mu_b(x | \bm{\theta}_1),
}
where \eqn{a_0 \in \mathbb{R}}, \eqn{a_1, c \geq 0} and \eqn{\bm{\theta}_1 \subset \bm{\theta}} are specified based on argument \code{model}, and with survival
\deqn{
S(x | \theta) = \exp\left[\frac{e^{a_0}}{a_1}\left(e^{a_1 x} - 1\right)-cx\right] S_b(x | \theta_1).
}

The most widely use \dQuote{\code{bathtub}} shaped model is the \bold{Siler mortality model} (Siler 1979), which provides considerably good fits to mammalian data. The arguments for the Siler model are:

\preformatted{CalcMort(..., model = "GO", shape = "bathtub", ...)}

\bold{2) Specifying theta parameters:}

Argument \code{theta} requires a numerical vector of parameters. For instance, for a Gompertz model with simple shape and with parameter vector \eqn{\bm{\theta}^{\top} = [b_0, b_1]} where \eqn{b_0 = -5} and \eqn{b_1 = 0.1}, the argument should be specified as 

\code{CalcMort(theta = c(b0 = -5, b1 = 0.1), ...)} 

Note that in this example the parameter names are specified directly, this is required when \code{checkTheta = FALSE}. Although assigning the names to each parameter is not necessary when \code{checkTheta = TRUE}, it is advisable to ensure that the right values are assigned to the right parameter.  

If argument \code{checkTheta} is set to \code{TRUE}, then the vector of \code{theta} parameters is verified for consistency with the requirements of the model and shape selected. 

}

\value{
\code{CalcMort} and \code{CalcSur} return a vector of class \dQuote{\code{numeric}} with the hazard rate or cumulative survival values at the specified ages, respectively.
}

\references{
Colchero, F. and J.S. Clark (2012) Bayesian inference on age-specific survival from capture-recapture data for censored and truncated data. \emph{Journal of Animal Ecology}. 81, 139-149.

Colchero, F., O.R. Jones and M. Rebke. (2012) BaSTA: an R package for Bayesian estimation of age-specific survival from incomplete mark-recapture/recovery data with covariates. \emph{Method in Ecology and Evolution}. 3, 466-470. 

Colchero, F., \emph{et al.} (2021) The long lives of primates and the "invariant rate of aging" hypothesis. \emph{Nature Communications} 12:3666

Cox, D. R., and Oakes D. (1984) \emph{Analysis of Survival Data}. Chapman and Hall, London.

Gompertz, B. (1825) On the nature of the function expressive of the law of human mortality, and on a new mode of determining the value of life contingencies. \emph{Philosophical Transactions of the Royal Society of London}, 115, 513-583.

Makeham, W. M. On the law of mortality (1866). \emph{Journal of the Institute of Actuaries} 13, 1-34.
  
Pinder III, J.E., Wiener, J.G. and Smith, M.H. (1978) The Weibull distribution: a new method of summarizing survivorship data. \emph{Ecology}, 59, 175-179.

Siler, W. A (1979) competing-risk model for animal mortality. \emph{Ecology} 60, 750-757.

}

\author{Fernando Colchero \email{fernando_colchero@eva.mpg.de}}

\seealso{\code{\link{CalcSurv}} to calculate age-specific survival, \code{\link{CalcFert}} to calculate age-specific fertility.}

\examples{
# Gompertz age specific mortality (default):
mort <- CalcMort(theta = c(b0 = -5, b1 = 0.1), x = 10)

# Siler age specific mortality:
mort <- CalcMort(theta = c(a0 = -1, a1 = 1, c = 0.0001, 
                          b0 = -6, b1 = 0.15), x = 10,
                  model = "GO", shape = "bathtub")

# Gompertz age specific survival (default):
surv <- CalcSurv(theta = c(b0 = -5, b1 = 0.1), x = 10)

# Siler age specific survival:
surv <- CalcSurv(theta = c(a0 = -1, a1 = 1, c = 0.0001, 
                          b0 = -6, b1 = 0.15), x = 10,
                  model = "GO", shape = "bathtub")


}

\keyword{distribution}

