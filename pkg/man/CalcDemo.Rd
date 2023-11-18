\name{CalcDemo}
\alias{CalcDemo}

\title{
 Calculating Parametric Age-Specific Mortality, Survival, and Fertility.
}

\description{
  \code{CalcDemo} is used to calculate age-specific demographic rates and summary statistics from parametric models of survival and fertility.
}

\usage{CalcDemo(theta = NULL, beta = NULL, x = NULL, dx = NULL, 
          model = "GO",  shape = "simple", modelFert = "quadratic", 
          type = "both", minSx = 0.01,  summarStats = TRUE, 
          ageMatur = 0, maxAge = NULL, agesAR = NULL,  
          SxValsAR = NULL)}

\arguments{
	\item{theta }{Numerical vector of age-specific mortality parameters (see details).}

	\item{beta }{Numerical vector of age-specific fertility parameters (see details).}

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

	\item{SxValsAR }{Numerical value for the levels of survival at which ageing rates should be calculated. If \code{NULL} (default) they are calculated at the ages when the survival function reaches 0.5, 0.2, and 0.05.}

}

\details{
\code{1) AGE-SPECIFIC MORTALITY AND SURVIVAL MODELS:}

The function \code{CalcSurv} uses parametric functions to calculate age-specific survival, defined as the probability of surviving to a given age. The function follows the conventions from package BaSTA (Colchero and Clark 2012, Colchero \emph{et al.} 2012, Colchero \emph{et al.} 2021). The mortality function describes how the risk of mortality changes with age, and is defined as 
\deqn{
\mu(x | \theta) = \lim_{\Delta x \rightarrow 0} \frac{\Pr[x < X < x + \Delta x | X > x]}{\Delta x},
}
where \eqn{X} is a random variable for ages at death, \eqn{x \geq 0} are ages and \eqn{\theta} is the vector of mortality parameters. From the mortality function, the survival function is then given by
\deqn{
S(x | \theta) = \exp[-\int_0^x \mu(t | \theta) dt].
}

\code{1.1)}\code{Argument }\dQuote{\code{model}}\code{:}

The \code{model} argument allows the user to choose between four basic mortality functions, namely 

(a) \code{model = }\dQuote{\code{EX}}: The exponential model (Cox and Oakes 1974), with constant mortality with age, specified as

\deqn{\mu_b(x | b) = b,}

where \eqn{b > 0}, with survival 
\deqn{
S_b(x | \theta) = \exp[-b x].
}

(b) \code{model = }\dQuote{\code{GO}}: The Gompertz mortality model (Gompertz 1925, Pletcher 1999), calculated as 

\deqn{\mu_b(x | b) = \exp(b_0 + b_1 x),}

where \eqn{-\infty < b_0, b_1 < \infty}, with survival
\deqn{
S_b(x | \theta) = \exp\left[\frac{e^{b_0}}{b_1}\left(1 - e^{b_1 x}\right)\right].
}

(c) \code{model = }\dQuote{\code{WE}}: The Weibull mortality model (Pinder III \emph{et al.} 1978) calculated as

\deqn{\mu_b(x | b) = b_0 b_1^{b_0} x^{b_0 -1},}

where \eqn{b_0, b_1 > 0}, with survival
\deqn{
S_b(x | \theta) = \exp\left[-(b_1 x)^{b_0}\right].
}

(d) \code{model = }\dQuote{\code{LO}}: The logistic mortality model (Pletcher 1999), calculated as 

\deqn{\mu_b(x | b) = \frac{\exp(b_0 + b_1 x)}{1 + b_2 \frac{e^{b_0}}{b_1} \left(e^{b_1 x}-1\right)},}
where \eqn{b_0, b_1, b_2 > 0}, with survival
\deqn{
S_b(x | \theta) = \left[1 + b_2 \frac{e^{b_0}}{b_1} \left(e^{b_1 x} - 1\right)\right]^{-1 / b_2}.
}

\code{1.2)}\code{Argument }\dQuote{\code{shape}}\code{:}

The \code{shape} argument allows the user to extend these models in order to explore more complex mortality shapes. 

(a) \code{shape = }\dQuote{\code{simple}}: (default) Leaves the model as defined above, with mortality given by
\deqn{
\mu(x | \theta) = \mu_b(x | \theta)
}
and survival
\deqn{
S(x | \theta) = S_b(x | \theta).
}

(b) \code{shape = }\dQuote{\code{Makeham}}: A constant is added to the mortality, such that the mortality is given by
\deqn{
\mu(x | \theta) = c + \mu_b(x | \theta_1),
}
where \eqn{\theta = [c, \theta_1]}, and with survival
\deqn{
S(x | \theta) = e^{-cx} S_b(x | \theta_1)
}

(c) \code{shape = }\dQuote{\code{bathtub}}: produces a concave shapes in mortality by adding a declining Gompertz term and a constant parameter to the basic mortality model, where the mortality function is
\deqn{
\mu(x | \theta) = \exp(a_0 - a_1 x) + c + \mu_b(x | \theta_1),
}
where \eqn{a_0 \in \mathbb{R}}, \eqn{a_1, c \geq 0} and \eqn{\theta_1 \in \theta} are specified based on argument \code{model}, and with survival
\deqn{
S(x | \theta) = \exp\left[\frac{e^{a_0}}{a_1}\left(e^{a_1 x} - 1\right)-cx\right] S_b(x | \theta_1).
}

\code{2) AGE-SPECIFIC FERTILITY MODELS:}

The age-specific fertility models correspond to the expected number of offspring produced by adults of a given age. Therefore, for a random variable \eqn{Y_{x}} with realizations \eqn{y_x} for the number of offspring produced by adults of age \eqn{x}, we have that \eqn{E(Y_x) = g(x | \beta)}, where \eqn{g: \mathbb{R}_{\ge 0} \rightarrow \mathbb{R}_{\ge 0}} is a non-negative smooth fertility function and \eqn{\beta} is a vector of parameters to be estimated. The functional forms of function \eqn{g} fall in two main categories, namely polynomial and distributional models.

2.a) Polynomial:

Of the models available in \code{paramDemo}, the \dQuote{\code{quadratic}}, \dQuote{\code{PeristeraKostaki}} (Peristera and Kostaki 2007), and \dQuote{\code{ColcheroMuller}} (Colchero  \emph{et al.} 2021, Muller \emph{et al.} 2020) fall within the polynomial category. Both, the \dQuote{\code{PeristeraKostaki}} and \dQuote{\code{ColcheroMuller}} are non-symmetric around the age at maximum fertility.

2.b) Distributional:

The distributional models are of the form

\deqn{g(x | \beta) = R f(x | \beta_1),}

where \eqn{f: \mathbb{R}_{\ge 0} \rightarrow \mathbb{R}_{\ge 0}} is a probability density function, \eqn{R > 0} is a parameter for the total fertility rate, and \eqn{\beta_1 \in \beta} is a vector of parameters. The
\dQuote{\code{Hadwiger}} (Hadwiger 1940), \dQuote{\code{gamma}} (Hoem \emph{et al.} 1981), \dQuote{\code{beta}}  (Hoem \emph{et al.} 1981), \dQuote{\code{skewNormal}} (Mazzuco and Scarpa 2011, 2015), \dQuote{\code{gammaMixture}} (Hoem \emph{et al.} 1981), \dQuote{\code{HadwigerMixture}} (Chandola \emph{et al.} 1991), \dQuote{\code{skewSymmetric}} (Mazzuco and Scarpa 2011, 2015), and \dQuote{\code{skewLogistic}} (Asili \emph{et al.} 2014) all fall in this category. Notably, the \dQuote{\code{gammaMixture}}, \dQuote{\code{HadwigerMixture}}, \dQuote{\code{skewSymmetric}}, \dQuote{\code{skewLogistic}} are appropriate when fertility might be bimodal (Hoem \emph{et al.} 1981, Chandola \emph{et al.} 1999, Mazzuco and Scarpa 2011, 2015, Asili \emph{et al.} 2014). 

\code{3) SPECIFYING THETA PARAMETERS:}

Argument \code{theta} requires a numerical vector of parameters. For instance, for a Gompertz model with simple shape and with parameter vector \eqn{\theta^{\top} = [-5, 0.1]}, the argument should be specified as \code{theta = c(b0 = -5, b1 = 0.1)}. Note that in this example the parameter names are specified directly, this is required when \code{checkTheta = FALSE}. Although assigning the names to each parameter is not necessary when \code{checkTheta = TRUE}, it is advisable to ensure that the right values are assigned to the right parameter.  

If argument \code{checkTheta} is set to \code{TRUE}, then the vector of \code{theta} parameters is verified for consistency with the requirements of the model and shape selected. 

\code{4) SPECIFYING BETA PARAMETERS:}

Argument \code{beta} requires a vector of parameters. For instance, for a quadratic model with parameter vector \eqn{\beta^{\top} = [0.5, 0.01, 10]}, the argument should be specified as \code{theta = c(b0 = 0.5, b1 = 0.01, b2 = 10)}. Note that in this example the parameter names are specified directly, this is required when \code{checkBeta = FALSE}. Although assigning the names to each parameter is not necessary when \code{checkBeta = TRUE}, it is advisable to ensure that the right values are assigned to the right parameter.  

If argument \code{checkBeta} is set to \code{TRUE}, then the vector of \code{beta} parameters is verified for consistency with the requirements of the model and shape selected. 

}

\value{
The output consist of two lists, one for survival and one for fertility. The survival list can be called as \code{object$surv} while the fertility list is in \code{object$fert}. Both include the following outputs:
  \item{functs }{data.frame including the follwing columns for survival \code{age}, \code{mort} for the mortality function, \code{surv} for the survival, \code{pdf} for the probability density function of ages at death, and \code{cumhaz} for the cumulative hazards. For fertility it includes columns for \code{age} and \code{fert} for the age-specific fertility.} 
  
  \item{summStats }{list with elements \code{calculated} indicating whether summary statistics were calculated, \code{ageingRates} with the resulting actuarial ageing rates, and \code{summStatsMort} with ages at which the survival function reaches different levels (default 0.5, 0.2, 0.05) (e.g., \code{AgeSx=0.5}), with the life expectancy \code{lifeExp}, \code{lifespIneq} for the lifespan inequality (see Demetrius 1974), \code{lifespEqual} with the lifespan equality as in Colchero \emph{et al.} (2021), \code{Gini} for the Gini coefficient, and \code{CoeffVar} for the coefficient of variation in the distribution of ages at death. For fertility, the list includes a logical \code{calculated} to indicate if summary statistics were calculated, and \code{summStatsFert} with the age at maximum fertility calculations (see \code{\link{CalcAgeMaxFert}}).} 
  
  \item{settings }{list with the details of the parameters and the models used.} 
  
  \item{analyzed }{logical indicating whether \code{surv} or \code{fert} were calculated.} 
}

\references{
Asili, S., Rezaei, S. & Najjar, L. (2014) Using Skew-Logistic Probability Density Function as a Model for Age-Specific Fertility Rate Pattern. \emph{BioMed Research International}, 2014, 790294.

Azzalini, A. (1985) A Class of Distributions Which Includes the Normal Ones. \emph{Scandinavian Journal of Statistics}, 2, 171-178.

Chandola, T., Coleman, D.A. & Hiorns, R.W. (1999) Recent European fertility patterns: Fitting curves to `distorted' distributions. \emph{Population Studies}, 53, 317-329.

Colchero, F. (\emph{In prep.}) Inference on age-specific fertility in ecology and evolution. Learning from other disciplines and improving the state of the art. 

Colchero, F., Eckardt, W. & Stoinski, T. (2021) Evidence of demographic buffering in an endangered great ape: Social buffering on immature survival and the role of refined sex-age classes on population growth rate. \emph{Journal of Animal Ecology}, 90, 1701-1713.

Colchero, F. and J.S. Clark (2012) Bayesian inference on age-specific survival from capture-recapture data for censored and truncated data. \emph{Journal of Animal Ecology}. 81, 139-149.

Colchero, F., O.R. Jones and M. Rebke. (2012) BaSTA: an R package for Bayesian estimation of age-specific survival from incomplete mark-recapture/recovery data with covariates. \emph{Method in Ecology and Evolution}. 3, 466-470. 

Colchero, F., \emph{et al.}  (2021) The long lives of primates and the "invariant rate of aging" hypothesis. \emph{Nature Communications} 12:3666

Cox, D. R., and Oakes D. (1984) \emph{Analysis of Survival Data}. Chapman and Hall, London.

Gompertz, B. (1825) On the nature of the function expressive of the law of human mortality, and on a new mode of determining the value of life contingencies. \emph{Philosophical Transactions of the Royal Society of London}, 115, 513-583.

Hadwiger, H. (1940) Eine analytische Reproduktionssunktion fur biologische Gesamtheiten. \emph{Scandinavian Actuarial Journal}, 1940, 101-113.

Hoem, J.M., Madien, D., Nielsen, J.L., Ohlsen, E.M., Hansen, H.O. & Rennermalm, B. (1981) Experiments in modelling recent Danish fertility curves. \emph{Demography}, 18, 231-244.

Mazzuco, S. & Scarpa, B. (2011) Fitting age-specific fertility rates by a skew- symmetric probability density function. (Working paper 10), University of Padua.

Mazzuco, S. & Scarpa, B. (2015) Fitting age-specific fertility rates by a flexible generalized skew normal probability density function. \emph{Journal of the Royal Statistical Society: Series A}, 178, 187-203.

Muller, M. N., Blurton Jones, N. G, Colchero, F., Thompson, M. E., Enigk, D. K. (2020) Sexual dimorphism in chimpanzee (\emph{Pan troglodytes schweinfurthii}) and human age-specific fertility. \emph{Journal of human evolution}, 144, 102795.

Peristera, P. & Kostaki, A. (2007) Modeling fertility in modern populations. \emph{Demographic Research}, 16, 141-194.

Pinder III, J.E., Wiener, J.G. and Smith, M.H. (1978) The Weibull distribution: a new method of summarizing survivorship data. \emph{Ecology}, 59, 175-179.

}
\author{Fernando Colchero \email{fernando_colchero@eva.mpg.de}}

\seealso{
\code{\link{CalcSurv}} to calculate age-specific survival, \code{\link{CalcMort}} to calculate age-specific mortality, \code{\link{CalcFert}} to calculate age-specific fertility.

\code{\link{CalcAgeingRate}} to calculate ageing rates from parametric models of age-specific mortality. \code{\link{CalcRemainLifeExp}} to calculate remaining life expectancy from parametric models of age-specific mortality. \code{\link{CalcAgeMaxFert}} to calculate the age at maximum fertility from parametric models of age-specific fertility.
}



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

