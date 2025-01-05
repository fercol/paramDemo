\name{CalcFert}
\alias{CalcFert}

\title{
 Calculating Parametric Age-Specific Fertility.
}

\description{
 \code{CalcFert} is used to calculate age-specific fertility from different parametric models.
}

\usage{
CalcFert(beta, x, modelFert = "quadratic", checkBeta = TRUE)
}

\arguments{
	\item{beta }{Numerical vector of age-specific fertility parameters (see details).}

	\item{x }{Numerical vector of ages at which to calculate fertility.}

	\item{modelFert }{Age-specific fertility model. Options are \dQuote{\code{quadratic}} (default), \dQuote{\code{PeristeraKostaki}}, \dQuote{\code{ColcheroMuller}}, \dQuote{\code{Hadwiger}}, \dQuote{\code{gamma}}, \dQuote{\code{beta}}, \dQuote{\code{skewNormal}}, \dQuote{\code{gammaMixture}}, \dQuote{\code{HadwigerMixture}}, \dQuote{\code{skewSymmetric}}, \dQuote{\code{skewLogistic}} (see details)}

	\item{checkBeta }{Logical to verify that the \code{beta} parameters conform with the model's specification (see details).}

}

\details{
\bold{1) Age-specific fertility models:}

The age-specific fertility models correspond to the expected number of offspring produced by adults of a given age. Therefore, for a random variable \eqn{Y_{x}} with realizations \eqn{y_x} for the number of offspring produced by adults of age \eqn{x}, we have that \eqn{E(Y_x) = g(x | \beta)}, where \eqn{g: \mathbb{R}_{\ge 0} \rightarrow \mathbb{R}_{\ge 0}} is a non-negative smooth fertility function and \eqn{\beta} is a vector of parameters to be estimated. The functional forms of function \eqn{g} fall in two main categories, namely polynomial and distributional models.

2.a) Polynomial:

Of the models available in \code{paramDemo}, the \dQuote{\code{quadratic}}, \dQuote{\code{PeristeraKostaki}} (Peristera and Kostaki 2007), and \dQuote{\code{ColcheroMuller}} (Colchero  \emph{et al.} 2021, Muller \emph{et al.} 2020) fall within the polynomial category. Both, the \dQuote{\code{PeristeraKostaki}} and \dQuote{\code{ColcheroMuller}} are non-symmetric around the age at maximum fertility.

2.b) Distributional:

The distributional models are of the form

\deqn{g(x | \beta) = R f(x | \beta_1),}

where \eqn{f: \mathbb{R}_{\ge 0} \rightarrow \mathbb{R}_{\ge 0}} is a probability density function, \eqn{R > 0} is a parameter for the total fertility rate, and \eqn{\beta_1 \in \beta} is a vector of parameters. The
\dQuote{\code{Hadwiger}} (Hadwiger 1940), \dQuote{\code{gamma}} (Hoem \emph{et al.} 1981), \dQuote{\code{beta}}  (Hoem \emph{et al.} 1981), \dQuote{\code{skewNormal}} (Mazzuco and Scarpa 2011, 2015), \dQuote{\code{gammaMixture}} (Hoem \emph{et al.} 1981), \dQuote{\code{HadwigerMixture}} (Chandola \emph{et al.} 1991), \dQuote{\code{skewSymmetric}} (Mazzuco and Scarpa 2011, 2015), and \dQuote{\code{skewLogistic}} (Asili \emph{et al.} 2014) all fall in this category. Notably, the \dQuote{\code{gammaMixture}}, \dQuote{\code{HadwigerMixture}}, \dQuote{\code{skewSymmetric}}, \dQuote{\code{skewLogistic}} are appropriate when fertility might be bimodal (Hoem \emph{et al.} 1981, Chandola \emph{et al.} 1999, Mazzuco and Scarpa 2011, 2015, Asili \emph{et al.} 2014). 

\bold{2) Specifying beta parameters:}

Argument \code{beta} requires a vector of parameters. For instance, for a quadratic model with parameter vector \eqn{\beta^{\top} = [0.5, 0.01, 10]}, the argument should be specified as \code{theta = c(b0 = 0.5, b1 = 0.01, b2 = 10)}. Note that in this example the parameter names are specified directly, this is required when \code{checkBeta = FALSE}. Although assigning the names to each parameter is not necessary when \code{checkBeta = TRUE}, it is advisable to ensure that the right values are assigned to the right parameter.  

If argument \code{checkBeta} is set to \code{TRUE}, then the vector of \code{beta} parameters is verified for consistency with the requirements of the model and shape selected. 
}

\value{
\code{CalcFert} returns a vector of class \dQuote{\code{numeric}} with the fertility values at the specified ages.
}

\references{
Asili, S., Rezaei, S. & Najjar, L. (2014) Using Skew-Logistic Probability Density Function as a Model for Age-Specific Fertility Rate Pattern. \emph{BioMed Research International}, 2014, 790294.

Azzalini, A. (1985) A Class of Distributions Which Includes the Normal Ones. \emph{Scandinavian Journal of Statistics}, 2, 171-178.

Chandola, T., Coleman, D.A. & Hiorns, R.W. (1999) Recent European fertility patterns: Fitting curves to `distorted' distributions. \emph{Population Studies}, 53, 317-329.

Colchero, F. (\emph{In prep.}) Inference on age-specific fertility in ecology and evolution. Learning from other disciplines and improving the state of the art. 

Colchero, F., Eckardt, W. & Stoinski, T. (2021) Evidence of demographic buffering in an endangered great ape: Social buffering on immature survival and the role of refined sex-age classes on population growth rate. \emph{Journal of Animal Ecology}, 90, 1701-1713.

Hadwiger, H. (1940) Eine analytische Reproduktionssunktion fur biologische Gesamtheiten. \emph{Scandinavian Actuarial Journal}, 1940, 101-113.

Hoem, J.M., Madien, D., Nielsen, J.L., Ohlsen, E.M., Hansen, H.O. & Rennermalm, B. (1981) Experiments in modelling recent Danish fertility curves. \emph{Demography}, 18, 231-244.

Mazzuco, S. & Scarpa, B. (2011) Fitting age-specific fertility rates by a skew- symmetric probability density function. (Working paper 10), University of Padua.

Mazzuco, S. & Scarpa, B. (2015) Fitting age-specific fertility rates by a flexible generalized skew normal probability density function. \emph{Journal of the Royal Statistical Society: Series A}, 178, 187-203.

Muller, M. N., Blurton Jones, N. G, Colchero, F., Thompson, M. E., Enigk, D. K. (2020) Sexual dimorphism in chimpanzee (\emph{Pan troglodytes schweinfurthii}) and human age-specific fertility. \emph{Journal of human evolution}, 144, 102795.

Peristera, P. & Kostaki, A. (2007) Modeling fertility in modern populations. \emph{Demographic Research}, 16, 141-194.

}
\author{Fernando Colchero \email{fernando_colchero@eva.mpg.de}}

\seealso{\code{\link{CalcMort}} to calculate age-specific mortality, \code{\link{CalcSurv}} to calculate age-specific survival}

\examples{
# Age specific fertility based on quadratic model (default):
fert <- CalcFert(beta = c(b0 = 0.5, b1 = 0.01, b2 = 10), x = 10)

# Age specific fertility based on gamma model:
fert <- CalcFert(beta = c(b0 = 13, b1 = 2, b2 = 0.15), x = 10,
                  modelFert = "gamma")
}

\keyword{models}

