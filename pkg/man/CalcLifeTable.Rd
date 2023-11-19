\name{CalcLifeTable}
\alias{CalcLifeTable}

\title{
 Calculating a Life Table from Data.
}

\description{
 \code{CalcLifeTable} uses non-parametric methods to calculate life tables.
}

\usage{CalcLifeTable (ageLast, ageFirst = NULL, departType, dx = 1)}

\arguments{
	\item{ageLast }{Numerical vector with the ages at last detection (i.e., death and censoring) (see \code{details})}

	\item{ageFirst }{Numerical vector of ages at first detection (i.e., truncation). If \code{NULL} then all values are set to 0 (see \code{details}).}

	\item{departType }{Character string vector for the type of departure (i.e., last detection), with values \dQuote{\code{D}} for death and \dQuote{\code{C}} for censoring (see \code{details}).}

	\item{dx }{Age interval size, default set at 1 (see \code{details})}

}

\details{
The function allows to construct life tables for data that includes the following types of records:

- \bold{Uncensored}: individuals with known ages at death;

- \bold{right-censored}: individuals last seen alive;

- \bold{left-truncated}: individuals born before the start of the study and are truncated at the age of entry. 

The data required are the ages at last detection (i.e., \emph{uncensored} or \emph{right-censored}) passed through argument \dQuote{\code{ageLast}}, the type of departure via argument \dQuote{\code{departType}}, which takes two values, namely \dQuote{\code{D}} for death, and \dQuote{\code{C}} for censored (i.e., right-censored). 

In addition, if there is \emph{left-truncation}, it takes the ages at entry to the study by means of argument \dQuote{\code{ageFirst}}. If all the individuals were born during the study, the value of \dQuote{\code{ageFirst}} can be left as \code{NULL}, which will make them all equal to 0. 
 
To calculate life tables, the function uses conventional formal demgraphic methods as depicted by Preston \emph{et al.} (2001). Argument \dQuote{\code{ageLast}} provides a vector of ages at last detection, \eqn{\bold{x}^{\top} = [x_1, x_2, \dots, x_n]}, while argument \dQuote{\code{ageFirst}} provides a vector of ages at first detection \eqn{\bold{y}^{\top} = [y_1, y_2, \dots, y_n]}. From argument \dQuote{\code{departType}} the function produces an indicator vector for censoring \eqn{\bold{v} = \{v_i\}_{i \in \mathbb{N}_n}} where \eqn{v_i = 1} if individual \eqn{i} is censored and 0 otherwise. 

}

\value{
The function produces a life table of class \code{matrix} with the following columns:
	\item{Ages }{Ages with increases given by \code{dx}}
	
	\item{Nx }{Number of individuals entering the interval, does not need to be an integer since it considers truncation}

	\item{Dx }{Number of individuals that died within the age interval}

	\item{lx }{Survival (i.e., cumulative survival)}

	\item{px }{Age-specific survival probability}

	\item{qx }{Age-specific mortality probability}

	\item{Lx }{Number of individual years lived within the interval}
	
	\item{Tx }{Number of individual years lived after age x}

	\item{ex }{Remaining life expectancy at each age}
}

\references{

Preston, S.H., Heuveline, P. and Guillot, M. (2001) Demography: Measuring and Modeling Population Processes. Blackwell, Oxford.
}
\author{Fernando Colchero \email{fernando_colchero@eva.mpg.de}}

\seealso{
\code{\link{CalcProductLimitEst}} to calculate product limit estimators.
}

\examples{
# Simulate age at death data from Gompertz model:
ages <- SampleRandAge(n = 100, theta = c(b0 = -5, b1 = 0.1))

# Calculate life table:
lt <- CalcLifeTable(ageLast = ages, departType = rep("D", 100))

}

\keyword{models}

