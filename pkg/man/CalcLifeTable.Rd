\name{CalcLifeTable}
\alias{CalcLifeTable}

\title{
 Calculating a Life Table from Data.
}

\description{
 \code{CalcLifeTable} uses non-parametric methods to calculate life tables and confidence intervals.
}

\usage{
CalcLifeTable(ageLast, ageFirst = NULL, departType, dx = 1, 
                          calcCIs = FALSE, nboot = 1000, alpha = 0.05)
}

\arguments{
	\item{ageLast }{Numerical vector with the ages at last detection (i.e., death and censoring) (see \code{details})}

	\item{ageFirst }{Numerical vector of ages at first detection (i.e., truncation). If \code{NULL} then all values are set to 0 (see \code{details}).}

	\item{departType }{Character string vector for the type of departure (i.e., last detection), with values \dQuote{\code{D}} for death and \dQuote{\code{C}} for censoring (see \code{details}).}

	\item{dx }{Age interval size, default set at 1 (see \code{details})}

  \item{calcCIs }{Logical indicating whether confidence intervals should be calculated}
  
	\item{nboot }{Number of bootstrap iterations}

	\item{alpha }{Alpha level. Default is 0.05 for 95\% CIs}

}

\details{

\bold{1) Data structure:}

\code{CalcLifeTable} allows to construct life tables for data that includes the following types of records:

\itemize{ 
\item\bold{Uncensored}: individuals with known ages at death;

\item\bold{right-censored}: individuals last seen alive;

\item\bold{left-truncated}: individuals born before the start of the study and are truncated at the age of entry. 
}

The data required are the ages at last detection (i.e., \emph{uncensored} or \emph{right-censored}) passed through argument \dQuote{\code{ageLast}}, the type of departure via argument \dQuote{\code{departType}}, which takes two values, namely \dQuote{\code{D}} for death, and \dQuote{\code{C}} for censored (i.e., right-censored). 

In addition, if there is \emph{left-truncation}, it takes the ages at entry to the study by means of argument \dQuote{\code{ageFirst}}. If all the individuals were born during the study, the value of \dQuote{\code{ageFirst}} can be left as \code{NULL}, which will make them all equal to 0. 
 
 
\bold{2) Computing life tables}

To calculate life tables, the function uses conventional formal demgraphic methods as depicted by Preston \emph{et al.} (2001). Argument \dQuote{\code{ageLast}} provides a vector of ages at last detection, \eqn{\bold{x}^{\top} = [x_1, x_2, \dots, x_n]}, while argument \dQuote{\code{ageFirst}} provides a vector of ages at first detection \eqn{\bold{y}^{\top} = [y_1, y_2, \dots, y_n]}. From argument \dQuote{\code{departType}} the function produces an indicator vector for censoring \eqn{\bold{v} = \{v_i\}_{i \in \mathbb{N}_n}} where \eqn{v_i = 1} if individual \eqn{i} is censored and 0 otherwise. 

The function creates a partition of the interval of ages between 0 and \eqn{\max(x)}, for age intervals \eqn{[x, x + \Delta x)} where \eqn{\Delta x} is specified by argument \dQuote{\code{dx}}. As default \code{dx = 1}. At each age interval, the function calculates the following variables:

\itemize{
\item \code{Nx}: which corresponds to number of individuals that entered the interval, but considering the proportion of time they were present within the interval as a function of left-truncation. It is given by
\deqn{
N_x = \sum_{i \in I_x} \lambda_{i,x},
}
where \eqn{I_x} is the subset of individuals recorded within the interval, and \eqn{\lambda_{i, x}} is the proportion of time during the age interval each individual was present in the study. For individuals that entered the study before \eqn{x} then \eqn{\lambda_{i,x} = 1}, while for those that were truncated within the interval \eqn{\lambda_{i,x} = (x + \Delta x - y_i) / \Delta x}. 

\item \code{Dx}: the number of individuals dying in the interval.

\item \code{qx}: the age-specific mortality probability, calculated as \eqn{q_x = D_x / N_x}. 

\item \code{px}: the age-specific survival probability, given by \eqn{p_x = 1 - q_x}.

\item \code{lx}: the life table survival calculated as 
\deqn{
l_x = \prod_{j = 0}^{x-1} p_j
}
where \eqn{l_0 = 1}. 

\item \code{ax}: Proportion of the interval lived by those that died in the interval, given by
\deqn{
a_x = \frac{\sum_{i\in J_{x}} \delta_{i,x}}{D_x}
}
where \eqn{J_{x}} is the subset of individuals that died within the interval, and \eqn{\delta_{i,x}} is the proportion lived by those individuals from the start of the interval to their deaths, this is \eqn{\delta_{i,x} = (x_i - x) / \Delta x}. 

\item \code{Lx}: The number of individual years lived within the interval, given by
\deqn{
L_x = l_x (1 - a_x q_x)
}

\item \code{Tx}: The total number of individual years lived after age \eqn{x}, given by \eqn{T_x = \sum_{j = x}^{\infty} L_j \Delta x}

\item \code{ex}: the remaining life expectancy at the beginning of each age interval, calculated as \eqn{e_x = T_x / l_x}.
}

\bold{3) Calculating confidence intervals}
If argument \dQuote{\code{calcCIs}} is set to \code{TRUE}, the function uses a non-parametric bootstrap by sampling with replacement the data. Argument \code{nboot} specifies the number of bootstrap steps, with default \code{nboot = 2000}. From each re-sampled dataset, it uses function \code{\link{CalcLifeTable}} to construct the corresponding life table and stores the values of $l_x$, $q_x$, $p_x$, and $e_x$ from each iteration. From these, it calculates quantiles at the given alpha level. 


}

\value{
\code{CalcLifeTable} returns an object of class \dQuote{\code{paramDemoLT}} with output consisting of a list with the life table and, if indicated by the user, with the confidence interval information. The life table in matrix format includes the following columns:
	\item{Ages }{Ages with increments given by \code{dx}}
	
	\item{Nx }{Number of individuals entering the interval, does not need to be an integer since it considers truncation}

	\item{Dx }{Number of individuals that died within the age interval}

	\item{lx }{Survival (i.e., cumulative survival)}

	\item{px }{Age-specific survival probability}

	\item{qx }{Age-specific mortality probability}

	\item{Lx }{Number of individual years lived within the interval}
	
	\item{Tx }{Number of individual years lived after age x}

	\item{ex }{Remaining life expectancy at each age}
	

If argument \code{calcCIs = TRUE}, the function also returns a list containing the following components:

	\item{lx }{Matrix with Ages and mean and upper and lower CIs for the life table survival}
	
	\item{qx }{Matrix with Ages and mean and upper and lower CIs for the age-specific mortality probability}
	
	\item{px }{Matrix with Ages and mean and upper and lower CIs for the age-specific survival probability}
	
	\item{ex }{Matrix with Ages and mean and upper and lower CIs for the remaining life expectancy}
	
	\item{Settings }{Numerical vector including whether CIs were calculated, the number of bootstrap iterations, \dQuote{\code{nboot}} and the alpha level \dQuote{\code{alpha}}}
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

# Calculate life table with 95\% CIs:
ltCIs <- CalcLifeTable(ageLast = ages, departType = rep("D", 100), 
                       calcCIs = TRUE, nboot = 100)

}

\keyword{survival}

