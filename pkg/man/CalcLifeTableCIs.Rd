\name{CalcLifeTableCIs}
\alias{CalcLifeTableCIs}

\title{
 Non-Parametric Bootstrap to Calculate Life Table Confidence Intervals.
}

\description{
 \code{CalcLifeTableCIs} uses non-parametric bootstrap to calculate confidence intervals for a life table.
}

\usage{CalcLifeTableCIs (ageLast, ageFirst = NULL, departType, 
                  nboot = 2000,  alpha = 0.05, dx = 1)}

\arguments{
	\item{ageLast }{Numerical vector with the ages at last detection (i.e., death and censoring)}

	\item{ageFirst }{Numerical vector of ages at first detection (i.e., truncation). If \code{NULL} then all values are set to 0.}

	\item{departType }{Character string vector for the type of departure (i.e., last detection), with values \dQuote{\code{D}} for death and \dQuote{\code{C}} for censoring.}

	\item{nboot }{Number of bootstrap iterations}

	\item{alpha }{Alpha level. Default is 0.05 for 95\% CIs}

	\item{dx }{Age interval size, default set at 1.}

}

\details{
\subsection{1) Data structure:}

The function allows to construct life tables for data that includes the following types of records:

- \bold{Uncensored}: individuals with known ages at death;

- \bold{right-censored}: individuals last seen alive;

- \bold{left-truncated}: individuals born before the start of the study and are truncated at the age of entry. 

The data required are the ages at last detection (i.e., \emph{uncensored} or \emph{right-censored}) passed through argument \dQuote{\code{ageLast}}, the type of departure via argument \dQuote{\code{departType}}, which takes two values, namely \dQuote{\code{D}} for death, and \dQuote{\code{C}} for censored (i.e., right-censored). 

In addition, if there is \emph{left-truncation}, it takes the ages at entry to the study by means of argument \dQuote{\code{ageFirst}}. If all the individuals were born during the study, the value of \dQuote{\code{ageFirst}} can be left as \code{NULL}, which will make them all equal to 0. 
 
 
\subsection{2) Computing life tables}

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

}

\value{
	\item{fill up }{FILL UP}
}

\author{Fernando Colchero \email{colchero@imada.sdu.dk}}

\seealso{FILL UP}

\examples{}

\keyword{models}

