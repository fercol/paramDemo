\name{CalcProductLimitEst}
\alias{CalcProductLimitEst}

\title{
 FILL UP
}

\description{
 FILL UP
}

\usage{CalcProductLimitEst (ageLast, ageFirst = NULL, departType)}

\arguments{
	\item{ageLast }{Numerical vector with the ages at last detection (i.e., death and censoring) (see \code{details})}

	\item{ageFirst }{Numerical vector of ages at first detection (i.e., truncation). If \code{NULL} then all values are set to 0 (see \code{details}).}

	\item{departType }{Character string vector for the type of departure (i.e., last detection), with values \dQuote{\code{D}} for death and \dQuote{\code{C}} for censoring (see \code{details}).}

}

\details{
\subsection{1) Data structure:}

The function allows to calculate product limit estimator (Wann \emph{et al.} 1987) for data that includes the following types of records:
\itemize{ 
\item\bold{Uncensored}: individuals with known ages at death;

\item\bold{right-censored}: individuals last seen alive;

\item\bold{left-truncated}: individuals born before the start of the study and are truncated at the age of entry. 
}

The data required are the ages at last detection (i.e., \emph{uncensored} or \emph{right-censored}) passed through argument \dQuote{\code{ageLast}}, the type of departure via argument \dQuote{\code{departType}}, which takes two values, namely \dQuote{\code{D}} for death, and \dQuote{\code{C}} for censored (i.e., right-censored). 

In addition, if there is \emph{left-truncation}, it takes the ages at entry to the study by means of argument \dQuote{\code{ageFirst}}. If all the individuals were born during the study, the value of \dQuote{\code{ageFirst}} can be left as \code{NULL}, which will make them all equal to 0. 
 }

\value{
	\item{fill up }{FILL UP}
}

\author{Fernando Colchero \email{colchero@imada.sdu.dk}}

\seealso{FILL UP}

\examples{}

\keyword{models}

