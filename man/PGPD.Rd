\name{PGPD}
\alias{PGPD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimation of the p-value using Pareto's distribution function
%%  ~~function to do ... ~~
}
\description{ Estimates the p-value of a data set y for a test statistic x0 with the Pareto's distribution function.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PGPD(x0,zexc,Nperm,seuil,estim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x0}{ Test statistic of the data set - real number
%%     ~~Describe \code{x} here~~
}
  \item{zexc}{ Data set which results from the initial one y - sorted list of real numbers
%%     ~~Describe \code{x} here~~
}
  \item{Nperm}{ Number of permutations - integer
%%     ~~Describe \code{x} here~~
}
\item{seuil}{ Threshold used to create the data set zexc composed of the values of y higher or equal to (x0 - seuil) - real number
%%     ~~Describe \code{x} here~~
}
\item{estim}{ Method to estimate the parameters of Pareto's function - String: either "EMV" for the method of maximum likelihood or "PWM" for the method of probability weighted moments
%%     ~~Describe \code{x} here~~
}
}

\value{ Returns a list composed of:
\item{p}{The value of the estimated p-value. -  real number}
  \item{k }{The estimated parameter k of Pareto's cumulative distribution function if it was used to compute the p-value - real number}
  \item{a }{The estimated parameter a of Pareto's cumulative distribution function if it was used to compute the p-value - real number}

}
\references{ Theo A. Knijnenburg, Lodewyk F. A. Wessels, Marcel J. T. Reinders
and Ilya Shmulevich, Fewer permutations, more accurate P-values, Bioinformatics.
%% ~put references to the literature/web site here ~
}
\author{ Marion
%%  ~~who you are~~
}

\examples{
PGPD(3,tail(sort(rnorm(1e6)),500),1e6,2,"PWM")

## The function is currently defined as
function (x0,zexc,Nperm,seuil,estim)
{
  M<-nb_exc(x0-seuil,zexc)

  if(M >= 10)
    return(list( p = M/Nperm,
                 k = NA,
                 a = NA ))
  else {
    result<-FGPD(x0-seuil,zexc,estim)
    return(list( p = length(zexc)/Nperm*(1-result$val),
                 k = result$k,
                 a = result$a ))
    }
}
}
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
