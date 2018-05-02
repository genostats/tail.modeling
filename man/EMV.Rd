\name{EMV}
\alias{EMV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of GPD's parameters by Maximum Likelihood
%%  ~~function to do ... ~~
}
\description{Estimates the two parameters k and a of Pareto's cumulative distribution function on a given data set z using the maximum of likelihood.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
EMV(z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{ Data set - sorted list of real numbers.
%%     ~~Describe \code{x} here~~
}
}

\value{Returns a list of the estimators k and a of the parameters, as real numbers.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{R. M. Hosking and J. R. Wallis, Parameter and Quantile Estimation for the Generalized Pareto Distribution, Technometrics. Available at http://www.jstor.org/stable/1269343.
%% ~put references to the literature/web site here ~
}
\author{ Marion
%%  ~~who you are~~
}

\examples{
EMV(tail(sort(rnorm(1e6)),500))

## The function is currently defined as
function (z)
{
  u <- optimize(function(u) Lp(u,z), c(-10, min(1/z)), maximum=TRUE )$maximum
  if(u == 0)
    return(list( a = mean(z), k = 0))
  k <- k_func(u, z)
  list(a = k/u, k = k)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
