\name{k_func}
\alias{k_func}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Useful function for the estimation of GPD's parameters with a given data set by Maximum Likelihood
%%  ~~function to do ... ~~
}
\description{ Performs a variable change needed to estimate Pareto's parameters with a given data set z by Maximum Likelihood.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
k_func(u,z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{ Real number on which the variable change will be done.
%%     ~~Describe \code{x} here~~
}
  \item{z}{ Data set - list of real numbers.
%%     ~~Describe \code{x} here~~
}
}

\value{ Returns the new variable k which results from the variable change of u.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{ Marion
%%  ~~who you are~~
}

\examples{
k_func(0,sort(rnorm(1e6)))

## The function is currently defined as
function (u,z)
{
  -mean(log(1-u*z))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
