\name{PBC_Z}
\alias{PBC_Z}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimation of the p-value using linear regression with Box-Cox transformation applied on test statistics
%%  ~~function to do ... ~~
}
\description{Estimates the p-value of a data set z for a test statistic Zobs with the Box-Cox function.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PBC_Z(z,Ntot,N,param,Zobs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{ Data set resulting from the initial one- sorted list of real numbers
%%     ~~Describe \code{x} here~~
}
\item{Ntot}{ Number of observations in the initial data set  - integer
%%     ~~Describe \code{x} here~~
}
\item{N}{ Number of observations in the data set z - integer
%%     ~~Describe \code{x} here~~
}
\item{param}{ Box-cox parameter lambda - real number: if param is missing, the parameter will be estimated with least squares; if it is a real number, this value will be used for lambda without performing any estimation
%%     ~~Describe \code{x} here~~
}
\item{Zobs}{ Test statistic of the data set  - real number
%%     ~~Describe \code{x} here~~
}
}
\details{ The method is to compute the coefficients of the linear regression between log(-log(p)) and BoxCox(log(z)) where p is a list of probabilities correspondings to the sorted list of test statistics z.
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns a list composed of:
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{p }{The value of the estimated p-value - real number}
  \item{interc }{The intercept of the linear regression used to estimate the p-value - real number}
  \item{pente }{The slope of the linear regression used to estimate the p-value - real number}
  \item{lbda }{ The estimated parameter lambda (or the lambda given if not estimated) - real number}

%% ...
}

\author{ Marion
%%  ~~who you are~~
}

\examples{
PBC_Z(z = tail(sort(rnorm(1e6)),500),Ntot = 1e6,N = 500,Zobs = 3)

## The function is currently defined as
function (x)
{
  p<-seq(N,1)/Ntot
  if (missing(param))
    lambda <- BoxCox_lm(log(-log(p)),log(z))
  else
    lambda <- param
  coeffs <- lm(  log(-log(p)) ~ BoxCox(log(z),lambda) )$coefficients
  return(list(p = exp(-exp(sum( coeffs*c(1, BoxCox(log(Zobs),lambda))))),
              interc = coeffs[[1]],
              pente = coeffs[[2]],
              lbda = lambda))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
