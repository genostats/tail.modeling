\name{PWM}
\alias{PWM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimation of GPD's parameters by Probability Weighted Moments
%%  ~~function to do ... ~~
}
\description{ Estimates the two parameters k and a of Pareto's cumulative distribution function on a given data set z using the probability weighted moments.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PWM(z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{ Data set - sorted list of real numbers.
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns a list of the estimators k and a of the parameters, as real numbers.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{ J. R. M. Hosking and J. R. Wallis, Parameter and Quantile Estimation for the Generalized Pareto Distribution, Technometrics. Available at http://www.jstor.org/stable/1269343.
%% ~put references to the literature/web site here ~
}
\author{ Marion
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
PWM(sort(rnorm(1e6)))

## The function is currently defined as
function (z)
{
  n <- length(z)
  p <- (seq.int(n)-0.35)/n
  t <- mean( (1-p)*z )
  mu <- mean(z)
  k <- mu/(mu-2*t)-2
  a <- 2*mu*t/(mu-2*t)
  list(a = a, k = k)

  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
