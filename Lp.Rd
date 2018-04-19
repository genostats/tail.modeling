\name{Lp}
\alias{Lp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Loglikelihood of Pareto's cumulative distribution function
%%  ~~function to do ... ~~
}
\description{ Computes the loglikelihood of Pareto's cumulative distribution function for a given data set z.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Lp(u,z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{ Real number
%%     ~~Describe \code{x} here~~
}
  \item{z}{ Data set - list of real numbers
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns the value of the loglikelihood.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
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
Lp(0,sort(rnorm(1e6)))

## The function is currently defined as
function (u,z)
{
  n <- length(z)
  if(u == 0)
    return( -n*(log(mean(z)) + 1) )
  k <- k_func(u,z)
  n*(log(u/k) + k - 1)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
