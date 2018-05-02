\name{BoxCox_lm}
\alias{BoxCox_lm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Box-Cox parameter lambda's estimation with Least Squares
%%  ~~function to do ... ~~
}
\description{Estimates the parameter lamnda of the Box-Cox function thanks to least squares.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BoxCox_lm(Y,X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{ List of real numbers
%%     ~~Describe \code{x} here~~
}
\item{X}{ List of real numbers
%%     ~~Describe \code{x} here~~
}
}
\details{ X and Y must have the same dimension.
%%  ~~ If necessary, more details than the description above ~~
}
\value{Returns the estimated parameter lambda with the method of least squares applied to the linear regression between Y and X.
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
BoxCox_lm(log(-log(seq(500,1)/1e6)),log(tail(sort(rlnorm(1e6)),500)))

## The function is currently defined as
function (x)
{
  f <- function(lambda) {
    if(lambda == 0)
      X <- log(X)
    else
      X <- (X**lambda - 1)/lambda

    b <-  cov(X,Y) / var(X)
    a <- mean(Y) - b*mean(X)
    sum( (Y - a - b*X)**2 )
  }
  optimize(f,c(0,1))$minimum
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
