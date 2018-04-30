\name{BoxCox}
\alias{BoxCox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Box-Cox function
%%  ~~function to do ... ~~
}
\description{ Computes the Box-Cox function evaluated in a real number z.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BoxCox(z,lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{ Real number - can also be a list of real numbers
%%     ~~Describe \code{x} here~~
}
\item{lambda}{ Box-Cox parameter - real number
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns the value of the function evaluated in z with the parameter lambda chosen.
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
BoxCox(56,1)

## The function is currently defined as
function (x)
{
  if (lambda != 0)
    return((z^lambda - 1)/lambda)
  else
    return(log(z))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
