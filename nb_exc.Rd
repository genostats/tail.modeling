\name{nb_exc}
\alias{nb_exc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Values higher or equal to Zobs
%%  ~~function to do ... ~~
}
\description{ Computes how many values of the data set z are higher or equal to the test statistic x0.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
nb_exc(x0,z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x0}{ the test statistic - real number.
%%     ~~Describe \code{x} here~~
}
  \item{z}{ the data set - list of real numbers.
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns an integer M which is the number of values of z higher or equal to x0.
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
nb_exc(5,1:10)

## The function is currently defined as
function (x0,z)
{
M<-length(z[z>=x0])
return(M)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
