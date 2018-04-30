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
PGPD(x0,zexc,y,seuil,estim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x0}{ Test statistic of the data set - real number
%%     ~~Describe \code{x} here~~
}
  \item{zexc}{ Data set which results from the initial one y - sorted list of real numbers
%%     ~~Describe \code{x} here~~
}
  \item{y}{ Data set - list of real numbers
%%     ~~Describe \code{x} here~~
}
\item{seuil}{ Threshold used to create the data set zexc composed of the values of y higher or equal to (x0 - seuil) - real number
%%     ~~Describe \code{x} here~~
}
\item{estim}{ Method to estimate the parameters of Pareto's function - String: either "EMV" for the method of maximum likelihood or "PWM" for the method of probability weighted moments
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns the value of the estimated p-value.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{ Theo A. Knijnenburg, Lodewyk F. A. Wessels, Marcel J. T. Reinders
and Ilya Shmulevich, Fewer permutations, more accurate P-values, Bioinformatics.
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
PGPD(3,tail(sort(rnorm(1e6)),500),rnorm(1e6),2,"PWM")

## The function is currently defined as
function (x0,zexc,y,seuil,estim)
{
  M<-nb_exc(x0,y)
  if(M >= 10)
    return(M/length(y))
  else
    return(length(zexc)/length(y)*(1-FGPD(x0-seuil,zexc,estim)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
