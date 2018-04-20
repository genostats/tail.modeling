\name{FGPD}
\alias{FGPD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Pareto's cumulative distribution function
%%  ~~function to do ... ~~
}
\description{ Computes the value of the distribution function of Pareto. The method of parameters' estimation can be chosen.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FGPD(z,zexc,estim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{ Real number in which the function will be evaluated - can also be a list of real numbers
%%     ~~Describe \code{x} here~~
}
  \item{zexc}{ Data set - sorted list of real numbers
%%     ~~Describe \code{x} here~~
}
  \item{estim}{ Method to estimate the parameters of Pareto's function - String: either "EMV" for the method of maximum likelihood or "PWM" for the method of probability weighted moments
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns the value of the cumulative distribution function evaluated in z using the method of estimation chosen.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Theo A. Knijnenburg, Lodewyk F. A. Wessels, Marcel J. T. Reinders
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
FGPD(3,sort(rnorm(1e6)),"PWM")

## The function is currently defined as
function (z, zexc, estim)
{
  if (estim == "EMV")
    coeff<-EMV(zexc)
  else
    coeff<-PWM(zexc)
  a<-coeff$a
  k<-coeff$k
  if (k!= 0)
    return(1-(1-k*z/a)^(1/k))
  else
    return(1-exp(-z/a))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
