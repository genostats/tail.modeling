\name{calcul_p}
\alias{calcul_p}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of the p-value with Pareto's function or Box-Cox function used on distribution's tail.
%%  ~~function to do ... ~~
}
\description{ Estimates the p-value of a given data set zsim with the test statistics Zobs thanks to Pareto's method or Box-Cox method with their different estimated parameters.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calcul_p(zsim,Ntail=500,estim=NA,Zobs,param=NA,method)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zsim}{ Data set -  list of real numbers
%%     ~~Describe \code{x} here~~
}
\item{Ntail}{ Length of the tail of the data set taken - integer
%%     ~~Describe \code{x} here~~
}
\item{estim}{ Method to estimate the parameters of Pareto's function - String: Default value is NA (when the selected method isn't Pareto's). If the selected method is the Pareto's one, estim takes either the value "EMV" for the method of maximum likelihood or the "PWM" for the method of probability weighted moments
%%     ~~Describe \code{x} here~~
}
\item{Zobs}{ Test statistic of the data set  - real number
%%     ~~Describe \code{x} here~~
}
\item{param}{ Box-cox parameter lambda - real number: Default value is NA that means that the parameter will be estimated with least squares; if it is remplaced by a real number, this value will be used for lambda without performing any estimation
%%     ~~Describe \code{x} here~~
}
\item{method}{ Method chosen to estimate the p-value  - String: either "GPD" for Pareto's method or "BC" for Box-Cox's method.
%%     ~~Describe \code{x} here~~
}
}
\details{ Both methods are applied on the distribution's tail of the data set.
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns a list composed of the estimated p-value and the parameter(s) of the selected method of estimation.
If the selected method is "BC", it returns a list composed of:
 \item{Pbc_z }{The estimated p-value with Box-Cox function -  real number}
 \item{lambda }{The estimated parameter lambda (or the lambda given if not estimated) - real number}

If the selected method is "GPD", it returns a list composed of:
  \item{Pgpd}{The estimated p-value with Pareto's function -  real number}
  \item{k }{The estimated parameter k of Pareto's cumulative distribution function - real number}
  \item{a }{The estimated parameter a of Pareto's cumulative distribution function - real number}

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
calcul_p(zsim=rnorm(1e6),Zobs=5,method="BC")

## The function is currently defined as
function (zsim,Ntail=500,estim=NA,Zobs,param=NA,method)
{
# les Ntail plus grandes valeurs (les dernieres)
  z1 <- tail( sort(zsim) , Ntail )
  #seuil pour la GDP
  t<-(z1[1] + head(sort(zsim,decreasing =TRUE),Ntail+1)[Ntail+1])/2
  #calcul des excedents de la GDP, ceux qui lui sont superieurs
  zgpd<-z1-t
  zgpd<-zgpd[zgpd>0] #uniquement ceux superieurs au seuil

  if (method == "BC")
    return(list(Pbc_z = PBC_Z(z1,length(zsim),Ntail,param,Zobs)$p,
                lbda = PBC_Z(z1,length(zsim),Ntail,param,Zobs)$lbda))

  if (method =="GPD"){
    if (estim == "PWM")
      return(list(Pgpd=PGPD(Zobs,zgpd,zsim,t,estim),
                  a=PWM(zgpd)$a,
                  k=PWM(zgpd)$k))
    else
      return(list(Pgpd=PGPD(Zobs,zgpd,zsim,t,estim),
                  a=EMV(zgpd)$a,
                  k=EMV(zgpd)$k))
}
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
