\name{BoxCox}
\alias{BoxCox}

\title{ Box-Cox function

}
\description{ Computes the Box-Cox function evaluated in a real number z.

}
\usage{
BoxCox(z,lambda)
}

\arguments{
  \item{z}{ Real number - can also be a list of real numbers

}
\item{lambda}{ Box-Cox parameter - real number

}
}

\value{ Returns the value of the function evaluated in z with the parameter lambda chosen.

}

\author{ Marion
%%  ~~who you are~~
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
