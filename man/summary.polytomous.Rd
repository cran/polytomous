\name{summary.polytomous}
\alias{summary.polytomous}
\alias{print.summary.polytomous}

\title{A summary of a Polytomous Logistic Regression model
}

\description{
  A summarization method for an object of the class \code{"polytomous"}.
}

\usage{
\method{summary}{polytomous}(object, \dots)

\method{print}{summary.polytomous}(x, digits = max(3, getOption("digits") - 3),
   parameter="odds", max.parameter=ifelse(parameter=="odds",10000,100),
   p.critical=.05, max.print=10, cycles=0, max.denominator=0, \dots)
}
\arguments{
\item{object}{
      An object of class \code{"polytomous"}, resulting from a call to
      \code{\link{polytomous}}.
     }
\item{x}{An object of class \code{"summary.polytomous"}, usually resulting from a
     call to \code{\link{summary.polytomous}}.
     }
\item{digits}{
    The number of significant digits to use when printing.
    }
\item{parameter}{
    The set of parameters to output in printing; by default
    \code{"odds"}, alternatively \code{"logodds"}.
    }
\item{max.parameter}{
    A value specifying an upper limit beyond which exceptionally large
    (and potentially unreliable and mostly nonsignificant) parameter
    values will be output as being beyond the scale, as \code{"Inf",
    "1/Inf", "-Inf"}; see \code{Details} for the default limit values
    and specific outputs for the different \code{parameter} types.
    }
\item{p.critical}{
   The critical P-value for considering a parameter value significant;
   by default set to \code{P=.05} as is common in the humanities.
   }
\item{max.print}{
    The maximum number of rows of the \code{parameter} to be output
    when printing with \code{print.summary.polytomous}; by default
    set to 10; if set to \code{NA} all rows will be output.
  }
\item{cycles}{
    A value to be passed on to the \code{fractions} function in the
    \code{MASS} package for representing \code{odds} as fractions;
    however, fractional representation will not result with the
    default value \code{=0}.
    }
\item{max.denominator}{
    A value to be passed on to the \code{fractions} function in the
    \code{MASS} package for representing \code{odds} as fractions;
    however, fractional representation will not result with the
    default value \code{=0}.
    }
\item{\dots}{further arguments passed to or from other methods.
    }
}

\details{
  Calculates descriptive statistics of a fitted Polytomous Logistic
  Regression model and prints a nice summary of the key results.

  Parameters for which the respective P-value is greater than the
  critical threshold value (set with \code{p.critical} by default as
  \code{P=.05}), i.e. which would be considered as not significant,
  are output within parentheses, e.g. \code{"(1.1)"}.

  For \code{parameter="odds"}, the default maximum output limit value is
  \code{max.parameter=10000}.  With this default value, \code{"Inf"}
  will be output if \code{odds > 10000} and \code{"1/Inf"} if
  \code{odds < 1/10000}.  For \code{parameter="logodds"}, the default
  maximum output limit value is \code{max.parameter=100}.  With this default
  value, \code{"Inf"} will be output if \code{logodds > 100} and
  \code{"-Inf"} if \code{logodds < -100}.
}

\value{
  \code{summary.polytomous} returns an object of the class
  \code{"summary.polytomous"}, a list with the following components:
\describe{
  \item{\code{formula}}{The formula specified for the \code{"polytomous"} object.
  }
  \item{\code{odds, logodds, p.values}}{The estimated odds, logodds, and their P-values
  }
  \item{\code{statistics}}{A range of descriptive statistics calculated with \code{\link{model.statistics}}.
  }
  }
}

\references{
Antti. A. (in prep.) Solutions for fixed and mixed effects modeling of
polytomous outcome settings.
}

\author{
Antti Arppe
}

\seealso{
\code{\link{polytomous}}, \code{\link{model.statistics}}
}

\examples{

data(think)
think.polytomous <- polytomous(Lexeme ~ Agent + Patient, data=think)
print(summary(think.polytomous), digits=2, parameter="odds")
print(summary(think.polytomous), digits=4, parameter="logodds")

## For more examples see examples(polytomous).

}

\keyword{ regression }
