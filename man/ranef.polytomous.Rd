\name{ranef.polytomous}
\alias{ranef.polytomous}
\alias{ranef}
\title{Extract the modes of the random effects for polytomous objects}
\description{

  A function to extract the conditional modes of the random effects
  from a fitted mixed-effects\code{"polytomous"} model object.  For
  linear mixed models the conditional modes of the random effects are
  also the conditional means.

}
\usage{
  \method{ranef}{polytomous}(object, \dots)
}
\arguments{
  \item{object}{

    an object of a class of fitted models with random effects,
    typically an \code{"polytomous"} object of the type
    \code{"mixed"}.

  }
  \item{\dots}{

    additional control arguments to be passed on to the underlying
    \code{\link[lme4]{ranef}} function; see \code{\link[lme4]{ranef}}
    (\code{postVar}, code{drop}, \code{whichel}).

  }
}
\value{

  A list of data frames, one for each grouping factor for the random
  effects.  The number of rows in the data frame is the number of
  levels of the grouping factor.  The number of columns is the
  dimension of the random effect associated with each level of the
  factor.

}
\details{

  If grouping factor \code{i} has \code{k} levels and \code{j} random
  effects per level the \emph{i'th} component of the list returned by
  \code{ranef} is a data frame with \code{k} rows and \code{j}
  columns. The \emph{k'th} face of this array is a positive definite
  symmetric \code{j} by \code{j} matrix.  If there is only one
  grouping factor in the model the variance-covariance matrix for the
  entire random effects vector, conditional on the estimates of the
  model parameters and on the data will be block diagonal and this
  \code{j} by \code{j} matrix is the \emph{k'th} diagonal block.

}

\author{
Antti Arppe
}
\seealso{
  \code{\link{polytomous.poisson.reformulation}}, \code{\link[lme4]{ranef}}
}
\examples{

\dontrun{
data(think)
think.lmer <- polytomous(Lexeme ~ Agent + Patient + (1|Section), data=think,
   heuristic="poisson.reformulation")
summary(think.lmer)
ranef(think.lmer)

think.lmer2 <- polytomous(Lexeme ~ Agent + Patient + (1|Section) + (1|Author),
   data=think, heuristic="poisson.reformulation")
summary(think.lmer)
ranef(think.lmer)
}

}
\keyword{ regression }
