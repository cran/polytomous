\name{wide2narrowcount}
\alias{wide2narrowcount}
\title{

Transformation of a count data table from \code{"wide"} to
\code{"narrow"} count format

}
\description{

Transforms a count data table in the \code{"wide"} into the
\code{"narrow"} format, so that a polytomous logistic regression model
can be fit with \code{heuristic="poisson.reformulation"} using
\code{glm} or \code{lmer} with \code{family=poisson} for count data.

}
\usage{
wide2narrowcount(data.table, variables, outcomes, outcome = "OUTCOME",
   variables.default = NULL, outcome.ordered = NULL)
}
\arguments{
  \item{data.table}{
  A data table in the \code{"wide"} format which contains a column for
  each of the indicated outcomes giving their frequency of occurrence
  for a combination of predictor values indicated on the same row in
  the \code{data.table}.

}
  \item{variables}{
  The predictor variables (columns in \code{data.table}) to be
  included in the \code{count.table}.

}
  \item{outcomes}{

  The outcome variables (columns in \code{data.table}) with the
  frequencies of the outcomes for the associated predictor variable
  value combinations.

}
  \item{outcome}{

  A character string designating a name for the outcome variable; by
  default \code{"OUTCOME"}.

}
  \item{variables.default}{

  a list indicating for selected categorical predictors the value(s)
  that should be designated as the default/reference levels; by
  default \code{NULL}, in which case the original default/reference
  levels as specified for predictors in the object referred to by the
  \code{data.table} argument will be used.

}
  \item{outcome.ordered}{

  a list specifying the order of the categories for the
  outcome/response variable; by default \code{NULL}, in which case the
  original order specified with \code{outcomes} argument will be used.

}
}
\details{

Transforms a count data table in the \code{"wide"} format into the
\code{"narrow"} format, so that a polytomous logistic regression model
can be fit with \code{heuristic="poisson.reformulation"} using
code{glm} or \code{lmer} with \code{family=poisson}.

}
\value{

A count data table with the frequency counts for each unique
combination of outcomes and predictor variable values. In addition to
columns with values for each included predictor, the count data table
has the following common columns:

\describe{

  \item{\code{"Proportion"}}{the relative proportion of the specific
  outcome in conjunction with the specific combination of selected
  predictor variables (in relation to the sum frequency of all the
  outcomes for the particular unique combination of predictor
  variables).}

  \item{\code{"Count"}}{the frequency count of the specific outcome
  value in conjunction with the specific combination of selected
  predictor variables.}

  \item{\code{outcome}}{the name of the response variable designated by the
  \code{outcome} argument; by default \code{"OUTCOME"}}

  \item{\code{"Observation"}}{the index number for each unique
  combination of values of selected predictor variables.}

  }
}

\references{

Arppe, A. (in prep.)

}
\author{
Antti Arppe
}
\seealso{
\code{\link{polytomous.poisson.reformulation}}, \code{\link{instance2narrowcount}}
}
\examples{

data(think)
think.Agent_Patient.counts <- instance2narrowcount(think, c("Agent","Patient"),
   "Lexeme")
think.Agent_Patient.wide <- cbind(matrix(think.Agent_Patient.counts$Count,,4,
   byrow=TRUE, dimnames=list(NULL,c("ajatella","harkita","miettia","pohtia"))),
   unique(think.Agent_Patient.counts[c("Agent","Patient")]))
think.Agent_Patient.wide

think.Agent_Patient.counts2 <- wide2narrowcount(think.Agent_Patient.wide,
   variables=c("Agent","Patient"),
   outcomes=c("ajatella","harkita","miettia","pohtia"), outcome="Lexeme")
think.Agent_Patient.counts2
identical(think.Agent_Patient.counts,think.Agent_Patient.counts2)

think.polytomous1 <- polytomous(ajatella|harkita|miettia|pohtia ~ Agent + Patient,
   data=think.Agent_Patient.wide, heuristic="one.vs.rest")
summary(think.polytomous1)

think.polytomous2 <- polytomous(ajatella|harkita|miettia|pohtia ~ Agent + Patient,
   data=think.Agent_Patient.wide, heuristic="poisson.reformulation",
   outcome="Lexeme")
summary(think.polytomous2)

}

\keyword{ regression }

