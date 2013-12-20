\name{instance2narrowcount}
\alias{instance2narrowcount}
\title{

Transformation of an uncounted instance-specific data table into a count data
table in the\code{"narrow"} format

}
\description{

Transforms a data table with uncounted instance-by-instance
information on the co-occurrences of individual outcomes with various
predictor values into a count table in the \code{"narrow"} format,
with frequency counts for the outcomes in conjunction with unique
combinations of predictor variable values.

}
\usage{
instance2narrowcount(data.table, variables, outcome = "OUTCOME",
   variables.default = NULL, outcome.ordered = NULL,
   numeric2discrete = function(x) cut2(x,levels.mean=TRUE,g=g.numeric),
   g.numeric = 2)
}
\arguments{
  \item{data.table}{

  a data table with instance-by-instance information on the occurrence
  of individual outcome variable variables in conjunction with
  specific values of predictor variables.

}
  \item{variables}{

  a list of the names of predictor variables (columns in
  \code{data.table}) to be included in the creation of the
  \code{count.table}.

}
  \item{outcome}{

  a name of the outcome variable (column in \code{data.table}); by
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
  original order specified in the object referred to by the
  \code{data} argument will be used.

}
  \item{numeric2discrete}{

  a function to transform a continuous numeric predictor into a
  discrete set of numeric values, by default \code{cut2} from the
  \code{Hmisc} package with the preset parameters
  \code{levels.mean=TRUE} and \code{g=g.numeric} (by default =2). If
  set to NULL, each value of each numeric predictor will be treated as
  a discrete value of its own.

}
  \item{g.numeric}{

  a parameter to be passed to the \code{numeric2discrete} function
  (parameter \code{g} for \code{Hmisc::cut2(..., g=g.numeric, ...)},
  or a user-defined function), determining the desired number of
  values for each numeric predictor; by default equal to
  2.

}
}
\details{

Transforms a data table with uncounted instance-by-instance
information on the co-occurrences of individual outcomes with various
predictor values into a count table in the \code{"narrow"} count
format, with frequency counts for the outcomes in conjunction with
unique combinations of predictor variable values.

Note that numeric variables will remain numeric despite the reduction
of their distinct values using the \code{numeric2discrete} function.

}
\value{
A count data table with the frequency counts for each unique
combination of outcomes and predictor variable values. In addition to
columns with values for each included predictor, the count data table
has the following common columns: 

\describe{

  \item{\code{"Proportion"}}{the relative proportion of the specific
  outcome in conjunction with the specific combination of selected
  predictor variables (in relation to all the outcomes for the
  particular unique combination of predictor variables.}

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

Antti. A. (in prep.) Solutions for fixed and mixed effects modeling of
polytomous outcome settings.

}
\author{
Antti Arppe
}
\seealso{
\code{\link{polytomous.poisson.reformulation}}, \code{\link{wide2narrowcount}}
}
\examples{

data(think)
think.counts <- instance2narrowcount(think, c("Agent","Patient"), "Lexeme")
think.counts

think.poisson <- glm(Count ~ Observation + Lexeme + Lexeme:Agent + Lexeme:Patient,
   data=think.counts, family=poisson)
summary(think.poisson)

think.polytomous.poisson1 <- polytomous(Lexeme ~ Agent + Patient, data=think.counts,
   frequency="Count", heuristic="poisson.reformulation")
summary(think.polytomous.poisson1)

think.polytomous.poisson2 <- polytomous(Lexeme ~ Agent + Patient, data=think,
   heuristic="poisson.reformulation")
summary(think.polytomous.poisson2)

\dontrun{

library(lme4)
think.counts2 <- instance2narrowcount(think, c("Agent","Patient","Section"), "Lexeme")
think.poisson.lmer <- lmer(Count ~ (1|Observation) + Lexeme + Lexeme:Agent +
   Lexeme:Patient + (1|Section), data=think.counts2, family=poisson)
summary(think.poisson.lmer)

think.polytomous.lmer <- polytomous(Lexeme ~ Agent + Patient + (1|Section), data=think)
summary(think.polytomous.lmer)

}

}

\keyword{ regression multivariate category }

