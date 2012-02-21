\name{multinomial2logical}
\alias{multinomial2logical}
\title{

Transformation of a data frame with multinomial variable columns into
a \code{"logical"} format

}
\description{

  Transforms a data frame (or an object that can be coerced in such)
  consisting of multinomial variables as columns into a logical
  format, where each unique value of each multinomial variable is
  represented by a corresponding logical variable.

}
\usage{multinomial2logical(data, outcome=NULL, variables=NULL)}
\arguments{ \item{data}{

  a data frame (or an object that can be coerced into such) with
  columns of multinomial variable.s

}
  \item{outcome}{

  the name of the outcome variable (column in \code{data}) that will
  be retained as a multinomial variable column in the resultant data
  frame; by default \code{NULL}.

}
  \item{variables}{

  a list of the names of predictor variables (columns in \code{data})
  to be included in the creation of the resultant data frame.

}
}
\details{

  Transforms a data frame (or an object that can be coerced in such)
  consisting of multinomial variables as columns into a logical
  format, where each unique value of each multinomial variable is
  represented by a corresponding logical variable. The logical format
  allows for an easy specification of individual variable values as
  predictors in the formula of a polytomous logistic regression model.

}
\value{

  A data frame with logical variables representing each unique value
  of all (or selected) multinomial variable columns in
  \code{data}. When such a value of a multinomial variable is present
  on an observation row in the original \code{data}, the value of the
  corresponding logical variable is \code{TRUE}, whereas otherwise it
  is \code{FALSE}. Each logical variable is named as a combination of
  the corresponding original multinomial variable (column) name and
  variable value, separated by a period (\sQuote{.}). For instance, the
  logical variable corresponding to any occurrence of \code{"Group"}
  class of the variable \code{"Agent"} in the \code{think} dataset is
  named \code{"Agent.Group"}.

  If an \code{outcome} column is specified, the resultant data frame
  will have as its first column the corresponding multinomial
  variable; this is necessary for using the resultant data frame as
  the input data for \code{\link{polytomous}}.

}

\references{

Arppe, A. (2008) Univariate, bivariate and multivariate methods in
corpus-based lexicography -- a study of synonymy. Publications of the
Department of General Linguistics, University of Helsinki,
No. 44. URN: http://urn.fi/URN:ISBN:978-952-10-5175-3.

Arppe, A. (in prep.)  Solutions for fixed and mixed effects modeling of
polytomous outcome settings.

}
\author{
Antti Arppe
}
\seealso{
\code{\link{polytomous}}, \code{\link{polytomous.one.vs.rest}},
\code{\link{polytomous.poisson.reformulation}},
\code{\link{instance2narrowcount}}, \code{\link{wide2narrowcount}}
}
\examples{

data(think)
think.logical <- multinomial2logical(think, outcome="Lexeme",
   variables=c("Agent","Patient"))
think.polytomous1 <- polytomous(Lexeme ~ Agent.Individual + Agent.Group +
  Patient.Abstraction + Patient.Activity, data=think.logical)
summary(think.polytomous1)

think.polytomous2 <- polytomous(Lexeme ~ Agent.Individual + Agent.Group +
   Patient.Abstraction + Patient.Activity, data=think.logical,
   heuristic="poisson.reformulation",
   outcome.ordered=c("harkita","miettia","pohtia","ajatella"))
summary(think.polytomous2)

}
\keyword{ regression }

