\name{polytomous.poisson.reformulation}
\alias{polytomous.poisson.reformulation}
\title{
  A function fitting a polytomous logistic regression model based on
  the Poisson-reformulation heuristic.
}
\description{

  A function fitting a polytomous logistic regression model based on
  the Poisson-reformulation heuristic. With the Poisson-reformulation
  heuristic, the polytomous setting is reformulated using one of the
  functions \code{\link{instance2narrowcount}} or
  \code{\link{wide2narrowcount}} as counts of outcome-predictor
  combinations, for which a logistic regression model can be fit using
  the \code{\link{glm}} or \code{\link[lme4]{lmer}} functions with
  setting \code{family=poisson}. See \code{Details} for the further
  specifics of the Poisson reformulation.

}
\usage{
polytomous.poisson.reformulation(formula, data, frequency = NA,
   variables.ordered = NULL, include.Observation = TRUE, \dots)
}
\arguments{
  \item{formula}{

   an object of class \code{"formula"} (or one that can be coerced to
   that class): a symbolic description of the model to be fitted.  The
   details of model specification are given under \code{Details} for
   the \code{\link{polytomous}} function.

}
  \item{data}{

   a data frame (or object coercible by \code{as.data.frame} to a data
   frame) containing the variables specified in the model. The data
   may be represented either in the \code{"narrow"} or the \code{"wide"}
   format; see \code{Details} for the function
   \code{\link{polytomous}}.

}
  \item{frequency}{ 

    a numeric vector (or the name of a column in the input data frame)
    with the frequencies of the instances.  If absent (set to
    \code{NA}), each instance is assigned a frequency equal to 1.
}
  \item{variables.ordered}{

   a list specifying the order of the predictor variables in the count
   data table resulting from the Poisson reformulation; by default
   \code{NULL}, in which case the order will be alphabetical; passed
   on to \code{\link{instance2narrowcount}} or
   \code{\link{wide2narrowcount}}.

}
  \item{include.Observation}{ 

   A logical whether a factor index specifying unique predictor value
   combinations should be included in the count data table resulting
   from the Poisson reformulation; by default \code{TRUE}.

}
  \item{\dots}{

  Control variables to be passed on to other functions, see
  \code{\link{instance2narrowcount}}, \code{\link{wide2narrowcount}}
  or \code{\link{model.statistics}}.

}
}
\details{

  With the Poisson reformulation heuristic, the polytomous setting is
  reformulated using one of the functions
  \code{\link{instance2narrowcount}} or \code{\link{wide2narrowcount}}
  as counts of outcome-predictor combinations, for which a logistic
  regression model can be fit using the \code{\link{glm}} or
  \code{\link[lme4]{lmer}} functions with setting
  \code{family=poisson}.

  With the Poisson reformulation, the original regression
  \code{formula} estimating probabilities of polytomous outcomes based
  on combinations of predictors is transformed into a formula
  estimating counts of polytomous outcomes in combination with
  (effectively discrete) predictor values. Consequently, the original
  outcome term is turned into a predictor and the original predictor
  terms are re-expressed as interaction terms between outcomes and
  predictors. For example, the original formula:

  \code{outcome ~ predictor1 + predictor2}

  is transformed into the following Poisson formula (provided as
  \code{formula.poisson} in the results):

  \code{Count ~ outcome + outcome:predictor1 + outcome:predictor2 +
  Observation}

  The auxiliary predictor variable \code{Observation} is an index
  designating unique combinations of predictor values specified in the
  original formula.

  When fitting a mixed model with random variables that have large
  numbers of values with frequencies close or equal to 1, it may be
  advisable in order to get proper estimation of their impact to set
  \code{include.Observation=FALSE}, as unique predictor value
  combinations indicated by each \code{Observation} value may
  effectively be convergent with the values of such random variables.

}
\value{

   \code{polytomous.poisson.reformulation} returns an object of class
   \code{c("polytomous","one.vs.rest", model.type)} (where
   \code{model.type} is either \code{"fixed"} or \code{"mixed"}), a
   list containing the following components:

  \describe{

  \item{\code{model}}{a list containing the underlying outcome-specific
   binary models fitted using \code{\link{glm}} in the case of a
   \code{"fixed"} effects model or \code{\link[lme4]{lmer}} in case of a mixed
  effects model, with \code{family=poisson} in either case.}
  \item{\code{data}}{the originally supplied \code{data} argument object.}
  \item{\code{data.poisson}}{the Poisson reformulated count data table used to fit the model}
  \item{\code{frequency}}{the originally supplied \code{frequency} argument.}
  \item{\code{fitted.poisson}}{a two-column matrix of the fitted counts and associated probability estimates corresponding to the input Poisson count data table}
  \item{\code{fitted}}{a matrix of the fitted outcome-specific probability estimates corresponding to the original data table in \code{"instance"}, \code{"narrow"} or \code{"wide"} format.}
  \item{\code{coefficients}}{the original coefficients of the model}
  \item{\code{logodds}}{a matrix of the outcome-by-predictor logodds estimated for the model.}
  \item{\code{odds}}{a matrix of the outcome-by-predictor odds (\code{exp(logodds)}) estimated for the model.}
  \item{\code{p.values}}{a matrix of the estimates of the significances of the outcome-by-predictor logodds/odds estimated for the model.}
  \item{\code{statistics}}{a range of descriptive statistics describing the goodness of fit and classification performance of the model; see \code{model.statistics} and \code{crosstable.statistics}. For mixed models with random effects, \code{statistics} will also include their standard deviations (\code{sd.ranef}) and variances (\code{var.ranef}).}
  \item{\code{formula}}{the formula specification used to fit the model.}
  \item{\code{formula.poisson}}{the reformulated formula with counts as the response, used to fit the model.}
  \item{\code{outcomes}}{the outcome categories.}
  \item{\code{heuristic}}{the heuristic used to fit the model, being \code{"poisson.reformulation"}.}
  \item{\code{data.format}}{the format type of the \code{data} argument object; having the value of either \code{"instance"}, \code{"narrow"}, or \code{"wide"}.}
  }
}
\references{

Arppe, A. (in prep.) Solutions for fixed and mixed effects modeling of
polytomous outcome settings.

Faraway, J. J. (2006). Extending the Linear Model with R. Boca Raton:
Chapman & Hall/CRC, 101-102.

Venables, W. N. and B. D. Ripley (2002). Modern Applied Statistics with S
(4th edition). New York: Springer, 199-202.

}
\author{
Antti Arppe, with ideas from R. Harald Baayen and Terrance M. Nearey
}

\seealso{
\code{\link{instance2narrowcount}}, \code{\link{wide2narrowcount}}, \code{\link{model.statistics}},
\code{\link{polytomous}}, \code{\link{glm}}, \code{\link[lme4]{lmer}}
}
\examples{

data(think)
think.polytomous <- polytomous.poisson.reformulation(Lexeme ~ Agent + Patient,
   data=think)
think.polytomous$statistics
think.polytomous$odds
think.polytomous$p.values

\dontrun{
library(nnet)
think.multinom <- multinom(Lexeme ~ Agent + Patient, data=think)
exp(coef(think.multinom))
}

think.counts <- instance2narrowcount(think, c("Agent","Patient"), "Lexeme")
think.poisson <- glm(Count ~ Observation + Lexeme + Lexeme:Agent +
   Lexeme:Patient, data=think.counts, family=poisson)
summary(think.poisson)

\dontrun{
think.lmer <- polytomous.poisson.reformulation(Lexeme ~ Agent + Patient +
   (1|Section), data=think)
summary(think.lmer)
ranef(think.lmer)
}

}

\keyword{ regression multivariate category classif }

