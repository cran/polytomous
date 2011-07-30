\name{polytomous}
\alias{polytomous}
\alias{print.polytomous}
\title{

  Fitting polytomous logistic regression models for fixed or mixed
  effects predictors.

}
\description{

  \code{polytomous} is a top-level function that is used to select one
  of the heuristics for fitting a polytomous logistic regression model
  for fixed or mixed effects predictors; the arguments supplied to
  \code{polytomous} are passed on to the appropriate function
  implementing the indicated heuristic.

}
\usage{
polytomous(formula, data, heuristic = "one.vs.rest", \dots)

\method{print}{polytomous}(x, max.print = 10, \dots)
}
\arguments{
  \item{formula}{

   an object of class \code{"formula"} (or one that can be coerced to
   that class): a symbolic description of the model to be
   fitted.  The details of model formula specification are given under
   \code{Details}.
}
  \item{data}{

   a data frame (or an object coercible by \code{as.data.frame} to a
   data frame) containing the variables specified in the model. The
   data may be represented either in the uncounted \code{"instance"}
   or the \code{"wide"} count format; see \code{Details}.

}
  \item{heuristic}{

   the heuristic to be used in fitting the model; by default
   \code{"one.vs.rest"}; currently the only alternative is
   \code{"poisson.reformulation"}. Fitting mixed-effects models will
   require using the \code{"poisson.reformulation"} heuristic.

}
  \item{x}{

    An object of the class \code{"polytomous"} fitted with
    \code{polytomous} to be printed with \code{print.polytomous}.

}

  \item{max.print}{

    The maximum number of rows of the \code{parameter} to be output
    when printing with \code{print.summary.polytomous}; by default
    equal to 10.

}
  \item{\dots}{

  further arguments passed to or from other methods.

}
}
\details{

     A typical predictor has the form \code{response ~ terms} where
     \code{response} is the name of the factorial response vector
     indicating the possible polytomous outcomes (a column in the
     \code{data} argument) and \code{terms} is a series of terms which
     specifies a linear predictor for \code{response}. In this case,
     the object referred to by the \code{data} argument is expected to
     be in the \code{"narrow"} format, with one row for each observation
     of an outcome and predictors. Then, the predictors may either be
     multinomial factors or binary variables (with the values
     \code{TRUE, FALSE}).
     
     Alternatively, \code{response} can consist of the names of the
     individual outcome values separated by \code{"|"},
     e.g. \code{response1|response2|response3 ~ \dots}.  In such a
     case, the object referred to by \code{data} argument is expected
     to be in the \code{"wide"} format so that is contains a column
     for each of the indicated outcomes giving their frequency of
     occurrence for a combination of predictor values indicated on the
     same row in the \code{data} argument.
     
     A terms specification of the form \code{first + second} indicates
     all the terms in \code{first} together with all the terms in
     \code{second} with any duplicates removed.

     A specification of the form \code{first:second} indicates the the
     set of terms obtained by taking the interactions of all terms in
     \code{first} with all terms in \code{second}.  The specification
     \code{first*second} indicates the _cross_ of \code{first} and
     \code{second}.  This is the same as \code{first + second +
     first:second}.

     A specification of the form \code{fixed|random} indicates a fixed
     effect predictor for which the impact of random effect (grouping
     factor) is to be evaluated, a formula with such terms specifying
     a mixed-effects model. Mixed-effects models with random effect
     terms can only be specified in the formula when using the
     \code{"poisson.reformulation"} heuristic. If \code{formula}
     contains random terms, \code{heuristic} will be automatically
     switched to \code{"poisson.reformulation"}.

     The actual fitting of the polytomous logistic regression model is
     undertaken by the function determined by the \code{heuristic}
     argument, i.e. \code{\link{polytomous.one.vs.rest}} or
     \code{\link{polytomous.poisson.reformulation}}. Models with only
     fixed-effect predictors will be fit using \code{\link{glm}};
     models with also random-effect predictors will be fit using
     \code{\link[lme4]{lmer}}.

}
\value{
   \code{polytomous} returns an object of class \code{"polytomous"},
   a list containing at least following components:

  \describe{
  \item{\code{model}}{the underlying model(s) fitted using \code{\link{glm}} or \code{\link{lmer}}.}
  \item{\code{data}}{the originally supplied \code{data} argument object.}
  \item{\code{frequency}}{the originally supplied \code{frequency} argument.}
  \item{\code{logodds}}{a matrix of the outcome-by-predictor logodds estimated for the model.}
  \item{\code{odds}}{a matrix of the outcome-by-predictor odds (\code{exp(logodds)}) estimated for the model.}
  \item{\code{p.values}}{a matrix of the estimates of the significances of the outcome-by-predictor logodds/odds estimated for the model.}
  \item{\code{fitted}}{a matrix of the fitted outcome-specific probability estimates corresponding to the original data table in uncounted \code{"instance"} or \code{"narrow"} count or \code{"wide"} count format.}
  \item{\code{statistics}}{a range of descriptive statistics describing the goodness of fit and classification performance of the model; see \code{\link{model.statistics}} and \code{\link{crosstable.statistics}}.}
  \item{\code{formula}}{the formula specification used to fit the model.}
  \item{\code{outcomes}}{the outcome categories.}
  \item{\code{heuristic}}{the heuristic used to fit the model.}
  \item{\code{data.format}}{the format type of the \code{data} argument object; having the value of either \code{"instance"}, \code{"narrow"}, or \code{"wide"}.}
  }

  The functions implementing the various heuristics may provide
  additional heuristic-specific results, see
  \code{\link{polytomous.one.vs.rest}} and
  \code{\link{polytomous.poisson.reformulation}}.

}
\note{

  In addition to polytomous logistic regression modeling and auxiliary
  functions, the \code{polytomous} package also contains two functions
  for the univariate analysis data with categorical (nominal), unordered
  variables, namely \code{\link{associations}} and
  \code{\link{posthoc.cells}}.

}
\references{

Arppe, A. 2008. Univariate, bivariate and multivariate methods in
corpus-based lexicography -- a study of synonymy. Publications of the
Department of General Linguistics, University of Helsinki,
No. 44. URN: http://urn.fi/URN:ISBN:978-952-10-5175-3.

Arppe, A. 2009. Linguistic choices vs. probabilities -- how much and
what can linguistic theory explain? In: Featherston, Sam & Winkler,
Susanne (eds.) The Fruits of Empirical Linguistics. Volume 1:
Process. Berlin: de Gruyter, pp. 1--24.

Antti. A. (in prep.)

}
\author{
Antti Arppe, with ideas and input from R. Harald Baayen
}

\seealso{
\code{\link{polytomous.one.vs.rest}}, \code{\link{polytomous.poisson.reformulation}}, 
\code{\link{summary.polytomous}}, \code{\link{anova.polytomous}},
\code{\link{predict.polytomous}}, \code{\link{plot.polytomous}},
\code{\link{posthoc.cells}}, \code{\link{associations}}
}
\examples{

data(think)
think.polytomous1 <- polytomous(Lexeme ~ Agent + Patient, data=think)
summary(think.polytomous1)

think.polytomous2 <- polytomous(Lexeme ~ Agent + Patient + Register, data=think)
summary(think.polytomous2)

think.polytomous.lmer1 <-  polytomous(Lexeme ~ Agent + Patient + (1|Register),
   data=think, heuristic="poisson.reformulation")
summary(think.polytomous.lmer1)

\dontrun{
think.polytomous3 <- polytomous(Lexeme ~ Agent + Patient + Section,
   data=think)
summary(think.polytomous2)

think.polytomous.lmer2 <-  polytomous(Lexeme ~ Agent + Patient + (1|Section),
   data=think, heuristic="poisson.reformulation")
summary(think.polytomous.lmer2)
}

}

\keyword{ regression }

