\name{polytomous.one.vs.rest}
\alias{polytomous.one.vs.rest}
\title{Fitting polytomous logistic regression models with the one-vs-rest heuristic
}
\description{
A function fitting a polytomous logistic regression model based on the
one-vs-rest heuristic. With the one-vs-rest heuristic, each individual
outcome is contrasted with all the other outcomes lumped together.
}
\usage{
polytomous.one.vs.rest(formula, data, frequency = NA, p.normalize = TRUE, \dots)
}
\arguments{
  \item{formula}{

   an object of class \code{"formula"} (or one that can be coerced to
   that class): a symbolic description of the model to be fitted.  The
   details of model specification are given under \code{Details} for
   the \code{polytomous} function.

}
  \item{data}{

   a data frame (or object coercible by \code{as.data.frame} to a data
   frame) containing the variables specified in the model. The data
   may be represented either in the uncounted \code{"instance"} or the
   \code{"wide"} count format; see \code{Details} for the function
   \code{\link{polytomous}}.

}
  \item{frequency}{ 

    A numeric vector (or the name of a column in the input data frame)
    with the frequencies of the instance.  If absent (set to
    \code{NA}), each exemplar is assigned a frequency equal to 1.

}
  \item{p.normalize}{

  a logical indicating whether outcome-specific probability estimates
  ought to be normalized so that sum(P(outcome|predictors))=1.

}
  \item{\dots}{

  further control arguments passed to or from other methods, see
  \code{\link{glm}} or \code{\link{model.statistics}}.

}
}
\details{

\code{polytomous.one.vs.rest} fits a polytomous logistic regression
model with fixed-effects predictors using the \code{one.vs.rest}
heuristic, where the occurrences of each individual outcome are
contrasted with all the other outcomes (i.e. the \dQuote{rest}) lumped
together. Consequently, a polytomous one-vs-rest model consits of a
set of binary logistic regression models fitted with \code{\link{glm}}, with
one such model for each possible outcome. Every such binary logistic
regression is fitted independently of each other.

The outcome-predictor \code{odds} (which are simply
\code{exp(logodds)}) indicate how much the occurrence of some predictor
increases or decreases the odds-ratio of the outcome in question to
occur, instead of any other outcome, all other predictors considered
equal.

With multiple predictors, due to the independent fitting of the
outcome-specific binary logistic regression models, the sums of
probability estimates for particular predictor combinations provided
by these binary models are not always exactly equal to one (though
mostly quite close). Thus, it is recommended that the probability
estimates are normalized by setting \code{p.normalize} as \code{TRUE}
(which is the default setting).

}
\value{
   \code{polytomous.one.vs.rest} returns an object of class \code{c("polytomous","one.vs.rest")},
   a list containing the following components:

  \describe{
  \item{\code{model}}{a list containing the underlying outcome-specific binary models fitted using \code{\link{glm}}.}
  \item{\code{data}}{the originally supplied \code{data} argument object.}
  \item{\code{frequency}}{the originally supplied \code{frequency} argument.}
  \item{\code{logodds}}{a matrix of the outcome-by-predictor logodds estimated for the model.}
  \item{\code{odds}}{a matrix of the outcome-by-predictor odds (\code{exp(logodds)}) estimated for the model.}
  \item{\code{p.values}}{a matrix of the estimates of the significances of the outcome-by-predictor logodds/odds estimated for the model.}
  \item{\code{fitted}}{a matrix of the fitted outcome-specific probability estimates corresponding to the original data table in uncounted \code{"instance"} or \code{"wide"} count format.}
  \item{\code{statistics}}{a range of descriptive statistics describing the goodness of fit and classification performance of the model; see \code{model.statistics} and \code{crosstable.statistics}.}
  \item{\code{formula}}{the formula specification used to fit the model.}
  \item{\code{outcomes}}{the outcome categories.}
  \item{\code{heuristic}}{the heuristic used to fit the model.}
  \item{\code{data.format}}{the format type of the \code{data} argument object; having the value of either \code{"instance"}, \code{"narrow"}, or \code{"wide"}.}
  }
}
\references{

Arppe, A. (2008) Univariate, bivariate and multivariate methods in
corpus-based lexicography -- a study of synonymy. Publications of the
Department of General Linguistics, University of Helsinki,
No. 44. URN: http://urn.fi/URN:ISBN:978-952-10-5175-3.

Arppe, A. (2009) Linguistic choices vs. probabilities -- how much and
what can linguistic theory explain? In: Featherston, S. and S. Winkler,
(eds.) The Fruits of Empirical Linguistics. Volume 1: Process. Berlin:
de Gruyter, pp. 1--24.

Rifkin, R. and A. Klautau (2004) In Defense of One-Vs-All Classification. 
Journal of Machine Learning Research, pp. 101--141.

}
\author{
Antti Arppe
}
\seealso{
\code{\link{polytomous}}, \code{\link{polytomous.one.vs.rest}}, \code{\link{anova.polytomous}},
\code{\link{predict.polytomous}}, \code{\link{model.statistics}}, \code{\link{glm}}
}
\examples{

data(think)
think.polytomous <- polytomous.one.vs.rest(Lexeme ~ Agent + Patient, data=think,
   heuristic="one.vs.rest")
summary(think.polytomous)

think.polytomous$statistics
think.polytomous$logodds
think.polytomous$odds
think.polytomous$p.values

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

think.polytomous2 <- polytomous(ajatella|harkita|miettia|pohtia ~ Agent + Patient,
   data=think.Agent_Patient.wide)
summary(think.polytomous2)

identical(think.polytomous$odds, think.polytomous2$odds)
identical(round(think.polytomous$odds,5),round(think.polytomous2$odds,5))

}
\keyword{ regression }

