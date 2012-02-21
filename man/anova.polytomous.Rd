\name{anova.polytomous}
\alias{anova.polytomous}
\alias{anova.polytomouslist}
\title{
Analysis of Model Fit for Polytomous Logistic Regression models
}
\description{

  Calculate an analysis of individual variable contributions or model
  comparisons for one or more Polytomous Logistic Regression models.

}
\usage{
\method{anova}{polytomous}(object, \dots, statistic = "deviance", test = "Chisq",
   outcome.specific = FALSE)
}

\arguments{
  \item{object, \dots}{

      objects of class \code{"polytomous"}, typically the result of a
      call to \code{polytomous}, or a list of \code{object}s for the
      \code{polytomouslist} method.

  }
  \item{statistic}{

      a character string, determining the statistic for evaluating
      model fit, by default \code{"deviance"}, alternatively
      \code{"AIC"} or \code{"BIC"}.

  }
  \item{test}{

      a character string, determining the statistical method by which
      the significance of the comparison are done, by default the
      Chi-squared test (\code{"Chisq"}); currently no other methods are
      implemented. If set to \code{NULL}, no significance testing will
      be undertaken.

  }
  \item{outcome.specific}{

      a logical, which, if set \code{TRUE} in the case of a single
      \code{"polytomous"} \code{object} fit using
      \code{heuristic="one.vs.rest"}, will result in the presentation
      of the application of \code{anova.glm} on the outcome-specific
      reduction of \code{deviance} of the constituent binary models;
      by default set \code{FALSE} resulting in a conventional ANOVA
      table.

  }
}

\details{

     Specifying a single object gives a table with sequential analysis
     of predictor impact with respect to the selected statistic of
     model fit. That is, the reductions in the residual statistic as
     each term of the formula is added in turn are given in as the
     rows of a table, plus the residual statistic values themselves.

     If more than one object is specified, the table has a row for the
     residual degrees of freedom and selected statistic for each
     model.  For all but the first model, the change in degrees of
     freedom and the statistic is also given. (This only makes
     statistical sense if the models are nested.)  It is conventional
     to list the models from smallest to largest, but this is up to
     the user.

     The table will optionally contain test statistics (and P values)
     comparing the reduction in deviance for the row to the residuals.
     Only a comparison of models or contributions of their components
     by the chi-squared test has been implemented, which is applicable
     only for the \code{"deviance"} statistic.

     The comparison between two or more models by
     \code{anova.polytomous}, redirected to the
     \code{anova.polytomouslist} method, will only be valid if they
     are fitted to the same dataset and with the same
     heuristic. \code{anova.polytomouslist} will look for such
     discrepancies, resulting in an error when detected.

     If \code{outcome.specific=TRUE}, the function will alternatively
     output the outcome-specific reductions of \code{deviance} for
     sequantially added predictors for each of the constituent binary
     models, using \code{anova.glm}. Thus,
     \code{outcome.specific=TRUE} is only applicable for a single
     \code{object} of the class \code{"polytomous"} fit with
     \code{heuristic="one.vs.rest"} and when
     \code{statistic="deviance"}.

}
\value{

     An object of class \code{"anova"} inheriting from class
     \code{"data.frame"}, when \code{outcome.specific=FALSE}.

     When \code{outcome.specific=TRUE}, the function will produce a
     list with the following components:

     \describe{
       \item{\code{model}}{A list of outcome-specific results procuded by
       \code{anova.glm} on the constituent binary models.}

       \item{\code{deviance}}{A table with outcome-specific reductions of
       \code{deviance} (columns) for each sequentially added predictor
       (rows).}

       \item{\code{df}}{A table with outcome-specific reductions of degrees
       of freedom (columns) for each sequentially added predictor
       (rows).}

       \item{\code{p.values}}{A table with outcome-specific evaluations
       (using the Chi-squared test) of the significance of reduction
       in \code{deviance} (columns) for each sequentially added
       predictor (rows).}

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
\code{\link{polytomous}}, \code{\link{anova.glm}}
}
\examples{

data(think)
think.polytomous1 <- polytomous(Lexeme ~ Agent * Patient, data=think)
anova(think.polytomous1)

\dontrun{
anova(think.polytomous1, statistic="AIC")
anova(think.polytomous1, statistic="BIC", test=NULL)
}

anova(think.polytomous1, outcome.specific=TRUE)

think.polytomous2 <- polytomous(Lexeme ~ Agent * Patient + Manner,
   data=think)
anova(think.polytomous1, think.polytomous2)

\dontrun{
anova(think.polytomous1, think.polytomous2, statistic="AIC", test=NULL)
}

}

\keyword{ regression }
