\name{model.statistics}
\alias{model.statistics}
\title{
Calculate statistics for goodness of fit and prediction accuracy for a model
}
\description{
  Calculate a range of goodness of fit measures for an model object fitted with some multivariate statistical method that yields probability estimates for outcomes.
}

\usage{
model.statistics(observed, predicted, p.values, frequency = NA,
   outcomes = NULL, p.normalize = TRUE, cross.tabulation = TRUE,
   p.zero.correction=1/(nrow(p.values)*ncol(p.values))^2, \dots)
}
\arguments{

  \item{observed}{

  observed values of the response variable

  }

  \item{predicted}{

  predicted values of the response variable; typically the outcome
  estimated to have the highest probability

}
  \item{p.values}{

  matrix of probabilities for all values of the response variable (i.e
  outcomes)

}

  \item{frequency}{ 

    A numeric vector (or the name of a column in the input data frame)
    with the frequencies of the instance.  If absent (set to
    \code{NA}), each exemplar is assigned a frequency equal to 1.

  }

  \item{outcomes}{

  the outcome categories

}
  \item{p.normalize}{

  if \code{TRUE}, probabilities are normalized so that \code{sum(P)}
  of all outcomes for each datapoint is equal to 1

}
  \item{cross.tabulation}{

  if \code{TRUE}, statistics on the crosstabulation of observed and
  predicted response values are calculated with
  \code{\link{crosstable.statistics}}

}
  \item{p.zero.correction}{

  a function to adjust slightly response/outcome-specific probability
  estimates which are exactly P=0; necessary for the proper
  calculation of pseudo-R-squared statistics; by default calculated on
  the basis of the dimensions of the matrix of probabilities
  \code{p.values} as \code{1/(nrow(p.values)*ncol(p.values))^2}

}
  \item{\dots}{
  further control arguments to be passed from and to other functions.
}

}

\value{ A list with the following components:
\describe{
  \item{\code{loglikelihood.null}}{Loglikelihood for null model}
  \item{\code{loglikelihood.model}}{Loglikelihood for fitted model}
  \item{\code{deviance.null}}{Null deviance}
  \item{\code{deviance.model}}{Model deviance}
  \item{\code{R2.likelihood}}{(McFadden's) R-squared}
  \item{\code{R2.nagelkerke}}{Nagelkerke's R-squared}
  \item{\code{crosstable}}{Crosstabulation of observed vs. predicted outcomes.}
  \item{\code{crosstable.statistics(crosstable)}}{Various statistics
   calculated on the crosstabulation (\code{crosstable}) with \code{\link{crosstable.statistics}}, if
   \code{cross.tabulation=TRUE}} 
  }
}
\references{

Arppe, A. (2008) Univariate, bivariate and multivariate methods in
corpus-based lexicography -- a study of synonymy. Publications of the
Department of General Linguistics, University of Helsinki,
No. 44. URN: http://urn.fi/URN:ISBN:978-952-10-5175-3.

Arppe, A. (in prep.) Solutions for fixed and mixed effects modeling of
polytomous outcome settings.

Hosmer, D. W., Jr., and S. Lemeshow (2000) Applied Regression Analysis 
(2nd edition). New York: Wiley.

}
\author{
Antti Arppe
}

\seealso{
\code{\link{crosstable.statistics}}, \code{\link{polytomous.one.vs.rest}}, \code{\link{polytomous.poisson.reformulation}}, \code{\link{polytomous}}
}
\examples{

## None for the time being

}

\keyword{ regression }

