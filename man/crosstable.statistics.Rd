\name{crosstable.statistics}
\alias{crosstable.statistics}
\title{
  Calculate prediction accuracy statistics for a contingency table 
}
\description{
  \code{crosstable.statistics} takes a contingency table of observed
  vs. predicted values for a binary or polytomous response variable as
  input, and calculates a range of statistics about prediction
  accuracy.
}
\usage{
crosstable.statistics(ctable)
}
\arguments{
  \item{ctable}{
    A contingency table cross-classifying observed and predicted values.
  }
}
\value{
  A list with the following components:
  \describe{
    \item{\code{accuracy}}{Overall prediction accuracy}
    \item{\code{recall.predicted}}{Recall of prediction for each outcome value}
    \item{\code{precision.predicted}}{Precision of prediction for each outcome value}
    \item{\code{lambda.prediction}}{lambda for prediction accuracy (improvement over baseline of always predicting mode)}
    \item{\code{tau.classification}}{tau for classification accuracy (improvement over baseline of homogeneous distribution of predicted outcomes)}
    \item{\code{d.lambda.prediction}}{d(lambda): used for calculating \code{P(lambda)}}
    \item{\code{d.tau.classification}}{d(tau): used for calculating \code{P(tau)}}
    \item{\code{p.lambda.prediction}}{P(lambda): probability of reaching \code{lambda} by chance}
    \item{\code{p.tau.classification}}{P(tau): probability of reaching \code{tau} by chance}
  }
}
\references{
  Arppe, Antti (in prep.).

  Menard, Scott (1995). Applied Logistic Regression Analysis. Sage University Paper 
  Series on Quantitative Applications in the Social Sciences 07-106. Thousand 
  Oaks: Sage Publications. 
}
\author{
Antti Arppe
}
\seealso{
See also \code{\link{model.statistics}}.
}
\examples{
ctable=matrix(c(30, 10, 5, 60),2,2)
crosstable.statistics(ctable)
}
\keyword{ classif }
