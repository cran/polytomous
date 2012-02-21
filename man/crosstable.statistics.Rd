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
    \item{\code{lambda.prediction}}{\eqn{\lambda_{prediction}}: improvement in prediction accuracy over baseline of always predicting mode}
    \item{\code{tau.classification}}{\eqn{\tau_{classification}}: improvement in classification accuracy over baseline of homogeneous distribution of predicted outcomes}
    \item{\code{d.lambda.prediction}}{\eqn{d(\lambda_{prediction})}: used for calculating \code{p.lambda.prediction}}
    \item{\code{d.tau.classification}}{\eqn{d(\tau_{classification})}: used for calculating \code{p.tau.classification}}
    \item{\code{p.lambda.prediction}}{\eqn{P(\lambda_{prediction})}: probability of reaching \eqn{\lambda_{prediction}} by chance}
    \item{\code{p.tau.classification}}{\eqn{P(\tau_{classification})}: probability of reaching \eqn{\tau_{classification}} by chance}
  }
}
\references{

  Arppe, A. (2008) Univariate, bivariate and multivariate methods in
  corpus-based lexicography -- a study of synonymy. Publications of
  the Department of General Linguistics, University of Helsinki,
  No. 44. URN: http://urn.fi/URN:ISBN:978-952-10-5175-3.

  Arppe, A. (in prep.) Solutions for fixed and mixed effects modeling of
  polytomous outcome settings.

  Menard, S. (1995). Applied Logistic Regression Analysis. Sage
  University Paper Series on Quantitative Applications in the Social
  Sciences 07-106. Thousand Oaks: Sage Publications.

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
