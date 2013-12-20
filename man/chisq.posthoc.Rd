\name{chisq.posthoc}
\alias{chisq.posthoc}
\title{
  Calculate cellwise posthoc analyses for a two-way contingency table 
}
\description{

  \code{chisq.posthoc} takes a contingency table crosstabulating two
  categorical, unordered variables, the overall
  independence/dependence of which has been evaluated with
  \code{\link{chisq.test}}, and calculates several variations of
  posthoc analyses concerning the impact of individual cells
  containing frequencies of value pairings of the two categorical
  variables, assessing the degree to which these individual cellwise
  observed values diverge (or not) from an overall hypothetical
  homogeneous distribution.

}
\usage{
chisq.posthoc(ctable, alpha = 0.05, reorder = "none",
   std.pearson.residual.min = 2, correct=FALSE)
}
\arguments{
    \item{ctable}{a two-way contigency table crosstabulating two, possible multiple-valued categorical, unordered variables}
    \item{alpha}{a numerical value between 0 and 1 specifying the critical P-value threshold for significance; by default set to 0.05}
    \item{reorder}{a character string specifying whether the rows or columns of \code{ctable}, or both, should be reordered according to their descending marginal frequencies; possible values \code{"none"} (default), \code{"both"}, \code{"rows"}, or \code{"cols"}; if none of these are provided, \code{ctable} will be left as it is.}
    \item{std.pearson.residual.min}{the minimum absolute value for considering a cellwise standardized Pearson residual in the contingency table to deviate significantly from the expected value (representing a homogeneous distribution); by default equal to 2.}
    \item{correct}{a logical indicating whether to apply Yates' continuity correction when computing the test statistic; by default set to \code{FALSE}}
  }
\value{
  A list with the following components:
  \describe{
    \item{\code{ctable}}{the two-way contingency table, reordered as specified by the \code{reorder} argument.}
    \item{\code{X2.df}}{the minimum value of the chi-squared test statistic for the degrees of freedom (\code{df=(nrow(ctable)-1)*(ncol(ctable)-1)}) of the contingency table so that the distribution of counts in the contingency table can be considered to overall diverge significantly (from the expected values representing a homogeneoues distribution), so that the probability of observing by chance such a distribution are at most \code{P=alpha}.}
    \item{\code{X2.df1}}{the minimum value of the chi-squared test statistic for the minimum \code{df=1}.}
    \item{\code{cells}}{a list of different assessments of the divergences of the cellwise values from expected values representing a homogeneous distribution, consisting of the following elements:}
    \describe{
      \item{\code{X2}}{the cellwise contributions to the chi-squared statistic, with the sign indicating whether the observed value is greater or less than the expected value.}
      \item{\code{X2.df.sign}}{the cellwise assessment of whether the chi-squared value of an individual cell by itself exceeds the overall minimum chi-squared value \code{X2.df} for a significantly non-homogeneous distribution of counts; having the values \sQuote{+} or \sQuote{-} when this is the case, with the sign indicating whether the cellwise observed value is greater or less than the expected value, or \sQuote{0} otherwise (indicating no significant divergence).}
      \item{\code{X2.df1.sign}}{the cellwise assessment of whether the chi-squared value of an individual cell exceeds the minimum chi-squared value \code{X2.df} when \code{df=1}; having the values \sQuote{+}, \sQuote{-}, or \sQuote{0}.}
      \item{\code{std.pearson.residuals}}{the cellwise standardized Pearson residuals.}
      \item{\code{std.pearson.residuals.sign}}{the cellwise assessment of whether the absolute value of a standardized Pearson residual is greater than \code{std.pearson.residual.min}, having the values \sQuote{+}, \sQuote{-}, or \sQuote{0}.}
    }
  }
}
\details{

 The cellwise posthoc analyses are variations based on the overall
 chi-squared test for homogeneity/heterogeneity of frequency
 distributions represented in contingency tables as implemented in
 \code{\link{chisq.test}}.

 Though for smaller contingency tables the suggested absolute minimum
 cellwise value for a standardized Pearson residual to be considered
 to signal a potentially significant cellwise divergence (in relation
 to an overall homogeneous distribution for the entire table) is 2 or
 more, this absolute threshold value should probably be increased to 3
 or more in the case of larger contingency tables (Agresti 2002).

 With respect to minimum cell-wise expected counts, in contrast to
 \code{\link{chisq.test}} no warnings will be output.

}
\references{

  Agresti, A. (2002) Categorical Data Analysis (2nd
  edition). Hoboken: John Wiley & Sons, Hoboken.

  Arppe, A. (2008) Univariate, bivariate and multivariate methods in
  corpus-based lexicography -- a study of synonymy. Publications of
  the Department of General Linguistics, University of Helsinki,
  No. 44. URN: http://urn.fi/URN:ISBN:978-952-10-5175-3.

  Liebetrau, A. M. (1983) Measures of Association. Sage University
  Paper series on Quantitative Applications in the Social Sciences,
  07-032. Beverly Hills and London: Sage Publications.

}
\author{
Antti Arppe
}
\seealso{
See also \code{\link{chisq.test}}, \code{\link{associations}}.
}
\examples{

data(think)
ctable <- table(think$Lexeme, think$Agent)
chisq.posthoc(ctable)

chisq.posthoc(table(think$Agent, think$Patient))
}
\keyword{ univar category }
