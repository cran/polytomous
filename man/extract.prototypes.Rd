\name{extract.prototypes}
\alias{extract.prototypes}
\title{Extract the prototypical features for a set of outcomes
}
\description{

A function that extracts for each outcome included in a polytomous
logistic regression model fitted with the function
\code{\link{polytomous}} the set of features/properties that together can
be interpreted to represent the prototypical characteristics of the
outcome in question.

}
\usage{
extract.prototypes(model.polytomous, p.critical=.05)

}
\arguments{
  \item{model.polytomous}{

   an object of class \code{"polytomous"} that has been fitted with
   the function \code{\link{polytomous}}.
}
  \item{p.critical}{

   a numeric value specifying the critical p-level (by default
   \code{=.05}) for the coefficient (i.e. odds/logodds) estimated by
   \code{\link{polytomous}} for each feature/property included
   in the regression function to be considered statistically significant
}
}
\details{

   This function in effect automatically selects and groups together
   in a convenient manner per each outcome in the polytomous logistic
   regression model the set of significant explanatory features/properties.
}
\value{

   \code{extract.prototypes} returns a named list (by outcomes)
   consisting per each outcome of a single-column matrix with the
   significant feature/properties as rownames and the associated
   estimated odds as values, sorted in decreasing order for each
   outcome.

}

\references{

Arppe, A. (2008) Univariate, bivariate and multivariate methods in
corpus-based lexicography -- a study of synonymy. Publications of the
Department of General Linguistics, University of Helsinki,
No. 44. URN: http://urn.fi/URN:ISBN:978-952-10-5175-3.

Divjak, D. & A. Arppe (2013). Extracting prototypes from
exemplars. What can corpus data tell us about concept representation?
Cognitive Linguistics, 24 (2): 221-274

}
\author{
Antti Arppe
}
\seealso{
\code{\link{extract.exemplars}}, code{\link{polytomous}}
}
\examples{

data(think)
think.polytomous <- polytomous.one.vs.rest(Lexeme ~ Agent + Patient, data=think,
   heuristic="one.vs.rest")
extract.prototypes(think.polytomous)

## For more details, see vignette.

}
\keyword{ regression multivariate category classif cluster }

  

