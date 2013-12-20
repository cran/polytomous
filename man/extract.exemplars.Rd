\name{extract.exemplars}
\alias{extract.exemplars}
\title{Extract a subset of exemplary exemplars from a dataset}

\description{

A function that extracts a subset of exemplary exemplars from a
dataset that has been used to fit a polytomous logistic regression
model. Hierarchical agglomeratative clustering (HAC) is used to divide
the dataset into distinct subsets in terms of their
features/properties; from each such cluster a single exemplar with an
outcome that has the highest expected probability for all the
exemplars within the cluster is selected. Consequently, the number of
exemplars that is extracted from the dataset equals the number of
clusters into which the dataset is divided.

}
\usage{
extract.exemplars(model.polytomous, model.hclust=NULL, n.clusters=10, p.bins=0,
   features=FALSE)
}
\arguments{
  \item{model.polytomous}{

   an object of class \code{"polytomous"} that has been fitted with
   the function \code{\link{polytomous}}.

}
  \item{model.hclust}{

   an object of class \code{"hclust"} that has been produced by
   applying the function \code{\link{hclust}} to the dataset. If none
   is provided \code{NULL} (default), one will be automatically
   created using the elements \code{data} (which must consist of
   logical features/properties) and \code{formula} included in the
   polytomous model inputted as the argument \code{model.polytomous};
   the lables in \code{formula} will specify the features that will be
   used in the clustering.

}
  \item{n.clusters}{ 

   a numeric argument specifying the number of clusters into which the
   dataset will be divided by applying the function
   \code{\link{cutree}} on the hierarchical agglomerative clustering
   model; by default \code{=10}.

}
  \item{p.bins}{

  a numeric argument specifying whether the exact probability
  estimates for the outcomes as provided by \code{model.polytomous}
  will be used in selecting the individual exemplars from the clusters
  (the default case determined by setting the value as \code{=0}), or
  the number of equal interval probability bins into which the exact
  probability estimates will be divided. In the latter case, the
  exemplar selected from the highest probability bin might not have
  received the absolutely highest probability estimate.

}
  \item{features}{

  a logical (by default \code{=FALSE}) indicating whether the number
  of features/properties evident in the individual exemplars in the
  dataset and belonging to the set of features specified in the
  \code{formula} will be used as a secondary ranking factor in
  addition to the outcome-specific probability estimates (or
  probability bins). When \code{features=TRUE} and \code{p.bins >= 2},
  the exemplar with the highest number of features/properties will be
  selected from the highest probability bin per each cluster.

}
}

\details{

The hierarchical agglomerative clustering, if done automatically
within the function, will be undertaken with the function
\code{\link{hclust}} using \code{method="ward"} on a distance matrix
created using the \code{\link{dist}} function with
\code{method="binary"}. Therefore, the \code{data} and
\code{formula}in \code{model.polytomous} must consist of logical
features/properties. The transformation of multinomial features to
logical ones can be undertaken with \code{\link{multinomial2logical}}.

}

\value{
  \code{extract.exemplars} returns a list with the following components:

  \describe{

  \item{\code{indices}}{A numeric vector of indices of the exemplar
  rows in the dataset that have been extracted.}

  \item{\code{outcomes}}{A factor containing the outcomes in each of
  the individual exemplars.}

  \item{\code{max.probs}}{A numeric vector of the exact probability
  estimates for the outcome apparent in each exemplar.}

  \item{\code{properties}}{A factor with the set of
  features/properties evident in each exemplar; the properties are
  separated by a semicolon and space.}
  }
}

\references{

Arppe, A. (2008) Univariate, bivariate and multivariate methods in
corpus-based lexicography -- a study of synonymy. Publications of the
Department of General Linguistics, University of Helsinki,
No. 44. URN: http://urn.fi/URN:ISBN:978-952-10-5175-3.

Divjak, D. and A. Arppe (2013). Extracting prototypes from
exemplars. What can corpus data tell us about concept representation?
Cognitive Linguistics, 24 (2): 221-274.

}

\author{
Antti Arppe
}

\seealso{
\code{\link{extract.prototypes}}, \code{\link{polytomous}}, \code{\link{hclust}},
   \code{\link{dist}}
}
\examples{

\dontrun{
data(think)
think.logical <- multinomial2logical(data=think, outcome="Lexeme",
   variables=names(think)[2:24])
think.formula <- as.formula(paste("Lexeme",
   paste(grep("(Other)|(None)|(FiniteVerbChain)|(Overt)",
   names(think.logical)[-1],value=TRUE,invert=TRUE),collapse=" + "),sep=" ~ "))
think.polytomous <- polytomous(think.formula, data=think.logical)
extract.exemplars(think.polytomous,n.clusters=50)
}

## For more details, see vignette.

}

\keyword{ regression multivariate category cluster classif }
