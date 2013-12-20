\name{shanghainese}
\alias{shanghainese}
\docType{data}
\title{
  Shanghainese topic markers.
}
\description{
  500 occurrences of the five most frequent topic markers \sQuote{ne, a, mo, zi, ma} in Shanghainese.
}
\usage{data(shanghainese)}
\format{
  A data frame with 500 observations on the following 6 variables.

  \describe{
    \item{\code{TOPIC_MARKER}}{A factor specifying one of the five topic markers}
    \item{\code{TOPIC_LENGTH}}{A numeric vector specifying the character length of the topic}
    \item{\code{TOPIC_POS}}{A factor specifying the part-of-speech of the topic}
    \item{\code{FUNCTION}}{A factor specifying function of the topic}
    \item{\code{COMMENT_TYPE}}{A factor specifying type of the comment}
    \item{\code{GENRE}}{A factor specifying the genre in which the topic marker had been used}
}}
\details{

The five most frequent topic markers \sQuote{ne, a, mo, zi, ma} in
Shanghainese were extracted SOURCES (REFERENCES). The
\code{shanghainese} dataset contains a selection of 5 contextual
features judged as most informative.

For extensive details of the data and its linguistic and statistical
analysis, see Han, Arppe & Newman (forthc.).

}

\references{

Han, W., A. Arppe & J. Newman (forthc.) Topic marking in a
Shanghainese corpus: From observation to prediction. Corpus
Linguistics and Linguistic Theory. DOI: 10.1515/cllt-2013-0014.

}
\examples{

## For examples see vignette

}
\keyword{ datasets }
