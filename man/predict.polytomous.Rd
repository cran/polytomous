\name{predict.polytomous}
\alias{predict.polytomous}
\title{
Predict method for polytomous objects
}
\description{

   Obtains predictions on the basis of a fitted \code{"polytomous"}
   object on data already incorporated in the object or on new
   data with the same predictors as the originally fitted model object.

}
\usage{
\method{predict}{polytomous}(object, newdata=NULL, type="response",
   p.normalize = TRUE, \dots)

}

\arguments{
  \item{object}{

      objects of class \code{"polytomous"}, typically the result of a
      call to \code{polytomous}.

  }
  \item{newdata}{

     optionally, a data frame in which to look for variables with
     which to predict.  If omitted (i.e. set to \code{NULL}), the
     fitted linear predictors of the \code{object} are used.

  }
  \item{type}{

    the type of prediction requested.  For the default
    \code{type="link"}, the predictions are cumulative log-odds
    (estimated probabilities on logit scale), while
    \code{type="response"} yields the distributions of predicted
    probabilities over the outcome responses.  The option
    \code{type="terms"} returns a matrix giving the fitted values of
    each term in the model formula on the linear predictor scale,
    whereas the option \code{type="choice"} produces the predicted
    individual discrete choices, given the selected predictors. The
    prediction type values \code{"link"} and \code{"terms"} are
    currently only implemented for the \code{"one.vs.rest"} heuristic,
    and are based on the function \code{\link{predict.glm}}.

  }
  \item{p.normalize}{

    a logical indicating whether outcome-specific probability
    estimates ought to be normalized so that
    sum(P(outcome|predictors))=1; only applicable and possibly
    desirable for \code{"polytomous"} objects fit with
    \code{heuristic="one.vs.rest"}, since the constituent binary
    models are fit independently of each other.

  }
  \item{\dots}{

    further arguments passed to and from other functions.

  }
}

\details{

    If \code{newdata} is omitted the predictions are based on the data
    used for the fit.

}
\value{

   a vector or matrix of predictions.

}
\references{
Arppe, A. (in prep.)
}
\author{
Antti Arppe
}

\seealso{
\code{\link{polytomous}}, \code{\link{polytomous.one.vs.rest}},
\code{\link{polytomous.poisson.reformulation}}, \code{\link{predict.glm}}
}
\examples{

data(think)
think.polytomous <- polytomous(Lexeme ~ Agent + Patient, data=think)
head(predict(think.polytomous, type="response"))
predict(think.polytomous, newdata=think[1:20,], type="choice")

}

\keyword{ regression }
