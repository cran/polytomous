\name{plot.polytomous}
\alias{plot.polytomous}
\alias{plot.polytomous.parameters}
\alias{plot.polytomous.probabilities}
\title{
  Plot function for selected results of \code{polytomous}.
}
\description{

  This function presents visually the estimated logodds or odds, or
  expected probabilities for a model fitted with
  \code{\link{polytomous}} and its auxiliary functions
  \code{\link{polytomous.one.vs.rest}} or
  \code{\link{polytomous.poisson.reformulation}}

}
\usage{
\method{plot}{polytomous}(x, values="probabilities", \dots)

\method{plot}{polytomous.parameters}(x, values="logodds",
type="density", predictors=NULL, outcomes=NULL, panes="single",
lty=NULL, col=NULL, mfrow=NULL, main=NULL,
legend.position="topright", \dots)

\method{plot}{polytomous.probabilities}(x, type="density",
select="all", panes="single", lty=NULL, col=NULL, pch=NULL,
mfrow=NULL, main=NULL, legend.position="topright", \dots)

}
\arguments{
  \item{x}{

     A object of the class \code{"polytomous"} produced by
     \code{\link{polytomous}} or its auxiliary functions
     \code{\link{polytomous.one.vs.rest}} or
     \code{\link{polytomous.poisson.reformulation}}, consisting of a
     list including estimated logodds or odds for predictors and
     estimated probabilities of outcomes for outcome-predictor
     combinations.

  }
  \item{values}{

     A character string specifiying whether expected
     \code{"probabilities"} (default) or estimated \code{"logodds"} or
     \code{"odds"} should be plotted.

  }
  \item{type}{ 

     A character string spefifying the type of plot to be drawn;
     \code{"density"} is available for both value types as default,
     while a histogram (\code{"hist"}) is available only for
     \code{plot.polytomous.parameters} and sorted values (\code{"sort"})
     only for \code{plot.polytomous.probabilities}.

  }
  \item{panes}{

     A character string specifying whether a \code{"single"} pane
     (default) integrating all component plots, or \code{"multiple"}
     panes for each individual component plot are to be plotted. If
     \code{"multiple"} panes are selected, the number or rows and
     columns is specified automatically. Alternatively, one can invoke
     the plotting of multiple panes by explicitly specifying the
     appropriate number of rows and columns with the parameter
     \code{mfrow} (N.B. this overrides \code{panes="single"}).

  }
  \item{predictors}{ 

     A regular expression specifying which predictors and their values
     should be included in the plot(s); by default \code{NULL} so that
     all predictors incorporated in the \code{"polytomous"} model will
     be included.

  }
  \item{outcomes}{

     A list of outcomes to be included in the plot; by default
     \code{NULL} so that all outcomes will be considered.

  }
  \item{select}{

    For the method \code{plot.polytomous.probabilities}, a character
    string specifying which instance-wise probability estimates should
    be plotted; by default \code{"all"}, other values are \code{"max"}
    for instance-wise maximum probabilities, \code{"min"} for
    instance-wise minimum probabilities, \code{"maxmin", "minmax"} for
    both maximum and minimum instance-wise
    probabilities. Alternatively, a numeric vector \code{c(1,2,\dots)}
    specifying selected ranks of the instance-wise probability
    estimates can be provided, with \code{1} corresponding to the
    instance-wise maximum probability estimates.

  }
  \item{lty, col, pch, mfrow, main, legend.position}{

    Specifications of various graphical parameters (see
    \code{\link{par}}) to be used in the plots; if any of these is set
    to \code{=NULL} default settings will be used (for
    \code{legend.position}, the default value is
    \code{topright}). Note that \code{lty} is relevant only to
    \code{plot.polytomous.parameters(\dots, type="density", \dots)}
    and \code{plot.polytomous.probabilities(\dots, type="density",
    \dots)}, and \code{pch} only to
    \code{plot.polytomous.probabilities(\dots, type="sort", \dots)}.

  }
  \item{\dots}{ 

     Arguments to be passed to methods, such as graphical parameters
     (see \code{\link{par}}).

  }
}
\value{
  A plot of the selected type is produced on the graphics device.
}
\references{

Arppe, A. (2008) Univariate, bivariate and multivariate methods in
corpus-based lexicography -- a study of synonymy. Publications of the
Department of General Linguistics, University of Helsinki,
No. 44. URN: http://urn.fi/URN:ISBN:978-952-10-5175-3.

Arppe, A. (2009) Linguistic choices vs. probabilities -- how much and
what can linguistic theory explain? In: Featherston, S. and
S. Winkler, (eds.) The Fruits of Empirical Linguistics. Volume 1:
Process. Berlin: de Gruyter, pp. 1--24.

Antti. A. (in prep.) Solutions for fixed and mixed effects modeling of
polytomous outcome settings.

}
\author{
Antti Arppe
}
\seealso{
\code{\link{polytomous}}, \code{\link{polytomous.one.vs.rest}},
\code{\link{polytomous.poisson.reformulation}}
}
\examples{

data(think)
think.polytomous <- polytomous(Lexeme ~ Agent + Patient + Section,
   data=think)

plot(think.polytomous, values="logodds")
plot(think.polytomous, values="logodds", type="hist", panes="multiple")
plot(think.polytomous, values="logodds", type="density", panes="multiple")
plot(think.polytomous, values="logodds", type="density", panes="multiple",
   predictors="Section*")
plot(think.polytomous, values="logodds", type="density", panes="multiple",
   predictors="Patient*")
plot(think.polytomous, values="logodds", type="hist", panes="multiple", col=1:4)
plot(think.polytomous, values="logodds", type="density", panes="single",
   outcomes=c("ajatella","miettia","pohtia","harkita"))

plot(think.polytomous, values="probabilities")
plot(think.polytomous, values="probabilities", panes="multiple")
plot(think.polytomous, values="probabilities", select="max")
plot(think.polytomous, values="probabilities", select=c(1:3))
plot(think.polytomous, values="probabilities", panes="multiple", select=c(1:3))
plot(think.polytomous, values="probabilities", type="sort", legend.position="topleft")
plot(think.polytomous, values="probabilities", type="sort", pch=".",
   legend.position="topleft")
plot(think.polytomous, values="probabilities", type="sort", pch=".", panes="multiple")

}
\keyword{ classif }
