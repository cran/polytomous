ranef.polytomous <- function(object, ...)
{ if(!"polytomous" %in% class(object))
    stop("'object' not of class 'polytomous'.")
  if(!"mixed" %in% class(object))
    stop("'object' does not contain a mixed-effects model.")

  return(ranef(object$model, ...))
}
