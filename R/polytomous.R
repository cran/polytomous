polytomous <-
function(formula, data, heuristic="one.vs.rest", ...)
{ 
  predictors <- attr(terms.formula(formula),"term.labels");
  if(length(grep("\\|",predictors))!=0 & heuristic=="one.vs.rest")
    { cat("Switching 'heuristic' to 'poisson.reformulation' due to random effect terms in 'formula'\n")
      heuristic="poisson.reformulation"
    }

  if(heuristic=="one.vs.rest")
    polytomous.one.vs.rest(formula, data, ...)
  else
  if(heuristic=="poisson.reformulation")
    polytomous.poisson.reformulation(formula, data, ...)
  else
    stop(paste(c("'heuristic': ",heuristic," not implemented/available.")));
}

print.polytomous <-
function(x, max.print=10, ...)
{
  digits=max(3,getOption("digits")-3)
  if(is.na(max.print))
    max.print <- nrow(x$odds)

  cat("\nFormula:\n")
  print.formula(x$formula)
  cat("\nHeuristic:\n")
  cat(x$heuristic)
  cat("\n")

  cat("\nOdds:\n")
  print.data.frame(as.data.frame(x$odds)[1:min(nrow(x$odds),max.print),], digits=digits)
  if(nrow(x$odds)>max.print)
    cat(paste("... [ omitted ",nrow(x$odds)-max.print," rows ] ...\n",sep=""))
  cat("\n")

  invisible(x)

}