summary.polytomous <-
function(object, ...)
{ 
  classes <- class(object)
  class(object) <- c("summary.polytomous",classes)
  return(object)
}

#summary.mer <-
#function(object, ...)
#{
#  require(lme4, quietly=TRUE)
#  summary(object, ...)
#}


print.summary.polytomous <- function(x, digits=max(3,getOption("digits")-3), parameter="odds", max.parameter=ifelse(parameter=="odds",10000,100), max.print=10, cycles=0, max.denominator=0, ...)
{ 
  require(MASS, quietly=TRUE)
  if(!is.null(x$max.print) & is.numeric(x$max.print))
    max.print=x$max.print;
  
  cat("\nFormula:\n")
  print(x$formula)

  cat("\nHeuristic:\n")
  cat(x$heuristic)
  cat("\n")

  if("mixed" %in% class(x))
    { cat("\nRandom effects:\n")
      sd.random <- unlist(x$model@ST)*lme4:::sigma(x$model) # summary(x$model)@sigma
      var.random <- sd.random^2
      names.random <- names(x$model@flist)
      summary.random <- data.frame(cbind(names.random,signif(var.random,digits),signif(sd.random,digits)))
      names(summary.random) <- c("Groups","Variance","Std.Dev.")
      print(data.frame(summary.random), row.names=FALSE, right=FALSE)
    }

  if(parameter=="odds")
    { cat("\nOdds:\n")
        parameters <- x$odds
    }
  if(parameter=="logodds")
   { cat("\nLog-odds:\n")
     parameters <- x$logodds
   }
  p.values <- x$p.values
  dim.parameters <- dim(parameters)
  if(cycles>0 & max.denominator>0 & parameter=="odds")
    char.parameters <- apply(parameters,c(1,2),function(x) as.character(fractions(x,cycles,max.denominator)))
  else
    char.parameters <- apply(parameters, c(1,2), function(x)
                            { if(is.na(x))
                                return("NA")
                              if(x > max.parameter)
                                return("Inf")
                              if(x < 1/max.parameter & parameter=="odds")
                                return("1/Inf")
                              if(x < -max.parameter & parameter=="logodds")
                                return("-Inf")
                              return(as.character(signif(x,digits)))
                            })
  for(i in 1:dim.parameters[1])
     for(j in 1:dim.parameters[2])
        if(!is.na(p.values[i,j]) & !is.na(char.parameters[i,j]))
          if(p.values[i,j]>=.05)
            char.parameters[i,j] <- paste("(",char.parameters[i,j],")",sep="")
  print.data.frame(as.data.frame(char.parameters)[1:min(nrow(char.parameters),max.print),], na.print="NA")

  if(nrow(char.parameters)>max.print)
    cat(paste("... [ omitted ",nrow(char.parameters)-max.print," rows ] ...\n",sep=""))
  cat("\n")
  
  statistics <- x$statistics
  deviances <- format(c(signif(statistics$deviance.null,digits),signif(statistics$deviance.model,digits)))
  DFs <- format(c(statistics$df.null,statistics$df.model))
  cat(c("Null deviance:             ", deviances[1], " on ", DFs[1], " degrees of freedom\n"))
  cat(c("Residual (model) deviance: ", deviances[2], " on ", DFs[2], " degrees of freedom\n"))
  cat(format(c("\nR2.likelihood: ", signif(statistics$R2.likelihood,digits), "\nAIC: ", signif(statistics$AIC.model,digits), "\nBIC: ", signif(statistics$BIC.model,digits))))
  cat("\n\n")

  invisible(x)
}
