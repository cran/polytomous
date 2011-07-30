anova.polytomous <-
function(object, ..., statistic = "deviance", test = "Chisq", outcome.specific = FALSE) 
{ 
    if(outcome.specific & statistic!="deviance")
      stop("'outcome.specific=TRUE' only applicable when statistic is 'deviance'.")
    if(!is.null(test))
      if(test!="Chisq")
        stop(paste("no test: ",test," implemented.",sep=""))
    if(statistic!="deviance" & !is.null(test))
      warning("no tests properly implemented for statistics other than 'deviance'.")

    if(statistic=="deviance")
      { statistic.model="deviance.model"; statistic.null="deviance.null"; }
    else   
    if(statistic=="AIC")
      { statistic.model="AIC.model"; statistic.null="deviance.null"; }
    else
    if(statistic=="BIC")
      { statistic.model="BIC.model"; statistic.null="deviance.null"; }
    else
      stop(paste("'statistic': ",statistic," not implemented.",sep=""))

    dotargs <- list(...)
    named <- if (is.null(names(dotargs))) 
        rep(FALSE, length(dotargs))
    else (names(dotargs) != "")
    if (any(named)) 
        warning("anova.polytomous: the following arguments to 'anova.polytomous' are invalid and dropped: ", 
            paste(deparse(dotargs[named]), collapse = ", "))
    dotargs <- dotargs[!named]
    is.polytomous <- unlist(lapply(dotargs, function(x) inherits(x, 
        "polytomous")))
    dotargs <- dotargs[is.polytomous]
    c(list(object), dotargs);
    if (length(dotargs) > 0) 
        return(anova.polytomouslist(c(list(object), dotargs), statistic = statistic, test = test))
    else
    if(!outcome.specific)
       { formula <- object$formula;
         heuristic <- object$heuristic;
         data <- object$data;
         frequency <- object$frequency;
         cues <- attr(terms(formula),"term.labels") # N.B. if interactions are determined in the formula, such terms are also extracted
         n.cues = length(cues)
         outcome <- as.character(formula)[2] 

         formula.model <- vector("character",n.cues+1); statistics <- statistic.difference <- significance <- vector("numeric",n.cues+1); df <- df.difference <- vector("integer",n.cues+1);
         formula.model[1] <- "NULL"; statistic.difference <- df.difference <- significance <- NA;

         for(i.cue in 1:n.cues)
           { formula.model[i.cue+1] <- paste(c(outcome, paste(cues[1:i.cue], collapse=" + ")), collapse=" ~ ")
             summary.model <- polytomous(as.formula(formula.model[i.cue+1]), data, frequency, heuristic=heuristic)$statistics
             statistics[i.cue+1] <- summary.model[[statistic.model]]
             df[i.cue+1] <- summary.model[["df.model"]];
             if(i.cue==1)
               { statistics[1] <- summary.model[[statistic.null]]; df[1] <- summary.model[["df.null"]]; }
                 statistic.difference[i.cue+1] = statistics[i.cue]-statistics[i.cue+1]; df.difference[i.cue+1] = df[i.cue]-df[i.cue+1];
                 significance[i.cue+1] = 1-pchisq(statistic.difference[i.cue+1], df.difference[i.cue+1])
               }     
         if(!is.null(test))
           if(test=="Chisq")
           { effect.cues <- data.frame(df, statistics, df.difference, statistic.difference, significance)
             colnames(effect.cues) <- c("Resid. Df",paste("Resid. ",statistic,sep=""),"Df",statistic,"P(>|Chi|)");
             rownames(effect.cues) <- formula.model;
             title <- paste("Analysis of ",statistic," Table", "\n\nModel: polytomous, heuristic: ", heuristic, 
                      "\n\nResponse: ", outcome, "\n\nTerms added sequentially (first to last)\n\n", sep = "")
           }
         if(is.null(test))
           { effect.cues <- data.frame(df, statistics, df.difference, statistic.difference)
             colnames(effect.cues) <- c("Resid. Df",paste("Resid. ",statistic,sep=""),"Df",statistic);
             rownames(effect.cues) <- formula.model;
             title <- paste("Analysis of ",statistic," Table", "\n\nModel: polytomous, heuristic: ", heuristic, 
                      "\n\nResponse: ", outcome, "\n\nTerms added sequentially (first to last)\n\n", sep = "")
           }

         return(structure(effect.cues, heading = title, class = c("anova", "data.frame")))
       }
     else
       { if(object$heuristic!="one.vs.rest")
           stop("'outcome.specific=TRUE' only applicable for a 'polytomous' object fit using the 'one.vs.rest' heuristic.")

         model <- lapply(object$model, function(x) anova.glm(x, dispersion = NULL, test = test));
         outcomes <- object$outcomes;
         n.outcomes <- length(outcomes);
         deviance <- df <- p.values <- vector("list", n.outcomes);
         names(deviance) <- names(df) <- names(p.values) <- outcomes;

         deviance <- sapply(model, function(x) deviance <- x$Deviance);
         df <- sapply(model, function(x) Df <- x$Df);
         if(!is.null(test))
           p.values <- sapply(model, function(x) p.values <- x[[5]])
         else
           p.values <- matrix(NA,NROW(model[[1]]),n.outcomes,dimnames=list(NULL,outcomes));

         rownames(deviance) <- rownames(df) <- rownames(p.values) <- rownames(model[[1]]);
         ANOVA <- list(model = model, deviance = deviance, df = df, p.values = p.values);
                  
         return(ANOVA) 
       }
}

