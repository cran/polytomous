model.statistics <-
function(observed, predicted, p.values, frequency=NA, outcomes=NULL, p.normalize=TRUE, cross.tabulation=TRUE, p.zero.correction=1/(nrow(p.values)*ncol(p.values))^2, ...)
{ if(p.zero.correction==0) warning("Loglikelihood and related statistics may be inestimable, if P=0 for any observed outcome.");
  dim.observed <- dim(observed)
  if(length(dim.observed)>1)
    { if(!identical(dim.observed,dim(p.values)))
        stop("Tables of observed outcome counts and estimated p-values do not match.")
      if(dim.observed[1]!=length(as.vector(predicted)))
        stop("Table of observed outcome counts and list of predicted outcomes do not match.")
      N <- sum(observed);
      if(is.null(outcomes))
        outcomes <- colnames(observed);
    }
  else
    { if(length(as.vector(observed))!=length(as.vector(predicted)))
        stop("Lists of observed and predicted outcomes do not match.")
      if(length(as.vector(observed))!=nrow(p.values))
        stop("List of observed outcomes and table of estimated p-values do not match.")
      if(all(!is.na(frequency)) & all(is.numeric(frequency)))
        if(length(as.vector(observed))!=length(as.vector(frequency)))
          stop("List of observed outcomes and frequencies do not match.")
      if(is.null(outcomes))
        outcomes <- levels(as.factor(observed));
      if(all(is.na(frequency)))
        N <- NROW(observed)
      else
        N <- sum(as.vector(frequency))
    }

  if(p.normalize)
    p.values <- p.values/apply(p.values,1,sum)
  p.outcomes <- colnames(p.values);

  if(length(dim.observed)>1)
    { observed <- observed[,outcomes];
      p.values <- apply(p.values[,outcomes], c(1,2), function(p) if(p==0) p.zero.correction else p);
      d <- as.vector(observed * log(p.values));
      n.outcomes <- apply(observed,2,sum);
    }
  else
    { d <- sapply(1:NROW(observed), function(i)
        { p <- p.values[i,which(observed[i]==p.outcomes)]
          if(p==0) p = p.zero.correction
          if(all(is.na(frequency)))
             return(log(p))
          else
             return(as.vector(frequency)[i]*log(p))
        })
      if(all(is.na(frequency)))
        n.outcomes <- sapply(outcomes, function(o) length(which(observed==o)))
      else
        n.outcomes <- sapply(outcomes, function(o) sum(as.vector(frequency)[which(observed==o)]))
    }

  loglikelihood.model <- sum(d);
  deviance.model <- -2 * loglikelihood.model;
  loglikelihood.null <- sum(sapply(outcomes, function(o) n.outcomes[o]*log(n.outcomes[o]/N)));
  deviance.null <- -2 * loglikelihood.null;
  R2.likelihood <- 1 - deviance.model/deviance.null;
#  R2.nagelkerke <- (1-(exp(loglikelihood.null)/exp(loglikelihood.model))^(2/N))/(1-(exp(loglikelihood.null)^(2/N)));
  R2.nagelkerke <- (1-exp(-2*(loglikelihood.model-loglikelihood.null)/N))/(1-exp(2*loglikelihood.null/N))
  statistics <- list(loglikelihood.null = loglikelihood.null, loglikelihood.model = loglikelihood.model, deviance.null = deviance.null, deviance.model = deviance.model, R2.likelihood = R2.likelihood, R2.nagelkerke = R2.nagelkerke);

  if(cross.tabulation)
     { if(length(dim.observed)>1)
         { crosstable <- matrix(0, length(outcomes), length(outcomes), dimnames=list(outcomes,outcomes));
           for(j in 1:nrow(observed))
              for(i in 1:ncol(observed))
                 crosstable[i,as.character(predicted[j])] = crosstable[i,as.character(predicted[j])] + observed[j,i];
         }
      else
        if(all(is.na(frequency)))
          crosstable <- table(factor(observed, levels=outcomes), factor(predicted, levels=outcomes))
        else
          { crosstable <- matrix(0, length(outcomes), length(outcomes), dimnames=list(outcomes,outcomes));
            for(i in 1:NROW(observed))
               crosstable[as.character(observed[i]),as.character(predicted[i])] = crosstable[as.character(observed[i]),as.character(predicted[i])] + as.vector(frequency)[i]
          }

       statistics <- c(statistics, list(crosstable = crosstable), crosstable.statistics(crosstable))
     }
  return(statistics);
}

