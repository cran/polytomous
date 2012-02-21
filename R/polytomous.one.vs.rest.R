polytomous.one.vs.rest <-
function(formula, data, frequency=NA, p.normalize=TRUE, ...)
{ equation <- gsub("[ ]+"," ",paste(deparse(formula[[3]],width.cutoff=500),collapse=""))
  response <- as.character(formula)[2];
  if(length(grep("\\|",response)!=0))
    { outcomes <- strsplit(response,"[ ]*\\|[ ]*")[[1]]; data.format="wide"; n.data=sum(data[outcomes]);
      if(length(outcomes)>length(unique(outcomes)))
        stop(paste(c("Duplicate response classes in 'formula': ",paste(outcomes,collapse=", ")),collapse=""))
    }
  else
  if(!is.na(frequency))
    { if((length(frequency)==1 & frequency %in% colnames(data)) | all(is.numeric(frequency)))
        return(polytomous.poisson.reformulation(formula, data, frequency))
    }
  else
    { outcomes <- levels(data[[response]]); data.format="instance"; outcome <- response; n.data=nrow(data); } 
  n.outcomes <- length(outcomes);
  predictors <- attr(terms.formula(formula),"term.labels");
  if(length(grep("\\|",predictors))!=0)
    stop("Mixed effects logistic regression not implemented for heuristic: one.vs.rest")
  variables <- levels(as.factor(unlist(lapply(attr(terms.formula(formula),"term.labels"),function(x) strsplit(x,"[ ]*([\\+]|[\\|]|[:])[ ]*")))))

  model <- vector("list", n.outcomes)
  for(m in 1:n.outcomes)
     { if(data.format=="instance")
         { fn <- as.formula(paste(c(outcomes[m],equation),collapse=" ~ "));
           model[[m]] <- glm(fn, data.frame(data[variables],matrix(data[outcome]==outcomes[m],,1,dimnames=list(NULL,outcomes[m]))), family=binomial);
         }
       if(data.format=="narrow")
         { fn <- as.formula(paste(c(outcomes[m],equation),collapse=" ~ "));
           model[[m]] <- glm(fn, data.frame(data[variables],matrix(data[outcome]==outcomes[m],,1,dimnames=list(NULL,outcomes[m]))), weights=data[,frequency], family=binomial);
         }
       if(data.format=="wide")
         { outcome.vs.others <- matrix(c(data[[outcomes[m]]],apply(data[setdiff(outcomes,outcomes[m])],1,sum)),,2,dimnames=list(NULL,c(outcomes[m],"OTHERS")));
           fn <- as.formula(paste(c("outcome.vs.others",equation),collapse=" ~ "));
           model[[m]] <- glm(fn, data[variables], family=binomial)
         }
     }
  names(model) <- outcomes;

 logodds <- p.values <- predictions <- NULL;

 for(m in 1:n.outcomes)
    { logodds.model <- as.matrix(coef(model[[m]])); colnames(logodds.model) <- outcomes[m];
      logodds <- merge(logodds,logodds.model,by=0,all=TRUE); rownames(logodds) <- logodds[,"Row.names"]; logodds <- logodds[-1];

      p.values.model <- as.matrix(summary(model[[m]])$coefficients[,"Pr(>|z|)"]); colnames(p.values.model) <- outcomes[m];
      p.values <- merge(p.values,p.values.model,by=0,all=TRUE); rownames(p.values) <- p.values[,"Row.names"]; p.values <- p.values[-1];

      predictions <- cbind(predictions, fitted(model[[m]]));
    }

# Alternative way to combine results of outcome-specific binary models - might not work properly if some coefficients cannot be estimated for some outcome models
#
#    for(m in 1:n.outcomes)
#       { logodds <- cbind(logodds,coef(model[[m]]));
#         p.values <- cbind(p.values, summary(model[[m]])$coefficients[,"Pr(>|z|)"]);
#         predictions <- cbind(predictions, fitted(model[[m]]));
#       }

  odds <- exp(logodds);

  colnames(logodds) <- colnames(odds) <- colnames(p.values) <- colnames(predictions) <- outcomes;

  if(p.normalize)
    { p.sums <- apply(predictions,1,sum); predictions <- predictions/p.sums; }

  outcomes.predicted <- apply(predictions,1,function(x) names(which.max(x)));

  if(data.format=="instance")  
    { outcomes.predicted <- apply(predictions,1,function(x) names(which.max(x)));
      outcomes.original <- data[[outcome]]
    }
  if(data.format=="wide")
    { outcomes.predicted <- apply(predictions,1,function(x) names(which.max(x)));
      outcomes.original <- data[outcomes];
    }

  statistics <- model.statistics(observed=outcomes.original, predicted=outcomes.predicted, p.values=predictions, frequency=NA, outcomes=outcomes, ...)
  df.null = n.outcomes * n.data;
  df.model <- df.null - length(which(as.vector(apply(logodds,c(1,2), function(x) !is.na(x)))));
  n.predictors <- length(which(as.vector(apply(logodds,c(1,2), function(x) !is.na(x)))))
  AIC.model <- 2 * n.predictors - 2 * statistics$loglikelihood.model;
  BIC.model <- n.predictors * log(n.data) - 2 * statistics$loglikelihood.model;
  statistics <- c(list(df.null=df.null, df.model=df.model, AIC.model=AIC.model, BIC.model=BIC.model), statistics)

  model <- list(model = model, data = data, frequency = frequency, logodds = logodds, odds = odds, p.values = p.values, fitted = predictions, statistics=statistics, formula = formula, outcomes = outcomes, heuristic = "one.vs.rest", data.format = data.format);

  class(model) <- c("polytomous","one.vs.rest")

  return(model);
}

