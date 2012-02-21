polytomous.poisson.reformulation <-
function(formula, data, frequency=NA, variables.ordered=NULL, include.Observation=TRUE, ...)
{ response <- as.character(formula)[2];
  if(length(grep("\\|",response)!=0))
    { outcomes <- strsplit(response,"[ ]*\\|[ ]*")[[1]];
      data.format="wide";
      formula.start <- "Count ~ "
    }
  else
  if(!all(is.na(frequency)))
    { outcomes <- levels(data[[response]])
      if(length(frequency)==1 & is.character(frequency))
        { if(frequency %in% colnames(data))
            formula.start <- paste(frequency," ~ ",sep="")
          else
            stop(paste("Column for 'frequency': ",frequency," not in 'data'.",sep=""))
        }
      else
      if(is.numeric(frequency) & length(as.vector(frequency))==NROW(data))
        { Count <- as.vector(frequency);
          data <- cbind(data, Count);
          formula.start <- "Count ~ ";
          }
      else
        stop("Numeric vector 'frequency' does not match 'data'.")
      data.format="narrow";
      outcome <- response; 
    }
  else
    { outcomes <- levels(data[[response]]);
      data.format="instance";
      outcome <- response;
      formula.start <- "Count ~ "
    } 
  n.outcomes <- length(outcomes);
  predictors <- attr(terms.formula(formula),"term.labels");
  if(length(grep("\\|",predictors))!=0)
    model.type="mixed"
  else
    model.type="fixed";
  variables <- levels(as.factor(unlist(lapply(attr(terms.formula(formula),"term.labels"),function(x) strsplit(x,"[ ]*([\\+]|[\\*]|[\\|]|[:])[ ]*")))))
  if(length(grep("^1$", variables))!=0) variables <- variables[-grep("^1$", variables)]

  if(!is.null(variables.ordered))
    { if(length(variables.ordered)>length(unique(variables.ordered)))
        stop("Duplicate terms in 'variables.ordered'")
      if(!all(sort(variables.ordered) %in% variables))
        stop("Not all 'variables.ordered' are variable terms in 'formula'")
      if(!all(variables %in% sort(variables.ordered)))
        stop("Variable terms in 'formula' are missing from 'variables.ordered'")
      variables <- variables.ordered;
    }
# The following is taken care of in 'instance2narrowcount':
#  if(!is.null(outcome.ordered))
#    { if(!(all(outcome.ordered %in% outcomes)))
#        stop("Classes in 'outcome.ordered' missing from outcome term in 'formula'")
#      if(!(all(outcomes %in% outcome.ordered)))
#        stop("Outcome classes of term in 'formula' missing from 'outcome.ordered'")
#      outcomes <- outcome.ordered;
#    }

  if(data.format=="instance")
     { data.poisson <- instance2narrowcount(data, variables, outcome, ...);
       n.data <- nrow(data);
       outcomes <- levels(data.poisson[[outcome]]);
     }
  if(data.format=="wide")
    { data.poisson <- wide2narrowcount(data, variables, outcomes, ...);
      n.data <- sum(data.poisson$Count); outcome=names(data.poisson)[3];
    }
  if(data.format=="narrow")
    { data.poisson <- data;
      if(is.character(frequency))
        n.data <- sum(data[,frequency])
      else
        n.data <- sum(data$Count)
    }


  if(model.type=="fixed")
    { if(include.Observation)
          formula.start <- paste(formula.start," Observation",sep="")
      predictors.fixed <- variables.fixed <- variables;

      if(data.format=="instance" | data.format=="narrow")
        formula.poisson <- paste(c(formula.start, outcome, unlist(lapply(predictors, function(x) paste(c(outcome,x),collapse=":")))),collapse=" + ");
      if(data.format=="wide")
        formula.poisson <- paste(c(formula.start, outcome, unlist(lapply(predictors, function(x) paste(c(outcome,x),collapse=":")))),collapse=" + ");

      model <- try(glm(as.formula(formula.poisson), data.poisson, family=poisson),silent=TRUE);
      if("try-error" %in% class(model))
        stop("Following error in fitting reformulated model using 'glm(..., family=poisson)':\n",model[1])

      if(include.Observation)
        { coefficients <- coef(summary(model))[-grep("^Observation",names(model$coefficients)),]
          coefficients.names <- names(model$coefficients)[-grep("^Observation",names(model$coefficients))];
        }
      else
        { coefficients <- coef(summary(model));
          coefficients.names <- names(model$coefficients);
        }
      predictors.classes <- c("(Intercept)", outcome, unique(gsub(paste(c("^",outcome,"(",paste(outcomes,collapse="|"),"):"),collapse=""),"",coefficients.names)[(n.outcomes+1):length(coefficients.names)]));
      logodds <- odds <- p.values <- matrix(NA, length(predictors.classes), n.outcomes, dimnames=list(predictors.classes, outcomes));
      logodds[1,outcomes[which(!(sapply(outcomes, function(oc) paste(c(outcome,oc),collapse="")) %in% rownames(coefficients)))]] <- coefficients["(Intercept)","Estimate"];
      p.values[1,outcomes[which(!(sapply(outcomes, function(oc) paste(c(outcome,oc),collapse="")) %in% rownames(coefficients)))]] <- coefficients["(Intercept)","Pr(>|z|)"];
      for(i in 1:n.outcomes)
         if(paste(c(outcome,outcomes[i]),collapse="") %in% coefficients.names)
           { logodds[2,outcomes[i]] <- coefficients[paste(c(outcome,outcomes[i]),collapse=""),"Estimate"];
             p.values[2,outcomes[i]] <- coefficients[paste(c(outcome,outcomes[i]),collapse=""),"Pr(>|z|)"];
           }
      for(i in 1:n.outcomes)
         for(j in 3:length(predictors.classes))
           if(paste(c(paste(c(outcome,outcomes[i]),collapse=""),predictors.classes[j]),collapse=":") %in% rownames(coefficients))
           { logodds[j,outcomes[i]] <- coefficients[paste(c(paste(c(outcome,outcomes[i]),collapse=""),predictors.classes[j]),collapse=":"),"Estimate"];
             p.values[j,outcomes[i]] <- coefficients[paste(c(paste(c(outcome,outcomes[i]),collapse=""),predictors.classes[j]),collapse=":"),"Pr(>|z|)"];
           }
      odds <- exp(logodds);
    }

  if(model.type=="mixed")
    { if(include.Observation)
        formula.start <- paste(formula.start,"(1|Observation)",sep="");
      predictors.fixed <- predictors[-grep("\\|",predictors)];
      variables.fixed <- unique(unlist(lapply(predictors.fixed,function(x) strsplit(x,":")[[1]])))
      predictors.random <- unlist(lapply(grep("\\|",predictors,value=TRUE), function(x) paste(c("(",x,")"),collapse="")));
      formula.poisson <- paste(c(formula.start, outcome, unlist(lapply(predictors.fixed, function(x) paste(c(outcome,x),collapse=":"))),predictors.random),collapse=" + ");

      require(lme4, quietly=TRUE);
      model <- try(lmer(formula.poisson, data.poisson, family=poisson),silent=TRUE);
      if("try-error" %in% class(model))
        stop("Following error(s) in fitting reformulated model using 'lmer(..., family=poisson)':\n",model[1])

      coefs.mer <- function(object)
      { 
        fcoef <- object@fixef
        dims <- object@dims
        vcov <- vcov(object)
        corF <- vcov@factors$correlation

        coefs <- cbind(Estimate = fcoef, `Std. Error` = corF@sd)
        if (nrow(coefs) > 0) {
           if (!dims[["useSc"]]) {
               coefs <- coefs[, 1:2, drop = FALSE]
               stat <- coefs[, 1]/coefs[, 2]
               pval <- 2 * pnorm(abs(stat), lower = FALSE)
               coefs <- cbind(coefs, `z value` = stat, `Pr(>|z|)` = pval)
           }
           else {
               stat <- coefs[, 1]/coefs[, 2]
               coefs <- cbind(coefs, `t value` = stat)
           }
        }
        return(coefs)
      }

      coefficients <- coefs.mer(model)
#      coefficients <- slot(summary(model),"coefs");

      predictors.classes <- c("(Intercept)", outcome, unique(gsub(paste(c("^",outcome,"(",paste(outcomes,collapse="|"),"):"),collapse=""),"",rownames(coefficients)[(n.outcomes+1):(nrow(coefficients))])));
      logodds <- odds <- p.values <- matrix(NA, (nrow(coefficients)/n.outcomes)+1, n.outcomes, dimnames=list(predictors.classes, outcomes));
      logodds[2:nrow(logodds),] <- matrix(c(NA,coefficients[2:nrow(coefficients),"Estimate"]),,n.outcomes,byrow=TRUE); logodds[1,1] <- coefficients[1,"Estimate"];
      odds <- exp(logodds);
      p.values[2:nrow(p.values),] <- matrix(c(NA,coefficients[2:nrow(coefficients),"Pr(>|z|)"]),,n.outcomes, byrow=TRUE); p.values[1,1] <- coefficients[1,"Pr(>|z|)"];
    }

  counts <- fitted(model);
  counts.matrix <- matrix(counts,,n.outcomes, byrow=TRUE); 
  probabilities.matrix <- counts.matrix/apply(counts.matrix,1,sum);
  colnames(counts.matrix) <- colnames(probabilities.matrix) <- outcomes;
  probabilities <- as.vector(t(probabilities.matrix));
  fitted.poisson <- cbind(counts, probabilities);

  variables.fixed.combinations <- apply(data.poisson[variables.fixed],1,function(x) paste(as.character(x), collapse="::"))
  variables2observations <- lapply(unique(variables.fixed.combinations), function(x) as.numeric(as.character(unique(data.poisson$Observation[which(x==variables.fixed.combinations)]))))
  names(variables2observations) <- unique(variables.fixed.combinations);
#  variables2observations <- rle(apply(data.poisson[variables],1,function(x) paste(x, collapse="::")))$values;

  transform.data <- function(data, numeric2discrete=function(x) cut2(x,levels.mean=TRUE,g=g.numeric), g.numeric=2, ...)
  { variables.numeric <- names(which(sapply(data, is.numeric)));
    if(length(variables.numeric)>0)
      for(i in 1:length(variables.numeric))
         if(is.null(numeric2discrete))
           data[[variables.numeric[i]]] <- data[[variables.numeric[i]]]
         else
           data[variables.numeric[i]] <- as.numeric(sapply(data[variables.numeric[i]],numeric2discrete))
#    variables.logical <- names(which(sapply(data, is.logical)))
#    if(length(variables.logical)>0)
#      for(i in 1:length(variables.logical))
#         data[[variables.logical[i]]] <- factor(as.character(as.integer(data[[variables.logical[i]]])));
    return(data);
  }

  predictions <- NULL;
  if(data.format=="instance")
    { data.transformed <- transform.data(data[variables.fixed], ...);
      for(i in 1:n.data)
#         predictions <- rbind(predictions, probabilities.matrix[which(variables2observations %in% apply(as.matrix(data[i,variables],1),1,function(x) paste(x, collapse="::"))),]);
         { indices <- variables2observations[[as.vector(apply(as.matrix(data.transformed[i,variables.fixed]),1,function(x) paste(as.character(x), collapse="::")))]];
           predictions <- rbind(predictions, apply(matrix(probabilities.matrix[indices,],nrow=length(indices)),2,mean))
         }
      colnames(predictions) <- outcomes
      outcomes.predicted <- apply(predictions,1,function(x) names(which.max(x)));
      outcomes.original <- data[[outcome]]
      frequencies <- NA
    }
  if(data.format=="wide")
    { predictions <- data.frame(probabilities.matrix)
      outcomes.predicted <- apply(predictions,1,function(x) names(which.max(x)))
      outcomes.original <- data[outcomes]
      frequencies <- NA
    }
  if(data.format=="narrow")
    { predictions <- data.frame(probabilities.matrix)
      predictions <- predictions[rep(1:nrow(predictions), rep(n.outcomes, nrow(predictions))),]
      rownames(predictions) <- seq(1:nrow(predictions))
      outcomes.predicted <- apply(predictions,1,function(x) names(which.max(x)))
      outcomes.original <- data.poisson[[response]]
      if(is.character(frequency))
        frequencies <- data.poisson[,frequency]
      if(is.numeric(frequency))
        frequencies <- frequency
    }

  statistics <- model.statistics(observed=outcomes.original, predicted=outcomes.predicted, p.values=predictions, frequency=frequencies, outcomes=outcomes, ...)
  df.null = n.outcomes * n.data;
  df.model <- df.null - length(which(as.vector(apply(logodds,c(1,2), function(x) !is.na(x)))));
  n.predictors <- length(which(as.vector(apply(logodds,c(1,2), function(x) !is.na(x)))))
  AIC.model <- 2 * n.predictors - 2 * statistics$loglikelihood.model;
  BIC.model <- n.predictors * log(n.data) - 2 * statistics$loglikelihood.model;
  statistics <- c(list(df.null=df.null, df.model=df.model, AIC.model=AIC.model, BIC.model=BIC.model), statistics)

  if(model.type=="mixed")
    { sd.ranef <- lapply(model@ST, function(x) x*lme4:::sigma(model)) # summary(model)@sigma)
      var.ranef <- lapply(sd.ranef, function(x) x^2)
      names(sd.ranef) <- names(var.ranef) <- names(model@flist)
      statistics <- c(statistics, sd.ranef=list(sd.ranef), var.ranef=list(var.ranef))
    }


  model <- list(model = model, data = data, data.poisson = data.poisson, frequency = frequency, fitted.poisson = fitted.poisson, fitted = predictions, coefficients = coefficients, logodds = logodds, odds = odds, p.values = p.values, statistics = statistics, formula = formula, formula.poisson = as.formula(formula.poisson), outcomes = outcomes, heuristic = "poisson.reformulation", data.format = data.format);

  class(model) <- c("polytomous", "poisson.reformulation", model.type);
  return(model);
}

