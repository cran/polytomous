predict.polytomous <-
function(object, newdata=NULL, type="response", p.normalize=TRUE, ...)
{ if(is.null(newdata))
     newdata <- object$data
  if(object$heuristic=="one.vs.rest")
  { outcomes <- object$outcomes;
    model <- object$model;
    predictors <- names(coef(model[[1]]));
    n.predictors <- length(predictors);
    if(type=="link")
      predictions <- sapply(outcomes, function(w) predict.glm(model[[w]],newdata,type="link"))
    if(type %in% c("response","choice"))
      { predictions <- sapply(outcomes, function(w) predict.glm(model[[w]],newdata,type="response"));
        if(is.null(dim(predictions)))
          { predictions <- data.frame(t(as.matrix(predictions))); colnames(predictions) <- outcomes; }
      }
    if(type=="terms")
      { predictions <- lapply(outcomes, function(w) predict.glm(model[[w]],newdata,type))
        names(predictions) <- outcomes
        predictions <- data.frame(predictions)
      }
    if(type %in% c("response","choice") & p.normalize)
      { p.sums <- apply(predictions,1,sum); predictions <- predictions/p.sums; }
    if(type=="choice")
      { outcomes.predicted <- apply(predictions,1,function(x) names(which.max(x)));
        names(outcomes.predicted) <- rownames(newdata)
        return(outcomes.predicted)
      }
    else
      return(predictions);
  }

  if(object$heuristic=="poisson.reformulation")
    { if(type %in% c("link","terms"))
        stop(paste("'type': ",type," not implemented for 'heuristic': ",object$heuristic,sep=""))
      outcomes <- object$outcomes;
      n.outcomes <- length(outcomes);
      data.poisson <- object$data.poisson;
      probabilities.matrix <- matrix(object$fitted.poisson[,"probabilities"],,n.outcomes,byrow=TRUE);
      colnames(probabilities.matrix) <- outcomes;
      predictors <- attr(terms.formula(object$formula),"term.labels");
      predictors.fixed <- grep("\\|",predictors,value=TRUE,invert=TRUE)
      n.newdata <- nrow(newdata);

      predictors.fixed.combinations <- apply(data.poisson[predictors.fixed],1,function(x) paste(x, collapse="::"))
      variables2observations <- lapply(unique(predictors.fixed.combinations), function(x) as.numeric(as.character(unique(data.poisson$Observation[which(x==predictors.fixed.combinations)]))))
      names(variables2observations) <- unique(predictors.fixed.combinations)
#      variables2observations <- rle(apply(data.poisson[variables],1,function(x) paste(x, collapse="::")))$values;

      predictions <- NULL;
      for(i in 1:n.newdata)
#         predictions <- rbind(predictions, probabilities.matrix[which(variables2observations %in% apply(as.matrix(newdata[i,variables],1),1,function(x) paste(x, collapse="::"))),]);
         { indices <- variables2observations[[apply(as.matrix(newdata[i,predictors.fixed]),1,function(x) paste(x, collapse="::"))]];
           predictions <- rbind(predictions, apply(matrix(probabilities.matrix[indices,],nrow=length(indices)),2,mean))
         }
      colnames(predictions) <- outcomes
      if(type=="choice")
        { outcomes.predicted <- apply(predictions,1,function(x) names(which.max(x)))
          names(outcomes.predicted) <- rownames(newdata)
          return(outcomes.predicted)
        }
      else
        return(data.frame(predictions))
    }
}

