anova.polytomouslist <-
function(object, ..., statistic = "deviance", test = "Chisq") 
{ 
    if(!is.null(names(object)))
      stop("'object' not a proper list of 'polytomous' models.")

    responses <- as.character(lapply(object, function(x) {
        deparse(x$formula[[2L]])
    }))
    sameresp <- responses == responses[1L]
    if (!all(sameresp)) {
        object <- object[sameresp]
        warning("models with response ", deparse(responses[!sameresp]), 
            " removed because response differs from model 1")
    }

    nmodels <- length(object)
    if (nmodels==1) 
        return(anova.polytomous(object[[1]], statistic = statistic, test = test))

    if(statistic=="deviance")
      statistic.model="deviance.model"
    else   
    if(statistic=="AIC")
      statistic.model="AIC.model"
    else
    if(statistic=="BIC")
      statistic.model="BIC.model"
    else
      stop(paste("'statistic': ",statistic," not implemented.",sep=""))

    heuristics <- lapply(object, function(x) x$heuristic)
    sameheuristic <- heuristics == heuristics[[1]]
    if(!all(sameheuristic))
      stop("models not fit with the same heuristic.")

    data.formats <- lapply(object, function(x) x$data.format)
    same.data.format <- data.formats == data.formats[[1]]
    if(!all(same.data.format))
      stop("models not fit with similarly formatted data.")

    if(data.formats[[1]]=="instance")
      { ns <- sapply(object, function(x) NROW(x$data))
        if(any(ns != ns[1])) 
          stop("models were not all fitted to the same size of dataset.")
      }
    if(data.formats[[1]]=="narrow")
      { frequencies <- lapply(object, function(x) x$frequency)
        if(length(frequencies[[1]])==1 & is.character(frequencies[[1]]))
          { samefrequencies <- frequencies == frequencies[[1]]
            if(!all(samefrequencies))
              stop("models were not all fitted with similarly formatted dataset.")
            same.cols <- sapply(object, function(x) frequencies[[1]] %in% colnames(x$data))
            if(!all(same.cols))
              stop("models were not all fitted with same datasets.")
            sums <- sapply(object, function(x) sum(x$data[,frequencies[[1]]]))
            same.sums <- sums == sums[1]
            if(!all(same.sums))
              stop("models were not all fitted to the same size of dataset.")
          }
        if(is.numeric(frequencies[[1]]))
          { sums <- sapply(frequencies, sum)
            same.sums <- sums == sums[1]
            if(!all(same.sums))
              stop("models were not all fitted to the same size of dataset.")
          }
      }
    if(data.formats[[1]]=="wide")
      { response.cols <- gsub("[ ]*","",strsplit(responses[[1]],"\\|")[[1]])
        all.responses <- sapply(object, function(x) all(response.cols %in% colnames(x$data)))
        if(!all(all.responses))
          stop("models were not all fitted with similarly formatted responses in datasets.")
        ns <- sapply(object, function(x) sum(x$data[,response.cols]))
        if(any(ns != ns[1]))
          stop("models were not all fitted to the same size of dataset")
      }

    df <- unlist(as.numeric(lapply(object, function(x) x$statistics$df.model)))
    statistics <- as.numeric(lapply(object, function(x) x$statistics[[statistic.model]]))
    formula.model <- unlist(lapply(object, function(x) { f <- as.character(x$formula); return(paste(c(f[2]," ",f[1]," ",f[3]),collapse=""))}))
    df.difference <- c(NA,df[1:(nmodels-1)]-df[2:nmodels])
    statistic.difference <- c(NA,statistics[1:(nmodels-1)]-statistics[2:nmodels])
    significance <- c(NA, sapply(1:(nmodels-1), function(i) 1-pchisq(statistic.difference[i+1], df.difference[i+1])))
 
    if(!is.null(test))
      if(test=="Chisq")
      { effects.models <- data.frame(df, statistics, df.difference, statistic.difference, significance);
        colnames(effects.models) <- c("Resid. Df","Resid. Dev.","Df","Deviance","P(>|Chi|)");
      }
    if(is.null(test))
      { effects.models <- data.frame(df, statistics, df.difference, statistic.difference);
        colnames(effects.models) <- c("Resid. Df",paste("Resid. ",statistic,sep=""),"Df",statistic);
      }
    nmodels <- length(object)
    rownames(effects.models) <- paste("Model ", seq(1:nmodels), sep="");

    title <- paste("Analysis of ",statistic," Table\n")
    topnote <- paste("Model ", format(1:nmodels), ": ", formula.model, sep = "", collapse = "\n")
     
    structure(effects.models, heading=c(title, topnote), class =c("anova", "data.frame"))
}

