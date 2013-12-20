nominal <- function(formula, data, sort.bivariate=NULL, std.pearson.residual.min=2, correct=FALSE, report.interval=100, factor2logical=FALSE)
{
  data <- as.data.frame(data)

  dependents = unique(as.vector(sapply(gsub("[ ]+"," ",paste(deparse(formula[[2]],width.cutoff=500),collapse="")), function(x) strsplit(x,"[ ]*([\\*]|[\\+]|[\\|]|[:])[ ]*")[[1]])))

  if(length(dependents)==1)
    { 
      if(dependents==".")
        only.independents=TRUE
      else
        { only.independents=FALSE
          dependents.internal = data[,dependents]
          dependents.values = sort(unique(as.character(data[,dependents])))
        }
    }
  else
    { 
      only.independents=FALSE
      dependents.internal <- unlist(apply(data[,dependents],1,function(x) names(which(x))))
      if(length(dependents.internal)!=NROW(data))
        stop("Dependent variables in 'formula' do not cover entire 'data'.")
      dependents.values = sort(dependents)
    }

  independents = unique(as.vector(sapply(gsub("[ ]+"," ",paste(deparse(formula[[3]],width.cutoff=500),collapse="")), function(x) strsplit(x,"[ ]*([\\*]|[\\+]|[\\|]|[:])[ ]*")[[1]])))
  if(length(independents)==1)
    if(independents==".")
      if(!only.independents)
        independents = colnames(data)[which(!colnames(data) %in% dependents)]
      else
        independents =  colnames(data)
#  independents = levels(as.factor(unlist(lapply(attr(terms.formula(formula, data=data),"term.labels"),function(x) strsplit(x,"[ ]*([\\+]|[\\|]|[:])[ ]*")))))

  independent.classes <- sapply(data[independents], class)
  if(!all(independent.classes=="logical"))
    if(!factor2logical)
      stop("Independent variables in 'data' not binary/logical.")
    else
      if(!all(independent.classes %in% c("factor","character","logical")))
        stop("Independent variables in 'data' to transformable to logical from factor.")
      else
        { data.new <- data[0]
          for(i in 1:ncol(data))
             { if(colnames(data)[i] %in% independents & class(data[[i]])!="logical")
                 data.new <- cbind(data.new, multinomial2logical(data[i]))
               if(colnames(data)[i] %in% independents & class(data[[i]])=="logical")
                 data.new <- cbind(data.new, data[i])
               if(colnames(data)[i] %in% dependents)
                 data.new <- cbind(data.new, data[i])
             }
          data <- data.new; rm(data.new)
          independents <- colnames(data)[!colnames(data) %in% dependents]
          #data <- cbind(data, multinomial2logical(data, variables=independents[which(independent.classes!="logical")]))
          # for(i in colnames(data)[which(colnames(data) %in% independents[which(independent.classes!="logical")])])
          #   data[,i] <- NULL
          #independents <- colnames(data)[!colnames(data) %in% dependents]
        }

  if(only.independents)
    {
      assocs <- NULL
      category1 <- category2 <- NULL
      n.independents = length(independents)
      n.rounds=.5*n.independents*(n.independents+1)
      if(n.rounds>report.interval)
        { cat(n.rounds,": ",sep="")
          ix=1
        }
      else
        ix=0

      if(!is.null(sort.bivariate))
         if(sort.bivariate %in% c("uc","lambda","tau"))
           { sort.key12 = paste(sort.bivariate,"RC",sep=".")
             sort.key21 = paste(sort.bivariate,"CR",sep=".")
           }
         else
           stop("Incorrect argument 'sort.bivariate': ", sort.bivariate)

      for(i in 1:(n.independents-1))
         for(j in (i+1):n.independents)
            {
              cat1 <- independents[i]
              cat2 <- independents[j]
              ctable <- table(data[,cat1], data[,cat2])
              if(!is.null(sort.bivariate))
                { assoc.tmp <- associations(ctable)
                  if(assoc.tmp[[sort.key21]]<assoc.tmp[[sort.key12]])
                    { cat1 <- independents[j]
                      cat2 <- independents[i]
                      ctable <- t(ctable)
                    }
                }
              assoc <- associations(ctable)
              names(assoc) <- gsub("RC","12",names(assoc))
              names(assoc) <- gsub("CR","21",names(assoc))
#              if(!is.null(sort.key))
#                if(sort.key %in% c("uc.RC","uc.CR","lambda.RC","lambda.CR","tau.RC","tau.CR"))
#                { if(length(grep("RC$",sort.key))==1)
#                    sort.key2 = gsub("RC","CR",sort.key)
#                  else
#                    sort.key2 = gsub("CR","RC",sort.key)
#                  assoc.tmp <- associations(ctable)
#                  if(assoc.tmp[[sort.key]]<assoc.tmp[[sort.key2]])
#                    { cat1 <- independents[j]
#                      cat2 <- independents[i]
#                      ctable <- t(ctable)
#                    }
#                 }
#                 else
#                  stop("Incorrect argument 'sort.key': ",sort.key)
              N1 <- length(which(data[,cat1]))
              N2 <- length(which(data[,cat2]))
              N12 <- ctable["TRUE","TRUE"]
              assocs <- rbind(assocs, c(N1=N1, N2=N2, N12=N12, assoc))
              category1 <- c(category1, cat1)
              category2 <- c(category2, cat2)

              if(ix!=0)
                { if((ix %% report.interval)==0)
                    cat("[",ix,"]", sep="")
                  ix=ix+1
                }
            }
    assocs <- data.frame(assocs)
    bivariate <- data.frame(cbind(category1, category2, assocs))

    results = list(bivariate=bivariate, independents=independents)
    class(results) = c("nominal","bivariate")
    return(results)
    }

  univariate <- NULL
  std.pearson.residuals.sign <- std.pearson.residuals <- X2.df.sign <-  X2.df1.sign <- X2 <- NULL
  assocs <- NULL

  data.internal <- data[,c(dependents, independents)]
  independents.internal <- colnames(data.internal)[which(!colnames(data.internal) %in% dependents)]

  independents.class <- sapply(data.internal[,independents.internal], function(x) class(x))
  if(!(all(independents.class=="logical") | all(independents.class!="logical")))
    stop("Independent variables in 'data' are not all of the same type (i.e. logical or NOT logical.")

  n.rounds = length(independents.internal)
  if(n.rounds>report.interval)
    { ix=1
      cat(n.rounds,": ",sep="")
    }
  else
    ix=0

  for(i in independents.internal)
     {
       ctable <- table(data.internal[,i],dependents.internal)
       N <- length(which(data.internal[,i]))
       posthoc <- chisq.posthoc(ctable, reorder="none", std.pearson.residual.min=std.pearson.residual.min, correct=correct)
       assoc <- c(N=N, associations(ctable))
       names(assoc) <- gsub("RC","12",names(assoc))
       names(assoc) <- gsub("CR","21",names(assoc))
       univariate <- c(univariate, list(c(list(posthoc=posthoc), list(assoc=assoc))))

       std.pearson.residuals.sign <- rbind(std.pearson.residuals.sign, apply(posthoc$cells$std.pearson.residuals.sign["TRUE",dependents.values],c(1,2),as.character))
       std.pearson.residuals <- rbind(std.pearson.residuals, posthoc$cells$std.pearson.residuals["TRUE",dependents.values])
       X2.df.sign <- rbind(X2.df.sign, apply(posthoc$cells$X2.df.sign["TRUE",dependents.values],c(1,2),as.character))
       X2.df1.sign <- rbind(X2.df1.sign, apply(posthoc$cells$X2.df1.sign["TRUE",dependents.values],c(1,2),as.character))
       X2 <- rbind(X2, posthoc$cells$X2["TRUE",dependents.values])
       assocs <- rbind(assocs, assoc)

       if(ix!=0)
         { if(ix %% report.interval == 0)
             cat("[",ix,"]",sep="")
           ix=ix+1
         }
     }
  names(univariate) <- independents.internal
  colnames(std.pearson.residuals.sign) <- colnames(std.pearson.residuals) <- colnames(X2.df.sign) <-colnames(X2.df1.sign) <- colnames(X2) <- dependents.values
  rownames(std.pearson.residuals.sign) <- rownames(std.pearson.residuals) <- rownames(X2.df.sign) <- rownames(X2.df1.sign) <- rownames(X2) <- rownames(assocs) <- independents.internal  

  std.pearson.residuals.sign <- data.frame(std.pearson.residuals.sign)
  X2.df.sign <- data.frame(X2.df.sign)
  X2.df1.sign <- data.frame(X2.df1.sign)
  assocs <- data.frame(assocs)

  results <- list(univariate=univariate, std.pearson.residuals.sign=std.pearson.residuals.sign, std.pearson.residuals=std.pearson.residuals, X2.df.sign=X2.df.sign, X2.df1.sign=X2.df1.sign, X2=X2, assocs=assocs, dependents=dependents, dependents.values=dependents.values, independents=independents)

  class(results) <- c("nominal","univariate")
  return(results)
}

print.nominal <- function(x, max.print=10, posthoc="std.pearson.residuals.sign", assoc=ifelse("univariate" %in% class(x), list(c("N","alpha.X2","uc.12","uc.21")), list(c("N1","N2","N12","uc.12","uc.21"))), sort.key=NULL, ...)
{
  sumry.table <- summary(x, posthoc=posthoc, assoc=assoc, sort.key=sort.key, ...)$sumry.table

  if("univariate" %in% class(x))
    {
      dependents = x$dependents
      dependents.values = x$dependents.values
      independents = x$independents

      cat("\n")
      cat("Univariate analysis of categorical variables:\n")
      cat("\n")
      cat("Dependents (2): ", dependents, "=" , paste(dependents.values,collapse=", "), "\n")
      cat("\n")
      cat("Independents (1): "); cat(independents, sep=", ", fill=TRUE)
      cat("\n\n")
    }

  if("bivariate" %in% class(x))
    {
      independents = x$independents
      cat("\n")
      cat("Bivariate analysis of categorical variables:\n")
      cat("\n")
      cat("Independents: "); cat(independents, sep=", ", fill=TRUE)
      cat("\n\n")
    }

   if(!is.na(max.print) & is.numeric(max.print) & nrow(sumry.table)>max.print)
     { print(sumry.table[1:max.print,])
       cat(paste("... [ omitted ",nrow(sumry.table)-max.print," rows ] ...\n",sep=""))
     }
   else
     print(sumry.table)

  if("univariate" %in% class(x) & any(c("uc.12","uc.21","lambda.12","lambda.21","tau.12","tau.21") %in% colnames(sumry.table)))
     {
       cat("\nAssociation measures:\n")
       cat("  statistic.12 ~ statistic(1|2) ~ statistic(independent|dependent) [ ~ statistic(R|C) ]\n")
       cat("  statistic.21 ~ statustic(2|1) ~ statistic(dependent|independent) [ ~ statistic(C|R) ]\n")
     }
  if("bivariate" %in% class(x) & any(c("uc.12","uc.21","lambda.12","lambda.21","tau.12","tau.21") %in% colnames(sumry.table)))
     {
       cat("\nAssociation measures:\n")
       cat("  statistic.12 ~ statistic(1|2) ~ statistic(category1|category2) [ ~ statistic(R|C) ]\n")
       cat("  statistic.21 ~ statustic(2|1) ~ statistic(category2|category1) [ ~ statistic(C|R) ]\n")
     }
  if("univariate" %in% class(x) & !is.null(posthoc))
     cat("\nPosthoc cellwise chi-squared analysis: ", posthoc, "\n")

  cat("\n")
  invisible(x)

}

summary.nominal <- function(object, posthoc="std.pearson.residuals.sign", assoc=ifelse("univariate" %in% class(object), list(c("N","alpha.X2","uc.12","uc.21")), list(c("N1","N2","N12","uc.12","uc.21"))), sort.key=NULL, ...)
{ 
  assoc <- unlist(assoc)
  if("bivariate" %in% class(object))
    assoc <- c("category1","category2",assoc)

  if("univariate" %in% class(object))
    { 
      if(!is.null(posthoc))
        cell.criterion = match.arg(posthoc,c("std.pearson.residuals.sign","std.pearson.residuals","X2.df.sign","X2.df1.sign","X2"))
      else
        cell.criterion = NULL

      if(!all(assoc %in% colnames(object$assocs)))
        stop("Statistic names in 'assoc' not valid results of functions 'chisq.posthoc' and 'associations'.")

      if(is.null(cell.criterion))
        sumry.table <- object$assocs[,assoc]
      else
      { if(cell.criterion=="std.pearson.residuals.sign")
          sumry.table <- cbind(object$assocs[,assoc], object$std.pearson.residuals.sign)
        if(cell.criterion=="std.pearson.residuals")
          sumry.table <- cbind(object$assocs[,assoc], object$std.pearson.residuals)
        if(cell.criterion=="X2.df.sign")
          sumry.table <- cbind(object$assocs[,assoc],object$X2.df.sign)
        if(cell.criterion=="X2.df1.sign")
          sumry.table <- cbind(object$assocs[,assoc], object$X2.df1.sign)
        if(cell.criterion=="X2")
          sumry.table <- cbind(object$assocs[,assoc], object$X2)
      }
      if(!is.null(sort.key))
        { if(!(sort.key %in% colnames(sumry.table)))
            stop("Invalid value of argument 'sort.key': ",sort.key)
          sumry.table <- sumry.table[order(unlist(object$assocs[,sort.key]),decreasing=TRUE),]
        }
    }

  if("bivariate" %in% class(object))
    {
      sumry.table <- object$bivariate[,assoc]
      if(!is.null(sort.key))
        {
          if(!(sort.key %in% assoc))
            stop("Invalid value of argument 'sort.key': ",sort.key)
          sumry.table <- sumry.table[order(unlist(object$bivariate[,sort.key]),decreasing=TRUE),]
        }
    }

   object$sumry.table <- sumry.table
   class(object) <- c("summary.nominal", class(object))
   attr(object, "posthoc") <- posthoc
   return(object)

}

print.summary.nominal <- function(x, max.print=10, ...)
{
  if("univariate" %in% class(x))
    {
      dependents = x$dependents
      dependents.values = x$dependents.values
      independents = x$independents

      cat("\n")
      cat("Univariate analysis of categorical variables:\n")
      cat("\n")
      cat("Dependents (2): ", dependents, "=" , paste(dependents.values,collapse=", "), "\n")
      cat("\n")
      cat("Independents (1): "); cat(independents, sep=", ", fill=TRUE)
      cat("\n\n")
    }

  if("bivariate" %in% class(x))
    {
      independents = x$independents
      cat("\n")
      cat("Bivariate analysis of categorical variables:\n")
      cat("\n")
      cat("Independents: "); cat(independents, sep=", ", fill=TRUE)
      cat("\n\n")
    }

   sumry.table <- x$sumry.table
   posthoc <- attr(x, "posthoc")

   if(!is.na(max.print) & is.numeric(max.print) & nrow(sumry.table)>max.print)
     { print(sumry.table[1:max.print,])
       cat(paste("... [ omitted ",nrow(sumry.table)-max.print," rows ] ...\n",sep=""))
     }
   else
     print(sumry.table)
  cat("\n")

  if("univariate" %in% class(x) & any(c("uc.12","uc.21","lambda.12","lambda.21","tau.12","tau.21") %in% colnames(sumry.table)))
     {
       cat("\nAssociation measures:\n")
       cat("  statistic.12 ~ statistic(1|2) ~ statistic(independent|dependent) [ ~ statistic(R|C) ]\n")
       cat("  statistic.21 ~ statustic(2|1) ~ statistic(dependent|independent) [ ~ statistic(C|R) ]\n")
     }
  if("bivariate" %in% class(x) & any(c("uc.12","uc.21","lambda.12","lambda.21","tau.12","tau.21") %in% colnames(sumry.table)))
     {
       cat("\nAssociation measures:\n")
       cat("  statistic.12 ~ statistic(1|2) ~ statistic(category1|category2) [ ~ statistic(R|C) ]\n")
       cat("  statistic.21 ~ statustic(2|1) ~ statistic(category2|category1) [ ~ statistic(C|R) ]\n")
     }
  if("univariate" %in% class(x) & !is.null(posthoc))
     cat("\nPosthoc cellwise chi-squared analysis: ", posthoc, "\n")

  cat("\n")
  invisible(x)

}
