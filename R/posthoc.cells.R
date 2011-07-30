posthoc.cells <- function(ctable, alpha=0.05, reorder="both", std.pearson.residual.min=2)
{
  if(reorder %in% c("rows","both"))
    ctable <- ctable[order(apply(ctable,1,sum),decreasing=TRUE),] 
  if(reorder %in% c("cols","both"))
    ctable <- ctable[,order(apply(ctable,2,sum),decreasing=TRUE)] 

  chisq.result <- suppressWarnings(chisq.test(ctable))
  X2.ctable <- as.vector(chisq.result$statistic)
  
  df=max((ncol(ctable)-1)*(nrow(ctable)-1),1)

  X2.df1 <- qchisq(alpha, df=1, lower=FALSE)

  X2.df <- qchisq(alpha, df=df, lower=FALSE)

  X2 <- chisq.result$residuals^2 * sign(chisq.result$residuals)

  cells.X2.df <- apply(X2, c(1,2), function(x) ifelse(abs(x)<=X2.df,"0",ifelse(x>X2.df,"+","-")))

  cells.X2.df1 <- apply(X2, c(1,2), function(x) ifelse(abs(x)<=X2.df1,"0",ifelse(x>X2.df1,"+","-")))

  std.pearson.residuals <- t(t(chisq.result$residuals/sqrt(1-apply(ctable,1,sum)/sum(ctable)))/sqrt(1-apply(ctable,2,sum)/sum(ctable)))

  std.pearson.residuals.sign <- apply(std.pearson.residuals,c(1,2),function(x) ifelse(abs(x)<=std.pearson.residual.min,"0",ifelse(x>std.pearson.residual.min,"+","-")))

  result <- list(ctable = ctable, X2.df1 = X2.df1, X2.df=X2.df, cells = list(X2 = X2, X2.df = data.frame(cells.X2.df), X2.df1 = data.frame(cells.X2.df1), std.pearson.residuals = std.pearson.residuals, std.pearson.residuals.sign = data.frame(std.pearson.residuals.sign)))

  return(result)
}