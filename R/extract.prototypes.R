extract.prototypes <- function(model.polytomous, p.critical=.05)
{ 
   if(!"polytomous" %in% class(model.polytomous))
     stop("Argument 'model.polytomous' not of class 'polytomous'.")

   odds <- model.polytomous$odds
   p.values <- model.polytomous$p.values
   outcomes <- colnames(odds)
   n.outcomes = NCOL(odds)

   prototypes <- sapply(outcomes,
     function(x)
     { properties <- odds[which(p.values[,x]<p.critical),x]
       names(properties) <- rownames(odds)[which(p.values[,x]<p.critical)]
       properties <- t(t(sort(properties,decreasing=TRUE)))
       colnames(properties) = "Odds"
       return(properties)
     })

   return(prototypes)

}   
