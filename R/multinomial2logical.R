multinomial2logical <- function(data, outcome=NULL, variables=NULL)
{ 
  data <- as.data.frame(data)
  if(!is.null(outcome))
    if(outcome %in% colnames(data))
      { outcome.multinomial <- data[outcome]
        data <- data[-which(colnames(data) %in% outcome)]
      }
    else
      stop(paste("'outcome': ",outcome," not a column in 'data'.",sep=""))
   
  if(!is.null(variables))
    if(all(variables %in% colnames(data)))
      data <- data[variables]
    else
      stop("One or more 'variables' not columns in 'data'.")
  
  if(length(data)==1)
    data.logical <- data.frame(lapply(data, function(y) sapply(levels(y), function(x) y==x)))
  else
    data.logical <- data.frame(sapply(data, function(y) sapply(levels(y), function(x) y==x)))

  if(!is.null(outcome))
     data.logical <- cbind(outcome.multinomial, data.logical)

  return(data.logical)

}