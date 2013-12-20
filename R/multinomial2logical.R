multinomial2logical <- function(data, outcome=NULL, variables=NULL, variable.value.separator="")
{ 
  data <- as.data.frame(data)

  if(!is.null(outcome) & !is.null(variables))
    if(outcome %in% variables)
      stop(paste("'outcome': ",outcome," also among 'variables'.",sep=""))

  if(!is.null(outcome))
    if(outcome %in% colnames(data))
      { if(class(data[[outcome]])!="factor")
          stop(paste("'outcome': ",outcome," in 'data' not a factor.",sep=""))
        outcome.multinomial <- data[outcome]
        data <- data[-which(colnames(data) %in% outcome)]
      }
    else
      stop(paste("'outcome': ",outcome," not a column in 'data'.",sep=""))
   
  if(!is.null(variables))
    if(all(variables %in% colnames(data)))
      data <- data[variables]
    else
      stop("One or more 'variables' not columns in 'data'.")

  if(is.null(variables))
    if(length(data)==0)
      stop("No variable columns left in 'data' after 'outcome' is excluded.")

  if(!all(sapply(data, class)=="factor"))
     stop("One or more 'variables' in 'data' not a factor.")

  data.logical <- data.frame(lapply(data, function(y) sapply(levels(y), function(x) y==x)))
    
#  if(length(data)==1)
    data.logical <- data.frame(lapply(data, function(y) sapply(levels(y), function(x) y==x)))
#  else
#    data.logical <- data.frame(sapply(data, function(y) sapply(levels(y), function(x) y==x),simplify=FALSE))

  names(data.logical) <- unlist(lapply(names(data), function(x) paste(x, levels(data[[x]]), sep=variable.value.separator)))

  if(!is.null(outcome))
     data.logical <- cbind(outcome.multinomial, data.logical)

  return(data.logical)

}