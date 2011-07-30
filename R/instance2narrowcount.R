instance2narrowcount <-
function(data.table, variables, outcome="OUTCOME", variables.default=NULL, outcome.ordered=NULL, numeric2discrete=function(x) cut2(x,levels.mean=TRUE,g=g.numeric), g.numeric=2)
{ require(Hmisc, quietly=TRUE)
  if(length(variables)>length(unique(variables)))
    stop(paste("Duplicate terms in 'variables': ",paste(variables,collapse=", "),sep=""));
  if(!all(variables %in% names(data.table)))
    stop("One or more 'variables' not in 'data.table'");
  if(!all(names(variables.default) %in%  names(data.table)))
    stop("One or more 'variables.default' not in 'data.table'");
  if(!(outcome %in% names(data.table)))
    stop(paste("Response variable 'outcome': ",outcome," not in 'data.table'",sep=""));
  if(!is.null(outcome.ordered))
    if(!all(outcome.ordered %in% levels(data.table[[outcome]])))
      stop(paste("One or more values of 'outcome.ordered' not a category of 'outcome': ",outcome,sep=""));
  if(!is.null(outcome.ordered))
    if(!all(levels(data.table[[outcome]]) %in% outcome.ordered))
      stop(paste("One or more values of 'outcome': ",outcome," missing from 'outcome.ordered'",sep=""));
  if(!is.null(variables.default))
    for(i in 1:length(variables.default))
       if(!(variables.default[[i]]) %in% data.table[[names(variables.default)[i]]])
         stop("One or more default classes in 'variables.default' not classes of the corresponding 'variables'");

  if(!is.null(outcome.ordered))
    outcomes <- outcome.ordered
  else
    outcomes <- levels(data.table[[outcome]]);

# Discretize numeric predictors

  variables.numeric <- names(which(sapply(variables, function(i) is.numeric(data.table[[i]]))));
  if(length(variables.numeric)>0)
    for(i in 1:length(variables.numeric))
       if(is.null(numeric2discrete))
         data.table[[variables.numeric[i]]] <- factor(data.table[[variables.numeric[i]]])
       else
         data.table[variables.numeric[i]] <- lapply(data.table[variables.numeric[i]],numeric2discrete)

  variables.logical <- names(which(sapply(data.table,is.logical)))
  data.table[variables.logical] <- lapply(variables.logical, function(x) factor(ifelse(data.table[[x]],x,NA)))

  variables.levels <- lapply(variables, function(i) levels(data.table[[i]]))
  variables.ordered <- sapply(variables, function(i) is.ordered(data.table[[i]]));
  n.outcomes <- length(outcomes);

  category.combinations.table <- data.frame(data.table[outcome],apply(data.table[variables],1,function(x) paste(x, collapse="::")));
  category.combinations <- levels(category.combinations.table[,2]);
  count.table <- data.frame(matrix(
    unlist(lapply(category.combinations,
      function(cc)
        lapply(outcomes,
           function(oc)
c(length(intersect(which(category.combinations.table[,1]==oc),which(category.combinations.table[,2]==cc)))/length(which(category.combinations.table[,2]==cc)),length(intersect(which(category.combinations.table[,1]==oc),which(category.combinations.table[,2]==cc))),oc,strsplit(cc,"::")[[1]],which(category.combinations==cc))))),,length(variables)+4,byrow=TRUE));

  colnames(count.table) <- c("Proportion","Count",outcome,variables,"Observation");
  count.table$Proportion <- as.numeric(as.character(count.table$Proportion));
  count.table$Count <- as.integer(as.character(count.table$Count));
  for(i in 1:length(variables))
     if(!variables[i] %in% variables.logical)
       count.table[[variables[i]]] <- factor(as.character(count.table[[variables[i]]]), levels=variables.levels[[i]], ordered=variables.ordered[[variables[[i]]]])
     else
       count.table[variables[i]] <- ifelse(count.table[[variables[i]]]==variables[i],TRUE,FALSE)
  if(!is.null(variables.default))
    for(i in 1:length(variables.default))
       count.table[[names(variables.default)[i]]] <- relevel(count.table[[names(variables.default)[i]]], variables.default[[i]]);
  if(length(variables.numeric)>0)
    for(i in 1:length(variables.numeric))
       count.table[[variables.numeric[i]]] <- as.numeric(as.character(count.table[[variables.numeric[i]]]));

  count.table[[outcome]] <- factor(as.character(count.table[[outcome]]), levels=outcomes);

  count.table;

}

