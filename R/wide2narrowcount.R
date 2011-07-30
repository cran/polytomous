wide2narrowcount <-
function(data.table, variables, outcomes, outcome="OUTCOME", variables.default=NULL, outcome.ordered=NULL)
{ if(length(variables)>length(unique(variables)))
    stop(paste("Duplicate terms in 'variables': ", paste(variables,collapse=", "),sep=""));
  if(!all(variables %in% names(data.table)))
    stop("One or more 'variables' not in 'data.table'");
  if(!all(names(variables.default) %in%  names(data.table)))
    stop("One or more 'variables.default' not in 'data.table'");
  if(!(all(outcomes %in% names(data.table))))
    stop("One or more values of 'outcomes' not in 'data.table'");
  if(!is.null(outcome.ordered))
    if(!all(outcome.ordered %in% outcomes))
      stop("One or more values of 'outcome.ordered' not a category of 'outcomes'");
  if(!is.null(outcome.ordered))
    if(!all(outcomes %in% outcome.ordered))
      stop("One or more values of 'outcome' missing from 'outcome.ordered'");
  if(!is.null(variables.default))
    for(i in 1:length(variables.default))
       if(!(variables.default[[i]]) %in% data.table[[names(variables.default)[i]]])
         stop("One or more default classes in 'variables.default' not classes of the corresponding 'variables'");

  if(!is.null(outcome.ordered))
    outcomes <- outcome.ordered;

  variables.numeric <- names(which(sapply(variables, function(i) is.numeric(data.table[[i]]))));

  if(length(variables.numeric)>0)
    for(i in 1:length(variables.numeric))
       data.table[[variables.numeric[i]]] <- factor(data.table[[variables.numeric[i]]]);

  variables.levels <- lapply(variables, function(i) levels(data.table[[i]]))
  variables.ordered <- sapply(variables, function(i) is.ordered(data.table[[i]]));
  n.variables <- length(variables);
  n.outcomes <- length(outcomes);
  n.observations <- nrow(data.table);
  
  data.table <- apply(data.table, c(1,2), as.character); count.table <- NULL; 
  for(obs in 1:n.observations)
     { sum.outcomes <- sum(as.numeric(data.table[obs,outcomes]));
       if(sum.outcomes>0)
         for(out in 1:n.outcomes)
            count.table <- rbind(count.table, c(as.character(as.numeric(data.table[obs,outcomes[out]])/sum.outcomes), as.character(data.table[obs,outcomes[out]]), outcomes[out], as.character(data.table[obs,variables]), as.character(obs)))
     }

  count.table <- data.frame(count.table);
  colnames(count.table) <- c("Proportion","Count",outcome,variables,"Observation");
  count.table$Proportion <- as.numeric(as.character(count.table$Proportion));
  count.table$Count <- as.integer(as.character(count.table$Count));
  for(i in 1:n.variables)
     count.table[[variables[i]]] <- factor(as.character(count.table[[variables[i]]]), levels=variables.levels[[i]], ordered=variables.ordered[[variables[[i]]]]);
  if(!is.null(variables.default))
    for(i in 1:length(variables.default))
       count.table[[names(variables.default)[i]]] <- relevel(count.table[[names(variables.default)[i]]], variables.default[[i]]);
  if(length(variables.numeric)>0)
    for(i in 1:length(variables.numeric))
       count.table[[variables.numeric[i]]] <- as.numeric(as.character(count.table[[variables.numeric[i]]]));
  count.table[[outcome]] <- factor(as.character(count.table[[outcome]]), levels=outcomes);

  count.table;

}

