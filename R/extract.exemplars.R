extract.exemplars <-
function(model.polytomous, model.hclust=NULL, n.clusters=10, p.bins=0, features=FALSE)
{
  if(!"polytomous" %in% class(model.polytomous))
    stop("Argument 'model.polytomous' not of class 'polytomous'.")

  response = gsub("[ ]+"," ",paste(deparse(model.polytomous$formula[[2]],width.cutoff=500),collapse=""))
  predictors = levels(as.factor(unlist(lapply(attr(terms.formula(model.polytomous$formula),"term.labels"),function(x) strsplit(x,"[ ]*([\\+]|[\\|]|[:])[ ]*")))))
  probs <- model.polytomous$fitted
  names.responses <- colnames(probs)

  if(is.null(model.hclust))
    model.hclust <- hclust(dist(model.polytomous$data[,predictors], method="binary"), method="ward")
  model.clusters <- cutree(model.hclust, k=n.clusters)

  n.features <- apply(model.polytomous$data[,predictors],1,function(x) length(which(x)))

  indices = NULL
  outcomes = NULL
  max.probs = NULL
  properties = NULL
  for(i in 1:n.clusters)
     { indices.cluster <- which(model.clusters==i)
       response.original <- as.character(model.polytomous$data[indices.cluster, response])
       response.predicted <- as.character(predict(model.polytomous, newdata=model.polytomous$data[indices.cluster,], type="choice"))
       indices.cluster.original <- indices.cluster[response.original==response.predicted]
       nrow.cluster.original = length(indices.cluster.original)
       if(nrow.cluster.original!=0)
       { if(p.bins>0)
            probs.cluster = round((probs[indices.cluster.original,,drop=FALSE] + 1/(p.bins * 2)) * p.bins, 0)
         else
            probs.cluster = probs[indices.cluster.original,,drop=FALSE]
         if(features)
            index.max = order(-as.vector(apply(probs.cluster,1,max)), -n.features[indices.cluster.original])[1]
         else
            index.max = order(-as.vector(apply(probs.cluster,1,max)))[1]
         response.exemplar.cluster = as.character(model.polytomous$data[indices.cluster.original,response])[index.max]
         index.exemplar.cluster = indices.cluster.original[index.max]
         # response.exemplar.cluster = names.responses[((index.max-1) %/% nrow.cluster.original + 1)]
         # index.exemplar.cluster = indices.cluster.original[(index.max-1) %% nrow.cluster.original + 1]

         indices <- c(indices, index.exemplar.cluster)
         outcomes <- c(outcomes, response.exemplar.cluster)
         max.probs <- c(max.probs, probs[index.exemplar.cluster, response.exemplar.cluster])
         properties <- c(properties, apply(model.polytomous$data[index.exemplar.cluster, predictors],1,function(x) paste(names(which(x)),collapse="; ")))
       }
     }
  exemplars <- data.frame(indices, outcomes, max.probs, properties)

  return(exemplars)
  
}
