#' @export plotDifferences
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
#' @title Plot a bipartite graph.
#' @description Plot a bipartite graph for each model to represent FBA flux changes between them.
plotDifferences <- function(model1,model2,...){
  S <- as.matrix(t(model1@S))
  colnames(S) <- model1@met_id
  rownames(S) <- model1@react_id
  S[S < 0] <- -1
  S[S > 0] <- 1
  S[model1@react_rev,][S[model1@react_rev,]!=0] <- 2
  fD <- as.data.frame.array(fluxDifferences(model1 = model1,model2 = model2,foldReport = 0)[getFluxDist(optimizeProb(model1))!=0 & getFluxDist(optimizeProb(model2))!=0,])
  S <- S[getFluxDist(optimizeProb(model1))!=0 & getFluxDist(optimizeProb(model2))!=0,]
  S <- S[,colSums(S!=0)!=0]
  types <- structure(c(rep(1,length(rownames(S))),rep(0,length(colnames(S)))),names=c(rownames(S),colnames(S)))
  edges <- NULL
  for(i in names(types)[types==1]){
    n_mets <- names(S[i,])
    if (length(grep(2,S[i,]))>0){
      edges <- c(edges,unlist(sapply(n_mets[S[i,]!=0],function(met){c(met,i)})))
      edges <- c(edges,unlist(sapply(n_mets[S[i,]!=0],function(met){c(i,met)})))
    } else{
      edges <- c(edges,unlist(sapply(n_mets[S[i,]<0],function(met){c(met,i)})))
      edges <- c(edges,unlist(sapply(n_mets[S[i,]>0],function(met){c(i,met)})))
    }
  }
  edges <- as.vector(unlist(as.vector(sapply(edges, function(edge){which(names(types)==edge)}))))
  g <- make_bipartite_graph(types = types,edges = edges,directed = TRUE)
  V(g)$name <- names(types)
  V(g)$color <- ifelse(types==1, "gray80", "gray90")
  V(g)$color[names(types)%in%rownames(fD)[fD$foldChange<0]] <- "red"
  V(g)$color[names(types)%in%rownames(fD)[fD$foldChange>0]] <- "green"
  V(g)$label.cex <- 0.7
  vSizes <- types
  vSizes[types==1] <- (abs(fD$foldChange)/max(abs(fD$foldChange)))*20
  vSizes[types==0] <- 10
  g <- simplify(g,count_multiple(g))
  topology <- layout_with_gem(g)
  plot.igraph(g,vertex.size=vSizes,vertex.label=names(types),edge.arrow.size=0.3)
  }
