#' @export plotDifferences
#' @importFrom "igraph" "make_bipartite_graph" "V" "V<-" "layout.davidson.harel" "plot.igraph"
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
#' @title Plot the fold change of fluxes between two models into a bipartite graph.
#' @description This functions calculates the fold change \code{"(fluxModel2/fluxModel1)-1"} for the fluxes of two given metabolic models and plot it into a bipartite graph. Vertex size is assigned proportional to the fold change; if fold change is positive, green color is assigned, in contrary case, red color is assigned.
#' @param model1 A valid model for the \code{'sybil'} package.
#' @param model2 A valid model for the \code{'sybil'} package. Must have the same reactions (reaction number and reaction identifiers) as \code{"model1"} with different restrictions.
#' @param ... Additional arguments affecting the plot
#' @examples
#' \dontrun{
#' # Loading a model
#' library("sybil")
#' library("Biobase")
#' data("Ec_core")
#' 
#' # Generating expressionSets
#' expressionData <- matrix(data = runif(3*length(Ec_core@allGenes),min = 1,max = 100),
#'                          nrow = length(Ec_core@allGenes),
#'                          dimnames = list(c(Ec_core@allGenes),c()))
#' expressionData <- ExpressionSet(assayData = expressionData)
#' 
#' # Applying exp2flux
#' Ec_coreGE <- exp2flux(model = Ec_core,
#'                       expression = expressionData,
#'                       missing = "mean")
#' 
#' # Plotting Differences
#' plotDifferences(model1 = Ec_core, 
#'                 model2 = Ec_coreGE)
#'                 }
#' @keywords plot flux differences between two model scenarios
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
  layout <- layout.davidson.harel(g)
  plot.igraph(g,vertex.size=vSizes,vertex.label=names(types),edge.arrow.size=0.3)
  }
