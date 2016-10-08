#' @export exp2flux
#' @importFrom "sybil" "findExchReact"
#' @importFrom "gage" "kegg.gsets"
#' @author Kelly Botero <kjboteroo@unal.edu.co> - Maintainer: Daniel Camilo Osorio <dcosorioh@unal.edu.co>
#' @title Convert Gene Expression Data to FBA fluxes
#' @param model A valid model for the \code{'sybil'} package
#' @param expression A valid ExpressionSet object (one by treatment)
#' @param missing A character string specifying the value to be used in missing cases; must be one of \code{'min'}, \code{'1q'}, \code{'mean'}, \code{'median'}, \code{'3q'}, or \code{'max'}
#' @param scale A boolean value to specify if data must be scaled
exp2flux <- function(model,expression,organism=NULL,typeID=NULL,missing="mean",scale=FALSE){
  if(!is.null(organism) && !is.null(typeID)){
    data <- try(kegg.gsets(species = organism, id.type = typeID)) 
    data <- matrix(gsub("[[:digit:]]+$","",names(unlist(data))),dimnames = list(as.vector(unlist(data)),c()))
  }
  gpr.expression <- function(gpr,expression,missing){
    gpr <- gsub("[()]","",gpr)
    gpr <- gsub("[[:space:]]","",gpr)
    complex <- lapply(gpr, function(gpr){unlist(strsplit(gpr,"or"))})
    genes <- lapply(complex, function(complex){strsplit(complex,"and")})
    genes[lengths(genes) == 0] <- NA
    min.complex <- lapply(genes, function(gene){
      lapply(gene, function(gene){
        gene <- unlist(gene)
        if(!all(gene%in%rownames(data))){
          gene <- gene[gene%in%rownames(data)]
        }
        if (length(gene)==0){
          minComplex <- 0
        } else {
          if(any(gene%in%rownames(expression@assayData$exprs))){
            minComplex <- min(rowMeans(expression@assayData$exprs,na.rm = TRUE)[gene],na.rm = TRUE)
          } else {
            if(!is.null(organism) && !is.null(typeID)){
              minComplex <- summary(rowMeans(expression@assayData$exprs,na.rm = TRUE)[names(data[data[,1]%in%names(sort(table(data[gene,]))[1]),])])[[match(missing,c("min","1q","median","mean","3q","max"))]]
            } else {
              minComplex <- 0
            }
          }
        }
        return(minComplex)
      })
    })
    exp <- unlist(lapply(min.complex, function(min.complex){sum(unlist(min.complex),na.rm = TRUE)}))
    exp[exp==0] <- summary(rowMeans(expression@assayData$exprs,na.rm = TRUE))[[match(missing,c("min","1q","median","mean","3q","max"))]]
    return(exp)
  }
  exp <- gpr.expression(model@gpr,expression,missing=missing)
  if(scale==TRUE){
    exp <- round((exp/max(exp,na.rm = TRUE)),6)*1000
  }
  lb <- model@lowbnd
  ub <- model@uppbnd
  model@lowbnd <- -1*exp
  model@lowbnd[!model@react_rev] <- 0
  model@uppbnd <- exp
  model@lowbnd[model@react_id%in%findExchReact(model)@react_id] <- lb[model@react_id%in%findExchReact(model)@react_id]
  model@uppbnd[model@react_id%in%findExchReact(model)@react_id] <- ub[model@react_id%in%findExchReact(model)@react_id]
  return(model)
}
