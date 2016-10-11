#' @export exp2flux
#' @importFrom "sybil" "findExchReact"
#' @importFrom "gage" "kegg.gsets"
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co> and Kelly Botero <kjboteroo@unal.edu.co>
#' @title Convert gene expression data to FBA fluxes
#' @description This function calculates the flux boundaries for each reaction based in their associated GPR. 
#' The value es obtained as follows: When two genes are associated by an \code{AND} operation according to the GPR rule, a \code{min} function is applied to their associated expression values. 
#' In the \code{AND} case, downregulated genes alter the reaction acting as enzyme formation limitant due two are required to complex formation. 
#' In turn, when the genes are associated by an \code{OR} rule, each one of then can code an entire enzyme to act as reaction catalyst. 
#' In this case, a \code{sum} function is applied for their associated expression values.To missing gene expression values, the function assigns one of: \code{'min'}, \code{'1q'}, \code{'mean'}, \code{'median'}, \code{'3q'}, or \code{'max'} expression value calculated from the genes associated to the same metabolic pathway.
#' In case of not possible pathway assignment to a gene, the value is calculated from all gene expression values. 
#' The fluxes boundaries of exchange reactions are not modified.
#' @param model A valid model for the \code{'sybil'} package.
#' @param expression A valid ExpressionSet object (one by treatment).
#' @param organism  A valid organism identifier for the KEGG database. List of valid organism identifiers are available in: http://rest.kegg.jp/list/organism.
#' @param typeID A string to define the type of ID used in GPR's. One of \code{"entrez"} or \code{"kegg"} must be given.
#' @param missing A character string specifying the value to be used in missing cases; must be one of \code{'min'}, \code{'1q'}, \code{'mean'}, \code{'median'}, \code{'3q'}, or \code{'max'}
#' @param scale A boolean value to specify if data must be scaled to assign a value of 1000 as max.
#' @examples 
#' \dontrun{
#' # Loading a model
#' library("sybil")
#' library("Biobase")
#' 
#' # Original model:
#' data("Ec_core")
#' 
#' # Original model evaluation:
#' optimizeProb(Ec_core)
#' 
#' # Generating simulated expressionSets
#' expressionData <- matrix(data = runif(3*length(Ec_core@allGenes),min = 1,max = 100),
#'                          nrow = length(Ec_core@allGenes),
#'                          dimnames = list(c(Ec_core@allGenes),c()))
#' expressionData <- ExpressionSet(assayData = expressionData)
#' 
#' # Applying exp2flux
#' Ec_coreGE <- exp2flux(model = Ec_core,
#'                       expression = expressionData,
#'                       missing = "mean")
#' # Evaluating exp2flux model
#' optimizeProb(Ec_coreGE)
#'                 }
#' @keywords integrate expression data metabolic network genome scale reconstruction tissue-specific
exp2flux <- function(model,expression,organism=NULL,typeID=NULL,missing="mean",scale=FALSE){
  if(!is.null(organism) && !is.null(typeID)){
    data <- try(kegg.gsets(species = organism, id.type = typeID)) 
    data <- matrix(gsub("[[:digit:]]+$","",names(unlist(data$kg.sets))),dimnames = list(as.vector(unlist(data$kg.sets)),c()))
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
        if(!is.null(organism) && !is.null(typeID)){
        if(!all(gene%in%rownames(data))){
          gene <- gene[gene%in%rownames(data)]
        }} else {
          gene <- gene[gene%in%rownames(expression@assayData$exprs)]
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
  exp <- gpr.expression(gpr = model@gpr,
                        expression = expression,
                        missing=missing)
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
