#' @export exp2flux
#' @importFrom "sybil" "findExchReact"
#' @importFrom "gage" "kegg.gsets"
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co> and Kelly Botero <kjboteroo@unal.edu.co>
#' @title Convert gene expression data to FBA flux boundaries.
#' @description Set the flux boundaries for each reaction based on the expression values of genes associated in their GPR rules.
#' @details This function computes the flux boundaries for each reaction based in their associated GPR. 
#' The value is obtained as follows: When two or more genes are associated by AND operators in the GPR rule, the minima value of their associated expression values is selected. 
#' The AND operator indicates that all asociated genes are required to the enzymatic complex formation. In that case, the gene with the smaller expression value acts as the enzyme complex formation limitant. 
#' In turn, when two or more genes are associated by OR operators in the GPR rule, the sum of all their associated expression values is computed. The OR operator indicates that each gene can code an entire isozyme to act as reaction catalyst. 
#' To missing gene expression values, the function assigns one of: \code{'min'}, \code{'1q'}, \code{'mean'}, \code{'median'}, \code{'3q'}, or \code{'max'} expression value computed from the genes associated to the same metabolic pathway. 
#' In case of not possible pathway assignment to a gene, the value is calculated from all gene expression values. 
#' The flux boundaries of exchange reactions are not modified.
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
#' }

exp2flux <- function(model,expression,organism=NULL,typeID=NULL,missing="mean",scale=FALSE){
  
  # Download valid organisms from the KEGG database
  keggDownload <- tempdir()
  download.file("rest.kegg.jp/list/organism",paste0(keggDownload,"organism.txt"),quiet = TRUE)
  keggOrganism <- as.vector(read.csv2(paste0(keggDownload,"organism.txt"),header = FALSE,sep ="\t")[,2])
  
  # Validate given organism
  if(isTRUE(organism %in% keggOrganism) & is.null(organism)){
    stop("Given organism is not included in the KEGG database")
  }
  
  # Validate typeID
  if(isTRUE(typeID %in% c("kegg","entrez")) & is.null(typeID)){
    stop("Given typeID is not allowed")
  }
  
  # Download pathways for given organism
  if(!is.null(organism) & !is.null(typeID)){
    organismPathways <- downloadPathways(organism = organism, typeID = typeID)
  } else {
    organismPathways <- NULL
  }

  # Expression Values for each GPR
  expressionValues <- gprExpression(gpr = model@gpr, 
                                    expression = expression,
                                    organismPathways = organismPathways, 
                                    missing = missing
                                    )
  
  # Original boundaries
  lb <- model@lowbnd
  ub <- model@uppbnd
  
  # Scale == TRUE
  expressionValues <- round((expressionValues / max(expressionValues))*1000)
  
  # New boundaries
  model@lowbnd <- -1 * expressionValues
  model@uppbnd <- expressionValues
  
  # Setting 0 as lower boundary for irreversible reactions
  model@lowbnd[!model@react_rev] <- 0
  
  # Exchange reactions
  model@lowbnd[model@react_id%in%findExchReact(model)@react_id] <- lb[model@react_id%in%findExchReact(model)@react_id]
  model@uppbnd[model@react_id%in%findExchReact(model)@react_id] <- ub[model@react_id%in%findExchReact(model)@react_id]
  
  # Return
  return(model)
}