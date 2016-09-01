#' @export exp2flux
#' @importFrom sybil findExchReact
#' @author Kelly Botero <kjboteroo@unal.edu.co> - 
#' Maintainer: Daniel Camilo Osorio Hurtado <dcosorioh@unal.edu.co>
#' @title Convert Gene Expression Data to FBA fluxes
#' @param model A valid model for the \code{'sybil'} package
#' @param expression A valid ExpressionSet object (one by treatment)
#' @param missing A character string specifying the value to be used in missing cases; must be one of \code{'gmean'}, \code{'min'}, \code{'1q'}, \code{'mean'}, \code{'median'}, \code{'3q'}, or \code{'max'}
#' @param scale A boolean value to specify if data must be scaled
#' @examples 
#' \dontrun{
#' # Required Libraries
#' library("sybil")
#' library("Biobase")
#' 
#' # Reference Data
#' data("Ec_core")
#' optimizeProb(Ec_core)
#' 
#' # Random Data
#' gE <- ExpressionSet(matrix(runif(50,1,100),dimnames = list(sample(Ec_core@allGenes,50),c())))
#' 
#' # Applying the function
#' modifiedModel <- exp2flux(model = Ec_core,expression = gE,missing = "gmean",scale = TRUE)
#' 
#' # Result
#' optimizeProb(modifiedModel)
#' }
exp2flux <- function(model,expression,missing="gmean",scale=FALSE){
  # This function was written by Paul McMurdie @ http://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  gpr.expression <- function(gpr,expression,missing){
    gpr <- gsub("[()]","",gpr)
    gpr <- gsub("[[:space:]]","",gpr)
    complex <- lapply(gpr, function(gpr){unlist(strsplit(gpr,"or"))})
    genes <- lapply(complex, function(complex){strsplit(complex,"and")})
    min.complex <- lapply(genes, function(gene){
      lapply(gene, function(gene){
        gExp <- rowMeans(expression@assayData$exprs,na.rm = TRUE)[unlist(gene)]
        if(all(is.na(gExp))){
          gExp <- 0
        }
        return(min(gExp,na.rm = TRUE))})
    })
    exp <- unlist(lapply(min.complex, function(min.complex){sum(unlist(min.complex),na.rm = TRUE)}))
    exp[exp==0] <- NA
    if(missing == "gmean"){
     exp[is.na(exp)] <- gm_mean(exp,na.rm = TRUE)
    } else {
      exp[is.na(exp)] <- summary(exp)[match(missing,c("min","1q","median","mean","3q","max"))]
    }
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



