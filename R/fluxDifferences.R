#' @export fluxDifferences
#' @importFrom "sybil" "getFluxDist" "optimizeProb"
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
#' @title Report the fold change of fluxes between two models
#' @description This functions calculates the fold change \code{"(fluxModel2/fluxModel1)-1"} for fluxes of two given metabolic models.
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
#' # Evaluating Differences
#' fluxDifferences(model1 = Ec_core, 
#'                 model2 = Ec_coreGE, 
#'                 foldReport = 0.5)
#'                 }

fluxDifferences <- function(model1,model2,foldReport=2){
  f_m1 <- getFluxDist(optimizeProb(model1))
  f_m2 <- getFluxDist(optimizeProb(model2))
  if(identical(length(f_m1),length(f_m2))){
    fold <- ((f_m2-f_m1)/f_m1)
    fold[f_m1==0] <- f_m2[f_m1==0]
    fold[is.na(fold)] <- 0
    fold[is.infinite(fold)] <- 0
    different <- (abs(fold)>=foldReport)
    differentFlux <- matrix(cbind(f_m1,f_m2,fold),nrow = model1@react_num,dimnames = list(model1@react_id,c("fluxModel1","fluxModel2","foldChange")))
    differentFlux <- differentFlux[different,]
    return(differentFlux)
  } 
  else {
    warning("Metabolic models must be the same with different restrictions")
  }
}