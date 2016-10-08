#' @export fluxDifferences
#' @importFrom "sybil" "getFluxDist" "optimizeProb"
#' @author Daniel Camilo Osorio <dcosorioh@unal.edu.co>
#' @title Report the fold change of fluxes between two models
#' 
#' @examples
#' \dontrun{
#' # Loading a model
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
#'                       missing = "median")
#' 
#' # Evaluating Differences
#' fluxDifferences(model1 = Ec_core, 
#'                 model2 = Ec_coreGE, 
#'                 foldReport = 0)
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
    different <-as.data.frame.array(cbind(model1@react_id[different],f_m1[different],f_m2[different],round(fold[different],2)))
    colnames(different) <- c("Reaction","fluxModel1","fluxModel2","foldChange")
    return(as.data.frame(different))
  } 
  else {
    warning("Models are not the same")
  }
}