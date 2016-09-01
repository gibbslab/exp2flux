#' @export fluxDifferences
#' @importFrom sybil getFluxDist optimizeProb
#' @author Daniel Camilo Osorio Hurtado <dcosorioh@unal.edu.co>
#' @title Report the fold change of fluxes between two models
fluxDifferences <- function(model1,model2,fold2report=2){
  f_m1 <- getFluxDist(optimizeProb(model1))
  f_m2 <- getFluxDist(optimizeProb(model2))
  if(identical(length(f_m1),length(f_m2))){
    fold <- (f_m2/f_m1)-1
    fold[is.na(fold)] <- 0
    fold[is.infinite(fold)] <- 0
    different <- (abs(fold)>fold2report)
    different <-as.data.frame.array(cbind(model1@react_id[different],f_m1[different],f_m2[different],round(fold[different],2)))
    colnames(different) <- c("Reaction","fluxModel1","fluxModel2","foldChange")
    return(as.data.frame(different))
  } 
  else {
    warning("Models are not the same")
  }
}