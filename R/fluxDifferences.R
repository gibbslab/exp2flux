#' @export fluxDifferences
fluxDifferences <- function(model1,model2,fold=2){
  f_m1 <- getFluxDist(optimizeProb(model1))
  f_m2 <- getFluxDist(optimizeProb(model2))
  if(identical(length(f_m1),length(f_m2))){
    fold <- (f_m2/f_m1)-1
    fold[is.na(fold)] <- 0
    fold[is.infinite(fold)] <- 0
    different <- abs(fold)>2
    different <-as.data.frame.array(cbind(model1@react_id[different],f_m1[different],f_m2[different],round(fold[different],2)))
    colnames(different) <- c("Reaction","FluxModel1","FluxModel2","Fold Change")
    return(as.data.frame(different))
  } 
  else {
    warning("Models are not the same")
  }
}