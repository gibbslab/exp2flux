getFluxes <- function(model){
  reaction <- model@react_id
  fluxes <- getFluxDist(optimizeProb(model))
  stdDev <- suppressMessages(apply(matrix(fluxVar(model)@lp_obj,nrow = 2,byrow = TRUE), 2, sd))
  stdErr <- stdDev/sqrt(2)
  output <- as.data.frame(cbind(flux = fluxes, sd = stdDev, se = stdErr))
  rownames(output) <- reaction
  attr(output,"reactionID") <- reaction
  return(output)
}