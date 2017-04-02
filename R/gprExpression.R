gprExpression <- function(gpr,expression,organismPathways,missing){
  # Split GPR
  gpr <- splitGPR(GPR = gpr)
  
  # Setting the min by complex
  complexExpression <- complexMin(gprList = gpr, expressionData = expression, organismPathways = organismPathways, missing = missing)
  
  # Sum complex values
  gprExpression <- unlist(lapply(complexExpression, function(complexExpression){sum(unlist(complexExpression),na.rm = TRUE)}))
  
  # Replacing empty using the missing option
  gprExpression[gprExpression == 0] <- summary(object = rowMeans(x = expression@assayData$exprs, na.rm = TRUE))[[match(missing,c("min","1q","median","mean","3q","max"))]]
  
  # Return
  return(gprExpression)
}