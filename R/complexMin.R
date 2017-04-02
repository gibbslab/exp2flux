complexMin <- function(gprList, expressionData, organismPathways = NULL, missing = "Mean"){
  expression <- rowMeans(expressionData@assayData$exprs)
  lapply(gprList, function(GPR){
    lapply(GPR, function(complex){
      if(all(is.na(complex))){
        return(NA)
      } else {
        if(any(complex %in% names(expression))){
          minValue <- min(expression[names(expression) %in% complex], na.rm = TRUE)
        } else {
          if(!is.null(organismPathways) & any(complex %in% rownames(organismPathways))){
            mPathway <- names(sort(table(organismPathways[rownames(organismPathways) %in% complex,"PATHWAY"]),decreasing = TRUE)[1])
            pathwaysGenes <- rownames(organismPathways)[organismPathways[,"PATHWAY"] %in% mPathway]
            summaryValues <- as.vector(summary(expression[names(expression) %in% pathwaysGenes])[match(missing,c("min","1q","median","mean","3q","max"))])
            return(summaryValues)
          } else {
            return(NA)
          }
        }
      }
    })
  })
  
}
