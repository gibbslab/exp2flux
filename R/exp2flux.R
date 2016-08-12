exp2flux <- function(model,expression,missing="max"){
  gpr.expression <- function(gpr,expression,missing){
    gpr <- gsub("[()]","",gpr)
    gpr <- gsub("[[:space:]]","",gpr)
    complex <- lapply(gpr, function(gpr){unlist(strsplit(gpr,"or"))})
    genes <- lapply(complex, function(complex){strsplit(complex,"and")})
    min.complex <- lapply(genes, function(gene){
      lapply(gene, function(gene){
        gExp <- rowMeans(expression@assayData$exprs,na.rm = TRUE)[unlist(gene)]
        if(identical(rep(TRUE,length(gExp)),as.vector(is.na(gExp)))){
          gExp <- 0
        }
        min(gExp,na.rm = TRUE)})
    })
    exp <- unlist(lapply(min.complex, function(min.complex){sum(unlist(min.complex),na.rm = TRUE)}))
    exp[exp==0]<-NA
    exp[is.na(exp)] <- as.vector(summary(exp)[match(missing,c("min","1q","median","mean","3q","max"))])
    exp <- round((exp/max(exp)),6)*1000
    return(exp)
  }
  exp <- gpr.expression(model@gpr,expression,missing=missing)
  model@lowbnd <- -1*exp
  model@lowbnd[!model@react_rev] <- 0
  model@uppbnd <- exp
  return(model)
}