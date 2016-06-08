map.gpr <- function(id, organism){
  id <- as.vector(id)
  gprs <- get.reference(organism)
  map <- function(id){
    return (as.vector(gprs[gprs[,"id"]%in%id,"gpr"]))
  }
  return(sapply(id, map))
}