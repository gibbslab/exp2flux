#' @importFrom "utils" "data" "download.file" "read.csv2"
 
# Download from KEGG database all pathway associations reported
downloadPathways <- function(organism, typeID){
  # Download pathways
  pathwaysGenes <- gage::kegg.gsets(species = organism, 
                                    id.type = typeID)
  # Removing non-required data
  pathwaysGenes <- unlist(pathwaysGenes$kg.sets)
  # Building a output matrix
  pathwaysGenes <- matrix(data = gsub(pattern = "[[:digit:]]+$",
                                      replacement = "",
                                      x = names(pathwaysGenes)
  ),
  dimnames = list(geneID = as.vector(x = pathwaysGenes), 
                  pathways = c("PATHWAY")))
  # Return
  return(pathwaysGenes)
}
