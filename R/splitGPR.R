# Split GPRs by complex
splitGPR <- function(GPR){
  # Remove brackets
  GPR <- gsub("[()]","",x = GPR)
  # Remove spaces
  GPR <- gsub("[[:space:]]+","", x = GPR)
  # Split complex
  complexGPR <- lapply(GPR, function(GPR){unlist(strsplit(GPR,"or"))})
  # Extract genes
  genesGPR <- lapply(complexGPR, function(complexGPR){strsplit(complexGPR,"and")})
  # Replacing empty by NA
  genesGPR[lengths(genesGPR) == 0] <- NA
  # Return
  return(genesGPR)
}
