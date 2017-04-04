#' @export getDifferences
getDifferences <- function(model1, model2, pLimit = 0.05) {
  # Getting fluxes by model
  fluxModel1 <- getFluxes(model1)
  idModel1 <- attr(fluxModel1, "reactionID")
  fluxModel2 <- getFluxes(model2)
  idModel2 <- attr(fluxModel2, "reactionID")
  # Identifying shared reactions
  fluxModel1 <- fluxModel1[idModel1 %in% idModel2, ]
  fluxModel2 <- fluxModel2[idModel2 %in% idModel1, ]
  attr(fluxModel1, "reactionID") <- idModel1[idModel1 %in% idModel2]
  attr(fluxModel2, "reactionID") <- idModel2[idModel2 %in% idModel1]
  # Set range
  iM1 <-
    matrix(
      data = c(
        fluxModel1[["flux"]],
        fluxModel1[["flux"]] + fluxModel1[["sd"]],
        fluxModel1[["flux"]] - fluxModel1[["sd"]]
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c(), c(attr(
        fluxModel1, "reactionID"
      )))
    )
  iM2 <-
    matrix(
      data = c(
        fluxModel2[["flux"]],
        fluxModel2[["flux"]] + fluxModel2[["sd"]],
        fluxModel2[["flux"]] - fluxModel2[["sd"]]
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c(), c(attr(
        fluxModel2, "reactionID"
      )))
    )
  # Identify difference
  pValue <-
    unlist(sapply(attr(fluxModel1, "reactionID"), function(reaction) {
      if (!isTRUE(all.equal(iM1[, reaction], iM2[, reaction]))) {
        t.test(iM1[, reaction], iM2[, reaction])$p.value
      } else {
        return(1)
      }
    }))
  differences <-
    (pValue < pLimit) &
    (abs(fluxModel1[["flux"]]) > fluxModel1[["se"]]) &
    (abs(fluxModel2[["flux"]]) > fluxModel2[["se"]])
  # Output Values
  pValue <- pValue[differences]
  log2foldChange <-
    log2(abs(fluxModel2[["flux"]][differences] / fluxModel1[["flux"]][differences]))
  valuesModel1 <- fluxModel1[["flux"]][differences]
  valuesModel2 <- fluxModel2[["flux"]][differences]
  errModel1 <- fluxModel1[["se"]][differences]
  errModel2 <- fluxModel2[["se"]][differences]
  output <-
    data.frame(
      cbind(
        fluxModel1 = valuesModel1,
        fluxModel2 = valuesModel2,
        fold2Change = log2foldChange,
        pValue = pValue
      )
    )
  attr(output, "seModel1") <- errModel1
  attr(output, "seModel2") <- errModel2
  # Return
  if (dim(output)[1] > 0) {
    return(output)
  } else {
    return(NA)
  }
}