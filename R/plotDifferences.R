#' @export plotDifferences
#' 
plotDifferences <- function(model1, model2) {
  if (!all(c(class(model1), class(model2)) == "modelorg")) {
    stop("class error")
  }
  mNames <- c(model1@mod_id, model2@mod_id)
  if (isTRUE(all.equal(mNames[1], mNames[2]))) {
    "names should be different"
  }
  differences <- getDifferences(model1 = model1, model2 = model2)
  valuesModel1 <- cbind(1,
                        differences[["fluxModel1"]],
                        attr(differences, "seModel1"))
  valuesModel2 <- cbind(2,
                        differences[["fluxModel2"]],
                        attr(differences, "seModel2"))
  diffValues <- as.data.frame.array(rbind(valuesModel1,
                                          valuesModel2))
  attr(diffValues, "reactionID") <- rep(rownames(differences), 2)
  colnames(diffValues) <- c("model", "fluxValue", "se")
  limits <- aes(
    ymax = diffValues$fluxValue + diffValues$se,
    ymin = diffValues$fluxValue - diffValues$se
  )
  p <-
    ggplot(data = diffValues, aes(
      x = factor(attr(diffValues, "reactionID")),
      y = fluxValue,
      fill = c(mNames)[factor(model)]
    )) + geom_bar(stat = "identity",
                  position = position_dodge(1)) +
    geom_errorbar(limits, position = position_dodge(1),
                  width = 0.25) +
    labs(x = "Reaction", y = "Flux Value") +
    ggtitle("Flux Differences") +
    scale_fill_discrete(name = "Model:") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  return(p)
}