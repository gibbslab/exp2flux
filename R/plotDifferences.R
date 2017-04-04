plotDifferences <- function(model1, model2) {
  differences <- getDifferences(model1 = model1, model2 = model2)
  valuesModel1 <-
    cbind(1, differences[["fluxModel1"]], attr(differences, "seModel1"))
  valuesModel2 <-
    cbind(2, differences[["fluxModel2"]], attr(differences, "seModel2"))
  diffValues <-
    as.data.frame.array(rbind(valuesModel1, valuesModel2))
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
      fill = factor(model)
    ))
  p + geom_bar(stat = "identity",
               position = position_dodge(1)) +
    geom_errorbar(limits, position = position_dodge(1),
                  width = 0.25) +
    labs(x = "Reaction", y = "Flux") +
    ggtitle("FDifference") +
    scale_fill_discrete(name = "Model") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
}