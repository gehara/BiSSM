plotPCS2 <- function (models, index, observed, subsample) 
  {
  labels <- unique(index)
  labels <- sort(c(labels, "observed"))
  sizes <- rep(2, length(labels))
  sizes[which(labels == "observed")] <- 10
  shapes <- rep(16, length(labels))
  shapes[which(labels == "observed")] <- 8
  data.PCA <- index[complete.cases(models)]
  models.PCA <- models[complete.cases(models), ]
  data.obs <- rep("observed", nrow(observed))
  x <- sample(1:length(data.PCA), length(data.PCA) * subsample)
  data.PCA <- data.PCA[x]
  models.PCA <- models.PCA[x, ]
  PCA <- prcomp(rbind(models.PCA, observed), center = T, scale. = T, 
                retx = T)
  scores <- data.frame(PCA$x[, 1:ncol(PCA$x)])
  PC <- colnames(scores)[1:10]
  plotPCA <- function(PCS) {
    PCS <- rlang::sym(PCS)
    p <- ggplot2::ggplot(scores, ggplot2::aes(x = PC1, y = !!PCS)) + 
      ggplot2::theme(legend.position = "none") + ggplot2::geom_point(ggplot2::aes(colour = c(data.PCA, 
                                                                                             data.obs), size = c(data.PCA, data.obs), shape = c(data.PCA, 
                                                                                                                                                data.obs))) + ggplot2::scale_shape_manual(values = shapes) + 
      ggplot2::scale_size_manual(values = sizes) + ggplot2::scale_color_brewer(palette = "Set1") + 
      if (PCS == "PC2") 
        ggplot2::theme(legend.position = "top", legend.direction = "horizontal", 
                       legend.title = ggplot2::element_blank())
    return(p)
  }
  P <- NULL
  for (i in 2:10) {
    P[[i]] <- plotPCA(PC[i])
  }
  gridExtra::grid.arrange(P[[2]], P[[3]], P[[4]], P[[5]], P[[6]], 
                          P[[7]], P[[8]], P[[9]], P[[10]], nrow = 3)
}
