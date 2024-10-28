gseaplot3<-function (x, geneSetID, title = "", color = "darkslateblue", base_size = 11, 
          rel_heights = c(1.5, 0.5, 1), subplots = 1:3, pvalue_table = FALSE, 
          ES_geom = "line") 
{
  ES_geom <- match.arg(ES_geom, c("line", "dot"))
  geneList <- position <- NULL
  if (length(geneSetID) == 1) {
    gsdata <- gsInfo(x, geneSetID)
  }
  else {
    gsdata <- do.call(rbind, lapply(geneSetID, gsInfo, object = x))
  }
  p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + theme_classic(base_size) + 
    theme(panel.grid.major = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(colour = "grey92"), 
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) + 
    scale_x_continuous(expand = c(0, 0))
  if (ES_geom == "line") {
    es_layer <- geom_line(aes_(y = ~runningScore, color = ~Description), 
                          size = 1)
  }
  else {
    es_layer <- geom_point(aes_(y = ~runningScore, color = ~Description), 
                           size = 1, data = subset(gsdata, position == 1))
  }
  p.res <- p + es_layer + theme(legend.position = c(0.8, 1), 
                                legend.title = element_blank(), legend.background = element_rect(fill = "transparent"))
  p.res <- p.res + ylab("Enrichment Score") + theme(axis.text.x = element_blank(), axis.title.y = element_text(size=16,face = "bold"),
                                                            axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
                                                            plot.margin = ggplot2::margin(t = 0.2, r = 0.2, b = 0, l = 0.2, 
                                                                                 unit = "cm"))
  i <- 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == 
                   term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }
  p2 <- ggplot(gsdata, aes_(x = ~x)) + geom_linerange(aes_(ymin = ~ymin, 
                                                           ymax = ~ymax, color = ~Description)) + xlab(NULL) + ylab(NULL) + 
    theme_classic(base_size) + theme(legend.position = "none", 
                                     plot.margin = ggplot2::margin(t = -0.1, b = 0, unit = "cm"), axis.ticks = element_blank(), 
                                     axis.text = element_blank(), axis.line.x = element_blank()) + 
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 
                                                                         0))
  if (length(geneSetID) == 1) {
    v <- seq(1, sum(gsdata$position), length.out = 9)
    inv <- findInterval(rev(cumsum(gsdata$position)), v)
    if (min(inv) == 0) 
      inv <- inv + 1
    col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))
    ymin <- min(p2$data$ymin)
    yy <- max(p2$data$ymax - p2$data$ymin) * 0.3
    xmin <- which(!duplicated(inv))
    xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
    d <- data.frame(ymin = ymin, ymax = yy, xmin = xmin, 
                    xmax = xmax, col = col[unique(inv)])
    p2 <- p2 + geom_rect(aes_(xmin = ~xmin, xmax = ~xmax, 
                              ymin = ~ymin, ymax = ~ymax, fill = ~I(col)), data = d, 
                         alpha = 0.9, inherit.aes = FALSE)
  }
  df2 <- p$data
  df2$y <- p$data$geneList[df2$x]
  p.pos <- p + geom_segment(data = df2, aes_(x = ~x, xend = ~x, 
                                             y = ~y, yend = 0), color = "darkslateblue",size=1.5)
  p.pos <- p.pos + ylab("Sum of samples") +  
    theme(plot.margin = ggplot2::margin(t = -0.1, r = 0.2, b = -0.6, 
                               l = 0.2, unit = "cm"),axis.text.x = element_blank())
  if (!is.null(title) && !is.na(title) && title != "") 
    p.res <- p.res + ggtitle(title)
  if (length(color) == length(geneSetID)) {
    p.res <- p.res + scale_color_manual(values = color)
    if (length(color) == 1) {
      p.res <- p.res + theme(legend.position = "none")
      p2 <- p2 + scale_color_manual(values = "black")
    }
    else {
      p2 <- p2 + scale_color_manual(values = color)
    }
  }
  if (pvalue_table) {
    pd <- x[geneSetID, c("Description", "pvalue", "NES")]
    rownames(pd) <- pd$Description
    pd <- pd[, -1]
    for (i in seq_len(ncol(pd))) {
      pd[, i] <- format(pd[, i], digits = 4)
    }
    tp <- tableGrob2(pd, p.res)
    p.res <- p.res + theme(legend.position = "none") + annotation_custom(tp, 
                                                                         xmin = quantile(p.res$data$x, 0.5), xmax = quantile(p.res$data$x, 
                                                                                                                             0.99), ymin = quantile(p.res$data$runningScore, 
                                                                                                                                                    0.85), ymax = quantile(p.res$data$runningScore, 
                                                                                                                                                                           0.99))
  }
  plotlist <- list(p.res, p2, p.pos)[subplots]
  n <- length(plotlist)
  plotlist[[n]] <- plotlist[[n]] + theme(axis.line.x = element_line(), 
                                         axis.ticks.x = element_line())
  plotlist[[3]]<-plotlist[[3]]+theme(axis.text.x = element_blank())
  if (length(subplots) == 1) 
    return(plotlist[[1]] + theme(plot.margin = ggplot2::margin(t = 0.2, 
                                                      r = 0.2, b = 0.2, l = 0.2, unit = "cm")))
  if (length(rel_heights) > length(subplots)) 
    rel_heights <- rel_heights[subplots]
  aplot::gglist(gglist = plotlist, ncol = 1, heights = rel_heights)
}


tableGrob2 <- function(d, p = NULL) {
  # has_package("gridExtra")
  d <- d[order(rownames(d)),]
  tp <- gridExtra::tableGrob(d)
  if (is.null(p)) {
    return(tp)
  }
  
  # Fix bug: The 'group' order of lines and dots/path is different
  p_data <- ggplot_build(p)$data[[1]]
  # pcol <- unique(ggplot_build(p)$data[[1]][["colour"]])
  p_data <- p_data[order(p_data[["group"]]), ]
  pcol <- unique(p_data[["colour"]])
  ## This is fine too
  ## pcol <- unique(p_data[["colour"]])[unique(p_data[["group"]])]  
  j <- which(tp$layout$name == "rowhead-fg")
  
  for (i in seq_along(pcol)) {
    tp$grobs[j][[i+1]][["gp"]] <- gpar(col = pcol[i])
  }
  return(tp)
}
