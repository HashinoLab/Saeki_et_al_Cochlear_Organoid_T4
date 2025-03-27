#Filters for (log-transformed) mitochondrial percentage genes
filterPlot <- function(obj) {
  obj[["mt.per"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  mt.per.vec <- obj$mt.per
  obj <- subset(obj, mt.per < 12.5)
  return(obj)
}


# Filtering for housekeeping gene expression
hk.filter <- 
  function(obj, hk.gene = NULL, subset = F, std.dev = 2) {
    # Sanitize input
    hk.gene <- as.character(hk.gene)
    # Creates a vector containing the reads found in each cell -> log normalize
    exp.vector <- obj@assays$RNA$counts[hk.gene, ] %>% log1p()
    # Add metadata column to the Seurat object to facilitate subsetting
    obj <- AddMetaData(obj, metadata = exp.vector, col.name = "hk")
    # Determines whether or not the object should be subset such that
    # cells greater or less than ~std.dev~ standard deviations are discarded
    if (subset){
      obj <- subset(obj, hk > mean(obj$hk) - std.dev*sd(obj$hk) &
                      hk < mean(obj$hk) + std.dev*sd(obj$hk))
    }
    return(obj)
  }


# Colors for grouped UMAP plots
cols.cond = c("#56A1FE", "#B656D4")
cols.day = c("#B656D4", "#C32500")


# Recurrently used ggplot2 themes

umap_theme <- function(){
  theme(plot.title = element_blank(), panel.background = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.line.x = element_blank(), axis.line.y = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        legend.position = "none",
        plot.margin=unit(c(0,0,0,0), "in"))
}


feat_theme_nolegend <- function(){
  theme(plot.title = element_text(size = 28, face = "italic"),
        panel.background = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.line.x = element_blank(), axis.line.y = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        legend.position = "none",
        plot.margin=unit(c(0,0,0,0), "in"))
}


feat_theme <- function(){
  theme(plot.title = element_text(size = 28, face = "italic"),
        panel.background = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.line.x = element_blank(), axis.line.y = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        plot.margin=unit(c(0,0,0,0), "cm"))
}


vln_theme <- function(){
  theme(
  plot.title = element_text(face = "bold.italic"),
  axis.line = element_line(color = "#888888"),
  axis.text = element_text(color = "#888888"),
  axis.ticks = element_line(color = "#888888"),
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.y = element_blank(),
  margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
  panel.background = element_blank(),
  legend.position = "none")
}
