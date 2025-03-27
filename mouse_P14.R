## Mouse P14 transcriptome data from Xu, 2022, Front Cell Neurosci
## www.frontiersin.org/journals/cellular-neuroscience/articles/10.3389/fncel.2022.962106/

setwd("C:/ ... /P14/")
library(Seurat)
library(dplyr)
library(ggplot2)

# Loading functions stored in the external file
source("./Script/functions.R")

# create v5 assays
cochlea_P14 <- Read10X(data.dir = "./")
cochlea_P14 <- CreateSeuratObject(cochlea_P14)
class(cochlea_P14[["RNA"]])

cochlea_P14[["mt.per"]] <- PercentageFeatureSet(cochlea_P14, pattern = "^mt-")

v_cochlea_P14 <- VlnPlot(cochlea_P14, features = c("nFeature_RNA", "nCount_RNA", "mt.per"), ncol = 3)
plot_grid(v_cochlea_P14)

#### ------------------------------ Filtering ------------------------------
# Filter each dataset according to mitochondrial percentage
counts.list.post <- counts.list.pre <- numeric()

  counts.list.pre <- length(Cells(cochlea_P14))
  cochlea_P14 <- filterPlot(cochlea_P14)
  counts.list.post <- length(Cells(cochlea_P14))

diff <- counts.list.pre - counts.list.post
diff

#same thing for Rpl27
counts.list.post.2 <- numeric()

  cochlea_P14 <- housekeeping.filter(cochlea_P14, housekeeping.gene = "Rpl27", subset = T, std.dev = 2)
  counts.list.post.2 <- length(Cells(cochlea_P14))

diff.2 <- counts.list.post - counts.list.post.2
diff.2

#### ------------------------------ Seurat ------------------------------
cochlea_P14 <- SCTransform(cochlea_P14) 
cochlea_P14 <- RunPCA(cochlea_P14)
cochlea_P14 <- FindNeighbors(cochlea_P14, dims = 1:20) %>%
  FindClusters(resolution = 0.3) %>% 
  RunUMAP(dims = 1:20)

saveRDS(cochlea_P14, "./cochlea_P14.rds")

## UMAP and dot plots
p <- DimPlot(cochlear_p14, pt.size = 0.4, label = T) + umap_theme()
ggsave("dimplot_cochlear_p14.tiff", plot = p, device = "tiff",
       path = "./Output/", width = 4, height = 4, units = "in", dpi = 600)

DotPlot(cochlear_p14, col.max = 2, dot.min = 0, dot.scale = 8,
        features = 
          c("Tmie","Pcp4","Gfi1","Myo7a","Slc26a5","Ikzf2",
            "Ocm","Kcnq4","Tmc1","Lmod3","Tbx2","Insm1",
            "Slc17a8","Fgf8","Otof","Calb1","Pvalb")) +
  theme(
    axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "italic"),
    axis.text.y = element_text(),
    axis.title = element_blank()
  ) 
ggsave("dot_cochlea_p14.tiff", device = "tiff", dpi = 600,
       path = "./Output/", width = 7, height = 4, units = "in")

## Outer hair cell markers 
cluster2.markers <- FindMarkers(cochlea_P14, ident.1 = 2, min.pct = 0.25, only.pos = T)
write.csv(cluster2.markers, "./cluster2_OHC.markers_v1.csv")
