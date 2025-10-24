
library(Seurat)
library(patchwork)
library(tidyverse)
library(cowplot)
options(future.globals.maxSize = 64000*1024^2)

# Load functions stored in the external file
source("./Script/functions_20250217.R")

# Load data for Seurat objects
d109_ctr_Ep_neg <- Read10X("./Data/d109_control_EPCAMneg/") %>% CreateSeuratObject()
d109_ctr_Ep_pos <- Read10X("./Data/d109_control_EPCAMpos/") %>% CreateSeuratObject()
d109_t4_Ep_neg <- Read10X("./Data/d109_T4_EPCAMneg/") %>% CreateSeuratObject()
d109_t4_Ep_pos <- Read10X("./Data/d109_T4_EPCAMpos/") %>% CreateSeuratObject()

# Add metadata
d109_ctr_Ep_neg$cond <- "CTR"
d109_ctr_Ep_pos$cond <- "CTR"
d109_t4_Ep_neg$cond <- "Thyroxine"
d109_t4_Ep_pos$cond <- "Thyroxine"

d109_ctr_Ep_neg$FACS <- "EPCAM-"
d109_ctr_Ep_pos$FACS <- "EPCAM+"
d109_t4_Ep_neg$FACS <- "EPCAM-"
d109_t4_Ep_pos$FACS <- "EPCAM+"


# List for filtering
list_all <- list(d109_t4_Ep_neg, d109_t4_Ep_pos,
                 d109_ctr_Ep_neg, d109_ctr_Ep_pos)

# Filter each dataset according to mitochondrial percentage
counts.list.post <- counts.list.pre <- numeric()
for (i in 1:length(list_all)) {
  counts.list.pre[i] <- length(Cells(list_all[[i]]))
  list_all[[i]] <- filterPlot(list_all[[i]])
  counts.list.post[i] <- length(Cells(list_all[[i]]))
}
diff <- counts.list.pre - counts.list.post; diff

# Same thing for RPL27
counts.list.post.2 <- numeric()
for (i in 1:length(list_all)) {
  list_all[[i]] <- hk.filter(list_all[[i]], hk.gene = "RPL27", subset = T, std.dev = 2)
  counts.list.post.2[i] <- length(Cells(list_all[[i]]))
}
diff.2 <- counts.list.post - counts.list.post.2; diff.2

length(Cells(list_all[[1]]))
length(Cells(list_all[[2]]))
length(Cells(list_all[[3]]))
length(Cells(list_all[[4]]))

list(d109_t4_Ep_neg, d109_t4_Ep_pos,
     d109_ctr_Ep_neg, d109_ctr_Ep_pos)

d109_Ep_neg <- merge(list_all[[1]], list_all[[3]])
d109_Ep_pos <- merge(list_all[[2]], list_all[[4]])

rm(d109_ctr_Ep_neg,
   d109_ctr_Ep_pos,d109_t4_Ep_neg,
   d109_t4_Ep_pos)


########### EPCAM- cells #################

d109_Ep_neg <- SCTransform(d109_Ep_neg)
d109_Ep_neg <- RunPCA(d109_Ep_neg)
d109_Ep_neg <- FindNeighbors(d109_Ep_neg, dims = 1:20) %>%
  FindClusters(resolution = 0.11) %>%  
  RunUMAP(n.neighbors = 20, dims = 1:20)

saveRDS(d109_Ep_neg, "./d109_EPCAM_neg/Data/d109_Ep_neg.rds")

Idents(d109_Ep_neg) <- factor(Idents(d109_Ep_neg),
                              levels = c("0", "1", "2", "3", "4", "5", "6", "7"))


p <- DimPlot(d109_Ep_neg, pt.size = 0.6, split.by = "cond",
             cols = c("#E41A1C", #0
                      "#377EB8", #1
                      "#984EA3", #2
                      "#4DAF4A", #3
                      "#BDBF33", #4
                      "#FF7F00", #5
                      "#999999", #6
                      "#A65628")) + umap_theme()
ggsave(paste0("dimplot_d109_Ep_neg_split.by.cond.tiff"), plot = p, device = "tiff",
       path = "./Output/", width = 8, height = 4, units = "in", dpi = 600)

p <- DimPlot(d109_Ep_neg, pt.size = 0.6, 
             cols = c("#E41A1C", #0
                      "#377EB8", #1
                      "#984EA3", #2
                      "#4DAF4A", #3
                      "#BDBF33", #4
                      "#FF7F00", #5
                      "#999999", #6
                      "#A65628")) + umap_theme()
ggsave(paste0("dimplot_d109_Ep_neg.tiff"), plot = p, device = "tiff",
       path = "./Output/", width = 4, height = 4, units = "in", dpi = 600)


DotPlot <- DotPlot(d109_Ep_neg, features = c("STMN2","TUBB3","NEUROD1","DCX", # neuron
                                             "S100B","PLP1","SOX10","GAP43","FABP7","GFAP", # glia
                                             "MSC","MYF5","TNNT1", # skeletal muscle
                                             "COL3A1","POSTN","COL1A1", # fibroblast
                                             "COL2A1",
                                             "DIO2","DIO3","THRA","THRB",
                                             "SLCO1A2","SLCO3A1","SLCO4A1","SLCO4C1",
                                             "SLC7A5","SLC7A8","SLC16A2","SLC16A10"),
                   split.by = "cond", cols = c("#401cfc", "#401cfc")) + theme(
                     axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1),
                     axis.text.x = element_text(face = "italic"),
                     axis.title = element_blank()
                   )
DotPlot
ggsave("Output/dot_d109_Ep_neg_split.tiff", device = "tiff", dpi = 600,
       width = 9, height = 5, units = "in")


########### EPCAM+ cells #################

d109_Ep_pos <- SCTransform(d109_Ep_pos)
d109_Ep_pos <- RunPCA(d109_Ep_pos)
d109_Ep_pos <- FindNeighbors(d109_Ep_pos, dims = 1:30) %>%
  FindClusters(resolution = 0.3) %>%   
  RunUMAP(n.neighbors = 20, dims = 1:30)

saveRDS(d109_Ep_pos, "./Data/d109_Ep_pos.rds")

Idents(d109_Ep_pos) <- factor(Idents(d109_Ep_pos),
                              levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))

p <- DimPlot(d109_Ep_pos, pt.size = 0.8, split.by = "cond",
             cols = c("#E41A1C", #0 "cochlear medial supporting cells",
                      "#377EB8", #1 "hair cells 1",
                      "#4DAF4A", #2 "cochlear epithelial cells",
                      "#984EA3", #3 "cochlear roof epithelium",
                      "#FF7F00", #4 "cochlear lateral supporting cells",
                      "#A65628", #5 "cochlear UBE2C+ cells",
                      "#B2DF8A", #6 "neurons",
                      "#FFD92F", #7 "fibrocytes", 
                      "#F781BF", #8 "spiral ganglion neurons",
                      "#66C2A5", #9 "hair cells 2",
                      "#999999" #10 "keratinocytes"
             )) + umap_theme()
ggsave(paste0("dimplot_d109_Ep_pos_split.by.cond.tiff"), plot = p, device = "tiff",
       path = "./Output/", width = 8, height = 4, units = "in", dpi = 600)


p <- DimPlot(d109_Ep_pos, pt.size = 0.8, 
             cols = c("#E41A1C", #0 "cochlear medial supporting cells",
                      "#377EB8", #1 "hair cells 1",
                      "#4DAF4A", #2 "cochlear epithelial cells",
                      "#984EA3", #3 "cochlear roof epithelium",
                      "#FF7F00", #4 "cochlear lateral supporting cells",
                      "#A65628", #5 "cochlear UBE2C+ cells",
                      "#B2DF8A", #6 "neurons",
                      "#FFD92F", #7 "fibrocytes", 
                      "#F781BF", #8 "spiral ganglion neurons",
                      "#66C2A5", #9 "hair cells 2",
                      "#999999" #10 "keratinocytes"
             )) + umap_theme()
ggsave(paste0("dimplot_d109_Ep_pos.tiff"), plot = p, device = "tiff",
       path = "./Output/", width = 4, height = 4, units = "in", dpi = 600)


p <- DimPlot(d109_Ep_pos, pt.size = 0.8, group.by = "cond",shuffle = T,
             cols = c("#56A1FE", "#B656D4")) + umap_theme()
ggsave(paste0("dimplot_d109_Ep_pos_merged.tiff"), plot = p, device = "tiff",
       path = "./Output/", width = 4, height = 4, units = "in", dpi = 600)


### dotplot split by cond ###
levels(Idents(d109_Ep_pos))
Idents(d109_Ep_pos) <- factor(Idents(d109_Ep_pos),
                              levels = c("1", "9", "0", "4", "2", "3", "5", "6", "8", "7", "10"))


DotPlot <- DotPlot(d109_Ep_pos, features = c("ATOH1", # hair cells
                                             "TECTB","LGR5","FGFR3","LFNG", # supporting cells
                                             "LMX1A","PAX2", #cochlear roof epithelium
                                             "UBE2C",
                                             "NEUROD1", # spiral ganglion
                                             "DCX", # neuron
                                             "TP63", #epidermis 
                                             "DIO2","DIO3","THRA","THRB",
                                             "SLCO1A2","SLCO3A1","SLCO4A1","SLCO4C1",
                                             "SLC7A5","SLC7A8","SLC16A2","SLC16A10"),
                   split.by = "cond", cols = c("#401cfc", "#401cfc")) + theme(
                     axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1),
                     axis.text.x = element_text(face = "italic"),
                     axis.title = element_blank()
                   )
DotPlot
ggsave("Output/dot_d109_Ep_pos_split_v2.tiff", device = "tiff", dpi = 600,
       width = 9, height = 7, units = "in")


#################### Cell Chat analysis ####################

library(CellChat)

d109_Ep_pos <- RenameIdents(d109_Ep_pos, 
                            `0` = "medial SCs",
                            `1` = "hair cells 1",
                            `2` = "cochlear epithelium",
                            `3` = "cochlear roof epithelium",
                            `4` = "lateral SCs",
                            `5` = "UBE2C+ cells",
                            `6` = "neurons",
                            `7` = "keratinocytes",
                            `8` = "spiral ganglion",
                            `9` = "hair cells 2",
                            `10` = "fibrocytes")


d109_Ep_pos_CTR <- subset(d109_Ep_pos, subset = cond == "CTR")
d109_EPCAM_pos_T4 <- subset(d109_Ep_pos, subset = cond == "Thyroxine")

DimPlot(d109_Ep_pos_CTR, reduction = "umap",label = T, repel = TRUE, label.size = 4) + NoLegend()
DimPlot(d109_EPCAM_pos_T4, reduction = "umap",label = T, repel = TRUE, label.size = 4) + NoLegend()

DefaultAssay(d109_Ep_pos_CTR) <- "RNA"
DefaultAssay(d109_EPCAM_pos_T4) <- "RNA"

Idents(d109_Ep_pos_CTR) <- "seurat_clusters"
Idents(d109_EPCAM_pos_T4) <- "seurat_clusters"

data.input.CTR <- GetAssayData(d109_Ep_pos_CTR, assay = "RNA", layer = "counts") 
data.input.T4 <- GetAssayData(d109_EPCAM_pos_T4, assay = "RNA", layer = "counts") 

labels.CTR <- Idents(d109_Ep_pos_CTR)
labels.T4 <- Idents(d109_EPCAM_pos_T4)

meta2.CTR <- data.frame(group = labels.CTR, row.names = names(labels.CTR)) 
meta2.T4 <- data.frame(group = labels.T4, row.names = names(labels.T4)) 

unique(meta2.CTR$group)
unique(meta2.T4$group)

cellchat.CTR <- createCellChat(object = data.input.CTR, meta = meta2.CTR, group.by = "group")
cellchat.T4 <- createCellChat(object = data.input.T4, meta = meta2.T4, group.by = "group")

groupSize <- as.numeric(table(cellchat.CTR@idents))
groupSize <- as.numeric(table(cellchat.T4@idents))

cellchat.CTR@DB <- CellChatDB.human
cellchat.T4@DB  <- CellChatDB.human

cellchat.CTR <- subsetData(cellchat.CTR)
cellchat.CTR <- identifyOverExpressedGenes(cellchat.CTR)
cellchat.CTR <- identifyOverExpressedInteractions(cellchat.CTR)
cellchat.CTR <- computeCommunProb(cellchat.CTR)
cellchat.CTR <- filterCommunication(cellchat.CTR)
cellchat.CTR <- computeCommunProbPathway(cellchat.CTR)
cellchat.CTR <- aggregateNet(cellchat.CTR)
cellchat.CTR <- netAnalysis_computeCentrality(cellchat.CTR)

cellchat.T4 <- subsetData(cellchat.T4)
cellchat.T4 <- identifyOverExpressedGenes(cellchat.T4)
cellchat.T4 <- identifyOverExpressedInteractions(cellchat.T4)
cellchat.T4 <- computeCommunProb(cellchat.T4)
cellchat.T4 <- filterCommunication(cellchat.T4)
cellchat.T4 <- computeCommunProbPathway(cellchat.T4)
cellchat.T4 <- aggregateNet(cellchat.T4)
cellchat.T4  <- netAnalysis_computeCentrality(cellchat.T4)

object.list <- list(CTR = cellchat.CTR, Thyroxine= cellchat.T4)
cellchat.merged <- mergeCellChat(object.list, add.names = names(object.list))

saveRDS(cellchat.CTR, "./Data/cellchat.CTR.rds")
saveRDS(cellchat.T4, "./Data/cellchat.T4.rds")
saveRDS(cellchat.merged, "./Data/cellchat.merged.rds")

gg1 <- compareInteractions(cellchat.merged, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat.merged, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

gg1 <- netVisual_heatmap(cellchat.merged)
gg2 <- netVisual_heatmap(cellchat.merged, measure = "weight")
gg1 + gg2

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

object.list <- lapply(object.list, function(x) {
  x <- netAnalysis_computeCentrality(x)
  return(x)
})

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)


gg1 <- netAnalysis_signalingChanges_scatter(cellchat.merged, idents.use = "hair cells 1", do.label = T ,xlims = c(-0.003, 0.003), ylims = c(-0.008, 0.008),label.size = 6, dot.size = 4)+ theme(
  axis.text = element_text(size = 12),
  axis.title = element_text(size = 14),
  plot.title = element_text(size = 16),
  text = element_text(size = 16))
p <- patchwork::wrap_plots(plots = gg1)

ggsave(paste0("d109_Ep_pos_HC_netAnalysis_signalingChanges_scatter_v3.tiff"), plot = p, device = "tiff",
       path = "./Output/", width = 8, height = 5, units = "in", dpi = 600)


##### dotplot #####
CellChatDB <- CellChatDB.human 
pairs <- CellChatDB$interaction %>% dplyr::filter(pathway_name == "ADGRL")
ligands <- unique(pairs$ligand)
receptors <- unique(pairs$receptor)
genes_ADGRL <- unique(c(ligands, receptors))
print(genes_ADGRL)

genes <- c("FLRT1", "FLRT2", "FLRT3", "NRXN1", "NRXN2", 
           "NRXN3", "TENM1", "TENM2", "TENM3", "TENM4", 
           "UNC5A", "ADGRL1", "ADGRL2", "ADGRL3")

d109_Ep_pos_HC_SC <- subset(d109_Ep_pos, idents = c("1","4","8"))

DotPlot(d109_Ep_pos_HC_SC, features = genes,
        col.min = -2, col.max = 2, dot.min = 0, dot.scale = 8,
        split.by = "cond", cols = c("#401cfc", "#401cfc")) +
  theme(
    axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text.x = element_text(face = "italic"),
    axis.title = element_blank()
  )
ggsave("Output/dot_d109_Ep_pos_ADGRL_pathway.tiff", device = "tiff", dpi = 600,
       width = 8, height = 2.5, units = "in")




