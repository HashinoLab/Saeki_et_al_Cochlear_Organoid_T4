setwd("C:/ ... /CochO_T4_project/")
library(Seurat)
library(SeuratData)
library(patchwork)
library(tidyverse)
options(future.globals.maxSize = 64000*1024^2)

# Load functions stored in the external file
source("./Script/functions.R")

# Load data for Seurat objects
d109_ctr <- Read10X("./Data/d109_ctr/") %>% CreateSeuratObject()
d109_t4 <- Read10X("./Data/d109_t4/") %>% CreateSeuratObject()
d140_t4 <- Read10X("./Data/d140_t4/") %>% CreateSeuratObject()

# Add metadata
d109_ctr$cond <- "CTR"
d109_t4$cond <- "T4"

d109_t4$day <- "109"
d140_t4$day <- "140"

# List for filtering
list_all <- list(d109_ctr, d109_t4, d140_t4)

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

# Merge data
MGD_d109 <- merge(list_all[[1]], list_all[[2]]) %>% JoinLayers()
MGD_d109[["RNA"]] <- split(MGD_d109[["RNA"]], f = MGD_d109$cond)
MGD_d109

MGD_t4 <- merge(list_all[[2]], list_all[[3]]) %>% JoinLayers()
MGD_t4[["RNA"]] <- split(MGD_t4[["RNA"]], f = MGD_t4$day)
MGD_t4

rm(d109_ctr, d109_t4, d140_t4, list_all); gc()

#### Seurat workflow ####
MGD_d109 <- SCTransform(MGD_d109)
MGD_d109 <- RunPCA(MGD_d109)

ElbowPlot(MGD_d109, ndims = 50, reduction = "pca")
DimHeatmap(MGD_d109, dims = 1:30, cells = 200, balanced = T)

MGD_t4 <- SCTransform(MGD_t4)
MGD_t4 <- RunPCA(MGD_t4)

ElbowPlot(MGD_t4, ndims = 50, reduction = "pca")
DimHeatmap(MGD_t4, dims = 1:30, cells = 200, balanced = T)

MGD_d109 <- FindNeighbors(MGD_d109, dims = 1:30) %>%
  FindClusters(resolution = 0.15) %>%
  RunUMAP(n.neighbors = 20, dims = 1:30)
p <- DimPlot(MGD_d109, pt.size = 0.4,
             cols = c(
               "#f26d43", #0 Immature cochlear HCs 1
               "#65c0a2", #1 Transitional
               "#5f4fa2", #2 Neuronal
               "#4292c3", #3 Vestibular HCs
               "#f02550", #4 OHCs
               "#f5b529", #5 Immature cochlear HCs 2
               "#07305d"  #6 Mesenchymal
             )) + umap_theme()
ggsave(paste0("dimplot_d109.tiff"), plot = p, device = "tiff",
       path = "./Output/", width = 4, height = 4, units = "in", dpi = 600)
p <- DimPlot(MGD_d109, pt.size = 0.4, group.by = "cond", cols = cols.cond) + umap_theme()
ggsave(paste0("dimplot_d109_grpd.tiff"), plot = p, device = "tiff",
       path = "./Output/", width = 4, height = 4, units = "in", dpi = 600)

MGD_t4 <- FindNeighbors(MGD_t4, dims = 1:30) %>%
  FindClusters(resolution = 0.15) %>%
  RunUMAP(n.neighbors = 20, dims = 1:30)
p <- DimPlot(MGD_t4, pt.size = 0.4,
             cols = c(
               "#f02550", #0 OHC
               "#4292c3", #1 Vestibular HCs
               "#a9dba1", #2 Transitional 2
               "#65c0a2", #3 Transitional 1
               "#5f4fa2", #4 Neuronal
               "#07305d", #5 Mesenchymal
               "#f26d43" #6 Immature cochlear HCs
             )) + umap_theme()
ggsave(paste0("dimplot_t4.tiff"), plot = p, device = "tiff",
       path = "./Output/", width = 4, height = 4, units = "in", dpi = 600)
p <- DimPlot(MGD_t4, pt.size = 0.4, group.by = "day", cols = cols.day) + umap_theme()
ggsave(paste0("dimplot_t4_grpd.tiff"), plot = p, device = "tiff",
       path = "./Output/", width = 4, height = 4, units = "in", dpi = 600)

#### Feature plots ####
genes <- c(
  # pan-HC
  "ATOH1", "GFI1", "MYO7A", "GRXCR1",
  # OHC
  "LMOD3", "SLC26A5",
  # IHC
  "OTOF", "LDB3", 
  # Cochlear HCs
  "GATA3", "NR2F1",
  # Vestibular HCs
  "PCDH20", "SKOR1", 
  # SC/Transitional
  "SOX2",  "LGR5",  "FGFR3", 
  # Neuronal
  "STMN2","ELAVL3"
)

for (i in genes) {
  p <- FeaturePlot(MGD_d109, i, pt.size = 0.4, max.cutoff = 2,
                   cols = c("lightgrey", "#002385")) + feat_theme_nolegend()
  ggsave(paste0(i, "_d109.tiff"), plot = p, device = "tiff",
         path = "./Output/", width = 4, height = 4, units = "in", dpi = 600)
}

for (i in genes) {
  p <- FeaturePlot(MGD_t4, i, pt.size = 0.4, max.cutoff = 2,
                   cols = c("lightgrey", "#002385")) + feat_theme_nolegend()
  ggsave(paste0(i, "_t4.tiff"), plot = p, device = "tiff",
         path = "./Output/", width = 4, height = 4, units = "in", dpi = 600)
}

#### DEGs ####
d109.deg <- FindAllMarkers(MGD_d109, logfc.threshold = 0.69, min.pct = 0.3, only.pos = T)
top.d109 <- d109.deg %>%
  select(cluster, gene, pct.1, pct.2, avg_log2FC, p_val_adj) %>%
  arrange(desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  slice(1:50)
write.csv(top.d109, "./Output/top.d109.csv")

t4.deg <- FindAllMarkers(MGD_t4, logfc.threshold = 0.69, min.pct = 0.3, only.pos = T)
top.t4 <- t4.deg %>%
  select(cluster, gene, pct.1, pct.2, avg_log2FC, p_val_adj) %>%
  arrange(desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  slice(1:50)
write.csv(top.t4, "./Output/top.t4.csv")

#### Dot plots ####
genes <- c(
  # pan-HC
  "ATOH1", "MYO7A", "MYO15A", "MYO6", "PCP4", "ANXA4", 
  "ESPN",  "CIB2", "XIRP2", "PCDH15", "USH2A", "BARHL1",
  "GRXCR1", "ESPNL", "JAG2", "GFI1", "LHX3", "LMO7", "ISL1",
  # OHC
   "LMOD3",  "LMOD1", "OCM", "IKZF2",
  "SLC26A5", "STRIP2", "SIX2", "KCNQ4", "INSM1", "BCL11B", "EFNA5", "KCNB1", "AQP11",
  # IHC
  "OTOF", "LDB3", "DNAJC5B", "TBX2", "HDAC3", "SLC17A8", 
  "BRIP1", "FGF8", "CALB2", 
  # Cochlear
  "GATA3", "NR2F1", "NR2F2", 
  # Vestibular
  "TEKT1", "PCDH20", "TCTEX1D1", "SKOR1", 
  "MEIS2", "VEPH1", "NDRG1", "NEUROD6",
  # SC/Transitional
  "SPARCL1", "SOX2", "COL9A2", 
  "LGR5", "KRT19", "GJB2", "GJB6", "FGFR3", "PROX1",
  # Neuronal
  "STMN2","ELAVL3","DCX",
  # Neural crest
  "PAX3","FABP7","NES",
  #Thyroid hormone receptors
  "THRA", "THRB"
)

DotPlot(MGD_d109, features = rev(genes),
        col.min = -2, col.max = 2, dot.min = 0, dot.scale = 8) +
  coord_flip() +
  scale_y_discrete(position = "right",
                   limits = c("4", "0", "5", "3", "1", "6", "2")) +
  theme(
    axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(face = "italic"),
    axis.title = element_blank()
    )
ggsave("output/dot_d109.tiff", device = "tiff", dpi = 600,
       width = 6, height = 10.5, units = "in")

DotPlot(MGD_t4, features = rev(genes),
        col.min = -2, col.max = 2, dot.min = 0, dot.scale = 8) +
  coord_flip() +
  scale_y_discrete(position = "right",
                   limits = c("0", "6", "2", "1", "3", "5", "4")) +
  theme(
    axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(face = "italic"),
    axis.title = element_blank()
  )
ggsave("output/dot_t4.tiff", device = "tiff", dpi = 600,
       width = 6, height = 10.5, units = "in")

#### Violin Plots ####
genes = c("SOX2", "ATOH1", "LMOD3", "SLC26A5")
for (i in genes){
  p <- VlnPlot(MGD_d109, i, idents = "4",
               split.plot = T, split.by = "cond",
               cols = cols.cond, pt.size = 0
               ) & vln_theme()
    p$layers[[1]]$aes_params$size <- 1
    ggsave(paste0("output/Violin_d109_cl4_" , i, ".tiff"), plot = p, device = "tiff",
           width = 3, height = 3, units = "in", dpi = 600)
}

for (i in genes){
  p <- VlnPlot(MGD_t4, i, idents = "0",
               split.plot = T, split.by = "day",
               cols = cols.day, pt.size = 0
               ) & vln_theme()
  p$layers[[1]]$aes_params$size <- 1
  ggsave(paste0("output/Violin_t4_cl0_" , i, ".tiff"), plot = p, device = "tiff",
         width = 3, height = 3, units = "in", dpi = 600)
}

#### Volcano plots ####
library(EnhancedVolcano)

MGD_d109 <- PrepSCTFindMarkers(MGD_d109)
MGD_t4 <- PrepSCTFindMarkers(MGD_t4)

DEG_immature_v_ohc <- FindMarkers(MGD_d109, ident.1 = "4", ident.2 = "0", verbose = F,
                                  logfc.threshold = 0, test.use = "wilcox", min.pct = 0.1,
                                  only.pos = F, return.thresh = 0.01) 

logfc.threshold <- 1
p.threshold = 1e-10
labels <- rownames(DEG_immature_v_ohc)
selectLab <- c("SOX2", "ATOH1", "INSM1", "ATOH7", "CORO2A", "CPM", "CHRNG", "GNG2",
               "SLC26A5", "IKZF2", "OCM", "LMOD3", 
               "TMC1", "TMC2", "CHRNA9", "PJVK", "P2RX2", "LOXHD1", "CLIC5", "TMIE",
               "SERPINB6", "GABRB2", "CDC14A", "MYO3A", "TJP2", "STRC", "MSRB3","CEMIP","CABP2")

# color for each gene
color.key <- ifelse(DEG_immature_v_ohc$avg_log2FC > logfc.threshold & DEG_immature_v_ohc$p_val_adj < p.threshold, "#B656D4",
                    ifelse(DEG_immature_v_ohc$avg_log2FC < -logfc.threshold & DEG_immature_v_ohc$p_val_adj < p.threshold, "#56A1FE", "lightgrey"))

# colors for highlighting
color.key[(rownames(DEG_immature_v_ohc) %in% selectLab) & DEG_immature_v_ohc$avg_log2FC > 1] <-"red4"
color.key[(rownames(DEG_immature_v_ohc) %in% selectLab) & DEG_immature_v_ohc$avg_log2FC < -1] <- "blue4"

# naming color.key
names(color.key)[color.key == "#B656D4"] <- "d109 T4"
names(color.key)[color.key == "#56A1FE"] <- "d109 CTR"
names(color.key)[color.key == "lightgrey"] <- "Below_Threshold"
names(color.key)[color.key == "red4"] <- "SelectedT4"
names(color.key)[color.key == "blue4"] <- "SelectedCTR"

p <-
  EnhancedVolcano(DEG_immature_v_ohc,
                  lab = labels,
                  x = "avg_log2FC",
                  y = "p_val_adj",
                  title = "Day 109 immature HCs vs. OHCs",
                  subtitle = NULL,
                  pCutoff = 1e-10,
                  FCcutoff = 1,
                  pointSize = 4,
                  colCustom = color.key,
                  selectLab = selectLab,
                  legendPosition = "right",
                  labSize = 5,
                  drawConnectors = T,
                  arrowheads = F,
                  maxoverlapsConnectors = Inf,
                  labFace = "italic") +
  xlim(-6, 6) + 
  NoLegend() + 
  theme(
    axis.title = element_text(color="#000000", size = 16),
    axis.text.x = element_text(color="#000000", size = 12),
    axis.text.y = element_text(color="#000000", size = 12),
    plot.margin = unit(c(0,0,0,0), "in")
  )

ggsave("d109_volcano_immaturevOHC.pdf", p, "pdf", "./Output/", width = 8, height = 5, unit = "in", dpi = 600)

DEG_ohc_109_v_140 <- FindMarkers(MGD_t4, ident.1 = "140", ident.2 = "109", verbose = F,
                                 group.by = "day", subset.ident = "0",
                                 logfc.threshold = 0, test.use = "wilcox", min.pct = 0.1,
                                 only.pos = F, return.thresh = 0.01,
                                 recorrect_umi = F)

logfc.threshold <- 1
p.threshold = 1e-10
labels <- rownames(DEG_ohc_109_v_140)
selectLab <- c("TMC1", "ATOH1", "SOX2", "SLC26A5","OCM","LMOD3", 
               "CALB1","SALL1","BHLHE40","RORB","PVALB","EGR1","FOS","MAFF","FZD4","P2RX2","GRM7","SMPX",
               "HES6","CCER2","RBFOX1","DLK2","PAX2")

# color for each gene
color.key <- ifelse(DEG_ohc_109_v_140$avg_log2FC > logfc.threshold & DEG_ohc_109_v_140$p_val_adj < p.threshold, "#C32500",
                    ifelse(DEG_ohc_109_v_140$avg_log2FC < -logfc.threshold & DEG_ohc_109_v_140$p_val_adj < p.threshold, "#9666D4", "lightgrey"))

# colors for highlighting
color.key[(rownames(DEG_ohc_109_v_140) %in% selectLab) & DEG_ohc_109_v_140$avg_log2FC > 1] <- "yellow1"
color.key[(rownames(DEG_ohc_109_v_140) %in% selectLab) & DEG_ohc_109_v_140$avg_log2FC < -1] <- "blue4"

# naming color.key
names(color.key)[color.key == "#C32500"] <- "d140 T4"
names(color.key)[color.key == "#9666D4"] <- "d109 T4"
names(color.key)[color.key == "lightgrey"] <- "Below_Threshold"
names(color.key)[color.key == "yellow1"] <- "Selected d140"
names(color.key)[color.key == "blue4"] <- "Selected d109"

p <-
  EnhancedVolcano(DEG_ohc_109_v_140,
                  lab = labels,
                  x = "avg_log2FC",
                  y = "p_val_adj",
                  title = "OHCs; T4 d109 vs. d140",
                  subtitle = NULL,
                  pCutoff = 1e-10,
                  FCcutoff = 1,
                  pointSize = 4,
                  colCustom = color.key,
                  selectLab = selectLab,
                  legendPosition = "right",
                  labSize = 5,
                  drawConnectors = T,
                  arrowheads = F,
                  maxoverlapsConnectors = Inf,
                  labFace = "italic") +
  xlim(-6, 6) + 
  NoLegend() + 
  theme(
    axis.title = element_text(color="#000000", size = 16),
    axis.text.x = element_text(color="#000000", size = 12),
    axis.text.y = element_text(color="#000000", size = 12),
    plot.margin = unit(c(0,0,0,0), "in")
  );p
ggsave("T4_volcano_OHC_D109vD140.pdf", p, "pdf", "./Output/", width = 8, height = 5, unit = "in", dpi = 600)

#### Bar chart ####
library(data.table)
library(magrittr)
library(reshape2)

N.CTR <- ncol(subset(MGD_d109, cond == "CTR"))
N.T4 <- ncol(subset(MGD_d109, cond == "T4"))

comp <- MGD_d109@meta.data %>% as.data.table()
comp <- comp[, .N, by = c("seurat_clusters", "cond")]
comp <- mutate(comp, rate = if_else(cond == "CTR", N/N.CTR, N/N.T4))

comp$cond <- factor(data$cond,
                    levels = c("CTR", "T4"))
comp$seurat_clusters <- factor(data$seurat_clusters,
                    levels = c("4", "0", "5", "3", "1", "6", "2"))

ggplot(comp) + 
  geom_bar(aes(x = cond, y = rate, fill = seurat_clusters),
           position="fill", stat = "identity",
           color = "#ffffff", size = 0.8, width = 0.8) + 
  scale_fill_manual(values=c(
                             "#f02550",#4 Outer HCs
                             "#f26d43",#0 Immature HCs 1
                             "#f5b529",#5 Immature HCs 2
                             "#4292c3",#3 Vestibular HCs
                             "#65c0a2",#1 Transitional
                             "#5f4fa2",#2 Neuronal
                             "#07305d" #6 Mesenchymal
  ))+ 
  theme(panel.background = element_rect(fill = "#ffffff"),
        legend.position = "none",
        axis.title = element_blank()
  )
ggsave("percentage_d109.tiff", device = "tiff",
       path = "./Output/",
       width = 4, height = 5, units = "in", dpi = 600)

N.4 <- ncol(subset(MGD_d109, idents = "4"))
N.0 <- ncol(subset(MGD_d109, idents = "0"))
N.5 <- ncol(subset(MGD_d109, idents = "5"))
N.3 <- ncol(subset(MGD_d109, idents = "3"))

comp2 <- comp %>% select(seurat_clusters:N) %>% filter(seurat_clusters == "4" |
                                                         seurat_clusters == "0" |
                                                         seurat_clusters == "5" |
                                                         seurat_clusters == "3")

comp2 <- mutate(comp2, rate = if_else(seurat_clusters == "4", N/N.4, 
                                      if_else(seurat_clusters == "0", N/N.0,
                                              if_else(seurat_clusters == "5", N/N.5, N/N.3))))

comp2$cond <- factor(data2$cond,
                     levels = c("T4","CTR"))
comp2$cluster <- factor(data2$cluster,
                               levels = rev(c("4", "0", "5", "3")))

ggplot(comp2) + 
  geom_bar(aes(x = seurat_clusters, y = rate, fill = cond),
           position="fill", stat = "identity",
           color = "#ffffff", size = 0.8, width = 0.8) + 
  coord_flip() +
  scale_fill_manual(values = rev(cols.cond))+ 
  theme(panel.background = element_rect(fill = "#ffffff"),
        legend.position = "none",
        axis.title = element_blank()
  )

ggsave("percentage_d109.tiff", device = "tiff",
       path = "./Output/",
       width = 5, height = 7, units = "in", dpi = 600)

### Heat maps ###
MGD_d109_s <- ScaleData(MGD_d109, features = rownames(MGD_d109))

MGD_d109_0_4 <- subset(MGD_d109_s, idents = c("0","4"))
MGD_d109_4 <- subset(MGD_d109_s, idents = c("4"))

ohc_gene_list <- read.csv("./OHC_Kolla_2020_NatCommun.csv")
ohc_gene_p14 <- read.csv("./OHC_Xu_2022_FrontCellNeurosci.csv")

genes_e16 <- toupper(ohc_gene_list$E16OHC)
genes_e16 <- genes_e16[genes_e16 != "POU4F3" & genes_e16 != ""]

genes_p7 <- toupper(ohc_gene_list$P7OHC)
genes_p7 <- genes_p7[genes_p7 != "POU4F3" & genes_p7 != ""]

genes_p14 <- toupper(ohc_gene_p14$P14OHC)
genes_p14 <- genes_p14[genes_p14 != "POU4F3" & genes_p14 != ""]

MGD_d109_OHC_s <- ScaleData(MGD_d109_OHC, features = rownames(MGD_d109_OHC))

# For the following plot, (168) genes were used for the plot, omitting the following (18) genes:
# 5033430I15RIK, 2510002D24RIK, MT-RNR2, ZFP667, 1110008P14RIK, DNAIC2, TOMT, GM12892, A730017C20RIK, KIF19A,
# ATP6V0E, BC030867, PLA2G16, DEFB25, GM30191, GM2694, CCL21A, GRP
p <- DoHeatmap(MGD_d109_0_4, genes_e16,  slot = "scale.data")+
  theme(panel.background = element_rect(fill = "#ffffff"),
        legend.position = "none",
        text = element_blank())
ggsave(paste0("heat_e16_s.tiff"), plot = p, device = "tiff",
       path = "./Output/", width = 6, height = 10, units = "in", dpi = 600)
  
# For the following plot, (418) genes were used for the plot, omitting the following (60) genes:
# TRP53INP1, TLDC1, GRCC10, ATP5E, MMP23, TXN1, FAM57B, 0610033M10RIK, 2210013O21RIK, AC154176.2,
# ATP5H, ARHGAP27OS2, 1700025G04RIK, ZFP503, ANK, GM42418, GM11942, AC153524.3, 1700001C02RIK, 1110008P14RIK,
# AC154305.3, ATP5J2, AC153524.4, 2010107E04RIK, 1810037I17RIK, E130308A19RIK, BMYC, ATP5K, B230219D22RIK, ZFP385A,
# GM29483, A330048O09RIK, ZFP423, 9530077C05RIK, FAM213A, 43717, ATP5J, SERPINB6A, DFNB59, GM2694,
# USMG5, MYO15, HIST3H2BA, MT-CYTB, 4930558C23RIK, MT-RNR1, 2210011C24RIK, 5730508B09RIK, 9530085L11RIK, GM4825,
# SERF1, CAR7, CAR2, URAH, DEFB25, TOMT, AI593442, GM30191, MT-RNR2, IGHM
p <- DoHeatmap(MGD_d109_0_4, genes_p7,  slot = "scale.data", disp.max = 2, disp.min = -2)+
  theme(panel.background = element_rect(fill = "#ffffff"),
        legend.position = "none",
        text = element_blank());p
ggsave(paste0("heat_p7_s.tiff"), plot = p, device = "tiff",
       path = "./Output/", width = 6, height = 10, units = "in", dpi = 600)

# For the following plot, (440) genes were used for the plot, omitting the following (60) genes:
# AHSA2, ZFP516, C9ORF72, ZFP385A, 2900026A02RIK, 4833439L19RIK, GM1673, ZFP101, 2510002D24RIK, B230219D22RIK,
# TMEM56, ZFP938, 2900093K20RIK, 9330020H09RIK, SERF1, 9530077C05RIK, PAKAP.1, FCOR, AKR1B3, CKMT1,
# 1700025G04RIK, HIST3H2BA, GM43672, ZFP618, 2810002D19RIK, 43717, IGHM, AI593442, GM2694, 1700001K23RIK,
# 2210011C24RIK, GM34653, DEFB25, TMEM254A, GM42702, FAM57B, TTC7, GM43260, A930033H14RIK, A530016L24RIK,
# TOMT, GM43259, CAR7, 4930426I24RIK, GM45716, GM42701, GM32742, 9530097N15RIK, GM4131, GM42696,
# MYO15, D7ERTD443E, GM43258, ZFP365, ARHGAP27OS2, 0610033M10RIK, MMP23, 4930558C23RIK, GM5111, 1600022D10RIK
p <- DoHeatmap(MGD_d109_0_4, head(genes_p14, n = 500),  slot = "scale.data")+
  theme(panel.background = element_rect(fill = "#ffffff"),
        legend.position = "none",
        text = element_blank())
ggsave(paste0("heat_p14_s.tiff"), plot = p, device = "tiff",
       path = "./Output/", width = 6, height = 10, units = "in", dpi = 600)

#### Saving Seurat objects ####
saveRDS(MGD_d109, "./Data/MGD_d109.rds")
saveRDS(MGD_t4, "./Data/MGD_t4.rds")

saveRDS(MGD_d109_0_4, "./MGD_d109_0_4.RDS")
saveRDS(MGD_t4_0_6, "./MGD_t4_0_6.RDS")
