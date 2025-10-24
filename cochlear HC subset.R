

### subclustering of cochlear hair cells

MGD_d109_HC <- subset(MGD_d109, idents = c(0,4,5))

MGD_d109_HC <- SCTransform(MGD_d109_HC)
MGD_d109_HC <- RunPCA(MGD_d109_HC)
MGD_d109_HC <- FindNeighbors(MGD_d109_HC, dims = 1:30) %>%
  FindClusters(resolution = 0.35) %>%   
  RunUMAP(n.neighbors = 20, dims = 1:30)

saveRDS(MGD_d109_HC, "./Data/MGD_d109_HC.rds")


MGD_t4_HC <- subset(MGD_t4, idents = c(0,6))

MGD_t4_HC <- SCTransform(MGD_t4_HC)
MGD_t4_HC <- RunPCA(MGD_t4_HC)
MGD_t4_HC <- FindNeighbors(MGD_t4_HC, dims = 1:30) %>%
  FindClusters(resolution = 0.3) %>%  
  RunUMAP(n.neighbors = 20, dims = 1:30)

saveRDS(MGD_t4_HC, "./Data/MGD_t4_cHC.rds")


############ feature plot ##########

genes <- c("OTOF", "LDB3", "DNAJC5B", "TBX2", "HDAC3", "SLC17A8", 
           "BRIP1", "FGF8", "CALB2", "SLC7A14", "ATP2A3", "SHTN1", "RPRM", 
           "SLC26A5","LMOD3","LBH","IKZF2","SIX2","SALL1","BHLHE40","PVALB",
           "BCL11B","INSM1",
           "ATOH1","SOX2", 
           "GATA3","RORB",
           "PCDH20", "SKOR1")


for (i in genes) {
  p <- FeaturePlot(MGD_t4_HC, i, pt.size = 0.6, max.cutoff = 2, order = T, 
                   cols = c("lightgrey", "#002385")) + feat_theme_nolegend()
  ggsave(paste0(i, "_t4.tiff"), plot = p, device = "tiff",
         path = "./Output_d109_d140_HC_v3/", width = 4, height = 4, units = "in", dpi = 600)
}


for (i in genes) {
  p <- FeaturePlot(MGD_d109_HC, i, pt.size = 0.8, max.cutoff = 2, order = T,
                   cols = c("lightgrey", "#002385")) + feat_theme_nolegend()
  ggsave(paste0(i, "_t4.tiff"), plot = p, device = "tiff",
         path = "./Output_d109_HC_v3/", width = 4, height = 4, units = "in", dpi = 600)
}

### featureplot split by day ###

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

for (i in genes) {
  p <- FeaturePlot(MGD_t4_HC, i, pt.size = 0.8, max.cutoff = 2, order = T, split.by = "day",
                   cols = c("lightgrey", "#002385")) + feat_theme_nolegend()
  ggsave(paste0(i, "_t4.tiff"), plot = p, device = "tiff",
         path = "./Output_d109_d140_HC_v3/featureplot_split.by_day/", width = 8, height = 4, units = "in", dpi = 600)
}


for (i in genes) {
  p <- FeaturePlot(MGD_d109_HC, i, pt.size = 0.8, max.cutoff = 2, order = T, split.by = "cond",
                   cols = c("lightgrey", "#002385")) + feat_theme_nolegend()
  ggsave(paste0(i, "_d109.tiff"), plot = p, device = "tiff",
         path = "./Output_d109_HC_v3/featureplot_split.by_cond/", width = 8, height = 4, units = "in", dpi = 600)
}

 
############ Dot plot ##########

genes <- c("ATOH1","PCP4","GFI1","MYO7A", #pan-hair cell 
           "GATA3","NR2F1","RORB", # cochlear hair cell
           "SLC26A5","IKZF2","OCM","KCNQ4","LMOD3","KCNB1","SIX2","EFNA5","STRIP2", # outer hair cell
           "INSM1","BCL11B", # immature outer hair cell
           "TBX2","BRIP1","FGF8","OTOF","CALB2","SLC7A14","SLC17A8","CASZ1" # inner hair cell
           )


DotPlot(MGD_d109_HC, features = rev(genes),
        col.min = -2, col.max = 2, dot.min = 0, dot.scale = 8) +
  coord_flip() +
  scale_y_discrete(position = "right",
                   limits = c("1", "3", "4", "5", "0", "2")) +
  theme(
    axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(face = "italic"),
    axis.title = element_blank()
  )
ggsave("Output_d109_HC_v3/dot_d109_HC_v2.tiff", device = "tiff", dpi = 600,
       width = 5.5, height = 7, units = "in")


DotPlot(MGD_t4_HC, features = rev(genes),
        col.min = -2, col.max = 2, dot.min = 0, dot.scale = 8) +
  coord_flip() +
  scale_y_discrete(position = "right",
                   limits = c("3", "2", "5", "1", "0", "4")) +
  theme(
    axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(face = "italic"),
    axis.title = element_blank()
  )
ggsave("Output_d109_d140_HC_v3/dot_MGD_t4_HC_v2.tiff", device = "tiff", dpi = 600,
       width = 5.5, height = 7, units = "in")


############ UMAP ##########

p <- DimPlot(MGD_d109_HC, pt.size = 0.6,
             cols = c(
               "#FF7F0E", #0 
               "#1F77B4", #1 
               "#D62728", #2 
               "#2CA02C", #3 
               "#9467DB", #4
               "#8C564B"  #5 
             )) + umap_theme()
ggsave(paste0("dimplot_d109.tiff"), plot = p, device = "tiff",
       path = "./Output_d109_HC_v3/", width = 4, height = 4, units = "in", dpi = 600)


p <- DimPlot(MGD_d109_HC, pt.size = 0.6, split.by = "cond",
             cols = c(
               "#FF7F0E", #0 
               "#1F77B4", #1 
               "#D62728", #2  
               "#2CA02C", #3 
               "#9467DB", #4
               "#8C564B"  #5 
             )) + umap_theme()
ggsave(paste0("dimplot_d109_split.by.cond.tiff"), plot = p, device = "tiff",
       path = "./Output_d109_HC_v3/", width = 8, height = 4, units = "in", dpi = 600)


p <- DimPlot(MGD_d109_HC, pt.size = 0.4, group.by = "cond",cols = cols.cond) + umap_theme()
ggsave(paste0("dimplot_d109_split.tiff"), plot = p, device = "tiff",
       path = "./Output_d109_HC_v3/", width = 4, height = 4, units = "in", dpi = 600)


p <- DimPlot(MGD_t4_HC, pt.size = 0.8,
             cols = c(
               "#8C564B", #0 
               "#FF7F0E", #1 
               "#1F77B4", #2 
               "#2CA02C", #3 
               "#D62728", #4 
               "#9467DB"  #5 
             )) + umap_theme()
ggsave(paste0("dimplot_t4.tiff"), plot = p, device = "tiff",
       path = "./Output_d109_d140_HC_v3/", width = 4, height = 4, units = "in", dpi = 600)

p <- DimPlot(MGD_t4_HC, pt.size = 0.8, split.by = "day",
             cols = c(
               "#8C564B", #0 
               "#FF7F0E", #1  
               "#1F77B4", #2  
               "#2CA02C", #3 
               "#D62728", #4 
               "#9467DB"  #5 
             )) + umap_theme()
ggsave(paste0("dimplot_t4_split.by.day.tiff"), plot = p, device = "tiff",
       path = "./Output_d109_d140_HC_v3/", width = 8, height = 4, units = "in", dpi = 600)


cols.day <- c("#9666d4","#c32500")

p <- DimPlot(MGD_t4_HC, pt.size = 0.8, group.by = "day",cols = cols.day) + umap_theme()
ggsave(paste0("dimplot_t4_split.tiff"), plot = p, device = "tiff",
       path = "./Output_d109_d140_HC_v3/", width = 4, height = 4, units = "in", dpi = 600)

             
#### Violin Plots ####

levels(Idents(MGD_d109_HC))
Idents(MGD_d109_HC) <- factor(Idents(MGD_d109_HC),
                              levels = c(1, 3, 4, 5, 0, 2))

genes <- c("ATOH1","PCP4","GFI1","MYO7A", #pan-hair cell
           "SOX2",
           "GATA3","NR2F1","RORB", # cochlear hair cell
           "SLC26A5","IKZF2","OCM","KCNQ4","LMOD3","KCNB1","SIX2","EFNA5","STRIP2", "LBH", # outer hair cell
           "INSM1","BCL11B", # immature outer hair cell
           "TBX2","BRIP1","FGF8","OTOF","CALB2","SLC7A14","SLC17A8", "CASZ1" # inner hair cell
)


for (i in genes){
  p <- VlnPlot(MGD_d109_HC, i, 
               cols = c(
                 "#1F77B4", #1
                 "#2CA02C", #3
                 "#9467DB", #4
                 "#8C564B", #5
                 "#FF7F0E", #0
                 "#D62728"), #2
               pt.size = 0
  )  & theme(
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
    legend.position = "none"
  )
  p$layers[[1]]$aes_params$size <- 1
  ggsave(paste0("Output_d109_HC_v3/Violin_d109_cl4_" , i, ".tiff"), plot = p, device = "tiff",
         width = 8, height = 4, units = "in", dpi = 600)
}

levels(Idents(MGD_t4_HC))
Idents(MGD_t4_HC) <- factor(Idents(MGD_t4_HC),
                              levels = c(3, 2, 5, 1, 0, 4))

for (i in genes){
  p <- VlnPlot(MGD_t4_HC, i, 
               cols =  c(
                 "#2CA02C", #3
                 "#1F77B4", #2
                 "#9467DB", #5
                 "#FF7F0E", #1
                 "#8C564B", #0
                 "#D62728"), #4
               pt.size = 0
  )+ scale_y_continuous(limits = c(1,5))
  & theme(
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
    legend.position = "none"
  )
  p$layers[[1]]$aes_params$size <- 1
  ggsave(paste0("Output_d109_d140_HC_v3/Violin_t4_cl0_" , i, ".tiff"), plot = p, device = "tiff",
         width = 8, height = 4, units = "in", dpi = 600)
}


#### Volcano plots ####

library(EnhancedVolcano)

DEG_immature_v_ohc <- FindMarkers(MGD_d109_HC, ident.1 = "2", ident.2 = "0", verbose = F,
                                  logfc.threshold = 0, test.use = "wilcox", min.pct = 0.1,
                                  only.pos = F, return.thresh = 0.01) 


write.csv(DEG_immature_v_ohc, "./DEG_ohc_109_markers.csv")

logfc.threshold <- 1
p.threshold = 1e-10
labels <- rownames(DEG_immature_v_ohc)
selectLab <- c("SOX2", "ATOH1", "INSM1", "ATOH7", "CORO2A", "CPM", "CHRNG", "GNG2","BCL11B",
               "SLC26A5", "IKZF2", "OCM", "LMOD3", 
               "TMC1", "TMC2", "CHRNA9", "PJVK", "P2RX2", "LOXHD1", "CLIC5", "TMIE",
               "SERPINB6", "GABRB2", "CDC14A", "MYO3A", "TJP2", "STRC", "MSRB3","CABP2")


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

ggsave("d109_volcano_immaturevOHC_v2.pdf", p, "pdf", "./Output_d109_HC_v3/", width = 8, height = 5, unit = "in", dpi = 600)
ggsave("d109_volcano_immaturevOHC_v2.tiff", p, "tiff", "./Output_d109_HC_v3/", width = 8, height = 5, unit = "in", dpi = 600)




DEG_ohc_109_v_140 <- FindMarkers(MGD_t4_HC, ident.1 = "4", ident.2 = "1", verbose = F,
                                  logfc.threshold = 0, test.use = "wilcox", min.pct = 0.1,
                                  only.pos = F, return.thresh = 0.01) 


write.csv(DEG_ohc_109_v_140, "./Output/DEG_ohc_4_marker.csv")


logfc.threshold <- 1
p.threshold = 1e-10
labels <- rownames(DEG_ohc_109_v_140)
selectLab <- c("TMC1", "ATOH1", "SOX2", "SLC26A5","OCM","LMOD3", "IKZF2", 
               "CALB1","SALL1","BHLHE40","RORB","PVALB","FOS","P2RX2","SMPX","CLIC5","CREB5",
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
ggsave("T4_volcano_OHC_D109vD140_v2.pdf", p, "pdf", "./Output_d109_d140_HC_v3/", width = 8, height = 5, unit = "in", dpi = 600)
ggsave("T4_volcano_OHC_D109vD140_v2.tiff", p, "tiff", "./Output_d109_d140_HC_v3/", width = 8, height = 5, unit = "in", dpi = 600)



