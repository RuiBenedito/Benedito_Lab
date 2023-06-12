########################################
# This file contains the plots produced in the Supplementary Figure 6
# 
# Author:
#   Alvaro Regano aregano@cnic.es
#
########################################

# Create Plot directory

dir.create("Plots/")
dir.create("Plots/Suppl_Figure_6/")

# IMPORTANT: Before starting this script, please execute all functions in the functions.R script

# rds (Produced from the Liver_non_ECs_preprocessing.R script)

Liver_Non_ECs_All_Cells <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Sample_Liver_Heart_20Jan22/rds/Liver_Non_ECs_All_Cells.rds")

# Color Palettes

Conditions_palette <- c("#606060", "#F94040")
Conditions_palette_clusters <- rep(Conditions_palette, 12)
custom_rdylblu_palette <- c("#005AD9", "#0070FF", "#FFE135", "#FFDA13", "#E60026", "#C60012")
BuRd <- c("#052F60", "#2066AC", "#4392C3", "#92C5DE", "#D1E4EF", "#F7F7F7", "#FCDBC6", "#F4A482", "#D55F4D", "#B1172B", "#67001F")
Bestholtz_palette <- c("#DEDAD6", "#FEE392", "#FEC44E", "#FE9929", "#ED6F11", "#CC4C17", "#993411", "#65260C")

# Seurat palette Final Clustering

# Load the "scales" package
require(scales)

# Create vector with levels of object@ident
identities <- levels(Liver_Non_ECs_All_Cells@meta.data$FinalClustering)

# Create vector of default ggplot2 colors
Final_Clustering_Seurat_palette <- hue_pal()(length(identities))

Condition_Clustering_Seurat_palette <- rep(Final_Clustering_Seurat_palette, each = 2)

######################################

# Suppl. Fig. 6B Heatmap

BenMoshe_markers <- c("Cd79b","Fcmr","Pax5","Ebf1","Ms4a1","Cd79a","Il2rb","Nkg7","Cd3d","Cd3g","Gzma","Smim5",
                      "Rpgrip1","Runx2","Ccr9","Siglech","Cox6a2","Gm21762","Klk1","H2-Ab1","Fscn1","Wdfy4","Tmem123",
                      "Tbc1d4","Ppt1","Cst3","Naaa","Chil3","S100a4","Ccr2","S100a6","S100a9","S100a8","Lpl","Pf4",
                      "Trem2","Fabp5","Gpnmb","Il18bp","Slc40a1","Folr2","Cd5l","Timd4","Marco","Clec4f","Vsig4","Hand2",
                      "Col6a1","Postn","Col5a1","Col3a1","Dcn","Col1a2","Col1a1","Clu","Fmo2","Tstd1","Sox9","Epcam",
                      "Pdzk1ip1","Pkhd1","Scara3","Fgb","Mup3","Slc10a1","Apoa2","Mup20","Apoa1","Apoa4","Hpx","Lyve1",
                      "Flt1","Pcdh17","Oit3","Jam2","Mmrn2","Aqp1","Ptprb")


BenMoshe_markers <- rev(BenMoshe_markers)

cairo_pdf("Plots/Suppl_Figure_6/Suppl_Figure_6B_Heatmap.pdf",  width = 18, height = 12, family = "Arial")
DoHeatmap(Liver_Non_ECs_All_Cells, features = BenMoshe_markers, size = 6, group.by = "FinalClustering", angle = 60)
dev.off()

######################################

# Suppl. Fig. 6C UMAP and Barplot

# UMAP Cell Types

Idents(Liver_Non_ECs_All_Cells) <- "FinalClustering"

p1 <- DimPlot(Liver_Non_ECs_All_Cells, pt.size = 1.2, 
              cols = Final_Clustering_Seurat_palette,
              split.by = "Condition", shuffle = F)

labels <- c("Control" = "Control", "Dll4 LOF" = "italic(Dll4)^iDEC")

p21 <- p1+facet_grid(~Condition, labeller = as_labeller(labels, default = label_parsed))+
  theme(strip.text.x = element_text(size = 24), plot.title = element_blank())+xlab("UMAP_1")+ylab("UMAP_2")

cairo_pdf("Plots/Suppl_Figure_6/Suppl_Figure_6C_UMAP_Cell_Types.pdf",  width = 12, height = 6, family = "Arial")
p21
dev.off()

#############################################

# BarPlot Cell Types

table(Liver_Non_ECs_All_Cells@meta.data$Condition)

p1 <- dittoBarPlot(Liver_Non_ECs_All_Cells, var = "FinalClustering", group.by = "Condition", 
                   main = "Liver non ECs All Cells", color.panel = Final_Clustering_Seurat_palette,
                   var.labels.reorder = c(4,9,8,1,11,2,10,3,5,6,7),
                   scale = "count", 
                   xlab = NULL)+
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 18), legend.text = element_text(size = 16), axis.text.x = element_text(angle = 25), plot.title = element_text(hjust = 0.5, face = "bold", size = 14))+
  scale_x_discrete(labels=c("Control", bquote(italic(Dll4)^iDEC)))+
  geom_col(width = 0.1)

cairo_pdf("Plots/Plots/Suppl_Figure_6/Suppl_Figure_6C_BarPlot_Celltypes.pdf",  width = 4, height = 6, family = "Arial")
p1
dev.off()

#############################################

# UMAP Conditions

new_labels <- c("Control" = "Control", "Dll4 LOF" = "italic(Dll4)^iDEC")
cols = Conditions_palette

cairo_pdf("Plots/Suppl_Figure_6/Suppl_Figure_6C_UMAP_Conditions.pdf",  width = 18, height = 6, family = "Arial")
UMAP_Conditions_Global(Liver_Non_ECs_All_Cells,2, new_labels,cols=cols, clustering="Condition")
dev.off()

######################################

# Suppl. Fig. 6D Violin Plots

Idents(Liver_Non_ECs_All_Cells) <- "FinalClustering"

levels(Liver_Non_ECs_All_Cells) <- c("Hep", "Chol", "HSC", "Macrophages", "KC", "Monocytes", "cDC", "pDC", "T+NK", "B", "Granulocytes")

genes <- c("Notch1", "Notch2", "Notch3", "Notch4", "Dll4", "Jag1", "Jag2", "Hes1", "Hes5", "Hey1", "Hey2")

Liver_Non_ECs_All_Cells@active.ident -> Liver_Non_ECs_All_Cells@meta.data$FinalClustering

genes <- as.data.frame(genes)

#  Script

myplots <- vector('list', nrow(genes))


for (i in 1:(nrow(genes)-1)) {
  
  p21 <- VlnPlot(Liver_Non_ECs_All_Cells, features = genes[i, 1], group.by = "FinalClustering", split.by = "Condition",
                 # idents = c("Hep", "KC", "HSC"),
                 cols = Conditions_palette)+ NoLegend() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), plot.title = element_text(face = "bold.italic", size = 16))
  
  myplots[[i]] <- local({
    i <- i
    p21
  })
  
}

p2 <- VlnPlot(Liver_Non_ECs_All_Cells, features = genes[nrow(genes), 1], group.by = "FinalClustering", split.by = "Condition",
              # idents = c("Hep", "KC", "HSC"),
              cols = Conditions_palette)+ NoLegend() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title = element_text(face = "bold.italic", size = 16))

myplots[[nrow(genes)]] <- local({
  p2
})


plegend <- VlnPlot(Liver_Non_ECs_All_Cells, features = genes[i, 1], group.by = "FinalClustering",
                   # idents = c("Hep", "KC", "HSC"),
                   split.by = "Condition", cols = Conditions_palette) +
  scale_fill_manual(values= Conditions_palette, labels=c("Control", expression(italic("Dll4"^"iDEC"))))

legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", legend.text = element_text(size = 18)))

p25 <- patchwork::wrap_plots(myplots, ncol = 1)

p26 <- plot_grid(p25, legend, nrow = 2, rel_heights = c(1, .03))

x.axis <- ggdraw() + draw_label("Expression Level", fontface='bold', angle = 90, size = 24)

p27 <- plot_grid(x.axis, p26, ncol = 2, rel_widths = c(.05, 1))

cairo_pdf("Plots/Suppl_Figure_6/Suppl_Figure_6D_VlnPlot.pdf", width = 8, height = 22, family = "Arial")
p27
dev.off()

######################################

# Suppl. Fig. 6E Dot Plot

genes <- c("Notch1", "Notch2", "Notch3", "Notch4", "Dll4", "Jag1", "Jag2", "Hes1", "Hes5", "Hey1", "Hey2")

genesR <- rev(genes)

Idents(Liver_Non_ECs_All_Cells) <- "Condition_Clustering"

levels(Liver_Non_ECs_All_Cells) <- c("Control_Hep", "Dll4 LOF_Hep", "Control_Chol", "Dll4 LOF_Chol",
                                     "Control_HSC", "Dll4 LOF_HSC", "Control_Macrophages", "Dll4 LOF_Macrophages", 
                                     "Control_KC", "Dll4 LOF_KC", "Control_Monocytes", "Dll4 LOF_Monocytes",
                                     "Control_cDC", "Dll4 LOF_cDC", "Control_pDC", "Dll4 LOF_pDC",
                                     "Control_T+NK", "Dll4 LOF_T+NK", "Control_B", "Dll4 LOF_B",
                                     "Control_Granulocytes", "Dll4 LOF_Granulocytes")


levels(Liver_Non_ECs_All_Cells)

Liver_Non_ECs_All_Cells@active.ident -> Liver_Non_ECs_All_Cells@meta.data$Clustering_Condition

Conditions <- rep(c("Control", "Dll4 LOF"), 11)

FinalClustering <- c("Hep","Chol","HSC","Macrophages","KC","Monocytes","cDC","pDC","T+NK","B","Granulocytes")

levels(Liver_Non_ECs_All_Cells) <- FinalClustering

Liver_Non_ECs_All_Cells@active.ident -> Liver_Non_ECs_All_Cells@meta.data$FinalClustering

Liver_Non_ECs_All_Cells@meta.data$FinalClustering <- droplevels(Liver_Non_ECs_All_Cells@meta.data$FinalClustering)

# Notch_targets

df <- data.frame(
  x = rep(c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5, 13.5, 15.5, 17.5, 19.5, 21.5)),
  y = rep(c(12), 1),
  cell_type = FinalClustering)

df$cell_type <- factor(df$cell_type, levels = c("Hep","Chol","HSC","Macrophages","KC","Monocytes","cDC","pDC","T+NK","B","Granulocytes"))

Idents(Liver_Non_ECs_All_Cells) <- "Condition_Clustering"

cairo_pdf("Plots/Suppl_Figure_6/Suppl_Figure_6E_Dotplot.pdf",
          width = 6, height = 3, family = "Arial")
DotPlot(Liver_Non_ECs_All_Cells, features = genes, col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 5), legend.text  = element_text(size = 4),
        legend.key.size = unit(0.125, "cm")) +
  scale_y_discrete(label = Conditions)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_raster(df, mapping = aes(y, x, fill = cell_type), inherit.aes = F)+
  scale_fill_manual(values = Final_Clustering_Seurat_palette)+
  geom_hline( 
    yintercept = c(2.5,4.5, 6.5,8.5, 10.5, 12.5, 14.5, 16.5, 18.5, 20.5, 22.5),
    linetype = 2 )+
  coord_flip()
dev.off()

######################################

# Suppl. Fig. 6F Heatmaps

# Full Transcriptome Heatmaps per cell types

# Separated by Celltype

Liver_Non_ECs_All_Cells@meta.data$Condition_Clustering <- as.factor(Liver_Non_ECs_All_Cells@meta.data$Condition_Clustering)

Liver_Non_ECs_All_Cells@meta.data$Condition_Clustering <- droplevels(Liver_Non_ECs_All_Cells@meta.data$Condition_Clustering)

Liver_Non_ECs_All_Cells@meta.data$FinalClustering <- droplevels(Liver_Non_ECs_All_Cells@meta.data$FinalClustering)

Idents(Liver_Non_ECs_All_Cells) <- "Condition_Clustering"

levels(Liver_Non_ECs_All_Cells) <- c("Control_Hep",  "Dll4 LOF_Hep", "Control_Chol",  "Dll4 LOF_Chol",
                                     "Control_HSC",  "Dll4 LOF_HSC", "Control_Macrophages",  "Dll4 LOF_Macrophages", "Control_KC",  "Dll4 LOF_KC",
                                     "Control_Monocytes",  "Dll4 LOF_Monocytes", "Control_cDC",  "Dll4 LOF_cDC", "Control_pDC",  "Dll4 LOF_pDC", 
                                     "Control_T+NK",  "Dll4 LOF_T+NK",  "Control_B",  "Dll4 LOF_B", "Control_Granulocytes",  "Dll4 LOF_Granulocytes",
                                     "Control_Endo",  "Dll4 LOF_Endo")

Liver_Non_ECs_All_Cells@meta.data$Condition_Clustering <- Liver_Non_ECs_All_Cells@active.ident 


cell_types <- levels(Liver_Non_ECs_All_Cells@meta.data$FinalClustering)
cell_types_df <- as.data.frame(cell_types)



myplots <- vector('list', nrow(cell_types_df))

Liver_Non_ECs_All_Cells <- ScaleData(Liver_Non_ECs_All_Cells)

Idents(Liver_Non_ECs_All_Cells) <- "FinalClustering"

for (i in 1:nrow(cell_types_df)) {
  
  selected_cells <- WhichCells(Liver_Non_ECs_All_Cells, idents = cell_types[i])  
  p24 <- DoHeatmap(Liver_Non_ECs_All_Cells, cells = selected_cells, group.by = "Condition_Clustering", label = F, group.colors = Conditions_palette_clusters)+NoLegend()+ggtitle(cell_types[i]," All genes CtlvDll4")
  
  myplots[[i]] <- local({
    i <- i
    p24
  })
  
  
}

p25 <- patchwork::wrap_plots(myplots, nrow = 1)

cairo_pdf("Plots/Suppl_Figure_6/Suppl_Figure_6F_Heatmaps.pdf", width = 90, height = 20, family = "Arial")
p25
dev.off()

######################################

# Suppl. Fig. 6J GSEA Heatmaps


# Trying GSEA with a for loop

table(Liver_Non_ECs_All_Cells@meta.data$Condition)

m_df<- msigdbr(species = "Mus musculus", category = 13)

h_df<- msigdbr(species = "Homo sapiens", category = 7)

m_gene_sets = msigdbr(species = "mouse", category = "H")

head(m_df)

fgsea_sets<- m_gene_sets %>% split(x = .$gene_symbol, f = .$gs_name)

fgsea_sets$HALLMARK_MYC_TARGETS_V1

condition_text <- c("Dll4vCtl")

clusters <- levels(Liver_Non_ECs_All_Cells@meta.data$FinalClustering)
clusters <- clusters[1:10]
clusters
clusters_text <- gsub('\\[', '\\(',
                      gsub('\\]', '\\)', clusters))

mywilcoxauc_DEG <- vector('list', length(clusters))
mywilcoxauc_DEG_names <- vector('list', length(clusters))


for (i in 1:length(clusters)) {
  
  Liver_Non_ECs_All_Cells_subset <- subset(Liver_Non_ECs_All_Cells, subset = Condition_Clustering == paste("Control", clusters[i], sep="_") | Condition_Clustering == paste("Dll4 LOF", clusters[i], sep="_"))
  
  
  liver_non_ecs.genes <- wilcoxauc(Liver_Non_ECs_All_Cells_subset, 'Condition_Clustering')
  
  mywilcoxauc_DEG[[i]] <- local({
    i <- i
    liver_non_ecs.genes
  })
  mywilcoxauc_DEG_names[[i]] <- local({
    i <- i
    paste(condition_text, clusters_text[i], sep = "_")
  })
}
names(mywilcoxauc_DEG) <- mywilcoxauc_DEG_names 

# All cells

write.xlsx(mywilcoxauc_DEG, "./Tables/GSEA_wilcoxauc_Liver_non_ECs_allcells_cluster&condition.xlsx", rowNames = T)

mywilcoxauc_DEG_stacked <- list.stack(mywilcoxauc_DEG)

dplyr::count(mywilcoxauc_DEG_stacked, group)

clusters <- levels(Liver_Non_ECs_All_Cells@meta.data$FinalClustering)
clusters <- clusters[1:10]
clusters
clusters_text <- gsub('\\[', '\\(',
                      gsub('\\]', '\\)', clusters))

myfgsea <- vector('list', length(clusters))
myfgsea_names <- vector('list', length(clusters))

for (i in 1:length(clusters)) {
  
  cluster.genes<- mywilcoxauc_DEG_stacked %>%
    dplyr::filter(group == paste("Dll4 LOF", clusters[i], sep="_")) %>%
    arrange(desc(auc)) %>% 
    dplyr::select(feature, auc)
  
  ranks<- deframe(cluster.genes)
  
  myfgsea[[i]] <- local({
    i <- i
    fgseaRes<- fgsea(fgsea_sets, stats = ranks)
  })
  myfgsea_names[[i]] <- local({
    i <- i
    paste(condition_text, clusters_text[i], sep = "_")
  })
  
}

names(myfgsea) <- myfgsea_names 

myfgsea[1]

write.xlsx(myfgsea, file = "./Tables/fgsea_Liver_non_ECs_per_cluster_and_Condition.xlsx")


# Curate Table Manually. It can be accessed in the github repository as well

cond.fgsea <- read.xlsx("./Tables/fgsea_Liver_non_ECs_per_cluster_and_Condition.xlsx", sheet = 12, sep.names = " ", rowNames = T)

liver_non_ecs.gs <- subset(Liver_Non_ECs_All_Cells, subset = Condition == "Dll4 LOF")

Idents(liver_non_ecs.gs) <- "FinalClustering"

avgexp <- AverageExpression(liver_non_ecs.gs, return.seurat = T)

avgexp@meta.data$FinalClustering <- avgexp@active.ident

avgexp <- Seurat::AddMetaData(avgexp, cond.fgsea)

avgexp@meta.data$FinalClustering <- avgexp@active.ident

levels(avgexp)

avgexp@meta.data$FinalClustering <- as.factor(avgexp@meta.data$FinalClustering)

names(cond.fgsea)

BuRd <- c("#052F60", "#2066AC", "#4392C3", "#92C5DE", "#D1E4EF", "#F7F7F7", "#F7F7F7", "#F4A482", "#D55F4D", "#B1172B", "#67001F")

cairo_pdf("Plots/Suppl_Figure_6/Suppl_Figure_6J_GSEA_Heatmaps.pdf", width = 8, height = 12, family = "Arial")
dittoHeatmap(avgexp, genes = NULL, metas = names(cond.fgsea), 
             annot.by = "FinalClustering",
             annot.colors = Final_Clustering_Seurat_palette,
             fontsize = 7, 
             cluster_cols = F,
             cluster_rows = F,
             heatmap.colors = BuRd,
             scaled.to.max = F,
             scale = "none",
             complex = T, data.out = F)
dev.off()
