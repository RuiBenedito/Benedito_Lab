########################################
# This file contains the plots produced in Figure 6
# 
# Author:
#   Alvaro Regano aregano@cnic.es
#
########################################

# Load Libraries from the Functions.R script

# Create Plot directory

dir.create("Plots/")
dir.create("Plots/Figure_6/")

# IMPORTANT: Before starting this script, please execute all functions in the functions.R script

# rds

Liver_all <- readRDS("")

Liver <- subset(Liver_all, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "Dll4/MycKO")

LiverDll4 <- subset(Liver, Condition == "Dll4KO")

Idents(Liver) <- Liver@meta.data$Condition

levels(Liver) <- c("Control", "Dll4KO", "Dll4/MycKO")

Liver@active.ident -> Liver@meta.data$Condition

table(Liver@meta.data$Condition)

# Color Palettes

my_palette_Rui_colors <- c("#E95A74", "#50B6EF", "#45FF8E", "#F4A753", "#A80519", "#880088", "#E28CF4", "#C1B80C","#FC0808",  "#0E47D8")

my_palette_Rui_colors_B <- c("#F4A753", "#E95A74", "#50B6EF", "#FC0808",  "#0E47D8", "#45FF8E", "#A80519", "#E28CF4", "#880088", "#C1B80C")

custom_rdylblu_palette <- c("#005AD9", "#0070FF", "#FFE135", "#FFDA13", "#E60026", "#C60012")

macarena_palette <- c("#FDDC85", "#6E2914")

Bestholtz_palette <- c("#DEDAD6", "#FEE392", "#FEC44E", "#FE9929", "#ED6F11", "#CC4C17", "#993411", "#65260C")

palette_Maca_Vln <- c("#606060", "#F94040", "#FF8080")

palette_Maca_GSEA <- c("#F94040", "mediumpurple1")

BuRd <- c("#052F60", "#2066AC", "#4392C3", "#92C5DE", "#D1E4EF", "#F7F7F7", "#FCDBC6", "#F4A482", "#D55F4D", "#B1172B", "#67001F")

###############################################################################################################################


## add the Arial font
font_add("Arial", regular = "arial.ttf",
         bold = "arialbd.ttf", italic = "ariali.ttf", bolditalic = "arialbi.ttf")


###############################################################################################################################

# Figure 6A GSEA

Liver <- subset(Liver_all, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "NOTCH1KO" | Condition == "RBPJKO")

liver.genes <- wilcoxauc(Liver, 'FigClustering')

dplyr::count(cluster.genes, group)

figclustering <- c("C0 - Unspecified quiescent capillaries", "C1a - Arterial capillaries", "C1v - Venous capillaries", "C2a - Large arteries", "C2v - Large veins",  "C3 - Activated capillaries", 
                   "C4 - Endothelial tip cells", "C5ip - Proliferating S-phase", "C5p - G2-M", "C6 - X")


liver.genes %>%
  dplyr::filter(group == figclustering[2]) %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 20)


myfgsea <- vector('list', 10)


for (i in 1:10) {
  
  cluster.genes<- liver.genes %>%
    dplyr::filter(group == figclustering[i]) %>%
    arrange(desc(auc)) %>% 
    dplyr::select(feature, auc)
  
  cluster.genes$feature <- tolower(cluster.genes$feature) %>% str_to_title()
  
  ranks<- deframe(cluster.genes)
  
  myfgsea[[i]] <- local({
    i <- i
    fgseaRes<- fgsea(fgsea_sets, stats = ranks)
  })
  
}

myfgsea[1]

class(myfgsea)

figclustering_short <- c("C0", "C1a", "C1v", "C2a", "C2v",  "C3", "C4", "C5ip", "C5p", "C6")

names(myfgsea) <- figclustering_short

write.xlsx(myfgsea, file = "./Tables/fgsea_Clusters_Ctrl_Dll4KO_RbpjKO_Notch1KO.xlsx")

# After Trimming and building a proper meta.data matrix in Excel. Load Table here below

cond.fgsea <- read.xlsx("Tables/fgsea_Clusters_Ctrl_Dll4KO_RbpjKO_Notch1KO.xlsx", sheet = 12, sep.names = " ", rowNames = T)

# Heatmap per cluster with AverageExpression

Idents(Liver) <- "FigClustering"

avgexp <- AverageExpression(Liver, return.seurat = T)

avgexp@meta.data$FigClustering <- avgexp@active.ident

avgexp <- Seurat::AddMetaData(avgexp, cond.fgsea)

levels(avgexp) <- c("C0 - Unspecified quiescent capillaries", "C1a - Arterial capillaries", "C1v - Venous capillaries", "C2a - Large arteries", "C2v - Large veins",  "C3 - Activated capillaries", 
                    "C4 - Endothelial tip cells", "C5ip - Proliferating S-phase", "C5p - G2-M", "C6 - X","C1")

avgexp@meta.data$FigClustering <- avgexp@active.ident

cairo_pdf("Plots/GSEA/Heatmaps/GSEA_Heatmap_Liver_Ctrl_Dll4_Notch1_Rbpj_unscaled.pdf", width = 14, height = 12, family = "Arial")
dittoHeatmap(avgexp, genes = NULL, metas = names(cond.fgsea), 
             annot.by = "FigClustering",
             annot.colors = my_palette_Rui_colors_B,
             fontsize = 7, 
             cluster_cols = F,
             heatmap.colors = BuRd,
             # heatmap.colors.max.scaled = BuRd,
             scaled.to.max = F,
             scale = "none",
             complex = F, data.out = F,
             # legend_breaks = c(1, 0.8, 0.6, 0.4, 0.2, 0),
             main = expression(Control~italic(Dll4)^iDEC~italic(Notch1)^iDEC~italic(Rbpj)^iDEC))
dev.off()

cond.fgsea.selected <- c('HALLMARK_E2F_TARGETS', 'HALLMARK_G2M_CHECKPOINT',  'HALLMARK_MITOTIC_SPINDLE', 'HALLMARK_MYC_TARGETS_V1', 'HALLMARK_MYC_TARGETS_V2',
                         'HALLMARK_PI3K_AKT_MTOR_SIGNALING', 'HALLMARK_MTORC1_SIGNALING', 'HALLMARK_GLYCOLYSIS', 'HALLMARK_FATTY_ACID_METABOLISM',
                         'HALLMARK_OXIDATIVE_PHOSPHORYLATION', 'HALLMARK_PROTEIN_SECRETION', 'HALLMARK_HYPOXIA', 'HALLMARK_P53_PATHWAY', 'HALLMARK_DNA_REPAIR',
                         'HALLMARK_APOPTOSIS', 'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY', 'HALLMARK_UNFOLDED_PROTEIN_RESPONSE', 'HALLMARK_INFLAMMATORY_RESPONSE')



cairo_pdf("Plots/GSEA/GSEA_Heatmap_Liver_Ctrl_Dll4_Notch1_Rbpj_unscaled.pdf", width = 14, height = 12, family = "Arial")
dittoHeatmap(avgexp, genes = NULL, metas = cond.fgsea.selected, 
             annot.by = "FigClustering",
             annot.colors = my_palette_Rui_colors_B,
             fontsize = 7, 
             cluster_cols = F,
             cluster_rows = F,
             heatmap.colors = BuRd,
             # heatmap.colors.max.scaled = BuRd,
             scaled.to.max = F,
             scale = "none",
             complex = T, data.out = F)
dev.off()



###############################################################################################################################

# Figure 6D BarPlot GSEA Dll4KO

# Load excel modified table for ease of plotting the BarChart

cluster.fgsea <- read.xlsx("./Tables/fgsea_Clusters_Ctrl_Dll4KO_RbpjKO_Notch1KO.xlsx", sheet = 13, sep.names = " ", rowNames = F)

cluster.fgsea[2:11] <- round(cluster.fgsea[2:11], digits = 2)

Hallmarks_quotes <- c("HALLMARK_E2F_TARGETS", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_OXIDATIVE_PHOSPHORYLATION")

#  Barplot of Gene Sets with positive and negative values stacked


cluster.fgsea.barchart <- expand_grid(
  HALMARK    = Hallmarks_quotes,  # Define all unique student names
  color = my_palette_Rui_colors_B,     # Define all unique HW assignments
  value      = NA)
cluster.fgsea.barchart

cluster.fgsea.barchart.long <- pivot_longer(cluster.fgsea, cols=c(2,3,5,6), names_to = "Hallmark", values_to = "Cluster")

cluster.fgsea.barchart <- cluster.fgsea.barchart.long[c(1, 8:10)] %>% arrange(Hallmark)

write.xlsx(cluster.fgsea.barchart, "./Tables/fgsea_Clusters_Ctrl_Dll4KO_RbpjKO_Notch1KO.Barchart.xlsx")


dat1 <- subset(cluster.fgsea.barchart.long,Cluster >= 0)
dat2 <- subset(cluster.fgsea.barchart.long,Cluster < 0)

p1 <- ggplot(cluster.fgsea.barchart.long, aes(x = Hallmark, y = Cluster, 
                                         fill = X1, label = Cluster)) +
  geom_bar(stat = "identity") + geom_text(
    size = 3, position = position_stack(vjust = 0.5),colour = "white")   +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 12), axis.title = element_blank(), legend.title = element_blank())+
  scale_fill_manual(values= cluster.fgsea.barchart$color)+
  scale_y_continuous(breaks=seq(-90,70,10))

cairo_pdf("Plots/Figure_6/Figure_6D_Barchart_HallMark_Gene_Sets_with_numbers.pdf", width = 5, height = 10, family = "Arial")
p1
dev.off()

###############################################################################################################################


# Figure 6J UMAP

Liver <- subset(Liver_all, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "Dll4/MycKO")

new_labels <- c("Control" = "Control", "Dll4KO" = "italic(Dll4)^iDEC", "Dll4/MycKO" = "italic(Dll4/Myc)^iDEC")
cols = my_palette_Rui_colors_B

cairo_pdf("Plots/Figure_6/Figure_6J_UMAP.pdf",  width = 28, height = 6, family = "Arial")
UMAP_Conditions_Global(Liver, 3, new_labels, cols)
dev.off()

jpeg("Plots/Figure_6/Figure_6J_UMAP.jpeg", width = 28, height = 6, units = 'in', res = 800)
UMAP_Conditions_Global(Liver, 3, new_labels, cols)
dev.off()

# BarPlot

table(Liver@meta.data$Condition)

p1 <- dittoBarPlot(Liver, var = "FigClustering", group.by = "Condition", x.reorder = c(1,3,2), color.panel = my_palette_Rui_colors_B, xlab = NULL, main = NULL, y.breaks = c(0,0.25,0.5, 0.75, 1), min = 0, max = 1)+
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 18), legend.text = element_text(size = 16), axis.text.x = element_text(angle = 25))+
  scale_x_discrete(labels=c("Control", bquote(italic(Dll4)^iDEC), bquote(italic(Dll4/Myc)^iDEC)))+
  geom_col(width = 0.1)+
  # scale_y_discrete(labels=c("0", "25", "50", "75", "100"))
  scale_y_continuous(labels = c("0", "25", "50", "75", "100"))
p1

cairo_pdf("Plots/Figure_6/Figure_6H.2_BarPlot.pdf",  width = 6.5, height = 6, family = "Arial")
p1
dev.off()

jpeg("Plots/Figure_6/Figure_6H.2_BarPlot.jpeg", width = 6.5, height = 6, units = 'in', res = 800)
p1
dev.off()


###############################################################################################################################


# Figure 6L GSEA Heatmap

# This function creates a Tables directory where it will produce a ranks table, a fgsea table and a wilcoxauc table

GSEA_loop_HallMark_GS(Liver, c("Dll4KO", "Dll4/MycKO"), "Homo sapiens")

# Edited GSEA Table output manually to match avgexp scRNASeq dataset (check table fgseaResD4&D4MycTidy.auc.xlsx)

cond.fgsea <- read.xlsx("./Tables/GSEA_Fig6/fgseaResD4&D4MycTidy.auc.xlsx", sheet = 3, sep.names = " ", rowNames = T)

Idents(Liver) <- "Condition"

liver.gs <- subset(Liver, subset = Condition == "Dll4KO" | Condition == "Dll4/MycKO")

avgexp <- AverageExpression(liver.gs, return.seurat = T)

avgexp <- Seurat::AddMetaData(avgexp, cond.fgsea)

avgexp@meta.data$Condition <- avgexp@active.ident

levels(avgexp) <- c("Dll4KO", "Dll4/MycKO")

avgexp@meta.data$Condition <- avgexp@active.ident

avgexp@meta.data$Condition <- as.factor(avgexp@meta.data$Condition)

avgexp@meta.data$Condition <- avgexp@active.ident

# All gene sets

cairo_pdf("Plots/Figure_6/Figure_6I_GSEA_Heatmap_Dll4vsControl_Dll4MyciDECvsControl_unscaled.pdf", width = 5, height = 12, family = "Arial")
dittoHeatmap(avgexp, genes = NULL, metas = names(cond.fgsea), 
             annot.by = "Condition",
             # order.by = c(2,1),
             annot.colors = palette_Maca_GSEA,
             fontsize = 7, 
             cluster_cols = F,
             cluster_rows = F,
             heatmap.colors = BuRd,
             # heatmap.colors.max.scaled = BuRd,
             scaled.to.max = F,
             scale = "none",
             complex = T, data.out = F)
# legend_breaks = c(1, 0.8, 0.6, 0.4, 0.2, 0),
# main = expression(Control~italic(Dll4)^iDEC~italic(Notch1)^iDEC~italic(Rbpj)^iDEC))
dev.off()

# Selected gene sets

cond.fgsea.selected <- c('HALLMARK_E2F_TARGETS', 'HALLMARK_G2M_CHECKPOINT',  'HALLMARK_MITOTIC_SPINDLE', 'HALLMARK_MYC_TARGETS_V1', 'HALLMARK_MYC_TARGETS_V2',
                         'HALLMARK_PI3K_AKT_MTOR_SIGNALING', 'HALLMARK_MTORC1_SIGNALING', 'HALLMARK_GLYCOLYSIS', 'HALLMARK_FATTY_ACID_METABOLISM',
                         'HALLMARK_OXIDATIVE_PHOSPHORYLATION', 'HALLMARK_PROTEIN_SECRETION', 'HALLMARK_HYPOXIA', 'HALLMARK_P53_PATHWAY', 'HALLMARK_DNA_REPAIR',
                         'HALLMARK_APOPTOSIS', 'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY', 'HALLMARK_UNFOLDED_PROTEIN_RESPONSE', 'HALLMARK_INFLAMMATORY_RESPONSE')

cairo_pdf("Plots/Figure_6/Figure_6I_GSEA_Heatmap_selected_Dll4vsControl_Dll4MyciDECvsControl_unscaled.pdf", width = 5, height = 6, family = "Arial")
dittoHeatmap(avgexp, genes = NULL, metas = cond.fgsea.selected, 
             annot.by = "Condition",
             order.by = c(2,1),
             annot.colors = palette_Maca_GSEA,
             fontsize = 7, 
             cluster_cols = F,
             cluster_rows = F,
             heatmap.colors = BuRd,
             # heatmap.colors.max.scaled = BuRd,
             scaled.to.max = F,
             scale = "none",
             complex = T, data.out = F
             # legend_breaks = c(1, 0.8, 0.6, 0.4, 0.2, 0),
)
dev.off()
