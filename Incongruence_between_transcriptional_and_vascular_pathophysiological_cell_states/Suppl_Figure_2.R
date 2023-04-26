########################################
# This file contains the plots produced in Suppl Figure 2
# 
# Author:
#   Alvaro Regano aregano@cnic.es
#
########################################

# Load Libraries from the Functions.R script

# Create Plot directory

dir.create("Plots/")
dir.create("Plots/Suppl_Figure_2/")

# IMPORTANT: Before starting this script, please execute all functions in the functions.R script

# rds

Liver_all <- readRDS("")

# Color Palettes

my_palette_Rui_colors <- c("#E95A74", "#50B6EF", "#45FF8E", "#F4A753", "#A80519", "#880088", "#E28CF4", "#C1B80C","#FC0808",  "#0E47D8")

my_palette_Rui_colors_B <- c("#F4A753", "#E95A74", "#50B6EF", "#FC0808",  "#0E47D8", "#45FF8E", "#A80519", "#E28CF4", "#880088", "#C1B80C")

custom_rdylblu_palette <- c("#005AD9", "#0070FF", "#FFE135", "#FFDA13", "#E60026", "#C60012")

macarena_palette <- c("#FDDC85", "#6E2914")

Bestholtz_palette <- c("#DEDAD6", "#FEE392", "#FEC44E", "#FE9929", "#ED6F11", "#CC4C17", "#993411", "#65260C")

palette_Maca_Vln <- c("#606060", "#F94040", "#BB005E")

palette_Maca_Heatmap <- c( "#BB005E", "#606060", "#F94040")

BuRd <- c("#052F60", "#2066AC", "#4392C3", "#92C5DE", "#D1E4EF", "#F7F7F7", "#FCDBC6", "#F4A482", "#D55F4D", "#B1172B", "#67001F")


palette_Maca_GSEA <- c("#F94040", "mediumpurple1", "#BB005E")

###################################################################

# Suppl. Fig. 2N Violin Plot

# orden : C2a  Arteries, C1a arterial capillaries, C0 unspecified quiescent capillaries,  C1v venous capillaries, C2v Veins (Central veins). 
# Only Ctl

Liver_ctl <- subset(Liver_all, subset = Condition == "Control")

table(Liver_ctl@meta.data$FigClustering)

Liver_ctl <- subset(Liver_ctl, subset = FigClustering == "C0 - Unspecified quiescent capillaries" |
                      FigClustering == "C1a - Arterial capillaries" |
                      FigClustering == "C1v - Venous capillaries" |
                      FigClustering == "C2a - Large arteries" |
                      FigClustering == "C2v - Large veins"
)

Idents(Liver_ctl) <- "FigClustering"

genes <- c("DLL4", "DLL1", "JAG1", "JAG2", "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "MFNG", "LFNG", "RFNG", "RBPJ", "HEY1", "HEY2", "HES1", "HEYL")
gene.names <- str_to_title(genes)
genes <- as.data.frame(genes)
gene.names <- as.data.frame(gene.names)
cols = my_palette_Rui_colors_B
Liver_ctl@meta.data$FigClustering <- as.factor(Liver_ctl@meta.data$FigClustering)
labels= levels(Liver_ctl@meta.data$FigClustering)

cairo_pdf("Plots/Suppl_Figure_2/Suppl_Figure_2N_VlnPlot.pdf", width = 6, height = 32, family = "Arial")
Violin_plot_stacked(Liver_ctl, genes, gene.names, cols = cols, labels = labels, legend.direction = "vertical")
dev.off()

jpeg("Plots/Suppl_Figure_2/Suppl_Figure_2N_VlnPlot.jpeg", width = 6, height = 32, units = 'in', res = 800)
Violin_plot_stacked(Liver_ctl, genes, gene.names, cols = cols, labels = labels, legend.direction = "vertical")
dev.off()

###################################################################

# Suppl. Fig. 2O Dot Plot

genes <- c("DLL4", "DLL1", "JAG1", "JAG2", "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "MFNG", "LFNG", "RFNG", "RBPJ", "HEY1", "HEY2", "HES1", "HEYL")

genes <- rev(genes)

gene.names <- str_to_title(genes)
  
cairo_pdf("Plots/Suppl_Figure_2/Suppl_Figure_2O_DotPlot.pdf", width = 5, height = 7, family = "Arial")
DotPlot(Liver_ctl, features = Notch_genes, col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label = gene.names)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  coord_flip()
dev.off()

###################################################################

# Suppl. Fig. 2O Dot Plot
