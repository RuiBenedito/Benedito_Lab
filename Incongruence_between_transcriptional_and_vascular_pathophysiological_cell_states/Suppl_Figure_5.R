########################################
# This file contains the plots produced in Suppl Figure 5
# 
# Author:
#   Alvaro Regano aregano@cnic.es
#
########################################

# Load Libraries from the Functions.R script

# Create Plot directory

dir.create("Plots/")
dir.create("Plots/Suppl_Figure_5/")

# IMPORTANT: Before starting this script, please execute all functions in the functions.R script

# Color Palettes

my_palette_Rui_colors <- c("#E95A74", "#50B6EF", "#45FF8E", "#F4A753", "#A80519", "#880088", "#E28CF4", "#C1B80C","#FC0808",  "#0E47D8")

my_palette_Rui_colors_B <- c("#F4A753", "#E95A74", "#50B6EF", "#FC0808",  "#0E47D8", "#45FF8E", "#A80519", "#E28CF4", "#880088", "#C1B80C")

custom_rdylblu_palette <- c("#005AD9", "#0070FF", "#FFE135", "#FFDA13", "#E60026", "#C60012")

macarena_palette <- c("#FDDC85", "#6E2914")

Bestholtz_palette <- c("#DEDAD6", "#FEE392", "#FEC44E", "#FE9929", "#ED6F11", "#CC4C17", "#993411", "#65260C")

palette_Maca_Vln <- c("#606060", "#BB005E", "#F94040")

palette_Maca_Heatmap <- c( "#BB005E", "#606060", "#F94040")

BuRd <- c("#052F60", "#2066AC", "#4392C3", "#92C5DE", "#D1E4EF", "#F7F7F7", "#FCDBC6", "#F4A482", "#D55F4D", "#B1172B", "#67001F")

palette_Maca_Vln_N1 <- c("#606060", "#05BE78")

palette_Maca_GSEA <- c("#F94040", "mediumpurple1", "#BB005E")

###################################################################

Liver_all <- readRDS("")

table(Liver_all@meta.data$Condition)

###################################################################

# Suppl Figure 5A Violin plot

Liver <- subset(Liver_all, subset = Condition == "Control" | Condition == "NOTCH1KO")

table(Liver@meta.data$Condition)

gene.names <- c("Notch4")
genes <- toupper(gene.names)
genes <- as.data.frame(genes)
gene.names <- as.data.frame(gene.names)
cols = palette_Maca_Vln_N1
labels = c("Control", expression(italic(Notch1)^"iDEC"))


cairo_pdf("Plots/Suppl_Figure_5/Suppl_Figure_5A_VlnPlot.pdf", width = 6, height = 3, family = "Arial")
Violin_plot_stacked(Liver, genes, gene.names, cols = cols, labels = labels)
dev.off()

jpeg("Plots/Suppl_Figure_5/Suppl_Figure_5A_VlnPlot.jpeg", width = 6, height = 3, units = 'in', res = 800)
Violin_plot_stacked(Liver_ctl, genes, gene.names, cols = cols, labels = labels)
dev.off()

###################################################################

# Suppl Figure 5H UMAP

Liver <- subset(Liver_all, subset = Condition == "Control" | Condition == "Dll4KO(4d)+Vehicle" | Condition == "Dll4KO")

Idents(Liver) <- "Condition"
levels(Liver) <- c("Control", "Dll4KO(4d)+Vehicle", "Dll4KO")
Liver@meta.data$Condition <- Liver@active.ident

new_labels <- c("Control" = "Control", "Dll4KO(4d)+Vehicle" = "italic(Dll4)^iDEC4d", "Dll4KO" = "italic(Dll4)^iDEC2w")
cols = my_palette_Rui_colors_B

cairo_pdf("Plots/Suppl_Figure_5/Suppl_Figure_5H_UMAP.pdf",  width = 22, height = 6, family = "Arial")
UMAP_Conditions_Global(Liver, 3, new_labels, cols)
dev.off()

jpeg("Plots/Suppl_Figure_5/Suppl_Figure_5H_UMAP.jpeg", width = 22, height = 6, units = 'in', res = 800)
UMAP_Conditions_Global(Liver, 3, new_labels, cols)
dev.off()

# Barplot

p1 <- dittoBarPlot(Liver, var = "FigClustering", group.by = "Condition", color.panel = my_palette_Rui_colors_B, xlab = NULL, main = NULL, y.breaks = c(0,0.25,0.5, 0.75, 1), min = 0, max = 1)+
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 18), legend.text = element_text(size = 16), axis.text.x = element_text(angle = 25))+
  scale_x_discrete(labels=c("Control", bquote(italic(Dll4)^iDEC), bquote(italic(Dll4)^iDEC+aVEGF)))+
  geom_col(width = 0.1)+
  # scale_y_discrete(labels=c("0", "25", "50", "75", "100"))
  scale_y_continuous(labels = c("0", "25", "50", "75", "100"))

cairo_pdf("Plots/Figure_5/Figure_5H_BarPlot.pdf",  width = 6, height = 6, family = "Arial")
p1
dev.off()


jpeg("Plots/Figure_5/Figure_5H_BarPlot.jpeg", width = 6, height = 6, units = 'in', res = 800)
p1
dev.off()


###################################################################

# Suppl Figure 5I Violin Plot

gene.names <- c("Kcne3", "Esm1", "Vegfa", "Apln",  "Odc1", "Stmn1")
genes <- toupper(gene.names)
genes <- as.data.frame(genes)
gene.names <- as.data.frame(gene.names)
cols = palette_Maca_Vln
labels = c("Control", expression(italic("Dll4"^"iDEC 4d")), expression(italic("Dll4"^"iDEC 2w")))


cairo_pdf("Plots/Suppl_Figure_5/Suppl_Figure_5I_VlnPlot.pdf", width = 6, height = 12, family = "Arial")
Violin_plot_stacked(Liver, genes, gene.names, cols = cols, labels = labels)
dev.off()

jpeg("Plots/Suppl_Figure_5/Suppl_Figure_5I_VlnPlot.jpeg", width = 6, height = 12, units = 'in', res = 800)
Violin_plot_stacked(Liver_ctl, genes, gene.names, cols = cols, labels = labels)
dev.off()

