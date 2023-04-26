########################################
# This file contains the plots produced in Figure 5
# 
# Author:
#   Alvaro Regano aregano@cnic.es
#
########################################

# Libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library("dittoSeq")
library("magick")
library("RColorBrewer")
library(cowplot)
library(dplyr)
library(extrafont)
loadfonts(device = "win")
library(openxlsx)
library(ComplexHeatmap)
library("stringr")
library("formattable")
library(readxl)
library(tidyverse)
library(showtext)
library(fgsea)

# Create Plot directory

dir.create("Plots/")
dir.create("Plots/Figure_5/")

# IMPORTANT: Before starting this script, please execute all functions in the functions.R script

# rds

Liver_all <- readRDS("")

LiverC <- subset(Liver_all, subset = Condition == "Control")

Liver <- subset(Liver_all, subset = Condition == "Control" | Condition == "Jag1/Jag2/Dll1KO(2w)")

# Color Palettes

my_palette_Rui_colors <- c("#E95A74", "#50B6EF", "#45FF8E", "#F4A753", "#A80519", "#880088", "#E28CF4", "#C1B80C","#FC0808",  "#0E47D8")

my_palette_Rui_colors_B <- c("#F4A753", "#E95A74", "#50B6EF", "#FC0808",  "#0E47D8", "#45FF8E", "#A80519", "#E28CF4", "#880088", "#C1B80C")

custom_rdylblu_palette <- c("#005AD9", "#0070FF", "#FFE135", "#FFDA13", "#E60026", "#C60012")

macarena_palette <- c("#FDDC85", "#6E2914")

Bestholtz_palette <- c("#DEDAD6", "#FEE392", "#FEC44E", "#FE9929", "#ED6F11", "#CC4C17", "#993411", "#65260C")

palette_Maca_Vln <- c("#606060", "#DBB49A")

palette_Maca_Vln_gene <- c("#A00000", "#FCB4B4", "#F2F0F0", "#F2F0F0")

BuRd <- c("#052F60", "#2066AC", "#4392C3", "#92C5DE", "#D1E4EF", "#F7F7F7", "#FCDBC6", "#F4A482", "#D55F4D", "#B1172B", "#67001F")

###############################################################################################################################

## add the Arial font
font_add("Arial", regular = "arial.ttf",
         bold = "arialbd.ttf", italic = "ariali.ttf", bolditalic = "arialbi.ttf")


###############################################################################################################################


# Figure 5A. Violin Plot

genes <- c("DLL4", "DLL1", "JAG1", "JAG2")
gene.names <- str_to_title(genes)
genes <- as.data.frame(genes)
gene.names <- as.data.frame(gene.names)
cols = palette_Maca_Vln
labels=c("Control")

cairo_pdf("Plots/Figure_5/Figure_5A_VlnPlot.pdf", width = 3, height = 12, family = "Arial")
Violin_plot_stacked(LiverC, genes, gene.names, cols = cols, labels = labels)
dev.off()

jpeg("Plots/Figure_4/Figure_4H_VlnPlot.jpeg", width = 6, height = 6, units = 'in', res = 800)
Violin_plot_stacked(LiverC, genes, gene.names, cols = cols, labels = labels)
dev.off()

###############################################################################################################################


# Figure 5H. Violin Plot

genes <- c("DLL4", "DLL1", "JAG1", "JAG2")
gene.names <- str_to_title(genes)
genes <- as.data.frame(genes)
gene.names <- as.data.frame(gene.names)
cols = palette_Maca_Vln
labels=c("Control", expression(italic(Jag1/2)^"iDEC"))

cairo_pdf("Plots/Figure_5/Figure_5H_VlnPlot.pdf", width = 5, height = 12, family = "Arial")
Violin_plot_stacked(Liver, genes, gene.names, cols = cols, labels = labels)
dev.off()

jpeg("Plots/Figure_5/Figure_5H_VlnPlot.jpeg", width = 5, height = 6, units = 'in', res = 800)
Violin_plot_stacked(Liver, genes, gene.names, cols = cols, labels = labels)
dev.off()

###############################################################################################################################

# Figure 5K. UMAP

new_labels <- c("Control" = "Control", "Jag1/Jag2/Dll1KO(2w)" = "italic(Jag1/2)^iDEC")
cols = my_palette_Rui_colors_B

cairo_pdf("Plots/Figure_5/Figure_5K_UMAP.pdf",  width = 22, height = 6, family = "Arial")
UMAP_Conditions_Global(Liver, 2, new_labels, cols)
dev.off()

jpeg("Plots/Figure_5/Figure_5K_UMAP.jpeg", width = 22, height = 6, units = 'in', res = 800)
UMAP_Conditions_Global(Liver, 2, new_labels, cols)
dev.off()

# Barplot

p1 <- dittoBarPlot(Liver, var = "FigClustering", group.by = "Condition", color.panel = my_palette_Rui_colors_B, xlab = NULL, main = NULL, y.breaks = c(0,0.25,0.5, 0.75, 1), min = 0, max = 1)+
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 18), legend.text = element_text(size = 16), axis.text.x = element_text(angle = 25))+
  scale_x_discrete(labels=c("Control", bquote(italic(Jag1/2)^iDEC)))+
  geom_col(width = 0.1)+
  # scale_y_discrete(labels=c("0", "25", "50", "75", "100"))
  scale_y_continuous(labels = c("0", "25", "50", "75", "100"))
p1

cairo_pdf("Plots/Figure_5/Figure_5G2_BarPlot.pdf",  width = 6, height = 6, family = "Arial")
p1
dev.off()


jpeg("Plots/Figure_5/Figure_5G2_BarPlot.jpeg", width = 6, height = 6, units = 'in', res = 800)
p1
dev.off()

###############################################################################################################################


# Figure 5L. Violin Plot

gene.names <- c("Esm1", "Kcne3", "Angpt2", "Stmn1", "Myc", "Odc1")
genes <- toupper(gene.names)
genes <- as.data.frame(genes)
gene.names <- as.data.frame(gene.names)
Liver <- subset(Liver_all, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "Jag1/Jag2/Dll1KO(2w)")
Idents(Liver) <- "Condition"
levels(Liver) <- c("Control", "Jag1/Jag2/Dll1KO(2w)", "Dll4KO")
Liver@active.ident -> Liver@meta.data$Condition

palette_Maca_Vln <- c("#606060", "#DBB49A", "#F94040")
cols = palette_Maca_Vln
labels=c("Control", expression(italic(Jag1/2)^"iDEC"), expression(italic(Dll4)^"iDEC"))

cairo_pdf("Plots/Figure_5/Figure_5L_VlnPlot.pdf", width = 6, height = 18, family = "Arial")
Violin_plot_stacked(Liver, genes, gene.names, cols = cols, labels = labels)
dev.off()

jpeg("Plots/Figure_5/Figure_5L_VlnPlot.jpeg", width = 6, height = 18, units = 'in', res = 800)
Violin_plot_stacked(Liver, genes, gene.names, cols = cols, labels = labels)
dev.off()

###############################################################################################################################

# Figure 5M. Violin Plot

gene.names <- c("Msr1", "Efnb2")
genes <- toupper(gene.names)
genes <- as.data.frame(genes)
gene.names <- as.data.frame(gene.names)
Liver <- subset(Liver_all, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "Jag1/Jag2/Dll1KO(2w)")
Idents(Liver) <- "Condition"
levels(Liver) <- c("Control", "Jag1/Jag2/Dll1KO(2w)", "Dll4KO")
Liver@active.ident -> Liver@meta.data$Condition

palette_Maca_Vln <- c("#606060", "#DBB49A", "#F94040")
cols = palette_Maca_Vln
labels=c("Control", expression(italic(Jag1/2)^"iDEC"), expression(italic(Dll4)^"iDEC"))

cairo_pdf("Plots/Figure_5/Figure_5M_VlnPlot.pdf", width = 6, height = 6, family = "Arial")
Violin_plot_stacked(Liver, genes, gene.names, cols = cols, labels = labels)
dev.off()

jpeg("Plots/Figure_5/Figure_5M_VlnPlot.jpeg", width = 6, height = 6, units = 'in', res = 800)
Violin_plot_stacked(Liver, genes, gene.names, cols = cols, labels = labels)
dev.off()

###############################################################################################################################

# Figure 5N. Violin Plot

gene.names <- c("Hes1", "Cd34", "Wnt2")
genes <- toupper(gene.names)
genes <- as.data.frame(genes)
gene.names <- as.data.frame(gene.names)
Liver <- subset(Liver_all, subset = Condition == "Control" | Condition == "Jag1/Jag2/Dll1KO(2w)")
Idents(Liver) <- "Condition"
levels(Liver) <- c("Control", "Jag1/Jag2/Dll1KO(2w)")
Liver@active.ident -> Liver@meta.data$Condition

palette_Maca_Vln <- c("#606060", "#DBB49A")
cols = palette_Maca_Vln
labels=c("Control", expression(italic(Jag1/2)^"iDEC"))

cairo_pdf("Plots/Figure_5/Figure_5N_VlnPlot.pdf", width = 5, height = 9, family = "Arial")
Violin_plot_stacked(Liver, genes, gene.names, cols = cols, labels = labels)
dev.off()

jpeg("Plots/Figure_5/Figure_5N_VlnPlot.jpeg", width = 5, height = 9, units = 'in', res = 800)
Violin_plot_stacked(Liver, genes, gene.names, cols = cols, labels = labels)
dev.off()
