########################################
# This file contains the plots produced in the Supplementary Figure 8
# 
# Author:
#   Alvaro Regano aregano@cnic.es
#
########################################

# Create Plot directory

dir.create("Plots/")
dir.create("Plots/Suppl_Figure_8/")

# Load Libraries from the Functions.R script

# IMPORTANT: Before starting this script, please execute all functions in the functions.R script

# rds

Liver_all <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.OldDll4Myc.rds")

table(Liver_all@meta.data$Condition)

Liver <- subset(Liver_all, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "Dll4/MycKO")

table(Liver@meta.data$Condition)
Idents(Liver) <- "Condition"
levels(Liver) <- c("Control", "Dll4KO", "Dll4/MycKO")
Liver@active.ident -> Liver@meta.data$Condition
# Color Palettes

Conditions_palette <- c("#606060", "#F94040", "#FF8080")
custom_rdylblu_palette <- c("#005AD9", "#0070FF", "#FFE135", "#FFDA13", "#E60026", "#C60012")
BuRd <- c("#052F60", "#2066AC", "#4392C3", "#92C5DE", "#D1E4EF", "#F7F7F7", "#FCDBC6", "#F4A482", "#D55F4D", "#B1172B", "#67001F")
Bestholtz_palette <- c("#DEDAD6", "#FEE392", "#FEC44E", "#FE9929", "#ED6F11", "#CC4C17", "#993411", "#65260C")
BuRd <- c("#052F60", "#2066AC", "#4392C3", "#92C5DE", "#D1E4EF", "#F7F7F7", "#FCDBC6", "#F4A482", "#D55F4D", "#B1172B", "#67001F")

########################################################################

# Suppl. Fig. 8A Violin Plot

genes <- c("MYC", "DLL4")
gene.names <- str_to_title(genes)
genes <- as.data.frame(genes)
gene.names <- as.data.frame(gene.names)
cols = Conditions_palette


labels = c("Control", expression(italic("Dll4")^"iDEC"), expression(italic("Dll4/Myc")^"iDEC"), expression(italic(Rbpj)^"iDEC"))
Idents(Liver) <- "Condition"

cairo_pdf("Plots/Suppl_Figure_8/Suppl_Figure_8A_VlnPlot.pdf",   width = 5, height = 6, family = "Arial")
Violin_plot_stacked(Liver, genes = genes, gene.names = gene.names, cols = Conditions_palette, labels = labels)
dev.off()

########################################################################

# Suppl. Fig. 8E GOBP GSEA Heatmap


# This function creates a Tables directory where it will produce a ranks table, a fgsea table and a wilcoxauc table

GSEA_loop_GOBP_GS(Liver, c("Dll4KO", "Dll4/MycKO"), "Homo sapiens")


# After manually curating the Table in order to select the GOBP of interest


cond.fgsea <- read.xlsx("./Tables/fgsea_GOBP_selected_Dll4_Dll4Myc.xlsx", sheet = 3, sep.names = " ", rowNames = T)

liver.gs <- subset(Liver, subset = Condition == "Dll4KO" | Condition == "Dll4/MycKO")
avgexp <- AverageExpression(liver.gs, return.seurat = T)
avgexp <- Seurat::AddMetaData(avgexp, cond.fgsea)
levels(avgexp) <- c("Dll4KO", "Dll4/MycKO")

avgexp@meta.data$Condition <- avgexp@active.ident
avgexp@meta.data$Condition <- as.factor(avgexp@meta.data$Condition)

# All gene sets

cairo_pdf("Plots/Suppl_Figure_8/Suppl_Figure_8E_GSEA_Heatmap.pdf", width = 5, height = 12, family = "Arial")
dittoHeatmap(avgexp, genes = NULL, metas = names(cond.fgsea), 
             annot.by = "Condition",
             annot.colors = Conditions_palette[2:3],
             fontsize = 7, 
             cluster_cols = F,
             cluster_rows = F,
             heatmap.colors = BuRd,
             # heatmap.colors.max.scaled = BuRd,
             scaled.to.max = F,
             scale = "none",
             complex = T, data.out = F)
dev.off()


