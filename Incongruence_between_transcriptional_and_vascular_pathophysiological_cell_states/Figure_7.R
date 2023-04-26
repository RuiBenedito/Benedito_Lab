########################################
# This file contains the plots produced in Figure 7
# 
# Author:
#   Alvaro Regano aregano@cnic.es
#
########################################

# Load Libraries from the Functions.R script

# Create Plot directory

dir.create("Plots/")
dir.create("Plots/Figure_7/")
dir.create("Tables")

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


###############################################################################################################################

# Figure 7E. UMAP

Liver <- subset(Liver_all, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "D4KO_aVEGF")

levels(Liver) <- c("Control", "Dll4KO", "Dll4/MycKO", "D4KO_aVEGF")

Liver@active.ident -> Liver@meta.data$Condition

new_labels <- c("Control" = "Control", "Dll4KO" = "italic(Dll4)^iDEC", "D4KO_aVEGF" = "italic(Dll4)^iDEC+aVEGF")
cols = my_palette_Rui_colors_B

cairo_pdf("Plots/Figure_7/Figure_7E_UMAP.pdf",  width = 28, height = 6, family = "Arial")
UMAP_Conditions_Global(Liver, 3, new_labels, cols)
dev.off()

jpeg("Plots/Figure_7/Figure_7E_UMAP.jpeg", width = 28, height = 6, units = 'in', res = 800)
UMAP_Conditions_Global(Liver, 3, new_labels, cols)
dev.off()

# Barplot

p1 <- dittoBarPlot(Liver, var = "FigClustering", group.by = "Condition", color.panel = my_palette_Rui_colors_B, xlab = NULL, main = NULL, y.breaks = c(0,0.25,0.5, 0.75, 1), min = 0, max = 1)+
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 18), legend.text = element_text(size = 16), axis.text.x = element_text(angle = 25))+
  scale_x_discrete(labels=c("Control", bquote(italic(Dll4)^iDEC), bquote(italic(Dll4)^iDEC+aVEGF)))+
  geom_col(width = 0.1)+
  # scale_y_discrete(labels=c("0", "25", "50", "75", "100"))
  scale_y_continuous(labels = c("0", "25", "50", "75", "100"))

cairo_pdf("Plots/Figure_7/Figure_7E_BarPlot.pdf",  width = 6, height = 6, family = "Arial")
p1
dev.off()


jpeg("Plots/Figure_7/Figure_7E_BarPlot.jpeg", width = 6, height = 6, units = 'in', res = 800)
p1
dev.off()

###############################################################################################################################

# Figure 7F. Heatmap

Liver <- subset(Liver_all, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "D4KO_aVEGF")

Idents(Liver) <- "Condition"
levels(Liver) <- c("Control", "Dll4KO", "D4KO_aVEGF")

avgexp@active.ident -> avgexp@meta.data$Condition

avgexp <- AverageExpression(Liver, return.seurat = T)

avgexp <- ScaleData(avgexp)

DoHeatmap(avgexp, features = "ESM1")

cairo_pdf("Plots/Figure_7/Figure_7F_Heatmap.pdf", width = 8, height = 16, family = "Arial")
dittoHeatmap(avgexp,
             annot.by = "Condition", 
             order.by = c(3,2,1),
             annot.colors = palette_Maca_Heatmap, 
             treeheight_row = 0,
             treeheight_col = 0,
             fontsize = 7, 
             cluster_cols = F,
             heatmap.colors = BuRd,
             highlight.features = c("CDKN1A", "MYC", "ODC1", "APLN", "STMN1", "KCNE3", "MCM2", "ANGPT2", "ESM1", "MSR1", "LTBP4", "WNT"),
             scaled.to.max = F,
             complex = T, data.out = F)
dev.off()

###############################################################################################################################

# Figure 7G Dotplot

Liver <- subset(Liver_all, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "Dll4/MycKO" | Condition == "D4KO_aVEGF")

Idents(Liver) <- "Condition"

levels(Liver) <- c("Control", "Dll4KO", "Dll4/MycKO", "D4KO_aVEGF")

# Some targets are taken from getting the ranks of specific gene sets, you can find the tables attached or produce them

GSEA_loop_HallMark_GS(Liver, condition = c("Dll4KO", "Dll4/MycKO", "D4KO_aVEGF"), "Homo sapiens")

# Get list of all top3 genes of each Condition, 27 total

E2F_genes_Condition <- read_excel_allsheets("Tables/GSEA_Fig7/E2F_ranks.Dll4KO.Dll4MycKO.Dll4KOaVEGF.xlsx")
Myc_genes_Condition <- read_excel_allsheets("Tables/GSEA_Fig7/Myc_ranks.Dll4KO.Dll4MycKO.Dll4KOaVEGF.xlsx")
OxPhos_genes_Condition <- read_excel_allsheets("Tables/GSEA_Fig7/OxPhos_ranks.Dll4KO.Dll4MycKO.Dll4KOaVEGF.xlsx")

Top3genes_GeneSets_Conditions <- c(E2F_genes_Condition[[1]][[1]][1:3], E2F_genes_Condition[[2]][[1]][c(2:4)], E2F_genes_Condition[[3]][[1]][c(2,3,5)],
                                   Myc_genes_Condition[[1]][[1]][1:3], Myc_genes_Condition[[2]][[1]][c(2, 4, 5)], Myc_genes_Condition[[3]][[1]][c(3,8,10)],
                                   OxPhos_genes_Condition[[1]][[1]][1:3], OxPhos_genes_Condition[[2]][[1]][c(1,5,8)], OxPhos_genes_Condition[[3]][[1]][c(2,3,5)])

unique(Top3genes_GeneSets_Conditions)
Top3genes_GeneSets_Conditions

C4_Tip_Cells <- c("Vegfa", "Kcne3", "Esm1", "Apln")
Cp_G2M <- c("Top2a", "Mki67", "Hmgb2", "Cenpf")
C1a_Arterial_capillaries <- c("Ntn4", "Ltbp4", "Msr1", "Ehd3", "Efnb2")

DotPlot_list <- c(Top3genes_GeneSets_Conditions, C4_Tip_Cells, Cp_G2M, C1a_Arterial_capillaries)

class(DotPlot_list)

DotPlot_list <- unique(DotPlot_list)

DotPlot_list_caps <- toupper(DotPlot_list)

DotPlot_list_caps

cairo_pdf("Plots/Figure_7/Figure_7E.3_DotPlot_GeneSets.pdf", width = 11, height = 2.2, family = "Arial")
DotPlot(Liver, features = DotPlot_list_caps, col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  scale_y_discrete(labels=c(bquote(italic(Dll4)^iDEC+aVEGF), bquote(italic(Dll4/Myc)^iDEC), bquote(italic(Dll4)^iDEC), "Control"))+
  scale_x_discrete(label = DotPlot_list)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_vline(xintercept = c(9.5, 18.5, 27.5, 31.5, 35.5),linetype = 2 )+
  theme(legend.box = "horizontal")
dev.off()

###############################################################################################################################

# Figure 7H Violin Plot

Idents(Liver) <- "Condition"
levels(Liver) <- c("Control", "Dll4KO", "D4KO_aVEGF")

genes <- c("MYC", "ODC1")
gene.names <- str_to_title(genes)
genes <- as.data.frame(genes)
gene.names <- as.data.frame(gene.names)
cols = palette_Maca_Vln
labels=c("Control", expression(italic(Dll4)^"iDEC"), expression(italic(Dll4)^"iDEC"+aVEGF))

cairo_pdf("Plots/Figure_7/Figure_7H_VlnPlot.pdf", width = 6, height = 12, family = "Arial")
Violin_plot_stacked(Liver, genes, gene.names, cols = cols, labels = labels)
dev.off()

jpeg("Plots/Figure_7/Figure_7H_VlnPlot.jpeg", width = 6, height = 6, units = 'in', res = 800)
Violin_plot_stacked(Liver, genes, gene.names, cols = cols, labels = labels)
dev.off()

###############################################################################################################################

# Figure 7J Dotplot

Liver <- subset(Liver_all, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "D4KO_aVEGF")

Idents(Liver) <- "Condition"

levels(Liver) <- c("Control", "Dll4KO", "D4KO_aVEGF")

Flow_genes <- c("KLF2", "KLF4")

genes <- rev(Flow_genes)

gene.names <- str_to_title(Flow_genes)

p1 <- DotPlot(Liver, features = Flow_genes, col.min = 0, dot.scale = 5, scale = F)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 5), legend.text  = element_text(size = 4),
        legend.key.size = unit(0.1, "cm")) +
  scale_y_discrete(labels=c("Control", bquote(italic(Dll4)^iDEC), bquote(italic(Dll4)^iDEC+aVEGF)))+
  scale_x_discrete(label = gene.names)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  coord_flip()

cairo_pdf("Plots/Figure_7/Figure_7J_DotPlot.pdf", width = 7, height = 2.5, family = "Arial")
p1
dev.off()

jpeg("Plots/Figure_7/Figure_7J_DotPlot.jpeg", width = 7, height = 2.5, res = 800, units = "in")
p1
dev.off()

###############################################################################################################################

# Figure 7K DEG Venn Diagram

Liver <- subset(Liver_all, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "D4KO_aVEGF" | Condition == "Dll4/MycKO")

Idents(Liver) <- "Condition"

levels(Liver) <- c("Control", "Dll4KO", "Dll4/MycKO", "D4KO_aVEGF")

Liver@active.ident -> Liver@meta.data$Condition

# Produce DEG Tables and tables of gene names upregulated and downregulated

conditions <- c("Dll4KO", "Dll4/MycKO", "D4KO_aVEGF")

DEG_Analysis_per_ConditionvCtl(Liver, conditions = conditions)

# You can find the tables in the "Tables" directory

DEGs_Upregulated <- read_excel("Tables/DEG_genes_up_per_condition_logfc0.25.xlsx")

DEGs_Upregulated <- DEGs_Upregulated[,2:4]

Venn_DEGs_Upregulated_3groups <- list(
  DEGs_Upregulated$Dll4KOvControl %>% na.omit(),
  DEGs_Upregulated$`Dll4/MycKOvControl` %>% na.omit(),
  DEGs_Upregulated$D4KO_aVEGFvControl %>% na.omit()
)

cairo_pdf("Plots/Figure_7/Supplementary/Venn_Diagram_DEGs_Upregulated_3groups.pdf", width = 5, height = 5, family = "Arial")
display_venn(Venn_DEGs_Upregulated_3groups,   category.names = c("CtrlvDll4" , "CtrlvDll4Myc" , "CtrlvDll4aVEGF"),
             fill = c("#F94040", "#C19EDD", "#BB005E"), main = "DEGs Upregulated 3groups")
dev.off()

###############################################################################################################################

# Figure 7L GSEA

# This function creates a Tables directory where it will produce a ranks table, a fgsea table and a wilcoxauc table

GSEA_loop_HallMark_GS(Liver, c("Dll4KO", "Dll4/MycKO", "D4KO_aVEGF"), "Homo sapiens")

# Edited GSEA Table output manually to match avgexp scRNASeq dataset (check table fgseaResD4&D4MycTidy.auc.xlsx)

cond.fgsea <- read.xlsx("./Tables/GSEA_Fig7/fgseaResD4&D4Myc&D4aVEGFTidy.auc.xlsx", sheet = 4, sep.names = " ", rowNames = T)

Idents(Liver) <- "Condition"

liver.gs <- subset(Liver, subset = Condition == "Dll4KO" | Condition == "Dll4/MycKO" | Condition == "D4KO_aVEGF")

avgexp <- AverageExpression(liver.gs, return.seurat = T)

avgexp <- Seurat::AddMetaData(avgexp, cond.fgsea)

avgexp@meta.data$Condition <- avgexp@active.ident

levels(avgexp) <- c("Dll4KO", "Dll4/MycKO", "D4KO_aVEGF")

avgexp@meta.data$Condition <- avgexp@active.ident

avgexp@meta.data$Condition <- as.factor(avgexp@meta.data$Condition)

avgexp@meta.data$Condition <- avgexp@active.ident

# All gene sets

cairo_pdf("Plots/Figure_6/Figure_6I_GSEA_Heatmap_Dll4vsControl_Dll4MyciDECvsControl_unscaled.pdf", width = 5, height = 12, family = "Arial")
dittoHeatmap(avgexp, genes = NULL, metas = names(cond.fgsea), 
             annot.by = "Condition",
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
             # order.by = c(2,1),
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

###############################################################################################################################

# Figure 7N UMAP

Liver4d <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver4d.Mapping.rds")

table(Liver4d@meta.data$Condition)

Liver4d@active.ident <- Liver4d@meta.data$Condition

new_labels <- c("Control(4d)" = "Control", "Dll4KO(4d)+Vehicle" = "italic(Dll4)^iDEC+Veh", "Dll4KO(4d)+SL327" = "italic(Dll4)^iDEC+SL327")
cols = my_palette_Rui_colors_B

cairo_pdf("Plots/Figure_7/Figure_7LNUMAP.pdf",  width = 28, height = 6, family = "Arial")
UMAP_Conditions_Global(Liver4d, 3, new_labels, cols)
dev.off()

jpeg("Plots/Figure_7/Figure_7N_UMAP.jpeg", width = 28, height = 6, units = 'in', res = 800)
UMAP_Conditions_Global(Liver4d, 3, new_labels, cols)
dev.off()

# Barplot

p1 <- dittoBarPlot(Liver4d, var = "FigClustering", group.by = "Condition", color.panel = my_palette_Rui_colors_B, xlab = NULL, main = NULL, y.breaks = c(0,0.25,0.5, 0.75, 1), min = 0, max = 1)+
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 18), legend.text = element_text(size = 16), axis.text.x = element_text(angle = 25))+
  scale_x_discrete(labels=c("Control", bquote(italic(Dll4)^iDEC), bquote(italic(Dll4)^iDEC+aVEGF)))+
  geom_col(width = 0.1)+
  # scale_y_discrete(labels=c("0", "25", "50", "75", "100"))
  scale_y_continuous(labels = c("0", "25", "50", "75", "100"))

cairo_pdf("Plots/Figure_7/Figure_7N_BarPlot.pdf",  width = 6, height = 6, family = "Arial")
p1
dev.off()


jpeg("Plots/Figure_7/Figure_7N_BarPlot.jpeg", width = 6, height = 6, units = 'in', res = 800)
p1
dev.off()

###############################################################################################################################

# Figure 7Q Violin Plot

genes <- c("KCNE3")
gene.names <- str_to_title(genes)
genes <- as.data.frame(genes)
gene.names <- as.data.frame(gene.names)
cols = palette_Maca_Vln
labels=c("Control", expression(italic(Dll4)^iDEC+Veh), expression(italic(Dll4)^iDEC+SL327))

cairo_pdf("Plots/Figure_5/Figure_5A_VlnPlot.pdf", width = 3, height = 12, family = "Arial")
Violin_plot_stacked(Liver4d, genes, gene.names, cols = cols, labels = labels)
dev.off()

jpeg("Plots/Figure_4/Figure_4H_VlnPlot.jpeg", width = 6, height = 6, units = 'in', res = 800)
Violin_plot_stacked(Liver4d, genes, gene.names, cols = cols, labels = labels)
dev.off()
