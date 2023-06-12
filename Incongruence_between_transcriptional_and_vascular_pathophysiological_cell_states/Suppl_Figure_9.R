########################################
# This file contains the plots produced in the Supplementary Figure 9
# 
# Author:
#   Alvaro Regano aregano@cnic.es
#
########################################

# Create Plot directory

dir.create("Plots/")
dir.create("Plots/Suppl_Figure_9/")

# Load Libraries from the Functions.R script

# IMPORTANT: Before starting this script, please execute all functions in the functions.R script

# rds

Liver_all <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.OldDll4Myc.rds")

table(Liver_all@meta.data$Condition)

Liver <- subset(Liver_all, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "Dll4/MycKO" | Condition == "D4KO_aVEGF")

table(Liver@meta.data$Condition)
Idents(Liver) <- "Condition"
levels(Liver) <- c("Control", "Dll4KO", "Dll4/MycKO", "D4KO_aVEGF")
Liver@active.ident -> Liver@meta.data$Condition
# Color Palettes

my_palette_Rui_colors_B <- c("#F4A753", "#E95A74", "#50B6EF", "#FC0808",  "#0E47D8", "#45FF8E", "#A80519", "#E28CF4", "#880088", "#C1B80C")
Conditions_palette <- c("#606060", "#F94040", "#FF8080", "#BB005E")
custom_rdylblu_palette <- c("#005AD9", "#0070FF", "#FFE135", "#FFDA13", "#E60026", "#C60012")
BuRd <- c("#052F60", "#2066AC", "#4392C3", "#92C5DE", "#D1E4EF", "#F7F7F7", "#FCDBC6", "#F4A482", "#D55F4D", "#B1172B", "#67001F")
Bestholtz_palette <- c("#DEDAD6", "#FEE392", "#FEC44E", "#FE9929", "#ED6F11", "#CC4C17", "#993411", "#65260C")
BuRd <- c("#052F60", "#2066AC", "#4392C3", "#92C5DE", "#D1E4EF", "#F7F7F7", "#FCDBC6", "#F4A482", "#D55F4D", "#B1172B", "#67001F")

########################################################################

# Suppl. Fig. 9A Violin Plot

genes <- c("VEGFA", "KDR")
gene.names <- str_to_title(genes)
genes <- as.data.frame(genes)
gene.names <- as.data.frame(gene.names)
cols = Conditions_palette

Idents(Liver) <- "FigClustering"


#  Script

myplots <- vector('list', nrow(genes))


for (i in 1:(nrow(genes)-1)) {
  
  p21 <- VlnPlot(Liver, features = genes[i, 1], group.by = "FigClustering",
                 cols = my_palette_Rui_colors_B)+ NoLegend() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), plot.title = element_text(face = "bold.italic", size = 16))
  
  myplots[[i]] <- local({
    i <- i
    p21
  })
  
}

p2 <- VlnPlot(Liver, features = genes[nrow(genes), 1], group.by = "FigClustering",
              cols = my_palette_Rui_colors_B)+ NoLegend() + theme(axis.title.y=element_blank(), plot.title = element_text(face = "bold.italic", size = 16))

myplots[[nrow(genes)]] <- local({
  p2
})


p2

plegend <- VlnPlot(Liver_ctl, features = genes[i, 1], group.by = "FigClustering",
                   cols = Conditions_palette) +
  scale_fill_manual(values= Conditions_palette, labels=c("Control", expression(italic("Dll4"^"iDEC"))))

legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", legend.text = element_text(size = 18)))

p25 <- patchwork::wrap_plots(myplots, ncol = 1)

p26 <- plot_grid(p25, legend, nrow = 2, rel_heights = c(1, .03))

x.axis <- ggdraw() + draw_label("Expression Level", fontface='bold', angle = 90, size = 24)

p27 <- plot_grid(x.axis, p26, ncol = 2, rel_widths = c(.05, 1))


cairo_pdf("Plots/Suppl_Figure_9/Suppl_Figure_9A_VlnPlot.pdf", width = 8, height = 22, family = "Arial")
p27
dev.off()


########################################################################

# Suppl. Figs. 9C,E, & G Violin Plots

Liver_vln <- subset(Liver, subset = Condition == "Dll4/MycKO", invert = T)
Idents(Liver_vln) <- "Condition"
levels(Liver_vln) <- c("Control", "Dll4KO", "D4KO_aVEGF")
Liver_vln@active.ident -> Liver_vln@meta.data$Condition
table(Liver_vln@meta.data$Condition)

genes <- c("DLL4", "MSR1", "LTBP4", "EFNB2", "KLF2", "KLF4")
gene.names <- str_to_title(genes)
genes <- as.data.frame(genes)
gene.names <- as.data.frame(gene.names)
cols = Conditions_palette[c(1,2,4)]
labels=c("Control", expression(italic(Dll4)^"iDEC"), expression(italic(Dll4)^"iDEC"+aVEGF))

cairo_pdf("Plots/Suppl_Figure_9/Suppl_Figure_9C_E_G_VlnPlot.pdf", width = 5, height = 12, family = "Arial")
Violin_plot_stacked(Liver_vln, genes, gene.names, cols = cols, labels = labels)
dev.off()

########################################################################

# Suppl. Fig. 9F DotPlot Cell Cycle 

Top10_Cell_Cycle_H <- c("Mcm6", "Mcm5", "Lig1", "Mcm2", "Mcm3", "Top2a", "Stmn1", "Mki67", "Hmgb2", "Cenpf")

Top10_Cell_Cycle <- Top10_Cell_Cycle_H %>% toupper()

cairo_pdf("Plots/Suppl_Figure_9/Suppl_Figure_9F_Dotplot.pdf", width = 7, height = 2, family = "Arial")
DotPlot(Liver, features = Top10_Cell_Cycle, col.min = 0, dot.scale = 10, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  scale_y_discrete(labels=c(bquote(italic(Dll4)^iDEC+aVEGF), bquote(italic(Dll4/Myc)^iDEC), bquote(italic(Dll4)^iDEC), "Control"))+
  scale_x_discrete(label = Top10_Cell_Cycle_H)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  theme(legend.box = "horizontal", title = element_text(size = 10))+
  ggtitle("GOBP: Cell Cycle")
dev.off()


########################################################################

# Suppl. Fig. 9H GOBP GSEA Heatmap


# This function creates a Tables directory where it will produce a ranks table, a fgsea table and a wilcoxauc table

GSEA_loop_GOBP_GS(Liver, c("Dll4KO", "Dll4/MycKO", "D4KO_aVEGF"), "Homo sapiens")


# After manually curating the Table in order to select the GOBP of interest


cond.fgsea <- read.xlsx("./Tables/fgsea_GOBP_selected_Dll4_Dll4Myc_Dll4aVEGF.xlsx", sheet = 4, sep.names = " ", rowNames = T)

liver.gs <- subset(Liver, subset = Condition == "Dll4KO" | Condition == "Dll4/MycKO" | Condition == "D4KO_aVEGF")
avgexp <- AverageExpression(liver.gs, return.seurat = T)
avgexp <- Seurat::AddMetaData(avgexp, cond.fgsea)
levels(avgexp) <- c("Dll4KO", "Dll4/MycKO", "D4KO_aVEGF")

avgexp@meta.data$Condition <- avgexp@active.ident
avgexp@meta.data$Condition <- as.factor(avgexp@meta.data$Condition)

# All gene sets

cairo_pdf("Plots/Suppl_Figure_9/Suppl_Figure_9H_GSEA_Heatmap.pdf", width = 5, height = 12, family = "Arial")
dittoHeatmap(avgexp, genes = NULL, metas = names(cond.fgsea), 
             annot.by = "Condition",
             annot.colors = Conditions_palette[2:4],
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

# Suppl. Fig. 9I DotPlot GOBP CELLULAR_MACROMOLECULE_BIOSYNTHETIC_PROCESS and RNA_PROCESSING

# Ranks produced in the GSEA_loop_GOBP_GS() function

ranksdll4aVEGF <- read_excel("./Tables/ranks_Condition.xlsx", sheet = 3 )

# Load fgsea dataset

msigdbr_species()

h_df<- msigdbr(species = "Homo sapiens", category = 7)

table(h_gene_sets$gs_cat)

# Using GOBP genes only

h_gene_sets = msigdbr(species = "Homo sapiens", subcategory = "GO:BP")

fgsea_sets<- h_gene_sets %>% split(x = .$gene_symbol, f = .$gs_name) %>% list.filter("Type" > 50)

fgsea_sets$GOBP_AGING

Cellular_Macromolecule_BS_Process_genes <- fgsea_sets$GOBP_CELLULAR_MACROMOLECULE_BIOSYNTHETIC_PROCESS

RNA_Processing_genes <- fgsea_sets$GOBP_RNA_PROCESSING

Cellular_Macromolecule_BS_Process_genes <- as.data.frame(Cellular_Macromolecule_BS_Process_genes)

RNA_Processing_genes <- as.data.frame(RNA_Processing_genes)

Cellular_Macromolecule_BS_Process_genes <- t(Cellular_Macromolecule_BS_Process_genes)

Cellular_Macromolecule_BS_Process_genes <- as.vector(Cellular_Macromolecule_BS_Process_genes)

Cellular_Macromolecule_rank <- ranksdll4aVEGF[ ranksdll4aVEGF$genes %in% Cellular_Macromolecule_BS_Process_genes, ]

RNA_Processing_genes <- t(RNA_Processing_genes)

RNA_Processing_genes <- as.vector(RNA_Processing_genes)

RNA_Processing_rank <- ranksdll4aVEGF[ ranksdll4aVEGF$genes %in% RNA_Processing_genes, ]

levels(Liver) <- c("D4KO_aVEGF", "Dll4/MycKO", "Dll4KO", "Control")


# Comparing repeated entries

Top10_GOBPs_check <- c(Top10_Cellular_Macromolecule, Top10_RNA_Processing, Top10_Cell_Cycle)

max_length <- max(c(length(Top10_GOBPs_check), length(Top10_GOBPs)))    # Find out maximum length
max_length 

Top10_GOBPs_check_1 <- data.frame(col1 = c(Top10_GOBPs_check),                 # Create data frame with unequal vectors
                                  col2 = c(Top10_GOBPs,
                                           rep(NA, max_length - length(Top10_GOBPs))))
Top10_GOBPs_check_1                                            # Print final da

# Repeated values are: 
# GOBP_RNA_Processing -> RPS8, RPL35A, RPS7
#  I will add to get those 10 per category

Top10_Cellular_Macromolecule <- Cellular_Macromolecule_rank$genes[c(1:10)]

Top10_Cellular_Macromolecule_H <- Cellular_Macromolecule_rank$genes[c(1:10)]

Top10_RNA_Processing <-RNA_Processing_rank$genes[c(1:13)]

Top10_RNA_Processing_H <- RNA_Processing_rank$genes[c(1:13)]


# Top10 Cellular Macromolecule Biosynthetic Process + RNA Processing


Top10_GOBPs <- c(Top10_Cellular_Macromolecule, Top10_RNA_Processing) %>% unique()

Top10_GOBPs_H <- c(Top10_Cellular_Macromolecule_H, Top10_RNA_Processing_H) %>% unique()

# Do the check

Top10_GOBPs_check <- c(Top10_Cellular_Macromolecule, Top10_RNA_Processing)

max_length <- max(c(length(Top10_GOBPs_check), length(Top10_GOBPs)))    # Find out maximum length
max_length 

Top10_GOBPs_check_1 <- data.frame(col1 = c(Top10_GOBPs_check),                 # Create data frame with unequal vectors
                                  col2 = c(Top10_GOBPs,
                                           rep(NA, max_length - length(Top10_GOBPs))))
Top10_GOBPs_check_1                                            # Print final da

cairo_pdf("Plots/Suppl_Figure_9/Suppl_Figure_9I_Dotplot.pdf", width = 8, height = 2, family = "Arial")
DotPlot(Liver, features = Top10_GOBPs, col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  scale_y_discrete(labels=c(bquote(italic(Dll4)^iDEC+aVEGF), bquote(italic(Dll4/Myc)^iDEC), bquote(italic(Dll4)^iDEC), "Control"))+
  scale_x_discrete(label = Top10_GOBPs_H)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_vline(xintercept = c(10.5),linetype = 2 )+
  theme(legend.box = "horizontal", title = element_text(size = 10))+
  ggtitle("GOBP: Cellular Macromolecule Biosynthetic Process; RNA Processing")
dev.off()



