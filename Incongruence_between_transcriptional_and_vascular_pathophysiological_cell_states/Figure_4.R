########################################
# This file contains the plots produced in Figure 4
# 
# Author:
#   Alvaro Regano aregano@cnic.es
#
########################################

# Load Libraries from the Functions.R script

# Create Plot directory

dir.create("Plots/")
dir.create("Plots/Figure_4/")

# IMPORTANT: Before starting this script, please execute all functions in the functions.R script

# rds

Liver_all <- readRDS("")

# Color Palettes

my_palette_Rui_colors <- c("#E95A74", "#50B6EF", "#45FF8E", "#F4A753", "#A80519", "#880088", "#E28CF4", "#C1B80C","#FC0808",  "#0E47D8")

my_palette_Rui_colors_B <- c("#F4A753", "#E95A74", "#50B6EF", "#FC0808",  "#0E47D8", "#45FF8E", "#A80519", "#E28CF4", "#880088", "#C1B80C")

custom_rdylblu_palette <- c("#005AD9", "#0070FF", "#FFE135", "#FFDA13", "#E60026", "#C60012")

macarena_palette <- c("#FDDC85", "#6E2914")

Bestholtz_palette <- c("#DEDAD6", "#FEE392", "#FEC44E", "#FE9929", "#ED6F11", "#CC4C17", "#993411", "#65260C")

palette_Maca_Vln <- c("#606060", "#F94040", "#05BE78", "#55A0FB")

BuRd <- c("#052F60", "#2066AC", "#4392C3", "#92C5DE", "#D1E4EF", "#F7F7F7", "#FCDBC6", "#F4A482", "#D55F4D", "#B1172B", "#67001F")

###############################################################################################################################

## add the Arial font
font_add("Arial", regular = "arial.ttf",
         bold = "arialbd.ttf", italic = "ariali.ttf", bolditalic = "arialbi.ttf")

###############################################################################################################################

Liver <- subset(Liver_all, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "NOTCH1KO" | Condition == "RBPJKO")

levels(Liver) <- c("Control", "Dll4KO", "NOTCH1KO", "RBPJKO")

table(Liver@meta.data$Condition)

Liver@active.ident -> Liver@meta.data$Condition


# Figure 4B. Violin Plot

genes <- list(c("DLL4", "NOTCH1", "RBPJ"))
genes <- as.data.frame(genes)
gene.names <- list(c("Dll4", "Notch1", "Rbpj"))
gene.names <- as.data.frame(gene.names)
cols = palette_Maca_Vln
labels = c("Control", expression(italic("Dll4"^"iDEC")), expression(italic(Notch1)^"iDEC"), expression(italic(Rbpj)^"iDEC"))


cairo_pdf("Plots/Figure_4/Figure_4B_VlnPlot.pdf", width = 9, height = 12, family = "Arial")
Violin_plot_stacked(Liver, genes, gene.names, cols, labels)
dev.off()

jpeg("Plots/Figure_4/Figure_4B_VlnPlot.jpeg", width = 8, height = 12, units = 'in', res = 800)
Violin_plot_stacked(Liver, genes, gene.names)
dev.off()

###############################################################################################################################

# Figure 4C. UMAP

new_labels <- c("Control" = "Control", "Dll4KO" = "italic(Dll4)^iDEC", "NOTCH1KO" = "italic(Notch1)^iDEC", "RBPJKO" = "italic(Rbpj)^iDEC")
cols = my_palette_Rui_colors_B

cairo_pdf("Plots/Figure_4/Figure_4C1_UMAP.pdf",  width = 34, height = 6, family = "Arial")
UMAP_Conditions_Global(Liver, 4, new_labels, cols)
dev.off()

jpeg("Plots/Figure_4/Figure_4C1_UMAP.jpeg", width = 34, height = 6, units = 'in', res = 800)
UMAP_Conditions_Global(Liver, 4, new_labels, cols)
dev.off()

###############################################################################################################################

# Figure 4D. Bar Plot

p1 <- dittoBarPlot(Liver, var = "FigClustering", group.by = "Condition", color.panel = my_palette_Rui_colors_B, xlab = NULL, main = NULL, y.breaks = c(0,0.25,0.5, 0.75, 1), min = 0, max = 1)+
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 18), legend.text = element_text(size = 16), axis.text.x = element_text(angle = 25))+
  scale_x_discrete(labels=c("Control", bquote(italic(Dll4)^iDEC), bquote(italic(Notch1)^iDEC), bquote(italic(Rbpj)^iDEC)))+
  geom_col(width = 0.1)+
  scale_y_continuous(labels = c("0", "25", "50", "75", "100"))
p1

cairo_pdf("Plots/Figure_4/Figure_4C2_BarPlot.pdf",  width = 7, height = 6, family = "Arial")
p1
dev.off()


jpeg("Plots/Figure_4/Figure_4C2_BarPlot.jpeg", width = 7, height = 6, units = 'in', res = 800)
p1
dev.off()

###############################################################################################################################

# Figure4E.  DotPlot with Top10 Upregulated genes per cluster.

DEG_FigClustering <- FindAllMarkers(Liver, logfc.threshold = 0.5, test.use = "wilcox", only.pos = T, verbose = T, min.diff.pct = -Inf)
DEG_FigClustering$order <- 1:nrow(DEG_FigClustering)
write.xlsx(DEG_FigClustering, "Tables/DEG_FigClustering.wilcox.xlsx", rowNames = T)
DEG_FigClustering.wilcox <- read_excel("Tables/DEG_FigClustering.wilcox.xlsx")

DEG_FigClustering <- as.data.frame(DEG_FigClustering.wilcox)
DEG_FigClustering$order <- 1:nrow(DEG_FigClustering)

Top10_logFC <- DEG_FigClustering %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
Top10 <- DEG_FigClustering %>% group_by(cluster) %>% top_n(n = 10, wt = -order)

write.xlsx(Top10, "Tables/Top10_per_cluster.wilcox.xlsx", rowNames = T)

# After comparison and curation of the DEG, specific ones were selected

#  DotPlot with Macarena selected genes

Dotplot_genes <- read.csv2("DotplotGenes_MFC.csv", sep = ",", header = F)

Dotplot_genes <- Dotplot_genes$V1

Idents(Liver) <- Liver@meta.data$FigClustering

Dotplot_genesR <- rev(Dotplot_genes) 

Dotplot_genesR

figclusteringR <- c("C0 - Unspecified quiescent capillaries", "C1a - Arterial capillaries", "C1v - Venous capillaries", "C2a - Large arteries", "C2v - Large veins",  "C3 - Activated capillaries", 
                    "C4 - Endothelial tip cells", "C5ip - Proliferating S-phase", "C5p - G2-M", "C6 - X") %>% rev()

figclusteringR

levels(Liver) <- figclusteringR

DotPlotGenesH <- unique(Dotplot_genes) %>% tolower() %>% str_to_title()

df <- data.frame(
  x = rep(c(5.5, 15.5, 25.5, 35.5, 45.5, 55.5, 65.5, 75.5, 85.5, 95.5)),
  y = rep(c(11), 1),
  z = factor(rep(1:10)),
  color = my_palette_Rui_colors)

jpeg("Plots/Figure_4/Figure_4E_DotPlot_Maca_genes.Final.jpeg", width = 22, height = 3, units = 'in', res = 800)
DotPlot(Liver, features = Dotplot_genes, col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  scale_x_discrete(label = DotPlotGenesH)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_vline(xintercept = c(10.5, 20.5, 30.5, 40.5, 50.5, 60.5, 70.5, 80.5, 90.5),linetype = 2 )+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = my_palette_Rui_colors_B)
dev.off()

cairo_pdf("Plots/Figure_4/Figure_4E_DotPlot_Maca_genes.Final.pdf", width = 22, height = 3, family = "Arial")
DotPlot(Liver, features = Dotplot_genes, col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  scale_x_discrete(label = DotPlotGenesH)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_vline(xintercept = c(10.5, 20.5, 30.5, 40.5, 50.5, 60.5, 70.5, 80.5, 90.5),linetype = 2 )+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = my_palette_Rui_colors_B)
dev.off()


###############################################################################################################################

# Figure 4F. Heatmap LSEC and CLEC

LSEC_genes_names <- c("Clec1b", "Col13a1", "Dpp4", "Wnt2", "Clec4g", "Maf", "Gpr182", "Bmp2", "Pde2a", "Hpse", "Fcgr2b", "Lyve1", 
                      "Stab2", "Tfec", "Ehd3", "Flt4", "Ctsk", "Ifitm1", "Il1a", "Rnd3", "Ccl24", "Cd5l", "Pianp", "Gata4")

LSEC_genes <- toupper(LSEC_genes_names)


CEC_genes_names <- c("Esm1", "Mcam", "Cxcr4", "Cd34", "Cav1", "Aplnr", "Vegfa", "Emcn", "Meox2", "Apln", "Gabre", "Acvrl1",
                     "Lama4", "Sox18", "Pecam1", "Fgfr1", "Nr2f2","Angptl2", "Rgs5", "Stc1", "Plod2", "Vwf", "Gata2", "Wif1")


all.genes <- rownames(Liver)
all.genes <- as.data.frame(all.genes)

CEC_genes <- toupper(CEC_genes_names)

genes_LSEC_CEC <- append(CEC_genes, LSEC_genes)

genes_names_LSEC_CEC <- str_to_title(genes_LSEC_CEC)

avgexpLiver <- AverageExpression(Liver, return.seurat = T)

avgexpLiver@active.ident -> avgexpLiver@meta.data$Condition

row.subsections <- c(24, 24)

p1 <- dittoHeatmap(avgexpLiver, genes = genes_LSEC_CEC, 
                   annot.by = "Condition",
                   annot.colors = palette_Maca_Vln,
                   fontsize = 7, 
                   cluster_cols = F,
                   cluster_rows = F,
                   heatmap.colors.max.scaled = BuRd,
                   heatmap.colors = BuRd,
                   complex = T, data.out = F,
                   scaled.to.max = T, 
                   row_split = data.frame(rep(c("CEC", "LSEC"), row.subsections)),
                   cluster_row_slices = FALSE,
                   row_gap = unit(3, "mm"),
                   row_labels = genes_names_LSEC_CEC,
                   heatmap_legend_param = list(title = NULL),
                   # row_names_gp = gpar(col = c("green", "orange"), fontsize = c(10, 14)),
                   # legend_breaks = c(-1,0,1), 
)

p1

cairo_pdf("Plots/Figure4/Figure_4F_Heatmap_LSEC_CEC_cluster_rows.pdf", width = 4, height = 7, family = "Arial")
p1
dev.off()

cairo_pdf("Plots/Final_queries/Figure_4F_Heatmap_LSEC_CEC.pdf", width = 3.5, height = 7, family = "Arial")
p1
dev.off()

###############################################################################################################################

# GSEA for LSEC and CLEC gene sets


condition_text <- c("Dll4vCtl")
table(Liver@meta.data$Condition)
class(Liver@meta.data$Condition)
levels(Liver@meta.data$Condition)
Liver@meta.data$Condition <- as.factor(Liver@meta.data$Condition)
Liver@meta.data$Condition <- droplevels(Liver@meta.data$Condition)
class(Liver@meta.data$Condition)
clusters <- levels(Liver@meta.data$Condition)
clusters <- clusters[2:4]
clusters_text <- gsub('\\[', '\\(',
                      gsub('\\]', '\\)', clusters))

mywilcoxauc_DEG <- vector('list', length(clusters))
mywilcoxauc_DEG_names <- vector('list', length(clusters))


for (i in 1:length(clusters)) {
  
  Liver_subset <- subset(Liver, subset = Condition == "Control" | Condition == clusters[i])
  
  
  liver.genes <- wilcoxauc(Liver_subset, 'Condition')
  
  mywilcoxauc_DEG[[i]] <- local({
    i <- i
    liver.genes
  })
  mywilcoxauc_DEG_names[[i]] <- local({
    i <- i
    paste(clusters_text[i])
  })
}
names(mywilcoxauc_DEG) <- mywilcoxauc_DEG_names 

# All cells

write.xlsx(mywilcoxauc_DEG, "./Tables/GSEA_wilcoxauc_Liver_ECs_condition_Dll4_Notch1_Rbpj.xlsx", rowNames = T)

mywilcoxauc_DEG <- read_excel_allsheets("./Tables/GSEA_wilcoxauc_Liver_ECs_condition_Dll4_Notch1_Rbpj.xlsx")

mywilcoxauc_DEG_stacked <- list.stack(mywilcoxauc_DEG)

dplyr::count(mywilcoxauc_DEG_stacked, group)

myfgsea <- vector('list', length(clusters))
myfgsea_names <- vector('list', length(clusters))
myranks <- vector('list', length(clusters))


fgsea_sets <- list(CEC_genes, LSEC_genes)
names(fgsea_sets) <- c("CEC gene set", "LSEC gene set")

for (i in 1:length(clusters)) {
  
  cluster.genes<- mywilcoxauc_DEG_stacked %>%
    dplyr::filter(group == clusters[i]) %>%
    arrange(desc(auc)) %>% 
    dplyr::select(feature, auc)
  
  ranks<- deframe(cluster.genes)
  
  myranks[[i]] <- local({
    i <- i
    ranks
  })
  myfgsea[[i]] <- local({
    i <- i
    fgseaRes<- fgsea(fgsea_sets, stats = ranks)
  })
  myfgsea_names[[i]] <- local({
    i <- i
    paste(clusters_text[i])
  })
  
}

class(ranks)


ranks
names(myranks) <- myfgsea_names
names(myfgsea) <- myfgsea_names 

write.xlsx(myranks, file = "./Tables/ranks_fgsea_Liver_ECs_per_Condition_Dll4_Rbpj_Notch1_LSEC_CEC.xlsx", rowNames =T)

write.xlsx(myfgsea, file = "./Tables/fgsea_Liver_ECs_per_Condition_Dll4_Rbpj_Notch1_LSEC_CEC.xlsx")


# Tidy Data with auc

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  head()

write.xlsx(fgseaResTidy, "Tables/Final_query_fgseaResTidy.auc.xlsx")

########################################################

# Figure 4G. GSEA LSEC and CLEC Dll4vCtl

cluster.genes<- mywilcoxauc_DEG_stacked %>%
  dplyr::filter(group == "Dll4KO") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

ranks<- deframe(cluster.genes)

# GSEA Style Plot
cairo_pdf("Plots/Figure_4/Figure_4G_CEC_GS_Enrichment.pdf", width = 4, height = 4, family = "Arial")
plotEnrichment(fgsea_sets[["CEC gene set"]],
               ranks, gseaParam = 1) + labs(title="CEC gene set Dll4vCtl")
dev.off()

cairo_pdf("Plots/Final_queries/Figure_4G_LSEC_GS_Enrichment.pdf", width = 4, height = 4, family = "Arial")
plotEnrichment(fgsea_sets[["LSEC gene set"]],
               ranks, gseaParam = 1) + labs(title="LSEC gene set Dll4vCtl")
dev.off()

########################################################################################################

# Figure 4H. Violin Plot

genes <- list(c("GATA4", "WNT2"))
genes <- as.data.frame(genes)
gene.names <- list(c("Gata4", "Wnt2"))
gene.names <- as.data.frame(gene.names)
cols = palette_Maca_Vln
labels = c("Control", expression(italic("Dll4"^"iDEC")), expression(italic(Notch1)^"iDEC"), expression(italic(Rbpj)^"iDEC"))


cairo_pdf("Plots/Figure_4/Figure_4H_VlnPlot.pdf", width = 6, height = 6, family = "Arial")
Violin_plot_stacked(Liver, genes, gene.names, cols, labels)
dev.off()

jpeg("Plots/Figure_4/Figure_4H_VlnPlot.jpeg", width = 6, height = 6, units = 'in', res = 800)
Violin_plot_stacked(Liver, genes, gene.names, cols, labels)
dev.off()



###############################################################################################################################

# Figure 4I. Violin Plot

genes <- c("LTBP4", "MSR1", "WNT2", "KCNE3", "ESM1", "ANGPT2", "MYC", "ODC1", "RPL32")
gene.names <- str_to_title(genes)
genes <- as.data.frame(genes)
gene.names <- as.data.frame(gene.names)
cols = palette_Maca_Vln
labels = c("Control", expression(italic("Dll4"^"iDEC")), expression(italic(Notch1)^"iDEC"), expression(italic(Rbpj)^"iDEC"))


cairo_pdf("Plots/Figure_4/Figure_4H_VlnPlot.pdf", width = 6, height = 6, family = "Arial")
Violin_plot_stacked(Liver, genes, gene.names)
dev.off()

jpeg("Plots/Figure_4/Figure_4H_VlnPlot.jpeg", width = 6, height = 6, units = 'in', res = 800)
Violin_plot_stacked(Liver, genes, gene.names)
dev.off()


###############################################################################################################################

# Figure 4N. Violin Plot

genes <- c("CDKN1A", "TRP53", "CDKN2A")
gene.names <- str_to_title(genes)
genes <- as.data.frame(genes)
gene.names <- as.data.frame(gene.names)
cols = palette_Maca_Vln[c(1,3:4)]
labels = c("Control", expression(italic(Notch1)^"iDEC"), expression(italic(Rbpj)^"iDEC"))
Liver <- subset(Liver_all, subset = Condition == "Control" | Condition == "NOTCH1KO" | Condition == "RBPJKO")

levels(Liver) <- c("Control", "NOTCH1KO", "RBPJKO")

cairo_pdf("Plots/Figure_4/Figure_4H_VlnPlot.pdf", width = 6, height = 6, family = "Arial")
Violin_plot_stacked(Liver, genes, gene.names, cols, labels)
dev.off()

jpeg("Plots/Figure_4/Figure_4H_VlnPlot.jpeg", width = 6, height = 6, units = 'in', res = 800)
Violin_plot_stacked(Liver, genes, gene.names)
dev.off()

###############################################################################################################################

# Subsetting and color palette for Liver "Control", "Dll4KO", "NOTCH1KO", "N1N4HET"

Liver <- subset(Liver_all, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "NOTCH1KO" | Condition == "N1N4HET")

Liver@meta.data$Condition <- droplevels(Liver@meta.data$Condition)

levels(Liver) <- c("Control", "Dll4KO", "NOTCH1KO", "N1N4HET")

palette_Maca_Vln <- c("#606060", "#F94040", "#05BE78", "#FFA040")

Liver@active.ident -> Liver@meta.data$Condition

###############################################################################################################################

# Figure 4O

genes <- c("NOTCH1", "NOTCH2", "NOTCH4")
gene.names <- str_to_title(genes)
genes <- as.data.frame(genes)
gene.names <- as.data.frame(gene.names)

cols = palette_Maca_Vln
labels = c("Control", expression(italic("Dll4"^"iDEC")), expression(italic(Notch1)^"iDEC"), expression(italic(Notch1/2/4)^"iDEC"))

cairo_pdf("Plots/Figure_4/Figure_4O_VlnPlot.pdf", width = 6, height = 9, family = "Arial")
Violin_plot_stacked(Liver, genes, cols, gene.name = gene.names, labels = labels)
dev.off()

jpeg("Plots/Figure_4/Figure_4O_VlnPlot.jpeg", width = 6, height = 9, units = 'in', res = 800)
Violin_plot_stacked(Liver, genes, gene.names)
dev.off()


# Figure 4O UMAP


new_labels <- c("Control" = "Control", "Dll4KO" = "italic(Dll4)^iDEC", "NOTCH1KO" = "italic(Notch1)^iDEC", "N1N4HET" = "italic(Notch1/2/4)^iDEC")
cols = my_palette_Rui_colors_B

cairo_pdf("Plots/Figure_4/Figure_4O_UMAP.pdf",  width = 34, height = 6, family = "Arial")
UMAP_Conditions_Global(Liver, 4, new_labels, cols)
dev.off()

jpeg("Plots/Figure_4/Figure_4O_UMAP.jpeg", width = 34, height = 6, units = 'in', res = 800)
UMAP_Conditions_Global(Liver, 4, new_labels, cols)
dev.off()

# Figure 4O BarPlot

p1 <- dittoBarPlot(Liver, var = "FigClustering", group.by = "Condition", x.reorder = c(1,2,4,3), color.panel = my_palette_Rui_colors_B, xlab = NULL, main = NULL, y.breaks = c(0,0.25,0.5, 0.75, 1), min = 0, max = 1)+
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 18), legend.text = element_text(size = 16), axis.text.x = element_text(angle = 25))+
  scale_x_discrete(labels=c("Control", bquote(italic(Dll4)^iDEC), bquote(italic(Notch1)^iDEC), bquote(italic(Notch1/2/4)^iDEC)))+
  geom_col(width = 0.1)+
  # scale_y_discrete(labels=c("0", "25", "50", "75", "100"))
  scale_y_continuous(labels = c("0", "25", "50", "75", "100"))
p1

cairo_pdf("Plots/Figure_4/Figure_4O_BarPlot.pdf",  width = 7, height = 6, family = "Arial")
p1
dev.off()


jpeg("Plots/Figure_4/Figure_4O_BarPlot.jpeg", width = 7, height = 6, units = 'in', res = 800)
p1
dev.off()

###############################################################################################################################

# Figure 4Q

genes <- c("HES1", "MSR1", "KCNE3", "ESM1", "MYC", "ODC1", "STMN1", "CDKN1A")
gene.names <- str_to_title(genes)
genes <- as.data.frame(genes)
gene.names <- as.data.frame(gene.names)

cols = palette_Maca_Vln
labels = c("Control", expression(italic("Dll4"^"iDEC")), expression(italic(Notch1)^"iDEC"), expression(italic(Notch1/2/4)^"iDEC"))

cairo_pdf("Plots/Figure_4/Figure_4Q_VlnPlot.pdf", width = 6, height = 16, family = "Arial")
Violin_plot_stacked(Liver, genes, cols, gene.name = gene.names, labels = labels)
dev.off()

jpeg("Plots/Figure_4/Figure_4Q_VlnPlot.jpeg", width = 6, height = 16, units = 'in', res = 800)
Violin_plot_stacked(Liver, genes, gene.names)
dev.off()

###############################################################################################################################

Liver <- subset(Liver_all, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "D1D4HET" | Condition == "WT_DBZ")

levels(Liver) <- c("Control", "Dll4KO", "D1D4HET", "WT_DBZ")

Liver@active.ident -> Liver@meta.data$Condition

new_labels <- c("Control" = "Control", "Dll4KO" = "italic(Dll4)^iDEC", "D1D4HET" = "italic(Dll4Het)^iDEC", "WT_DBZ" = "italic(DBZ)^iDEC")
cols <- my_palette_Rui_colors_B

cairo_pdf("Plots/Figure_4/Figure_4S_UMAP.pdf",  width = 34, height = 6, family = "Arial")
UMAP_Conditions_Global(Liver, 4, new_labels, cols)
dev.off()

jpeg("Plots/Figure_4/Figure_4O_UMAP.jpeg", width = 34, height = 6, units = 'in', res = 800)
UMAP_Conditions_Global(Liver, 4, new_labels, cols)
dev.off()

###############################################################################################################################

# Figure 4S

genes <- c("HES1")
gene.names <- str_to_title(genes)
genes <- as.data.frame(genes)
gene.names <- as.data.frame(gene.names)

cols = c("#606060", "#F94040", "#9E1F63", "#ACC2D9")

labels = c("Control", expression(italic("Dll4"^"iDEC")), expression(italic(Dll4Het)^"iDEC"), "DBZ")

cairo_pdf("Plots/Figure_4/Figure_4Q_VlnPlot.pdf", width = 6, height = 3, family = "Arial")
Violin_plot_stacked(Liver, genes, cols, gene.name = gene.names, labels = labels)
dev.off()

jpeg("Plots/Figure_4/Figure_4Q_VlnPlot.jpeg", width = 6, height = 3, units = 'in', res = 800)
Violin_plot_stacked(Liver, genes, gene.names)
dev.off()


