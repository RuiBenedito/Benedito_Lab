########################################
# This file contains the plots produced in Suppl Figure 4
# 
# Author:
#   Alvaro Regano aregano@cnic.es
#
########################################

# Load Libraries from the Functions.R script

# Create Plot directory

dir.create("Plots/")
dir.create("Plots/Suppl_Figure_4/")

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



###################################################################

# Suppl Figure 4B Feature plots

Liver <- subset(Liver_all, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "NOTCH1KO" | Condition == "RBPJKO")

levels(Liver) <- c("Control", "Dll4KO", "NOTCH1KO", "RBPJKO")

table(Liver@meta.data$Condition)

Liver@active.ident -> Liver@meta.data$Condition

DEG.names.text <- c("Gja5","Ltbp4", "Msr1", "Rspo3", "Eif3f", "Wnt2", 
                     "Kcne3","Esm1", "Cd34",  "Odc1", "Myc",
                    "Mki67", "Stmn1", "Mcm2", "Ctnnb1", "Malat1",
                    "Eif3f", "Meis1", "Rpl10a", "Rab3b")

DEG.names <- toupper(DEG.names.text)

myplots <- vector('list', length(DEG.names))

for (i in 1:length(DEG.names)) {
  
  p3 <- FeaturePlot(Liver, features = DEG.names[i], order = T, pt.size = 1, slot = "data", combine = F)
  
  p4 <- lapply(X = p3, FUN = function(p) p + scale_colour_gradientn(colors = Bestholtz_palette) + ggtitle(bquote(~italic(.(DEG.names.text[i])))))
  
  p4 <- Reduce( `+`, p4 )+patchwork::plot_layout( ncol = 1 )
  
  myplots[[i]] <- local({
    i <- i
    p4
  })
  
  # plegend <- VlnPlot(Liver, features = genes[i, 1], group.by = "Condition")
}

myplots[10]

p25 <- patchwork::wrap_plots(myplots, ncol = 5)

p25

cairo_pdf("Plots/Plots/Suppl_Figure_4/Suppl_Figure_4B_FeaturePlots.pdf", width = 25, height = 20, family = "Arial")
p25
dev.off()

###################################################################

# Suppl Figure 4C Violin Plots

gene.names <- c("Ltbp4", "Msr1", "Rspo3", "Rab3b",  "Rpl10a", 
           "Esm1", "Cd34",  "Odc1", "Mcm2")
genes <- toupper(gene.names)
genes <- as.data.frame(genes)
gene.names <- as.data.frame(gene.names)
cols = palette_Maca_Vln
labels = c("Control", expression(italic("Dll4"^"iDEC")), expression(italic(Notch1)^"iDEC"), expression(italic(Rbpj)^"iDEC"))


cairo_pdf("Plots/Suppl_Figure_4/Suppl_Figure_4C_VlnPlot.pdf", width = 6, height = 32, family = "Arial")
Violin_plot_stacked(Liver, genes, gene.names, cols = cols, labels = labels)
dev.off()

jpeg("Plots/Suppl_Figure_4/Suppl_Figure_4C_VlnPlot.jpeg", width = 6, height = 32, units = 'in', res = 800)
Violin_plot_stacked(Liver_ctl, genes, gene.names, cols = cols, labels = labels)
dev.off()
