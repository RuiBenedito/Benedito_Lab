########################################
# This file contains the functions used for plotting
# 
# Author:
#   Alvaro Regano aregano@cnic.es
#
########################################

# Libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(tibble)
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
library(SingleCellExperiment)
library(escape)
library(reshape2)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(presto)
library(msigdbr)
library("tidyselect")
library(magrittr)
library(tidyr)
library(forcats)
library("ggVennDiagram")
library("gridExtra") 


# Functions used in the writing of the scRNASeq data

read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

##################################################################################

cbind.all <- function (...) 
{
  nm <- list(...)
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow))
  do.call(cbind, lapply(nm, function(x) rbind(x, matrix(, n - nrow(x), ncol(x)))))
}

##################################################################################

Violin_plot_stacked <- function(seurat.object, genes, gene.names, cols, labels, legend.direction = "horizontal"){
  
  require(Seurat)
  require(ggplot2)
  require(patchwork)
  require(cowplot)
  require(stringr)
  
  myplots <- vector('list', nrow(genes))  
  
  for (i in 1:nrow(genes)) {
    
    p21 <- VlnPlot(seurat.object, features = genes[i, 1], cols = cols)+ NoLegend() + theme(text = element_text(family="Arial"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), plot.title = element_text(face = "bold.italic"))+ggtitle(gene.names[i, 1])
    
    myplots[[i]] <- local({
      i <- i
      p21
    })
    
  }
  
  plegend <- VlnPlot(seurat.object, features = genes[i, 1]) +
    scale_fill_manual(values= cols, labels= labels)
  plegend
  legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center", legend.direction = legend.direction, legend.text = element_text(size = 14)))
  
  p25 <- patchwork::wrap_plots(myplots, ncol = 1)
  
  p26 <- plot_grid(p25, legend, ncol = 1, rel_heights = c(1, .1))
  x.axis <- ggdraw() + draw_label("Expression Level", fontface='bold', angle = 90, size = 16+2*length(genes))
  
  p27 <- plot_grid(x.axis, p26, ncol = 2, rel_widths = c(.05, 1))
  
  p27}

UMAP_Conditions_Global <- function(seurat.object, ncond, labels, cols){
  
  require(Seurat)
  require(ggplot2)
  require(patchwork)
  
  p1 <- DimPlot(Liver, cols = cols, pt.size = 1.2, split.by = "Condition", group.by = "FigClustering", combine = T)
  
  p21 <- p1+facet_grid(~Condition, labeller = as_labeller(labels, default = label_parsed))+
    theme(strip.text.x = element_text(size = 24), plot.title = element_blank())+NoLegend()+xlab("UMAP_1")+ylab("UMAP_2")
  
  p22 <- DimPlot(Liver, cols = cols, group.by = "FigClustering", pt.size = 1.2)+theme(legend.text = element_text(size = 14), axis.title.x = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), plot.title = element_blank())
  
  p23 <- list(p21, p22)
  
  design <- c(patchwork::area(1, 1, 1, ncond), patchwork::area(1, ncond+1, 1, ncond+1.5))
  
  p24 <- Reduce( `+`,  p23)+patchwork::plot_layout(design = design)
  
  p24}

##########################################################################

GSEA_loop_HallMark_GS <- function(seurat.object, condition, species){
  
  require(fgsea)
  require(presto)
  require(openxlsx)
  require(msigdbr)
  require(rlist)
  require(dplyr)
  require(tibble)
  
  clusters_text <- gsub('\\[', '\\(',
                        gsub('\\]', '\\)', 
                             gsub('\\/', '\\_', condition)))
  
  if(!dir.exists("Tables/")){
    dir.create("Tables/")
  }else{
    
  }
  
  m_gene_sets = msigdbr(species = species, category = "H")
  
  fgsea_sets<- m_gene_sets %>% split(x = .$gene_symbol, f = .$gs_name)
  
  condition_text <- c("Ctlv")
  
  mywilcoxauc_DEG <- vector('list', length(condition))
  mywilcoxauc_DEG_names <- vector('list', length(condition))
  
  
  for (i in 1:length(condition)) {
    
    seurat.object.subset <- subset(seurat.object, subset = Condition == "Control" | Condition == condition[i])
    
    seurat.object.DEG <- wilcoxauc(seurat.object.subset, 'Condition')
    
    mywilcoxauc_DEG[[i]] <- local({
      i <- i
      seurat.object.DEG
    })
    mywilcoxauc_DEG_names[[i]] <- local({
      i <- i
      paste(condition_text, clusters_text[i], sep = "_")
    })
  }
  names(mywilcoxauc_DEG) <- mywilcoxauc_DEG_names 
  
  # All cells
  
  write.xlsx(mywilcoxauc_DEG, "./Tables/GSEA_wilcoxauc_condition.xlsx", rowNames = T)
  
  mywilcoxauc_DEG_stacked <- list.stack(mywilcoxauc_DEG)
  
  dplyr::count(mywilcoxauc_DEG_stacked, group)
  
  myfgsea <- vector('list', length(condition))
  myfgsea_names <- vector('list', length(condition))
  myranks <- vector('list', length(condition))
  
  for (i in 1:length(condition)) {
    
    condition.genes<- mywilcoxauc_DEG_stacked %>%
      dplyr::filter(group == condition[i]) %>%
      arrange(desc(auc)) %>% 
      dplyr::select(feature, auc)
    
    ranks<- deframe(condition.genes)
    
    myfgsea[[i]] <- local({
      i <- i
      fgseaRes<- fgsea(fgsea_sets, stats = ranks)
    })
    myfgsea_names[[i]] <- local({
      i <- i
      paste(condition_text, clusters_text[i], sep = "")
    })
    
    ranks<- as.data.frame(ranks)
    ranks$genes <- rownames(ranks)
    myranks[[i]] <- local({
      i <- i
      ranks
    })
  }
  
  names(myranks) <- myfgsea_names 
  names(myfgsea) <- myfgsea_names 
  
  write.xlsx(myfgsea, file = "./Tables/fgsea_Condition.xlsx")
  write.xlsx(myranks, file = "./Tables/ranks_Condition.xlsx")
}

##########################################################################


DEG_Analysis_per_ConditionvCtl <- function(seurat.object, control.name = "Control", conditions){
  
  require(Seurat)
  require(dplyr)
  require(openxlsx)
  
  condition_text <- gsub('\\[', '\\(',
                        gsub('\\]', '\\)', 
                             gsub('\\/', '\\_', conditions)))
  
  myDEGs <- vector("list", length = length(conditions))
  DEG_names <- vector("list", length = length(conditions))
  # up = data.frame(matrix(nrow = 0, ncol = length(conditions))) 
  # colnames(up) = paste(conditions, control.name, sep ="v")
  # down = data.frame(matrix(nrow = 0, ncol = length(conditions))) 
  # colnames(down) = paste(conditions, control.name, sep ="v")
  
  up = data.frame()
  down = data.frame()
  
  for (i in 1:length(conditions)) {
  
  
  DEG <- FindMarkers(seurat.object, ident.1 = conditions[i],  ident.2 = control.name, test.use = "wilcox", only.pos = F, verbose = T, min.diff.pct = -Inf)
  filter_DEG <- filter(DEG, p_val_adj < 0.05)
  up_DEG <- filter(filter_DEG, avg_log2FC > 0.25)
  down_DEG <- filter(filter_DEG, avg_log2FC < -0.25)
  
  up <- cbind.all(up, rownames(up_DEG))
  down <- cbind.all(down, rownames(down_DEG))
  myDEGs[[i]] <- local({
      i <- i
      DEG
    })
    DEG_names[[i]] <- local({
      i <- i
      paste(condition_text[i], control.name, sep = "v")
    })

  }

up <- as.data.frame(up)
colnames(up) <- paste(conditions, control.name, sep = "v")
down <- as.data.frame(down)
colnames(down) <- paste(conditions, control.name, sep = "v")

names(myDEGs) <- DEG_names

write.xlsx(myDEGs, "./Tables/DEG_per_condition_logfc0.25.xlsx", rowNames = T)
write.xlsx(up, "./Tables/DEG_genes_up_per_condition_logfc0.25.xlsx", rowNames = T)
write.xlsx(down, "./Tables/DEG_genes_down_per_condition_logfc0.25.xlsx", rowNames = T)

}

##########################################################################

display_venn <- function(x, ...){
  require(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}
