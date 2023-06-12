########################################
# This file contains the pre-processing steps for assembling the Liver non ECs scRNASeq dataset
# 
# Author:
#   Alvaro Regano aregano@cnic.es
#
########################################


library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(BiocStyle)
library(dittoSeq)
library(dplyr)
library(mcmcplots)
library(celldex)
library(SingleR)
library(scRNAseq)
library(SingleCellExperiment)
library(escape)
library(reshape2)
library(presto)
library(tidyverse)
library(scuttle)
library(scater)
library(scran)
library(pheatmap)

# Create Plot directory

dir.create("Plots/")
dir.create("Plots/Suppl_Figure_6/")

# Produce the raw data objects for Control and Dll4LOF separately and then merge them together

Liver_Non_ECs.data <- Read10X(data.dir = "raw_data/Liver_Non_ECs/Control/")
Liver_Non_ECs.data <- Read10X(data.dir = "raw_data/Liver_Non_ECs/Dll4LOF/")

Liver_Non_ECs <- CreateSeuratObject(counts = Liver_Non_ECs.data$`Gene Expression`, 
                                    project = 
                                      # "Liver_Dll4LOF_NonECs",
                                      "Liver_Control_NonECs",
                                    assay = "RNA",
                                    min.cells = 3 , min.features = 200)

Liver_Non_ECs[["HTO"]] <- CreateAssayObject(counts = Liver_Non_ECs.data$`Antibody Capture`[, colnames(x = Liver_Non_ECs)])

Liver_Non_ECs[["CMO"]] <- CreateAssayObject(counts = Liver_Non_ECs.data$`Multiplexing Capture`[, colnames(x = Liver_Non_ECs)])

Liver_Non_ECs -> Liver_Control_Non_ECs

Liver_Non_ECs -> Liver_Dll4KO_Non_ECs

rm(Liver_Non_ECs.data)

rm(Liver_Non_ECs)

Liver_Non_ECs <- merge(Liver_Control_Non_ECs, Liver_Dll4KO_Non_ECs)

table(Liver_Non_ECs@meta.data$orig.ident)

# there should be 1356 cells for Liver_Control_NonECs and 3647 cells for Liver_Dll4LOF_NonECs

saveRDS(Liver_Non_ECs, "../Sample_Liver_Heart_20Jan22/rds/Liver_Non_ECs_raw.rds")

VlnPlot(Liver_Non_ECs, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

# QC Analysis
DefaultAssay(Liver_Non_ECs) <- "RNA"
Liver_Non_ECs[["percent.mt"]] <- PercentageFeatureSet(Liver_Non_ECs, pattern = "^mt-")
head(Liver_Non_ECs@meta.data, 50)

VlnPlot(Liver_Non_ECs, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "orig.ident", ncol = 3, pt.size = 0)

plot1 <- FeatureScatter(Liver_Non_ECs, feature1 = "nFeature_RNA",
                        group.by = "orig.ident", feature2 = "percent.mt") +
  geom_vline(xintercept = c(200,5000),linetype = 2 ) +
  geom_hline(yintercept = 15 ,linetype = 2)
plot2 <- FeatureScatter(Liver_Non_ECs, feature1 = "nCount_RNA",
                        group.by = "orig.ident",feature2 = "nFeature_RNA")+
  geom_hline(yintercept = c(200,5000),linetype = 2 )
plot1 / plot2

nFeature_RNA <- Liver_Non_ECs@meta.data$nFeature_RNA

nFeature_RNA <- as.data.frame(nFeature_RNA)

# Produce Density plot

nFeature_RNA %>%
  ggplot(aes(x=nFeature_RNA)) +
  geom_density(fill='#C8DCE7', size = 1, alpha = 0.7) +
  scale_x_continuous() +
  ggtitle("Genes Detected Density Plot")+
  ylab("density") +
  xlab("detected")+
  geom_vline(xintercept = 600, color = "red", size = 1)

# Produce Cummulative plot

ggplot(nFeature_RNA, aes(nFeature_RNA)) + stat_ecdf(geom = "point", color = "#F8756C")+
  scale_x_continuous() +
  ggtitle("Genes Detected Cummulative Plot")+
  ylab("CumSum_Features") +
  xlab("total_Features")+
  geom_vline(xintercept = 600, color = "red", size = 1)

# Going over mitochondrial genes

percent_mt <- Liver_Non_ECs@meta.data$percent.mt
percent_mt <- as.data.frame(percent_mt)

percent_mt %>%
  ggplot(aes(x=percent_mt)) +
  geom_density(fill='#C8DCE7', size = 1, alpha = 0.7) +
  scale_x_continuous() +
  scale_y_discrete()+
  ggtitle("MT Content Density Plot")+
  ylab("density") +
  xlab("subset_MT_percent")+
  ylim(-0.4, 0.4)+
  geom_vline(xintercept = 25, color = "red", size = 1)+
  geom_jitter(aes(x = percent_mt, y = 0), height = 0.1, alpha = 0.01)

# Ribosomal genes

Liver_Non_ECs <- PercentageFeatureSet(Liver_Non_ECs, "^Rp[sl]", col.name = "percent_ribo")

# Hemoglobin genes

Liver_Non_ECs <- PercentageFeatureSet(Liver_Non_ECs, "^Hb[^(p)]", col.name = "percent_hb")

# Violin Plots

feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent_ribo", "percent_hb")
VlnPlot(Liver_Non_ECs, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) +
  NoLegend()

# Clear up the plots

Liver_Non_ECs <- subset(Liver_Non_ECs,
                        subset = nFeature_RNA > 500 & nFeature_RNA < 6000 &
                          nCount_RNA > 2000 & nCount_RNA < 30000  & percent.mt < 5 & percent_ribo < 35 & percent_hb < 1)


table(Liver_Non_ECs@meta.data$orig.ident) 

# there should be 986 cells for Liver_Control_NonECs and 2644 cells for Liver_Dll4LOF_NonECs

Liver_Non_ECs <- RenameIdents(Liver_Non_ECs, 'Liver_Control_NonECs' = 'Control', 'Liver_Dll4LOF_NonECs' = 'Dll4 LOF')

Liver_Non_ECs@meta.data$Condition <- Liver_Non_ECs@active.ident

Idents(Liver_Non_ECs) <- "Condition"

#  Liver_Non_ECs <- subset(Liver_Non_ECs, downsample = 1000)

table(Liver_Non_ECs@meta.data$Condition)


# It seems I need to go to the specific folders for each CMO section in the cellranger multi output

Liver_Non_ECs <- NormalizeData(Liver_Non_ECs)
Liver_Non_ECs <- FindVariableFeatures(Liver_Non_ECs)
Liver_Non_ECs <- ScaleData(Liver_Non_ECs)
Liver_Non_ECs <- RunPCA(Liver_Non_ECs, verbose = FALSE)
Liver_Non_ECs <- FindNeighbors(Liver_Non_ECs, dims = 1:30)
Liver_Non_ECs <- FindClusters(Liver_Non_ECs, resolution = 0.3, verbose = FALSE)
Liver_Non_ECs <- RunUMAP(Liver_Non_ECs, dims = 1:30)
DimPlot(Liver_Non_ECs, label = TRUE, split.by = "Condition") 

# Producing a UMAP

UMAP_Conditions_Global(Liver_Non_ECs, 2)

DimPlot(Liver_Non_ECs, label = T)

Idents(Liver_Non_ECs) <- Liver_Non_ECs@meta.data$RNA_snn_res.0.8

p1 <- DimPlot(Liver_Non_ECs, pt.size = 1.2, group.by = "RNA_snn_res.0.3", split.by = "Condition", combine = T)

p1

p21 <- p1+
  # facet_grid(~Condition, labeller = as_labeller(new_labels, default = label_parsed))+
  theme(strip.text.x = element_text(size = 24), plot.title = element_blank())+NoLegend()+xlab("UMAP_1")+ylab("UMAP_2")

p22 <- DimPlot(Liver_Non_ECs, group.by = "RNA_snn_res.0.3", pt.size = 1.2, label = T, label.box = T, repel = T)+theme(legend.text = element_text(size = 14), axis.title.x = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), plot.title = element_blank())

p21
p22

p23 <- list(p21, p22)

table(Liver_Non_ECs@meta.data$Condition)

design <- c(patchwork::area(1, 1, 1, 2), patchwork::area(1, 3, 1, 3.5))

p24 <- Reduce( `+`,  p23)+patchwork::plot_layout(design = design)


cairo_pdf("Plots/Liver_Non_ECS/UMAP_Condition_Split.pdf",  width = 20, height = 6, family = "Arial")
p24
dev.off()


##############################Label Transfer#############################################



sce <- as.SingleCellExperiment(Liver_Non_ECs, assay = "RNA")

# Using built-in references

hpca.se <- celldex::HumanPrimaryCellAtlasData()

hpca.se@NAMES <- tolower(hpca.se@NAMES) %>% str_to_title()

mmdb <- celldex::MouseRNAseqData()

mmdb$label.main

bpen <- celldex::BlueprintEncodeData()

bpen@NAMES <- tolower(bpen@NAMES) %>% str_to_title()

pred.Liver_Heart.hs <- SingleR(test = sce, ref = hpca.se, assay.type.test=1, 
                               labels = hpca.se$label.main)

table(pred.Liver_Heart.hs@listData$first.labels)
table(pred.Liver_Heart.hs@listData$labels)
table(pred.Liver_Heart.hs@listData$pruned.labels)


pred.Liver_Heart.mmdb <- SingleR(test = sce, ref = mmdb, assay.type.test=1, 
                                 labels = mmdb$label.main)

pred.Liver_Heart.bpen <- SingleR(test = sce, ref = bpen, assay.type.test=1, 
                                 labels = bpen$label.main)

pred.Liver_Heart <- pred.Liver_Heart.mmdb

pred.Liver_Heart.hs@listData$labels -> pred.Liver_Heart@listData$labels.hpca
pred.Liver_Heart.bpen@listData$labels -> pred.Liver_Heart@listData$labels.bpen

########################################################################################

# Annotations diagnostics

plotScoreHeatmap(pred.Liver_Heart,grid.vars = list())

plotScoreHeatmap(pred.Liver_Heart.bpen)

# plotHeatmap(pred.Liver_Heart, silent=TRUE, 
# order_columns_by="labels", features="labels")

plotDeltaDistribution(pred.Liver_Heart.mmdb, ncol = 5)

summary(is.na(pred.Embryo$pruned.labels))

all.markers <- metadata(pred.Liver_Heart.mmdb)$de.genes
sce$labels <- pred.Liver_Heart.mmdb$labels

# Add metadata EASY

Liver_Non_ECs@meta.data$mmdb <- pred.Liver_Heart.mmdb$labels

Liver_Non_ECs@meta.data$hpca <- pred.Liver_Heart.hs$labels

Liver_Non_ECs@meta.data$bpen <- pred.Liver_Heart.bpen$labels

DimPlot(Liver_Non_ECs, group.by = "hpca", pt.size = 1.2, label = T, repel = T, label.box = T)+ NoLegend()

DimPlot(Liver_Non_ECs, group.by = "bpen", pt.size = 1.2, label = T, repel = T, label.box = T)+ NoLegend()

DimPlot(Liver_Non_ECs, group.by = "mmdb", pt.size = 1.2, label = T, repel = T, label.box = T)+ NoLegend()


#####################################  Mapping using Ben Moshe et al dataset  ############################

# Paper link: https://www.sciencedirect.com/science/article/pii/S1934590922001618?via%3Dihub

BenMoshe <- readRDS("rds/BenMoshe.rds")

BenMoshe.Anchors <- FindTransferAnchors(reference = BenMoshe, query = Liver_Non_ECs,
                                        dims = 1:30, reference.reduction = "pca")

predictions <- TransferData(anchorset = BenMoshe.Anchors, refdata = BenMoshe@meta.data$cell_type,
                            dims = 1:30)

Liver_Non_ECs <- AddMetaData(Liver_Non_ECs, metadata = predictions)

table(Liver_Non_ECs$predicted.id)

DimPlot(Liver_Non_ECs, group.by = "predicted.id", label = T, label.box = T)

###############################################################################################################################

# UMAP of al labels and conditions

p1 <- DimPlot(Liver_Non_ECs, pt.size = 1.2, group.by = "hpca", split.by = "Condition")+ NoLegend()

p21 <- p1+
  theme(strip.text.x = element_text(size = 24), plot.title = element_blank())+NoLegend()+xlab("UMAP_1")+ylab("UMAP_2") 

p22 <- DimPlot(Liver_Non_ECs, pt.size = 1.2, group.by = "hpca", label = T, repel = T, label.box = T)+theme(legend.text = element_text(size = 14), axis.title.x = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust = 0.5, size = 24))+ ggtitle("HPCA")

p21
p22

p23 <- DimPlot(Liver_Non_ECs, pt.size = 1.2, group.by = "bpen", split.by = "Condition")+ NoLegend()+
  theme(strip.text.x = element_text(size = 24), plot.title = element_blank())+NoLegend()+xlab("UMAP_1")+ylab("UMAP_2")

p24 <- DimPlot(Liver_Non_ECs, pt.size = 1.2, group.by = "bpen", label = T, repel = T, label.box = T)+theme(legend.text = element_text(size = 14), axis.title.x = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust = 0.5, size = 24))+ ggtitle("BPEN")


p25 <- DimPlot(Liver_Non_ECs, pt.size = 1.2, group.by = "mmdb", split.by = "Condition")+ NoLegend()+
  theme(strip.text.x = element_text(size = 24), plot.title = element_blank())+NoLegend()+xlab("UMAP_1")+ylab("UMAP_2")

p26 <- DimPlot(Liver_Non_ECs, pt.size = 1.2, group.by = "mmdb", label = T, repel = T, label.box = T)+theme(legend.text = element_text(size = 14), axis.title.x = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust = 0.5, size = 24))+ ggtitle("Mouse Database")

p27 <- DimPlot(Liver_Non_ECs, pt.size = 1.2, group.by = "predicted.id", split.by = "Condition")+ NoLegend()+
  theme(strip.text.x = element_text(size = 24), plot.title = element_blank())+NoLegend()+xlab("UMAP_1")+ylab("UMAP_2")

p28 <- DimPlot(Liver_Non_ECs, pt.size = 1.2, group.by = "predicted.id", label = T, repel = T, label.box = T)+theme(legend.text = element_text(size = 14), axis.title.x = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust = 0.5, size = 24))+ ggtitle("Ben Moshe et al dataset")



p29 <- list(p21, p22, p23, p24, p25, p26, p27, p28)

design <- c(patchwork::area(1, 1, 1, 2), patchwork::area(1, 3, 1, 3.5), patchwork::area(2, 1, 2, 2), patchwork::area(2, 3, 2, 3.5), patchwork::area(3, 1, 3, 2), patchwork::area(3, 3, 3, 3.5), patchwork::area(4, 1, 4, 2), patchwork::area(4, 3, 4, 3.5))

p30 <- Reduce( `+`,  p29)+patchwork::plot_layout(design = design)

p30

cairo_pdf("Plots/Suppl_Figure_6/UMAP_LabelTransfer.All_Cells.pdf",  width = 24, height = 24, family = "Arial")
p30
dev.off()

###############################################################################################################################

# Trimmming clusters

# Finishing touches on the clustering, making everything less disgregated. This step is Manually curated

p1 <- DimPlot(Liver_Non_ECs_All_Cells, pt.size = 1.2, label = T, label.box = T, repel = T, group.by = "FinalClustering")+theme(legend.text = element_text(size = 14), axis.title.x = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), plot.title = element_blank())

Idents(Liver_Non_ECs_All_Cells) <- "FinalClustering"

p1

new_B <- CellSelector(p1)

Idents(Liver_Non_ECs_All_Cells, cells = new_B) <- "B"

new_KC <- CellSelector(p1)

Idents(Liver_Non_ECs_All_Cells, cells = new_KC) <- "KC"

new_Chol <- CellSelector(p1)

Idents(Liver_Non_ECs_All_Cells, cells = new_Chol) <- "Chol"

new_pDC <- CellSelector(p1)

Idents(Liver_Non_ECs_All_Cells, cells = new_pDC) <- "pDC"

new_TNK <- CellSelector(p1)

Idents(Liver_Non_ECs_All_Cells, cells = new_TNK) <- "T+NK"

new_HSC <- CellSelector(p1)

Idents(Liver_Non_ECs_All_Cells, cells = new_HSC) <- "HSC"

new_Hep <- CellSelector(p1)

Idents(Liver_Non_ECs_All_Cells, cells = new_Hep) <- "Hep"

new_Monocytes <- CellSelector(p1)

Idents(Liver_Non_ECs_All_Cells, cells = new_Monocytes) <- "Monocytes"

new_Endo <- CellSelector(p1)

Idents(Liver_Non_ECs_All_Cells, cells = new_Endo) <- "Endo"

neutrophils <- CellSelector(p1)

Idents(Liver_Non_ECs_All_Cells, cells = neutrophils) <- "Granulocytes"

table(Liver_Non_ECs_All_Cells@active.ident)

table(Liver_Non_ECs_All_Cells@meta.data$ZhuClustering)

Liver_Non_ECs_All_Cells@active.ident -> Liver_Non_ECs_All_Cells@meta.data$FinalClustering

levels(Liver_Non_ECs_All_Cells) <- c("Endo", "Hep", "Chol", "HSC", "Macrophages", "KC", "Monocytes", "cDC", "pDC", "T+NK", "B", "Granulocytes")

table(Liver_Non_ECs_All_Cells@meta.data$FinalClustering)

DimPlot(Liver_Non_ECs_All_Cells)

DimPlot(Liver_Non_ECs_All_Cells, split.by = "FinalClustering")

saveRDS(Liver_Non_ECs_All_Cells, "rds/Liver_Non_ECs_All_Cells.rds")
