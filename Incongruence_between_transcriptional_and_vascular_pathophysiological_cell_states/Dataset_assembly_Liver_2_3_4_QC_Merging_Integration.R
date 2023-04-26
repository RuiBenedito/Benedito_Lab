########################################
# This file contains the QC and Integration of the Liver_2, Liver_3 and Liver_4 datasets using Liver_1 as the query dataset
# 
# Author:
#   Alvaro Regano aregano@cnic.es
#
########################################

# Starting analysis on the Liver sample

library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dittoSeq)
library(patchwork)

###############################

# EXAMPLE OF QC: used in all samples besides Liver_1

# Normalize data

Liver <- NormalizeData(Liver, 
                       normalization.method = "LogNormalize",
                       scale.factor = 10000)

# Identification of highly variable features (feature selection)

Liver <- FindVariableFeatures(Liver, 
                              selection.method = "vst",
                              nfeatures = 2000)
top10 <- head(VariableFeatures(Liver), 10)
top10

plot1 <- VariableFeaturePlot(Liver)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# Scaling the data

all.genes <- rownames(Liver)
Liver <- ScaleData(Liver, features = all.genes)


## Perform linear dimensional reduction

Liver <- RunPCA(Liver,
                features = VariableFeatures(object = Liver))


# Examine and visualize PCA results a few different ways

print(Liver[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Liver, dims = 1:2, reduction = "pca")
DimPlot(Liver,group.by = "orig.ident", reduction = "pca")
DimHeatmap(Liver, dims = 1:9, cells = 500, balanced = TRUE)

ElbowPlot(Liver, ndims = 30)


Liver <- FindNeighbors(Liver, dims = 1:30)
Liver <- FindClusters(Liver, resolution = 0.1)
Liver <- FindClusters(Liver, resolution = 0.2)
Liver <- FindClusters(Liver, resolution = 0.35)
Liver <- FindClusters(Liver, resolution = 0.5)


DimPlot(Liver)

#Run UMAP

Liver <- RunUMAP(Liver, dims = 1:30)

DimPlot(Liver, group.by = "RNA_snn_res.0.35")
DimPlot(Liver, group.by = "RNA_snn_res.0.35", split.by = "Condition")

table(Liver@meta.data$hash.ID)

# Subsetting by their different hashtags, each liver sample will slightly differ in this step
# Liver_1: BC303=Control	BC304=Dll4iDEC(Dll4KO)	BC305=RbpjiDEC(RbpjKO)	BC302=Notch1iDEC(NOTCH1KO)
# Liver_2: BC301=Control(WT)	BC302=DBZ(WT+DBZ)	BC305=Dll4iDEC+antiVEGF (D4KO_aVEGF)	BC306=Notch1/2/4iDEC	BC307=Dll4HETiDEC
# Liver_3: BC304=Jag1/Jag2/Dll1iDEC (Jag1/Jag2/Dll1KO(2w))
# Liver_4: BC305=Control + Vehicle	BC302=Dll4iDEC+Vehicle(D4KO_Veh)	BC303=Dll4iDEC+SL327(D4KO_SL327)	BC304=Dll4/MYCiDEC (Dll4/MycKO)
# Liver_5: CMO302=Control	CMO303=Dll4iDEC (Dll4LOF)

Liver_singlets <- subset(x = Liver, subset = (hash.ID == "HTO302" | hash.ID == "HTO307" | hash.ID == "HTO306" | hash.ID == "HTO301" | hash.ID == "HTO305"))


#  Now to take out Blood cells

FeaturePlot(Liver_singlets,
            features = c("PTPRC", # Leukocyte lineage
                         "CD14", # Macrophage
                         "CD3E", # T cell                        
                         "CD79A", # B cell
                         "PECAM1", # ECs
                         "CDH5") # ECs
)

Liver_singlets <- FindNeighbors(Liver_singlets, dims = 1:30)
Liver_singlets <- FindClusters(Liver_singlets, resolution = 1.5)

DimPlot(Liver_singlets, label = T, label.box = T, group.by = "RNA_snn_res.1.5")

# Detect clusters with high amounts of markers related to blood cells 


# Extract said clusters

VlnPlot(Liver_singlets, features =  "PTPRC",
        group.by =  "RNA_snn_res.1.5", 
        pt.size = 0.05 ) + theme(legend.position="none")

Liver <- Liver_singlets@meta.data$RNA_snn_res.1.5  %in%  as.character(c(0:17, 19, 24:26, 28))

p1 <- DimPlot(Liver_singlets, reduction = "umap", cells = Liver,
              group.by = "RNA_snn_res.1.5",
              label = T, pt.size = 0.5)
p1

x <- Liver_singlets@meta.data$RNA_snn_res.1.5  %in%  as.character(c(0:17, 19, 24:26, 28))
Liver_nonblood <- subset(Liver_singlets, cells = colnames(Liver_singlets)[x])
DefaultAssay(Liver_nonblood) <- "RNA"

Liver_nonblood <- FindVariableFeatures(Liver_nonblood, 
                                       selection.method = "vst",
                                       nfeatures = 2000)
Liver_nonblood <- RunPCA(Liver_nonblood,
                         features = VariableFeatures(object = Liver_nonblood))
Liver_nonblood <- FindNeighbors(Liver_nonblood, dims = 1:20)
Liver_nonblood <- FindClusters(Liver_nonblood, resolution = 0.1)
Liver_nonblood <- FindClusters(Liver_nonblood, resolution = 0.3)
Liver_nonblood <- FindClusters(Liver_nonblood, resolution = 0.05)

Liver_nonblood <- RunUMAP(Liver_nonblood, dims = 1:20)

DimPlot(Liver_nonblood)

####################################################

# DATASET MERGING

# Here I will gather all rds, with all conditions, and make 3 final objects:
# 1. All Conditions
# 2. Liver 2weeks deletion
# 3. Liver 4 days deletion
# In order to group them correctly, I will also downsample to a max # of cells of 1000, using a set seed, 42.

library(Seurat)

Liver_1 <- readRDS("")
Liver_2 <- readRDS("")
Liver3 <- readRDS("")
Liver_4 <- readRDS("")


table(Liver_1@meta.data$Condition)
table(Liver_4@meta.data$Condition)
table(Liver_2@meta.data$Condition)
table(Liver_3@meta.data$Condition)

# Liver_2

Idents(Liver_2) <- "Condition"

Liver_2 <- subset(Liver_2, downsample = 1000, seed = 42)

# Liver_3

Idents(Liver_3) <- "hash.ID"

Liver_3 <- subset(Liver_3, subset = hash.ID == "BC304", downsample = 1000, seed = 42)

condition <- c("Jag1/Jag2/Dll1KO")

names(condition) <- levels(Liver_3)

Liver_3 <- RenameIdents(Liver_3, condition)

Liver_3@active.ident -> Liver_3@meta.data$Condition

# Liver_4

Idents(Liver_4) <- "Condition"

Liver_4 <- subset(Liver_4, downsample = 1000, seed = 42)

# Liver_MultiplGroups

Idents(Liver_2) <- Liver_2@meta.data$Condition

Liver_2 <- subset(Liver_2, downsample = 1000, seed = 42)

# Merging datasets

Liver_A <- merge(Liver_2, Liver_4)
Liver_B <- merge(Liver_A, Liver_3)
Liver <- merge(Liver_B, Liver_1)

table(Liver@meta.data$Condition)

saveRDS(Liver, "rds/Groups/Liver.Figures.merge.All.Conditions.rds")

# Mapping them using the Liver_1 dataset produced by Carlos Torroja

# Read rds object with Ctl, Dll4KO, Notch1KO and RbpjKO where initial DimRed analysis was performed

Liver_1 <- readRDS("")

# Read rds object where all conditions have been merged to a normalized count of 1000 cells/condition

Liver <- readRDS("")

# Perform Intergration using Liver_1 as reference from the UMAP assembly performed by Carlos

Liver.query <- Liver
Liver.anchors <- FindTransferAnchors(reference = Liver_1, query = Liver.query,
                                     dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = Liver.anchors, refdata = Liver_1$NewClustering,
                            dims = 1:30)
Liver.query <- AddMetaData(Liver.query, metadata = predictions)

table(Liver.query@mata.data$predicted.id)

Liver.query@mata.data$predicted.id -> Liver.query@meta.data$FigClustering
