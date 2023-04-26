########################################
# This file contains the assembly, QC and UMAP Dimensionality reduction of the Liver_1 sample
# 
# Author:
#   Carlos Torroja, carlos.torroja@cnic.es
#
########################################


library(tidyverse)
library(scater)
library(Seurat)
library(biomaRt)
library(SingleR)

analysisPath <- "D:/LABS/LAB_RB/AReganoScript"
dataPath <- file.path(analysisPath,"data")

sampleMetadata <- data.frame(Sample=rep("Liver",4)
                             ,Condition=c("Control","Dll4KO","RBPJKO","NOTCH1KO")
                             ,Index=rep("All",4)
                             ,Genotype=c("WT","Dll4KO","RBPJKO","NOTCH1KO")
                             ,BarcodeID=c("BC303","BC304","BC305","BC302")
                             ,BarcodeSeq=c("CTTGCCGCATGTCAT","AAAGCATTCTTCACG","CTTTGTCTTTGTGAG","GGTCGAGAGCATTCA")
                             ,Tissue=rep("Liver",4)
                             ,Strain=rep("C57",4)
                             ,Assay_Type=rep("3_UMI",4)
                             ,Instrument=rep("10X",4)
                             ,Organism=rep("Mouse",4)
                             ,Platform=rep("ILLUMINA",2)
                             )

rownames(sampleMetadata) <- sampleMetadata$BarcodeID

myGeneDB <- useEnsembl(biomart = 'genes', 
                       dataset = "mmusculus_gene_ensembl",
                       version = "98")

genesMetadata <- getBM(attributes = c("external_gene_name"
                                      ,"mgi_symbol"
                                      ,"ensembl_gene_id","gene_biotype"
                                      ,"chromosome_name","start_position","end_position","strand"
                                      ,"description")
                       ,uniqueRows=T, mart = myGeneDB )

genesMetadata <- genesMetadata %>%
  mutate(external_gene_name = case_when(
    external_gene_name == "" ~ ensembl_gene_id,
    TRUE ~ external_gene_name
  )) %>%
  filter(!duplicated(genesMetadata$ensembl_gene_id)) %>%
  mutate(uniq_name = toupper(make.names(external_gene_name,unique = T)))

transGenes <- read_delim(file.path(dataPath,"mTmt_CRE.Final.gtf"),delim = "\t",col_names = F)
transGenes <- transGenes %>% filter(X3=="gene")
transGenes <- transGenes %>%
  mutate(X9=gsub("\\w+\\s","",X9)) %>%
  mutate(X9=gsub("\"","",X9)) %>%
  separate(X9,into = c("gene_id","gene_name","gene_source","gene_type"),sep = "; ") %>%
  mutate(description = "transgene", gene_type = sub(";","",gene_type)) %>% 
  dplyr::rename(chromosome_name = "X1", start_position = "X4", end_position = "X5", strand = "X7", external_gene_name = "gene_name", ensembl_gene_id = "gene_id", gene_biotype = "gene_type") %>% 
  mutate(mgi_symbol=external_gene_name,uniq_name=toupper(external_gene_name),strand=ifelse(strand=="+",1,-1)) %>%
  dplyr::select(external_gene_name,mgi_symbol,ensembl_gene_id,gene_biotype,chromosome_name,start_position,end_position,strand,description,uniq_name)

genesMetadata <- rbind(genesMetadata,transGenes)

rawSC <- Read10X(file.path(dataPath,"filtered_feature_bc_matrix"),gene.column = 1)

sce <- SingleCellExperiment(
  assays=list(counts=rawSC$`Gene Expression`),
  mainExpName = "RNA"
)

htoSCE <- SingleCellExperiment(
  assays=list(counts=rawSC$`Antibody Capture`)
)

altExps(sce) <- list(HTO=htoSCE)

genesMetadata <- genesMetadata %>%
  filter(ensembl_gene_id %in% rownames(sce))

rownames(sce) <- genesMetadata[match(rownames(sce),genesMetadata$ensembl_gene_id),"uniq_name"]

HBBGenes <- genesMetadata %>%
  filter(grepl("hemoglobin",description,ignore.case = T))

MTGenes <- genesMetadata %>%
  filter(chromosome_name == "MT")

TmtGenes <- genesMetadata %>%
  filter(ensembl_gene_id %in% c("mTmt_CRE_G1"))

phiYFPGenes <- genesMetadata %>%
  filter(ensembl_gene_id %in% c("PhiYFP_G1"))

sce <- logNormCounts(sce)

# Centered Log Ratio
clr <- function(x) {
  return(log1p(x = x / (expm1(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x)))))
}

assay(altExp(sce,"HTO"),"logcounts") <- clr(assay(altExp(sce,"HTO"),slot = "counts"))

sce <- addPerCellQC(
  sce
  ,subsets = list(
    MT = MTGenes$uniq_name
    ,HBB = HBBGenes$uniq_name
    ,Tmt = TmtGenes$uniq_name
    ,YFP = phiYFPGenes$uniq_name
  )
  ,percent.top = c(50,500)
  ,use_altexps = T
)

sce <- addPerFeatureQC(
  sce
)

altExp(sce,"HTO") <- addPerFeatureQC(
  altExp(sce,"HTO")
)
altExp(sce,"HTO") <- addPerCellQC(
  altExp(sce,"HTO")
)

colData(sce)$Cell_Fraction_Counts <- sce$sum/median(sce$sum)

colData(sce)$mTmt_CRE_counts <- assay(sce,"counts")["MTMT_CRE",]
colData(sce)$PhiYFP <- assay(sce,"counts")["PHIYFP",]

minDepth <- 1500
maxDepth <- 40000
minGenesDetected <- 600
maxMT <- 25
minCF <- 0.1
maxPct <- 65
maxHBB <- 0.1
minHTO <- 100

fcells <- colData(sce) %>% as.data.frame() %>%
  filter(sum > minDepth & sum < maxDepth) %>%
  filter(detected > minGenesDetected) %>%
  filter(subsets_MT_percent < maxMT) %>%
  filter(Cell_Fraction_Counts > minCF) %>%
  filter(percent.top_50 < maxPct) %>%
  filter(subsets_HBB_percent < maxHBB) %>%
  filter(altexps_HTO_sum > minHTO)

myGenes <- !rownames(sce) %in% c(MTGenes$uniq_name,
                                 HBBGenes$uniq_name,
                                 phiYFPGenes$uniq_name,
                                 TmtGenes$uniq_name)

seuratSC <- as.Seurat(sce[myGenes,rownames(fcells)]
                      , counts = "counts"
                      , data = "logcounts")

DefaultAssay(seuratSC) <- "HTO"

seuratSC <- seuratSC %>%
  NormalizeData(
    assay = "HTO",
    normalization.method = "CLR",
    margin = 1) %>%
  HTODemux(
    assay = "HTO",
    positive.quantile = 0.99)

seuratSC@meta.data <- seuratSC@meta.data %>%
  left_join(sampleMetadata,by = c("hash.ID"="BarcodeID"))
rownames(seuratSC@meta.data) <- colnames(seuratSC)

DefaultAssay(seuratSC) <- "RNA"

seuratSC <- seuratSC %>%
  NormalizeData(
    assay = "RNA",
    normalization.method = "LogNormalize"
  ) %>%
  FindVariableFeatures(
    assay = "RNA",
    nfeatures = 1000,
    selection.method = "vst") %>%
  ScaleData(
    features = rownames(.)
  ) %>%
  RunPCA() %>%
  FindNeighbors(
    dims = seq(10),
    reduction = "pca") %>%
  FindClusters(
    assay = "RNA",
    resolution = c(0.5)
    ) %>%
  RunUMAP(
    dims = seq(10),
    return.model = T
  )


hpca <- celldex::HumanPrimaryCellAtlasData()
bpen <- celldex::BlueprintEncodeData()

cellannot <- SingleR(
  clusters=seuratSC$seurat_clusters
  , test = Seurat::Assays(seuratSC,slot = "RNA")@data
  , ref = hpca
  , labels = hpca$label.fine
  , genes = "de"
  , quantile = 0.8
  , fine.tune = T
  , tune.thresh = 0.05
  , sd.thresh = 1
)
cellannot <- cellannot %>% as.data.frame %>%
  mutate(seurat_clusters=rownames(.)) %>%
  dplyr::select(seurat_clusters,labels) %>%
  dplyr::rename("hpca_labels" = labels)

seuratSC@meta.data <- seuratSC@meta.data %>%
  left_join(cellannot,by=c("seurat_clusters"))
rownames(seuratSC@meta.data) <- colnames(seuratSC)


cellannot <- SingleR(
  clusters=seuratSC$seurat_clusters
  , test = Seurat::Assays(seuratSC,slot = "RNA")@data
  , ref = bpen
  , labels = bpen$label.fine
  , genes = "de"
  , quantile = 0.8
  , fine.tune = T
  , tune.thresh = 0.05
  , sd.thresh = 1
)
cellannot <- cellannot %>% as.data.frame %>%
  mutate(seurat_clusters=rownames(.)) %>%
  dplyr::select(seurat_clusters,labels) %>%
  dplyr::rename("bpen_labels" = labels)

seuratSC@meta.data <- seuratSC@meta.data %>%
  left_join(cellannot,by=c("seurat_clusters"))
rownames(seuratSC@meta.data) <- colnames(seuratSC)

seuratSCEndoSinglet <- subset(seuratSC,HTO_classification.global == "Singlet" & bpen_labels == "Endothelial cells")

seuratSCEndoSinglet <- seuratSCEndoSinglet %>%
  DietSeurat(counts = TRUE, data = TRUE, scale.data = TRUE)

seuratSCEndoSinglet <- seuratSCEndoSinglet %>%
  FindVariableFeatures(
  assay = "RNA",
  selection.method = "vst") %>%
  RunPCA() %>%
  FindNeighbors(
    dims = seq(7),
    reduction = "pca") %>%
  FindClusters(
    assay = "RNA",
    resolution = c(0.3)
  ) %>%
  RunUMAP(
    dims = seq(7),
    seed=123456,
    return.model = T
  )
