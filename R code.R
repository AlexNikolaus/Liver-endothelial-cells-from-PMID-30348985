library(Seurat)
library(tidyverse)
library(SCINA)

#I have separately preprocessed every dataset to remove the ENSG gene names
load('D:/Endothelial scRNA-Seq project/PMID 30348985 liver/GSE115469/SRA716608_SRS3391633.sparse.RData')
Liver5 <- sm
Liver5.rownames <- rownames(Liver5)
Liver5.rownames <- gsub('_ENSG.*', '', Liver5.rownames)
rownames(Liver5) <- Liver5.rownames
rm(sm, Liver5.rownames)

load('D:/Endothelial scRNA-Seq project/PMID 30348985 liver/GSE115469/SRA716608_SRS3391629.sparse.RData')
Liver1 <- sm
Liver1.rownames <- rownames(Liver1)
Liver1.rownames <- gsub('_ENSG.*', '', Liver1.rownames)
rownames(Liver1) <- Liver1.rownames
rm(sm, Liver1.rownames)

load('D:/Endothelial scRNA-Seq project/PMID 30348985 liver/GSE115469/SRA716608_SRS3391630.sparse.RData')
Liver2 <- sm
Liver2.rownames <- rownames(Liver2)
Liver2.rownames <- gsub('_ENSG.*', '', Liver2.rownames)
rownames(Liver2) <- Liver2.rownames
rm(sm, Liver2.rownames)

load('D:/Endothelial scRNA-Seq project/PMID 30348985 liver/GSE115469/SRA716608_SRS3391631.sparse.RData')
Liver3 <- sm
Liver3.rownames <- rownames(Liver3)
Liver3.rownames <- gsub('_ENSG.*', '', Liver3.rownames)
rownames(Liver3) <- Liver3.rownames
rm(sm, Liver3.rownames)

load('D:/Endothelial scRNA-Seq project/PMID 30348985 liver/GSE115469/SRA716608_SRS3391632.sparse.RData')
Liver4 <- sm
Liver4.rownames <- rownames(Liver4)
Liver4.rownames <- gsub('_ENSG.*', '', Liver4.rownames)
rownames(Liver4) <- Liver4.rownames
rm(sm, Liver4.rownames)

#Create Seurat objects
Liver1 <- CreateSeuratObject(Liver1, project = 'Liver1', min.cells = 3, min.features = 200)
Liver2 <- CreateSeuratObject(Liver2, project = 'Liver2', min.cells = 3, min.features = 200)
Liver3 <- CreateSeuratObject(Liver3, project = 'Liver3', min.cells = 3, min.features = 200)
Liver4 <- CreateSeuratObject(Liver4, project = 'Liver4', min.cells = 3, min.features = 200)
Liver5 <- CreateSeuratObject(Liver5, project = 'Liver5', min.cells = 3, min.features = 200)

#Standard Seurat Workflow
Liver1[["percent.mito"]] <- PercentageFeatureSet(Liver1, pattern = "^MT-")
Liver2[["percent.mito"]] <- PercentageFeatureSet(Liver2, pattern = "^MT-")
Liver3[["percent.mito"]] <- PercentageFeatureSet(Liver3, pattern = "^MT-")
Liver4[["percent.mito"]] <- PercentageFeatureSet(Liver4, pattern = "^MT-")
Liver5[["percent.mito"]] <- PercentageFeatureSet(Liver5, pattern = "^MT-")

Liver1 <- subset(Liver1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mito < 25)
Liver2 <- subset(Liver2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mito < 25)
Liver3 <- subset(Liver3, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mito < 25)
Liver4 <- subset(Liver4, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mito < 50)
Liver5 <- subset(Liver5, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mito < 30)

Liver <- merge(Liver1, y = c(Liver2, Liver3, Liver4, Liver5), add.cell.ids = c('GSM3178782','GSM3178783', 'GSM3178784', 'GSM3178785', 'GSM3178786'))

Liver <- NormalizeData(Liver, normalization.method = "LogNormalize", scale.factor = 10000)
Liver <- FindVariableFeatures(Liver, selection.method = "vst", nfeatures = 2000)
Liver <- ScaleData(Liver)
Liver <- RunPCA(Liver, features = VariableFeatures(object = Liver))
Liver <- FindNeighbors(Liver, dims = 1:10)
Liver <- FindClusters(Liver, resolution = 0.01)
Liver <- RunUMAP(Liver, dims = 1:10)
DimPlot(Liver, reduction = "umap", label = TRUE)
VlnPlot(Liver, features = c('CLDN5', 'VWF', 'FLT1', 'PECAM1'), ncol = 2)
Liver <- subset(Liver, idents = c('3'))

#List of endothelial marker genes
Endothelial_Markers <- list(c('CLDN5', 'VWF', 'FLT1', 'PECAM1'))
names(Endothelial_Markers) <- 'Endothelial cells'

#Run SCINA to annotate endothelial cells
SCINA_results <- SCINA(Liver@assays$RNA@data,
                       Endothelial_Markers,
                       max_iter = 2000, 
                       convergence_n = 100, 
                       convergence_rate = 0.999, 
                       sensitivity_cutoff = 0.9, 
                       rm_overlap=FALSE, 
                       allow_unknown=TRUE)
Liver$cell_labels <- SCINA_results$cell_labels
DimPlot(Liver,reduction = "umap", pt.size = 1, label = TRUE, group.by = 'cell_labels')

#Subset for endothelial and write object
Liver <- subset(Liver, cell_labels == 'Endothelial cells')
Liver$organ <- 'Liver'
write_rds(Liver, 'Liver.rds')
