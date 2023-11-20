source("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data/utils.R")
library(Seurat)
library(SeuratData)
library(scater)
library(stringr)
options(stringsAsFactors=FALSE)
library(reshape2)
library(plyr)

load("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data/MMNeu_2023.RData")

# extract raw_count from seurat object
# subset 1. new_cluster == Other_MatNeu
met_1 <- Seurat::FetchData(object = mm, vars = "new_cluster")
seu_obj <- mm[, which(met_1 == "Other_MatNeu")]
seu_obj_m <- as.matrix(Seurat::GetAssay(seu_obj, "SCT")@data)
# seu_obj_m <- as.matrix(Seurat::GetAssay(seu_obj, "RNA")@data)

# column data
met_cell_type <- as.matrix(Seurat::FetchData(object = seu_obj, vars = "hpca_fine", clean = FALSE))
met_cell_type[is.na(met_cell_type)] <- "Unknown"
cell_type <- as.character(met_cell_type[,1])
col_data <- data.frame(cellType=cell_type, row.names=colnames(seu_obj_m))

# marker the metabolic genes
pathways <- gmtPathways("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data/KEGG_metabolism.gmt")
metabolics <- unique(as.vector(unname(unlist(pathways))))
row_data <- data.frame(metabolic=rep(FALSE,nrow(seu_obj_m)),row.names = rownames(seu_obj_m))
row_data[rownames(row_data)%in%metabolics,"metabolic"]=TRUE

# build scater object
raw_tpm <- 2^ seu_obj_m - 1

sce <- SingleCellExperiment(
  assays = list(tpm=raw_tpm, exprs=seu_obj_m),
  colData = col_data,
  rowData = row_data
)

saveRDS(sce,file.path("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data","sce_Other_MatNeu.rds"))





