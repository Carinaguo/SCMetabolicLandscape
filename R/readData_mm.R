source("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data/utils.R")
library(Seurat)
library(SeuratData)
library(scater)
library(stringr)
options(stringsAsFactors=FALSE)
library(reshape2)
library(plyr)

load("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data/MMNeu_2023.RData")

## extract raw_count from seurat object
seu_obj <- mm
seu_obj_m <- as.matrix(Seurat::GetAssay(seu_obj, "RNA")@data)

met_cell_type <- as.matrix(Seurat::FetchData(object = seu_obj, vars = "new_cluster", clean = FALSE))
met_cell_type[is.na(met_cell_type)] <- "Unknown"
cell_type <- as.character(met_cell_type[,1])
col_data <- data.frame(cellType=cell_type, row.names=colnames(seu_obj_m))

# marker the metabolic genes
pathways <- gmtPathways("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data/KEGG_metabolism.gmt")
metabolics <- unique(as.vector(unname(unlist(pathways))))
row_data <- data.frame(metabolic=rep(FALSE,nrow(seu_obj_m)),row.names = rownames(seu_obj_m))
row_data[rownames(row_data)%in%metabolics,"metabolic"]=TRUE # Total number of metabolic genes: 1201

# build scater object
raw_tpm <- seu_obj_m

sce <- SingleCellExperiment(
  assays = list(tpm=raw_tpm, exprs=seu_obj_m),
  colData = col_data,
  rowData = row_data
)

saveRDS(sce,file.path("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data","sce_all.rds"))

