library(scater)
library(reshape2)
options(stringsAsFactors=FALSE)
library(scran)
library(uwot)

load("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data/MMNeu_2023.RData")

# Loading the data
# selected_sce <- readRDS(file.path("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data","sce_all.rds"))

# Use imputed data
selected_sce <- readRDS(file.path("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data","selected_impute_sce.rds"))
selected_metabolic_sce <- selected_sce[rowData(selected_sce)$metabolic,]

# tsne
set.seed(12345)
library("Rtsne")
tsne_out <- Rtsne(t(assay(selected_metabolic_sce,"exprs")),initial_dims=20,perplexity=10,theta=0.0)
tmp <- data.frame(x=tsne_out$Y[,1],y=tsne_out$Y[,2],group = colData(selected_metabolic_sce)$cellType)
g <- ggplot(tmp) + geom_point(aes(x, y, colour = group), size = 1) +
  labs(x = "tSNE1",y = "tSNE2") +theme_bw() + ggtitle("Rtsne")
ggsave("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data/metabolic_tsne_cellTypeColor_all_imputed.pdf",g,width=10,height=10)

# UMAP
set.seed(12345) 
pca_out <- runPCA(selected_metabolic_sce, ncomponents = 50, exprs_values = "exprs")
umap_out <- runUMAP(pca_out, dimred = "PCA")
umap_coords <- reducedDims(umap_out)$UMAP
cell_types <- colData(umap_out)$cellType
umap_df <- data.frame(UMAP1 = umap_coords[, 1], UMAP2 = umap_coords[, 2], CellType = cell_types)
g_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = CellType)) +
  geom_point(alpha = 0.8) +  # Adjust alpha for point transparency, if desired
  labs(x = "UMAP 1", y = "UMAP 2", color = "Cell Type") +
  theme_bw() + 
  ggtitle("UMAP")
ggsave("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data/metabolic_umap_cellTypeColor_all.pdf",g_umap,width=10,height=10)
