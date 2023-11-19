library(scater)
library(reshape2)
options(stringsAsFactors=FALSE)

# Loading the data
selected_sce <- readRDS(file.path("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data","sce_Other_MatNeu.rds"))
selected_metabolic_sce <- selected_sce[rowData(selected_sce)$metabolic,]

# tsne
set.seed(12345)
library("Rtsne")
tsne_out <- Rtsne(t(assay(selected_metabolic_sce,"exprs")),initial_dims=20,perplexity=30,theta=0.0)
tmp <- data.frame(x=tsne_out$Y[,1],y=tsne_out$Y[,2],group = colData(selected_metabolic_sce)$cellType)
g <- ggplot(tmp) + geom_point(aes(x, y, colour = group), size = 1) +
  labs(x = "tSNE1",y = "tSNE2") +theme_bw() + ggtitle("Rtsne")
ggsave("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data/metabolic_tsne_cellTypeColor.pdf",g,width=4,height=3)