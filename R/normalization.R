library(scater)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(RColorBrewer)
library(viridis)
library(edgeR)
library(scran)
options(stringsAsFactors=FALSE)
library(scuttle)

#1 load the data
selected_impute_sce <- readRDS(file.path("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data","selected_impute_sce.rds"))
selected_impute_sce <- selected_impute_sce[rowData(selected_impute_sce)$metabolic,]
selected_sce <- readRDS(file.path("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data","sce_all.rds"))
cell_types <- unique(selected_sce$cellType)

## choose the genes with higher expression and low-dropout rate 
dropout_cutoff <- 0.6
gene_select_mat <- matrix(FALSE,nrow=nrow(selected_impute_sce),
                          ncol=length(cell_types),
                          dimnames = list(rownames(selected_impute_sce),cell_types))
for(c in cell_types){
  each_sce <- selected_impute_sce[,selected_impute_sce$cellType == c]
  each_exp <- assay(each_sce,"exprs")
  dropout_rate <- apply(each_exp,1, function(x) sum(x>0)/ncol(each_exp))
  select <- dropout_rate >= dropout_cutoff
  gene_select_mat[select,c] <- TRUE
}
print("the number of genes selected:")
print(sum(rowSums(gene_select_mat) >= length(cell_types)))
low_dropout_genes <- rownames(gene_select_mat)[rowSums(gene_select_mat) >= length(cell_types)]

##########################################################################################
## different normalization methods:
##########################################################################################

# prepare the gene length file
# subset the tpm data
all_gene_lengths <- read.table("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data/gene_length.txt",sep="\t",header=F,row.names=1)
common_genes <- intersect(rownames(selected_impute_sce), rownames(all_gene_lengths))
selected_impute_sce <- selected_impute_sce[common_genes, ]

# subset genelength data
genelen <- all_gene_lengths[common_genes, ]
genelen <- as.numeric(as.vector(genelen))

#1. up-quantile normalization
selected_impute_tpm <- tpm(selected_impute_sce)
selected_impute_counts <- sweep(selected_impute_tpm, 1, genelen, FUN = "*")

median.sf <- calcNormFactors(selected_impute_counts[low_dropout_genes,],method="upperquartile",p=0.75)
selected_impute_tpm_norm <- t(t(selected_impute_tpm) / median.sf)
selected_impute_exp_norm <- log2(selected_impute_tpm_norm + 1)
#save
saveRDS(selected_impute_tpm_norm,file.path("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data","UpperQuartile_tpm.rds"))


#2.DESeq2 normalization
# get the function from "https://github.com/mikelove/DESeq2/blob/master/R/core.R"
source("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data/utils.R")
deseq2_sf <- estimateSizeFactorsForMatrix(selected_impute_counts[low_dropout_genes,])
selected_impute_tpm_norm <- t(t(selected_impute_tpm) / deseq2_sf)
selected_impute_exp_norm <- log2(selected_impute_tpm_norm + 1)
#save
saveRDS(selected_impute_tpm_norm,file.path("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data","RLE_tpm.rds"))

#3. edgeR (TMM)
edgeR_sf <- calcNormFactors(selected_impute_counts[low_dropout_genes,],method="TMM")
selected_impute_tpm_norm <- t(t(selected_impute_tpm) / edgeR_sf)
selected_impute_exp_norm <- log2(selected_impute_tpm_norm + 1)
#save
saveRDS(selected_impute_tpm_norm,file.path("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data","TMM_tpm.rds"))


#4. scran, deconvolution
scran.sf <- scuttle::pooledSizeFactors(selected_impute_counts[low_dropout_genes,],clusters=selected_impute_sce$cellType)
summary(scran.sf)
selected_impute_tpm_norm <- t(t(selected_impute_tpm) / scran.sf)
selected_impute_exp_norm <- log2(selected_impute_tpm_norm+1)
#save
saveRDS(selected_impute_tpm_norm,file.path("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data","Deconvolution_tpm.rds"))


###Evaluation of the normalization methods:
all_cell_type <- as.vector(selected_sce$cellType)
for(m in c("RLE","TMM","UpperQuartile","Deconvolution"))
{
  rds_file <- file.path("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data",paste0(m,"_tpm.rds"))
  norm_tpm <- readRDS(rds_file) 
  
  ##### check the ratio distribution
  low_dropout_genes_tpm <- norm_tpm[low_dropout_genes, ]
  low_dropout_genes_tpm_mean <- apply(low_dropout_genes_tpm, 1, function(x) by(x, all_cell_type, mean))
  low_dropout_genes_tpm_ratio <- t(low_dropout_genes_tpm_mean) / colMeans(low_dropout_genes_tpm_mean)
  dat <- reshape2::melt(low_dropout_genes_tpm_ratio)
  p <- ggplot(dat,aes(x=Var2,y=value)) +
    geom_boxplot(outlier.alpha=0.1)+ theme_classic() + 
    ylab("expression ratio") + xlab("") +    
    theme(axis.text.x = element_text(angle=45,hjust=1))
  ggsave(file.path("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data",paste0(m,"_ratio_distribution.pdf")),p,width=3.5,height=2.5)
}