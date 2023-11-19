library(scImpute)
library(scater)

num_cores <- 4 #for windows it should be 1

selected_sce <- readRDS(file.path("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data","sce_Other_MatNeu.rds"))

#write the tpm matrix
selected_tpm <- tpm(selected_sce) 
labels <- selected_sce$cellType

write.csv(selected_tpm,file.path("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data","tumor.tpm"))

##prepare the gene length file
all_gene_lengths <- read.table("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data/gene_length.txt",sep="\t",header=F,row.names=1)
tmp <- intersect(rownames(all_gene_lengths),rownames(selected_tpm))
if (length(tmp) != nrow(selected_tpm)){
  warning("check the length file")
  print(setdiff(rownames(selected_tpm),rownames(all_gene_lengths)))
  q()
}
genelen <- all_gene_lengths[rownames(selected_tpm),]
genelen <- as.numeric(as.vector(genelen))
scimpute(file.path("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data","tumor.tpm"),infile="csv",outfile="csv",out_dir=file.path("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data","malignant_"),
         labeled=TRUE,labels=as.vector(labels),
         type="TPM",genelen=genelen,drop_thre=0.5,ncores=num_cores)

imputed_tpm <- read.csv(file.path("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data","malignant_scimpute_count.csv"),header=T,row.names=1)
tpm(selected_sce) <- data.matrix(imputed_tpm) 
assay(selected_sce,"exprs") <- data.matrix(log2(imputed_tpm + 1))

#save as sce
impute_tpm <- tpm(selected_sce)
impute_exprs <- assay(selected_sce, "exprs")
impute_tpm <- impute_tpm[,colnames(selected_sce)]
impute_exprs <- impute_exprs[,colnames(selected_sce)]

selected_impute_sce <- SingleCellExperiment(
  assays = list(tpm = impute_tpm, exprs=impute_exprs),
  colData = colData(selected_sce),
  rowData = rowData(selected_sce)
)
saveRDS(selected_impute_sce,file.path("/Users/zhiruiguo/UW-Madison/23fall/DinhLab/metabolic_lab_data","selected_impute_sce.rds"))