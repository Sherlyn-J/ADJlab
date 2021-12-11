## Mapping ENSEMBL gene IDs to HGNC gene names

library(biomaRt)
library(tibble)

df <- read.table("D:/Bryce/Everolimus-sensitivity/GSE155923_SCLC_BETi_mTORi_combo_Chen2020_normed_counts.tsv",header=TRUE,sep="\t")
df$ensgene <- as.character(df$ensgene)
genes <- df$ensgene

mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl",mart)

gene_IDs <- getBM(filters= "ensembl_gene_id",
                  attributes=c("ensembl_gene_id","hgnc_symbol","hgnc_id"),
                  values = genes, 
                  mart= mart )

#df <- add_column(df, HGNC = gene_IDs$hgnc_symbol, .after = 1)
write.table(gene_IDs,"D:/Bryce/Everolimus-sensitivity/HGNC_GS.txt", sep="\t",row.names=FALSE,quote=FALSE)
