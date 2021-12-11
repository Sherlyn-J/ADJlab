library(DESeq2)
library(pheatmap)
library(RColorBrewer)

#combined raw count matrix for all 4 cell lines
cts <- as.matrix(read.table("D:/Asia_Jan2021/ctmatrix-male.txt",header=TRUE,sep="\t",row.names="X"))
#sample information
coldata <- as.matrix(read.table("D:/Asia_Jan2021/coldata-male.txt",header=TRUE,sep="\t",row.names=1))

#initialize DESeq obj
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design=formula(~Condition))
dds <- DESeq(dds)

#results
res <- results(dds, name="Condition_Parent_vs_Knockout", alpha=0.05)
res <- results(dds, contrast=c("Condition","Parent","Knockout"), alpha=0.05)

write.csv(res, "D:/Asia_Jan2021/condition_male.csv")

#significant genes
res <- na.omit(res)
res <- res[ res$pvalue < 0.05, ]
fc2 <- res[ res$log2FoldChange > 2 | res$log2FoldChange < -2, ]
top <- fc2[order(fc2$pvalue),]
genes <- rownames(top)

#get heatmap
dds <- estimateSizeFactors(dds)
nrmd <- counts(dds, normalized=TRUE)
nrmd <- nrmd[ genes, ]
nrmd <- na.omit(nrmd)
mat <- log2(nrmd+1)
mat <- as(mat, "matrix")

write.csv(fc2, "D:/Asia_Jan2021/condition_male.csv")

hplot <- pheatmap(mat, color=colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")))(100), kmeans_k=NA, breaks=NA, 
                  border_color="grey60", cellwidth=NA, cellheight=10, scale="none", legend=TRUE, legend_breaks=NA, 
                  legend_labels=NA, annotation_row=NA, annotation_col=NA, annotation=NA, annotation_colors=NA, 
                  annotation_legend=TRUE, annotation_names_row=TRUE, annotation_names_col=TRUE, drop_levels=TRUE, 
                  show_rownames=T, show_colnames=T, main=NA, fontsize=10, angle_col=c("270", "0", "45", "90", "315"), 
                  display_numbers=F, number_format="%.2f", number_color="grey30", gaps_row=NULL, gaps_col=NULL, 
                  labels_row=NULL, labels_col=NULL, filename="D:/Asia_Jan2021/heatmap_bycondition-males.png", width=NA, 
                  height=NA, silent=FALSE, na_col="#DDDDDD" )
