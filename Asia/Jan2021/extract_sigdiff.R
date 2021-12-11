library(DESeq2)
library(pheatmap)
library(RColorBrewer)

#combined raw count matrix for all 4 cell lines
cts <- as.matrix(read.table("D:/Asia_Jan2021/ctmatrix.txt",header=TRUE,sep="\t",row.names="X"))
#sample information
coldata <- as.matrix(read.table("D:/Asia_Jan2021/coldata.txt",header=TRUE,sep="\t",row.names=1))

#initialize DESeq obj
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design=formula(~Sex+Condition))
dds <- DESeq(dds)

#dds$Condition <- relevel(dds$Condition,ref="Parent")

#results
res <- results(dds, name="Sex_Male_vs_Female", alpha=0.05)
res <- results(dds, contrast=c("Sex","Male","Female"), alpha=0.05)

write.csv(res, "D:/Asia_Jan2021/sexdiff_all.csv")

#significant genes
res <- na.omit(res)
res <- res[ res$pvalue < 0.05, ]
res <- res[ res$padj < 0.05, ]
top <- rbind(res[res$log2FoldChange < -2,], res[res$log2FoldChange > 2,] )
top <- top[order(top$log2FoldChange),]
top <- rbind(head(top,30),tail(top,30))
genes <- rownames(top)
#genes <- c('ADGRB1', 'APOBEC3G', 'AXDND1', 'BSN', 'CCDC38', 'CHST3', 'CNN3', 'CYTOR', 'DNAJC22', 'EFCAB6', 'EFNA1', 'EPHA3', 'FNDC10', 'HAAO', 'HCG22', 'KATNBL1P6', 'LGI4', 'LINC00539', 'LINC00654', 'LINC00862', 'LINC01366', 'LOC105370401', 'LRRC73', 'MICB-DT', 'MIR1273C', 'MIR3127', 'MIR640', 'MIR6732', 'MIR6767', 'MIR6797', 'MMRN1', 'MTUS2-AS1', 'NRG4', 'PEG13', 'POU4F1', 'PSTPIP2', 'RAB34', 'SCARNA6', 'SERPINH1', 'SHC3', 'SIRPA', 'SLC7A2', 'SLC8A1', 'SNORA53', 'SNORD102', 'SNORD49A', 'SNORD7', 'SNORD89', 'SOAT2', 'SPATA32', 'SPEF1', 'SUMO1P1', 'TNFRSF14-AS1', 'WDR17', 'ZNF236-DT', 'ZNF396')

#get heatmap
dds <- estimateSizeFactors(dds)
nrmd <- counts(dds, normalized=TRUE)
nrmd <- nrmd[ genes, ]
nrmd <- na.omit(nrmd)
mat <- log2(nrmd+1)
mat <- as(mat, "matrix")

hplot <- pheatmap(mat, color=colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")))(100), kmeans_k=NA, breaks=NA, 
                  border_color="grey60", cellwidth=NA, cellheight=10, scale="none", legend=TRUE, legend_breaks=NA, 
                  legend_labels=NA, annotation_row=NA, annotation_col=NA, annotation=NA, annotation_colors=NA, 
                  annotation_legend=TRUE, annotation_names_row=TRUE, annotation_names_col=TRUE, drop_levels=TRUE, 
                  show_rownames=T, show_colnames=T, main=NA, fontsize=10, angle_col=c("270", "0", "45", "90", "315"), 
                  display_numbers=F, number_format="%.2f", number_color="grey30", gaps_row=NULL, gaps_col=NULL, 
                  labels_row=NULL, labels_col=NULL, filename="D:/Asia_Jan2021/heatmap_bysex.png", width=NA, height=NA, 
                  silent=FALSE, na_col="#DDDDDD" )
                  
#write sig results
write.csv(rbind(res[res$log2FoldChange < -2,], res[res$log2FoldChange > 2,]), "D:/Asia_Jan2021/sexdiff.csv")
