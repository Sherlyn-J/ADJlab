## load libraries
library(DEGreport)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)

library(org.Hs.eg.db)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(AnnotationHub)
library(ensembldb)
library(tidyverse)
library(annotables)

library(Cairo)

## load count data and sample information

#  files
count_file <- "D:/Ito-data-plots-Michal-paper-Feb21/A_HDT_vs_control.csv"
# smple_info <- "D:/Ito-data-plots-Michal-paper-Feb21/B_coldata.txt"

## count data
cts <- read.csv(count_file,header=TRUE,sep=",")
cts <- cts[!duplicated(cts$NAME),]
row.names(cts) <- cts$NAME
## sample information/metadata
# coldata <- read.table(smple_info,sep="\t",row.names=1)

# immgenes <- as.vector(scan("D:/Ito-data-plots-Michal-paper-Feb21/Mouse_immune_genes.txt",what=character(),sep="\n"))
# cts <- cts %>%
#   dplyr::filter(NAME %in% immgenes)

for(i in 2:nrow(cts))
{
  
  group1 <- cts[i,grep('^Vehicle',names(cts))]
  group2 <- cts[i,grep('^HDT',names(cts))]
  
  group1_avg <- sum(group1)/length(group1)
  group2_avg <- sum(group2)/length(group2)
  
  cts[i, "C_avg"] <- group1_avg
  cts[i, "D_avg"] <- group2_avg
  
  logfc <- log2(group1_avg/group2_avg)
  
  
  cts[i, "LogFC"] <-  logfc
  
  
  temp_data1 <- data.frame(values=c(group1,group2),vars = rep(c("C","D"), times = c(length(group1),length(group2))))
  vars1 = c(rep(c("C","D"), times = c(length(group1),length(group2))))
  values1=c(t(group1),t(group2))
  
  x <- t.test(values1 ~ vars1, data = temp_data1)
  pv <- c(x$p.value)
  
  y <- wilcox.test(values1 ~ vars1, data = temp_data1)
  
  wilcox <- c(y$p.value)
  
  cts[i, "TTEST"] <-  pv
  cts[i, "WilcoxTest"] <-  wilcox
  
  
}

cts <- na.omit(cts)
cts <- cts %>%
  dplyr::filter( LogFC < -1.5 | LogFC > 1.5)

write.table(cts,"D:/Ito-data-plots-Michal-paper-Feb21/A_DEG.txt", sep="\t", quote=FALSE, row.name=FALSE)
