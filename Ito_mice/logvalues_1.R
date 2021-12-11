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
count_file <- "D:/Ito-data-plots-Michal-paper-Feb21/C_Kras_vs_WT.csv"

## count data
cts <- read.table(count_file,header=TRUE,sep=",")
cts <- cts[!duplicated(cts$NAME),]
row.names(cts) <- cts$NAME

## filter 0s

cts <- cts %>%
  dplyr::filter(cts$Tumor1!=0)

cts <- cts %>%
  dplyr::filter(cts$Tumor2!=0)

cts <- cts %>%
  dplyr::filter(cts$Tumor3!=0)

cts <- cts %>%
  dplyr::filter(cts$Control1!=0)

cts <- cts %>%
  dplyr::filter(cts$Control2!=0)

cts <- cts %>%
  dplyr::filter(cts$Control3!=0)

## Fold change calculation (FPKM to log2FC)
cts$log2FC <- log2(((cts$Tumor1+cts$Tumor2+cts$Tumor3)/3)/((cts$Control1+cts$Control2+cts$Control3)/3))
write.table(cts,"D:/Ito-data-plots-Michal-paper-Feb21/C_diff.csv",sep=",")

