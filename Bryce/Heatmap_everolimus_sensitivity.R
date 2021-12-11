library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(pheatmap)
Metabolomics <- read.csv("C:/Users/csislyn/Desktop/CCLE_metabolomics_20190502.csv", header=TRUE)
# Gene_dependency <- read.csv("D:/ADJ_Lab Data/Everolimus + Palbociclib Study/Depmap Analyses/Everolimus IC50/gene_effect (CERES).csv", header=TRUE)
Everolimus_sensitivity <- read.csv("C:/Users/csislyn/Desktop/primary-screen-replicate-collapsed-logfold-change (Everolimus only).csv", header=TRUE)
  
m1 <- merge(Metabolomics, Everolimus_sensitivity, by.x = "line", by.y = "line")

#drop line column
m1 <- m1[!(names(m1) %in% c('line')) ]

# arrange in ascending order by Everolimus_sensitivity
m1 <- m1[ order(m1$Everolimus.Sensitivity),]
m1 <- t(m1)
# scaled the matrix to check if the heamap captures any differences better?
m1 <- scale(m1)

### Set a color palette
heat_colors <- brewer.pal(9, "YlGnBu")
hmap<- pheatmap(m1,
                show_rownames = T,
                show_colnames = F,
                cluster_rows = F,
                cluster_cols = F,
                fontsize = 10, 
                scale = "row", 
                fontsize_row = 8, 
                height = 20)

CairoPNG(file="C:/Users/csislyn/Desktop/Ev-sensitivity-heatmap.png", units="in", width=30, height=40, pointsize=12, dpi=400)
print(hmap)
dev.off()
