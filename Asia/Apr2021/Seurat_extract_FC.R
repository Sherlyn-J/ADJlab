library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Cairo)

dlbc.d1 <- Read10X(data.dir="D:/Michal-Asia/raw/DLBCL1/")
d1 <- CreateSeuratObject(counts = dlbc.d1, project = "D1", min.cells = 3, min.features = 200)
d1@meta.data[['batch']] <- 1
d1@meta.data[['sex']] <- 1

dlbc.d2 <- Read10X(data.dir="D:/Michal-Asia/raw/DLBCL2/")
d2 <- CreateSeuratObject(counts = dlbc.d2, project = "D2", min.cells = 3, min.features = 200)
d2@meta.data[['batch']] <- 1
d2@meta.data[['sex']] <- 0

dlbc.d3 <- Read10X(data.dir="D:/Michal-Asia/raw/DLBCL3/")
d3 <- CreateSeuratObject(counts = dlbc.d3, project = "D3", min.cells = 3, min.features = 200)
d3@meta.data[['batch']] <- 2
d3@meta.data[['sex']] <- 1

dlbc.l1 <- Read10X(data.dir="D:/Michal-Asia/raw/rLN1/")
l1 <- CreateSeuratObject(counts = dlbc.l1, project = "L1", min.cells = 3, min.features = 200)
l1@meta.data[['batch']] <- 1
l1@meta.data[['sex']] <- 0

dlbc.l2 <- Read10X(data.dir="D:/Michal-Asia/raw/rLN2/")
l2 <- CreateSeuratObject(counts = dlbc.l2, project = "L2", min.cells = 3, min.features = 200)
l2@meta.data[['batch']] <- 2
l2@meta.data[['sex']] <- 0

dlbc.l3 <- Read10X(data.dir="D:/Michal-Asia/raw/rLN3/")
l3 <- CreateSeuratObject(counts = dlbc.l3, project = "L3", min.cells = 3, min.features = 200)
l3@meta.data[['batch']] <- 2
l3@meta.data[['sex']] <- 1

dlbc.f1 <- Read10X(data.dir="D:/Michal-Asia/raw/FL1/")
f1 <- CreateSeuratObject(counts = dlbc.f1, project = "F1", min.cells = 3, min.features = 200)
f1@meta.data[['batch']] <- 1
f1@meta.data[['sex']] <- 0

dlbc.f2 <- Read10X(data.dir="D:/Michal-Asia/raw/FL2/")
f2 <- CreateSeuratObject(counts = dlbc.f2, project = "F2", min.cells = 3, min.features = 200)
f2@meta.data[['batch']] <- 1
f2@meta.data[['sex']] <- 0

dlbc.f3 <- Read10X(data.dir="D:/Michal-Asia/raw/FL3/")
f3 <- CreateSeuratObject(counts = dlbc.f3, project = "F3", min.cells = 3, min.features = 200)
f3@meta.data[['batch']] <- 1
f3@meta.data[['sex']] <- 0

dlbc.f4 <- Read10X(data.dir="D:/Michal-Asia/raw/FL4/")
f4 <- CreateSeuratObject(counts = dlbc.f4, project = "F4", min.cells = 3, min.features = 200)
f4@meta.data[['batch']] <- 2
f4@meta.data[['sex']] <- 0

# batch 1 normalization
b1 <- merge(f1, y=f2, add.cell.ids = c("F1","F2"), project = "batch1")
b1 <- merge(b1, y=f3, add.cell.ids = c( "" ,"F3"), project = "batch1")
b1 <- merge(b1, y=d1, add.cell.ids = c( "" ,"D1"), project = "batch1")
b1 <- merge(b1, y=d2, add.cell.ids = c( "" ,"D2"), project = "batch1")
b1 <- merge(b1, y=l1, add.cell.ids = c( "" ,"L1"), project = "batch1")

# batch 2 normalization
b2 <- merge(d3, y=f4, add.cell.ids = c("D3","F4"), project = "batch2")
b2 <- merge(b2, y=l2, add.cell.ids = c( "" ,"L2"), project = "batch2")
b2 <- merge(b2, y=l3, add.cell.ids = c( "" ,"L3"), project = "batch2")

# female DLBC merge
fem_dlbc <- merge(d1, y=d3, add.cell.ids = c("D1","D3"), project = "female-DLBC")

# female only
fem_only <- merge(d1, y=d3, add.cell.ids = c("D1","D3"), project = "female-only")
fem_only <- merge(fem_only, y=l3, add.cell.ids = c("","L3"), project = "female-only")

# male only
male_only <- merge(d2, y=l1, add.cell.ids = c("D2","L1"), project = "male-only")
male_only <- merge(male_only, y=l2, add.cell.ids = c("","L2"), project = "male-only")

dlbc <- merge(b1,y=b2,add.cell.ids = c("B1","B2"), project = "DLBC" )

# QC
dlbc[["percent.mt"]] <- PercentageFeatureSet(dlbc, pattern = "^MT-")
dlbc <- subset(dlbc, subset = nFeature_RNA > 200 & percent.mt < 5)
dlbc <- NormalizeData(dlbc)

a3gdata <- GetAssayData(object = dlbc, assay = "RNA", slot = "data")["APOBEC3G",]
dlbc@meta.data['A3G_quartiles'] <- ntile(a3gdata,4)

# find markers
# for each subset
Idents(dlbc) <- dlbc@meta.data$orig.ident
f <- subset(dlbc,idents = "D1")
Idents(f) <- f@meta.data$A3G_quartiles
m <- FindMarkers(f,ident.1 = 2)
write.table(m,"D:/Michal-Asia/D1_A3G-Q1-vs-others.txt",quote=F,sep="\t")

Idents(dlbc) <- dlbc@meta.data$orig.ident
f <- subset(dlbc,idents = "D2")
Idents(f) <- f@meta.data$A3G_quartiles
m <- FindMarkers(f,ident.1 = 1)
write.table(m,"D:/Michal-Asia/D2_A3G-Q1-vs-others.txt",quote=F,sep="\t")

Idents(dlbc) <- dlbc@meta.data$orig.ident
f <- subset(dlbc,idents = "D3")
Idents(f) <- f@meta.data$A3G_quartiles
m <- FindMarkers(f,ident.1 = 1)
write.table(m,"D:/Michal-Asia/D3_A3G-Q1-vs-others.txt",quote=F,sep="\t")

Idents(dlbc) <- dlbc@meta.data$orig.ident
f <- subset(dlbc,idents = "L1")
Idents(f) <- f@meta.data$A3G_quartiles
m <- FindMarkers(f,ident.1 = 1)
write.table(m,"D:/Michal-Asia/L1_A3G-Q1-vs-others.txt",quote=F,sep="\t")

Idents(dlbc) <- dlbc@meta.data$orig.ident
f <- subset(dlbc,idents = "L2")
Idents(f) <- f@meta.data$A3G_quartiles
m <- FindMarkers(f,ident.1 = 1)
write.table(m,"D:/Michal-Asia/L2_A3G-Q1-vs-others.txt",quote=F,sep="\t")

Idents(dlbc) <- dlbc@meta.data$orig.ident
f <- subset(dlbc,idents = "L3")
Idents(f) <- f@meta.data$A3G_quartiles
m <- FindMarkers(f,ident.1 = 1)
write.table(m,"D:/Michal-Asia/L3_A3G-Q1-vs-others.txt",quote=F,sep="\t")

Idents(dlbc) <- dlbc@meta.data$orig.ident
f <- subset(dlbc,idents = "F1")
Idents(f) <- f@meta.data$A3G_quartiles
m <- FindMarkers(f,ident.1 = 1)
write.table(m,"D:/Michal-Asia/F1_A3G-Q1-vs-others.txt",quote=F,sep="\t")

Idents(dlbc) <- dlbc@meta.data$orig.ident
f <- subset(dlbc,idents = "F2")
Idents(f) <- f@meta.data$A3G_quartiles
m <- FindMarkers(f,ident.1 = 1)
write.table(m,"D:/Michal-Asia/F2_A3G-Q1-vs-others.txt",quote=F,sep="\t")

Idents(dlbc) <- dlbc@meta.data$orig.ident
f <- subset(dlbc,idents = "F3")
Idents(f) <- f@meta.data$A3G_quartiles
m <- FindMarkers(f,ident.1 = 1)
write.table(m,"D:/Michal-Asia/F3_A3G-Q1-vs-others.txt",quote=F,sep="\t")

Idents(dlbc) <- dlbc@meta.data$orig.ident
f <- subset(dlbc,idents = "F4")
Idents(f) <- f@meta.data$A3G_quartiles
m <- FindMarkers(f,ident.1 = 1)
write.table(m,"D:/Michal-Asia/F4_A3G-Q1-vs-others.txt",quote=F,sep="\t")
