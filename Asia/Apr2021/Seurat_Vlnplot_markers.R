library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Cairo)

dlbc.d1 <- Read10X(data.dir="D:/Michal-Asia/raw/DLBCL1/")
d1 <- CreateSeuratObject(counts = dlbc.d1, project = "D1-F", min.cells = 3, min.features = 200)
dlbc.d2 <- Read10X(data.dir="D:/Michal-Asia/raw/DLBCL2/")
d2 <- CreateSeuratObject(counts = dlbc.d2, project = "D2-M", min.cells = 3, min.features = 200)
dlbc.d3 <- Read10X(data.dir="D:/Michal-Asia/raw/DLBCL3/")
d3 <- CreateSeuratObject(counts = dlbc.d3, project = "D3-F", min.cells = 3, min.features = 200)
dlbc.l1 <- Read10X(data.dir="D:/Michal-Asia/raw/rLN1/")
l1 <- CreateSeuratObject(counts = dlbc.l1, project = "L1-M", min.cells = 3, min.features = 200)
dlbc.l2 <- Read10X(data.dir="D:/Michal-Asia/raw/rLN2/")
l2 <- CreateSeuratObject(counts = dlbc.l2, project = "L2-M", min.cells = 3, min.features = 200)
dlbc.l3 <- Read10X(data.dir="D:/Michal-Asia/raw/rLN3/")
l3 <- CreateSeuratObject(counts = dlbc.l3, project = "L3-F", min.cells = 3, min.features = 200)
dlbc.f1 <- Read10X(data.dir="D:/Michal-Asia/raw/FL1/")
f1 <- CreateSeuratObject(counts = dlbc.f1, project = "F1-M", min.cells = 3, min.features = 200)
dlbc.f2 <- Read10X(data.dir="D:/Michal-Asia/raw/FL2/")
f2 <- CreateSeuratObject(counts = dlbc.f2, project = "F2-M", min.cells = 3, min.features = 200)
dlbc.f3 <- Read10X(data.dir="D:/Michal-Asia/raw/FL3/")
f3 <- CreateSeuratObject(counts = dlbc.f3, project = "F3-M", min.cells = 3, min.features = 200)
dlbc.f4 <- Read10X(data.dir="D:/Michal-Asia/raw/FL4/")
f4 <- CreateSeuratObject(counts = dlbc.f4, project = "F4-M", min.cells = 3, min.features = 200)


x1 <- merge(d1, y = d2, add.cell.ids = c("D1", "D2"), project = "combined1")
x1 <- merge(x1, y = d3, add.cell.ids = c(""  , "D3"), project = "combined1")
x1 <- merge(x1, y = l1, add.cell.ids = c(""  , "L1"), project = "combined1")
x1 <- merge(x1, y = l2, add.cell.ids = c(""  , "L2"), project = "combined1")
x1 <- merge(x1, y = l3, add.cell.ids = c(""  , "L3"), project = "combined1")
x1 <- merge(x1, y = f1, add.cell.ids = c(""  , "F1"), project = "combined1")
x1 <- merge(x1, y = f2, add.cell.ids = c(""  , "F2"), project = "combined1")
x1 <- merge(x1, y = f3, add.cell.ids = c(""  , "F3"), project = "combined1")

dlbc <- merge(x1, y = f4, add.cell.ids = c(""  , "F4"), project = "DLBC")
dlbc[["percent.mt"]] <- PercentageFeatureSet(dlbc, pattern = "^MT-")
dlbc <- subset(dlbc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

dlbc <- NormalizeData(dlbc, normalization.method = "RC", scale.factor = 10000)

# check APOBEC family genes expression; RNA editing and hypoxia factor
plot <- VlnPlot(dlbc,
                features = c("APOBEC2" , "APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D",
                             "APOBEC3F", "APOBEC3G", "APOBEC3H", "AICDA"),
                #features = c("HIF1A","TM7SF3","EIF3I","RFX7"),  
                #features = c("EIF3I","RSL24D1","RPL10A","RPL18A","RPL23A","RPS16","RPS2","PABPC4"),
                same.y.lims=T)
# 
# 
CairoPNG(file="D:/Michal-Asia/APOBEC_family-exp.png",
         units="in",
         width=14,
         height=10,
         pointsize=8,
         dpi=300)
print(plot)
dev.off()
