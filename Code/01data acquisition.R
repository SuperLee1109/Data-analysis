setwd("~/Desktop/designer/xiaohan")
library(Seurat)
library(Seurat)
library(harmony)
library(tidyverse)
library(patchwork)
rm(list = ls())
## 批量读取数据
### 设置数据路径与样本名称
assays <- dir("GSE254762_RAW/")
dir <- paste0("GSE254762_RAW/", assays)
dir
# 按文件顺序给样本命名，名称不要以数字开头，中间不能有空格 
samples_name = c("M_P01","M_L01","M_P02",
                 "M_P03","P_G04","P_G05",
                 "P_G06","P_C07","M_L07")

scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  #不设置min.cells过滤基因会导致CellCycleScoring报错：
  #Insufficient data values to produce 24 bins.  
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=samples_name[i],
                                       min.cells=3, min.features = 200)
  #给细胞barcode加个前缀，防止合并后barcode重名
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = samples_name[i])   
  #计算线粒体基因比例
  if(T){    
    scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-") 
  }
  #计算核糖体基因比例
  if(T){
    scRNAlist[[i]][["percent.rb"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^RP[SL]")
  }
  #计算红细胞基因比例
  if(T){
    HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
    HB.genes <- CaseMatch(HB.genes, rownames(scRNAlist[[i]]))
    scRNAlist[[i]][["percent.HB"]]<-PercentageFeatureSet(scRNAlist[[i]], features=HB.genes) 
  }
}

### 给列表命名并保存数据
dir.create("Integrate")
setwd("~/Desktop/designer/xiaohan/Integrate")

names(scRNAlist) <- samples_name
########
scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
scRNA
# An object of class Seurat 
table(scRNA$orig.ident)
scRNA@meta.data
saveRDS(scRNA, file = "scRNA_orig.rds")
dior::write_h5(scRNA, file='./my.h5',object.type = 'seurat') #singlecellexperiment
#scRNAlist <- SplitObject(scRNA, split.by = "orig.ident") #分割Seurat对象
# pbmc33k  pbmc3k  pbmc4k  pbmc6k  pbmc8k 
# 8381   33141    2700    4340    5418 
