library(Seurat)
library(dior)
library(reticulate)
library(SingleCellExperiment)
library(ggplot2)
###
scRNA <- read_h5(file = "scRNA.h5",
                 assay.name = 'RNA', 
                 target.object = 'seurat')
###创建分组
sample_names <- scRNA$orig.ident
# 创建分组信息向量
group <- ifelse(sample_names %in% c("M_P01", "M_L01", "M_P03"), "R",
                ifelse(sample_names %in% c("M_P02", "P_G04", "P_G05", "P_C07", "M_L07"), "S", NA))
# 将分组信息添加到 Seurat 对象的元数据中
scRNA$group <- group
##提取细胞进行二次聚类
Idents(scRNA) <- "cell_type"
scRNA <- subset(scRNA,cell_type==c("tumor"))
DimPlot(scRNA, group= "cell_type",reduction = "umap",label = T)
######标准化基因
scRNA <- NormalizeData(object = scRNA, normalization.method = "LogNormalize", 
                       scale.factor = 10000)
###找出2000个高变基因
scRNA <- FindVariableFeatures(object = scRNA, 
                              selection.method = "vst", nfeatures = 2000)
####对所有基因进行标准化
all.genes <- rownames(scRNA)
scRNA <- ScaleData(scRNA, features = VariableFeatures(object = scRNA))   
####进行pca分析，选择高变基因2000个进行整合
scRNA=RunPCA(object= scRNA,npcs = 50,pc.genes=VariableFeatures(object = scRNA))

# Determine percent of variation associated with each PC
pct <- scRNA [["pca"]]@stdev / sum(scRNA [["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2
# Minimum of the two calculation
pcs <- min(co1, co2)
pcs
# Create a dataframe with values
plot_df <- data.frame(pct = pct,   cumu = cumu,   rank = 1:length(pct))
# Elbow plot to visualize 
library(ggplot2)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()

####确定几个pc.num
ElbowPlot(scRNA, ndims = 50)
pc.num=1:13
library(harmony)
scRNA <- RunHarmony(scRNA, group.by.vars="orig.ident", 
                    max.iter.harmony = 20)
####寻找最佳区分度
library(clustree)
library(patchwork)
###重新载入
scRNA <- FindNeighbors(scRNA,reduction="harmony",dims = pc.num) 
scRNA <- FindClusters(scRNA,reduction="harmony",resolution = c(0.1,0.2,0.3,0.4,0.6,0.8,1.0,1.2))
scRNA@meta.data
p_clu = clustree(scRNA@meta.data, prefix = "RNA_snn_res.")
p_clu
scRNA <- FindClusters(scRNA,reduction="harmony",resolution = 0.3)
#scRNA = RunTSNE(scRNA, reduction="harmony",dims = pc.num)
scRNA <- RunUMAP(scRNA,reduction="harmony", dims = pc.num)
p<- DimPlot(scRNA, reduction = "umap",label = T)
p
table(scRNA$group)
P1=FeaturePlot(scRNA, features = "CTGF",split.by = "group",pt.size = 1.5, cols = c("lightgrey", "blue"))


P1=DimPlot(scRNA, group.by = "seurat_clusters",split.by = "group",reduction = "umap",label = T)
ggsave("细注释图.pdf", P1, width=10 ,height=4)
dior::write_h5(scRNA,file="scRNA_sub.h5",object.type = "seurat")
###开启多核
library(future)
library(tidyverse)
options(future.globals.maxSize = 20 * 1024^3)
plan(multisession, workers = 4) #开启多核运算 (8个核)
logFCfilter=0.5
adjPvalFilter=0.05
scRNA.markers <- FindAllMarkers(object = scRNA,
                                only.pos = FALSE,
                                min.pct = 0.25,
                                logfc.threshold = logFCfilter)
sig.markers=scRNA.markers[(abs(as.numeric(as.vector(scRNA.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(scRNA.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.markers,file="markers1.xls",sep="\t",row.names=F,quote=F)
top10 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
plan('sequential') #终止多核运算
saveRDS(scRNA,file = "scRNA_sub.rds")







