library(Seurat)
library(DT)
library(ggplot2)
library(patchwork)
library(reshape2)
stringsAsFactors = FALSE
set.seed(12345)
fdir <- "D:/M-GSGC0246918/result/ExpressionData/123/"
(dir.files <- list.files(pattern="*.tsv.gz", path=fdir)) #�鿴����
(dir.list <- paste0(fdir, dir.files))
mice.files <- lapply(dir.list, read.table) #������ȡ��������
length(mice.files)
mice.files[[1]][1:5, 1:5]
names(mice.files) <- sub("\\D-1.counts.txt", "", dir.files) 
names(mice.files)
d1 <- read.table("D:/M-GSGC0246918/result/ExpressionData/123/D-1.counts.tsv.gz", header=T,row.names = 1)
d1[1:5, 1:5]
d1.obj <- CreateSeuratObject(counts = d1, project = "d1")
d1.obj
d1.obj@assays$RNA@data[1:10, 1:10]
datatable(d1.obj@meta.data, options = list(pageLength = 10))
d3 <- read.table("D:/M-GSGC0246918/result/ExpressionData/123/D-3.counts.tsv.gz", header=T,row.names = 1)
d3[1:5, 1:5]
d3.obj <- CreateSeuratObject(counts = d3, project = "d3")
d3.obj
d3.obj@assays$RNA@data[1:10, 1:10]
datatable(d1.obj@meta.data, options = list(pageLength = 10))
#sample merge
d1d3.combind <- merge(d1.obj, y = d3.obj,add.cell.ids = c("d1", "d3"), project = "d1d3")
head(d1d3.combind@meta.data)
table(d1d3.combind$orig.ident)
#���ν�������
mice.obj.list <- SplitObject(d1d3.combind, split.by = "orig.ident")
mice.obj.list
mice.obj.list <- lapply(X = mice.obj.list, FUN = function(x) { 
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
}) #�ֱ��ÿ��������׼�����ҳ��߱����

features <- SelectIntegrationFeatures(object.list = mice.obj.list) #�ҳ�ÿ���������ڼ���ϸ�������Եĸ߱����
dat.anchors <- FindIntegrationAnchors(object.list = mice.obj.list, anchor.features = features, normalization.method = "LogNormalize", dims = 1:15, reduction="cca") #����������ġ�ê�㡱
cca.merged.obj <- IntegrateData(anchorset = dat.anchors) #ͨ����ê�㡱�ϲ�������ݼ�
save.image("G:/run.Rdata")
#����������
d1d3.combind$Group <- rep("D1", ncol(d1d3.combind)) 
D3 <- c("d3")
d1d3.combind$Group[d1d3.combind$orig.ident %in% d3] <- "D3"
addmargins(table(d1d3.combind$orig.ident, d1d3.combind$Group), 1)
#������������
d1d3.combind$sample <- as.vector(sapply(d1d3.combind$orig.ident, function(x){strsplit(x, "_")[[1]][2]}))
# add sample - ����������Ϊ����
cca.merged.obj$sample <- as.vector(sapply(cca.merged.obj$orig.ident, function(x){strsplit(x, "_")[[1]][2]}))
all.genes <- rownames(d1d3.combind)
#???
LabelPoints(plot = VariableFeaturePlot(d1d3.combind), 
            points = head(VariableFeatures(d1d3.combind), 10), 
            repel = TRUE)
#cca�ϲ����������
set.seed(12345)
DefaultAssay(cca.merged.obj) <- "integrated" #�ϲ������ݼ����ڷ���
cca.merged.obj <- ScaleData(cca.merged.obj, verbose = FALSE)
cca.merged.obj <- RunPCA(cca.merged.obj, verbose = FALSE)
DimPlot(cca.merged.obj, reduction = "pca")
ElbowPlot(cca.merged.obj, ndims=50)

cca.merged.obj <- FindNeighbors(cca.merged.obj, reduction = "pca", dims = 1:15)
cca.merged.obj <- FindClusters(cca.merged.obj, resolution = 0.8)
cca.merged.obj <- RunUMAP(cca.merged.obj, reduction = "pca", dims = 1:15)
cca.merged.obj <- RunTSNE(cca.merged.obj, reduction = "pca", dims = 1:15)
#��ͼ
p1 <- DimPlot(cca.merged.obj, reduction = "umap", label=T) + NoLegend()
p2 <- DimPlot(cca.merged.obj, reduction = "tsne", label=T) + NoLegend()
p1+p2
p3 <- DimPlot(cca.merged.obj, reduction = "tsne", label=T,  split.by = "orig.ident") + NoLegend()
p3
ggsave("G:Seurat_output/cca/tsne_splitbySample0402.png", p3,  width = 10, height = 5, units = "in")
ggsave("G:Seurat_output/cca/tsne0402.png", p2,  width = 5, height = 5, units = "in")
# Look at cluster IDs of the first 5 cells
head(cca.merged.obj@meta.data)
table(cca.merged.obj$seurat_clusters)
table(Idents(cca.merged.obj))
#Ѱ�Ҳ�����������
library(limma)
library(dplyr)
library(magrittr)
logFCfilter=0.25
adjPvalfilter=0.05
cluster.markers <- FindAllMarkers(obj=cca.merged.obj,
                                only.pos=FALSE,
                                min.pct=0.25,
                                logfc.thrshold=logFCfilter)
sig.markers <- cluster.markers[(abs(as.numeric(as.vector(cluster.markers$avg_log2FC))))>logFCfilter & as.numeric((as.vector(cluster.markers$p_val_adj))<adjPvalfilter),]
write.xls(sig.markers,file="cluster.markers.xls",sep="",row.names = F,quote = F)
top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
#����marker�ڸ�clusterd����ͼ
pdf(file="cluster.markers.pdf",width=48,height=36)
DoHeatmap(obj = cca.merged.obj, features = top10$gene) + NoLegend()
dev.off()
#����marker��CLUSTER��С����ͼ
write.csv(top10,file="G:/top10.csv")
pdf(file="cluster.markerViolin.pdf",width=10,height=6)
VlnPlot(obj = cca.merged.obj, features = c("S100a9","S100a8","Ifitm1","Il1f9","Cxcr2","Fpr1","Rdh12","Stfa2l1","Asprv1"),pt.size=0)
dev.off()
#����marker��cluster��ɢ��ͼ
pdf(file="cluser0markerScatter0402.pdf",width = 10,height = 6)
FeaturePlot(obj = cca.merged.obj, features = c("S100a9","S100a8","Ifitm1","Il1f9","Cxcr2","Fpr1","Rdh12","Stfa2l1","Asprv1"),pt.size=0.1)
dev.off()
#����marker�ڸ�cluster������ͼ
pdf(file="cluster0markerBubble.pdf",width=20,height=6)
cluster0Markert= c("S100a9","S100a8","Ifitm1","Il1f9","Cxcr2","Fpr1","Rdh12","Stfa2l1","Asprv1")
DotPlot(obj = cca.merged.obj,features = cluster0Markert )
dev.off()                                           
#SingleRע��ϸ������
library(celldex)
immgen <-ImmGenData()
library(SingleR)
library(BiocParallel)
#�����������ϸ�����ͣ���ѡ��method ����Ϊ��single��������ÿ��ϸ���ļ��������clusters�����ĸ���������
pred.cca.merged.obj <- SingleR(test = cca.merged.obj@assays$RNA@data,ref = immgen,labels = immgen$label.fine,clusters = cca.merged.obj@active.ident,fine.tune = TRUE,BPPARAM = SnowParam(40) )
pred.cca.merged.obj$pruned.labels
#�鿴ע��׼ȷ��
pdf(file="celllable0403",width=10,height=10)
plotScoreHeatmap(pred.cca.merged.obj,clusters = pred.cca.merged.obj@rownames,fontsize.row=9,show_colnames = F )

il1 <- subset(cluster.markers,p_val=0)
il1