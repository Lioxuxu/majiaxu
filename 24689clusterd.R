library(Seurat)
library(DT)
library(ggplot2)
library(patchwork)
library(reshape2)
library(dplyr)
stringsAsFactors = FALSE
set.seed(12345)
clustered <- readRDS("/Users/lioxuxu/Desktop/M-GSGC0246918/clustered.rds")
DimPlot(clustered)
table(clustered$seurat_clusters)
table(Idents(clustered))
Idents(object = clustered)
levels(x = clustered)
head(x = rownames(x = clustered))
head(x = colnames(x = clustered))
names(clustered)
#提取特定细胞群构建新的seurat对象
mclustered <- subset(x=clustered,idents=c("2","4","6","8","9"))
#PCA
mclustered=RunPCA(object= mclustered,npcs = 50,pc.genes=VariableFeatures(object = mclustered)) 
ElbowPlot(mclustered)
pdf(file="pca.elbowplot.pdf",width=10,height=8)
ElbowPlot(mclustered)
dev.off()

mclustered <- FindNeighbors(object = mclustered, dims = 1:20)
#????ͬresolution?Է?Ⱥ??Ӱ?? 
mclustered <- FindClusters(object = mclustered, resolution = c(seq(0,1.6,.2))) #????resolutionΪ??0??1.6????Ϊ0.2?ĸ???ֵʱ?ķ?Ⱥ????
#???ϲ?????????????mclustered@metadata,ÿ???ֱ??ʶ??е???һ?У?ʹ??active.identҲ???Բ鿴??active.ident??level???Ǹ?????Ⱥ
library(clustree)
pdf(file="06.clustree.pdf",width=10,height=14)  #???ӻ???ͬresolution?Է?Ⱥ??Ӱ??
clustree(mclustered@meta.data, prefix = "RNA_snn_res.")
dev.off()
#?Է????Խ?ά?Ľ??????ӻ?????ͨ??idents()????��ָ???ֱ???
Idents(object = mclustered) <- "RNA_snn_res.1.0"
levels(mclustered)

#Ҳ??ֱ??ָ??resolution
mclustered <- FindClusters(object = mclustered, resolution = 0.6)
#TSNE????
mclustered <- RunTSNE(object = mclustered, dims = 1:20)                      
pdf(file="07.0.6TSNE.pdf",width=6.5,height=6)
TSNEPlot(object = mclustered, pt.size = 0.5, label = TRUE)   
dev.off()
write.table(mclustered$seurat_clusters,file="07.0.6tsneCluster.txt",quote=F,sep="\t",col.names=F)

#umap???? 
mclustered <- RunUMAP(mclustered, dims = 1:20)
pdf(file="08.0.6uMAP.pdf",width=6.5,height=6)
DimPlot(mclustered, reduction = "umap",pt.size=0.5,label = TRUE)
dev.off()
#????????????
logFCfilter=0.5
adjPvalFilter=0.05
mclustered.markers <- FindAllMarkers(object = mclustered,
                                only.pos = FALSE,
                                min.pct = 0.25,
                                logfc.threshold = logFCfilter)
sig.markers=mclustered.markers[(abs(as.numeric(as.vector(mclustered.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(mclustered.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.markers,file="09.0.6markers.xls",sep="\t",row.names=F,quote=F)
top10 <- mclustered.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10,file="09.0.6top10markers.xls",sep="\t",row.names=F,quote=F)
#????marker?ڸ???cluster????ͼ
pdf(file="09.06tsneHeatmap.pdf",width=48,height=36)
DoHeatmap(object = mclustered, features = top10$gene) + NoLegend()
dev.off()  #?????ٶȽ???

#????marker??С????ͼ
pdf(file="02.markerViolin.pdf",width=10,height=6)
VlnPlot(object = mclustered, features = c("Cd83", "Cxcl2"))
dev.off()

#????marker?ڸ???cluster??ɢ??ͼ
pdf(file="02.markerScatter.pdf",width=10,height=6)
FeaturePlot(object = mclustered, features = c("Cd83", "Cxcl2"))
dev.off()

#????marker?ڸ???cluster??????ͼ
pdf(file="06.markerBubble.pdf",width=12,height=6)
cluster0Marker=c("Cbr2", "Lyve1", "Selenop", "Folr2", "Ednrb", "F13a1", "Mrc1", "Igf1", "Slc40a1
", "Cd163")
DotPlot(object = mclustered, features = cluster0Marker)
dev.off()
#ϸ???ֲ?????
head(mclustered@meta.data)
table(mclustered$seurat_clusters)
table(Idents(mclustered))
table(mclustered$seurat_clusters, mclustered$orig.ident) 
 

#ϸ??ע??
library(celldex)
ref <- ImmGenData()
library(SingleR)
library(BiocParallel)
pred.mclustered <- SingleR(test = mclustered@assays$RNA@data, ref = ref,labels = ref$label.main, clusters = mclustered@active.ident, fine.tune = TRUE, BPPARAM = MulticoreParam(40))
#????????????ϸ?????ͣ???ѡ??method ????Ϊ??single????????ÿ??ϸ???ļ?????????clusters????ָ???ĸ???????????
pred.mclustered$pruned.labels

#?鿴ע??׼ȷ?? 
pdf(file="10.0.6celllabel.pdf",width=10,height=10)   
plotScoreHeatmap(pred.mclustered, clusters=pred.mclustered@rownames, fontsize.row = 9,show_colnames = T)
dev.off()

#???ƴ?cell label??tsne??umapͼ
new.cluster.ids <- pred.mclustered$pruned.labels
names(new.cluster.ids) <- levels(mclustered)
levels(mclustered)
# 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
mclustered <- RenameIdents(mclustered,new.cluster.ids)
levels(mclustered)

pdf(file="10.0.6TSNE_lable.pdf",width=6.5,height=6)
TSNEPlot(object = mclustered, pt.size = 0.5, label = TRUE)  
dev.off()

pdf(file="10.uMAP_label.pdf",width=6.5,height=6)
DimPlot(mclustered, reduction = "umap",pt.size=0.5,label = TRUE)
dev.off()
  