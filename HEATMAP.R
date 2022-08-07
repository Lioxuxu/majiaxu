#heatmap
library(Seurat)
stringsAsFactors = FALSE
set.seed(12345)
clustered <- readRDS("D:/gk/clustered.rds")
# 构建比较标签
Idents(clustered) <- clustered$seurat_clusters
clustered$seurat_clusters_day <- paste(Idents(clustered), clustered$orig.ident, sep = "_")
table(clustered$seurat_clusters_day)
Idents(clustered) <- clustered$seurat_clusters_day
table(Idents(clustered)) 
name <- as.vector(unique(clustered$seurat_clusters))
D_7 <- paste0(name, "_D-7")
W_7 <- paste0(name, "_W-7")
D_5 <- paste0(name, "_D-5")
W_5 <- paste0(name, "_W-5")
D_3 <- paste0(name, "_D-3")
W_3 <- paste0(name, "_D-3")
D_1 <- paste0(name, "_D-1")
W_1 <- paste0(name, "_D-1")
table(D_7)
table(name)

D1_W1_C0_df <- FindMarkers(clustered[0], ident.1 = W_7, ident.2 = W_1, logfc.threshold = 0.25, min.pct = 0.1)
table(Idents(clustered))
Sclustered <- ScaleData(clustered, verbose = FALSE)
x=D1_W1_C0_df$avg_log2FC #deg取logFC这列并将其重新赋值给x
names(x)=rownames(D1_W1_C0_df) #deg取probe_id这列，并将其作为名字给x
cg=c(names(head(sort(x),20)),#对x进行从小到大排列，取前100及后100，并取其对应的探针名，作为向量赋值给cg
     names(tail(sort(x),20)))
DoHeatmap(subset(clustered, idents=c("0","1")),  features = cg)
p6 <- DoHeatmap(cca.merged.obj, features = as.vector(marker_known$Marker), cells = 1:10000)
ggsave("./Seurat_output/cca/marker_Heatmap.png", p6,  width = 42, height = 8, units = "in")

rownames(D1_W1_C0_df)
