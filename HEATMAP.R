#heatmap
library(Seurat)
stringsAsFactors = FALSE
set.seed(12345)
clustered <- readRDS("D:/gk/clustered.rds")
# �����Ƚϱ�ǩ
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
x=D1_W1_C0_df$avg_log2FC #degȡlogFC���в��������¸�ֵ��x
names(x)=rownames(D1_W1_C0_df) #degȡprobe_id���У���������Ϊ���ָ�x
cg=c(names(head(sort(x),20)),#��x���д�С�������У�ȡǰ100����100����ȡ���Ӧ��̽��������Ϊ������ֵ��cg
     names(tail(sort(x),20)))
DoHeatmap(subset(clustered, idents=c("0","1")),  features = cg)
p6 <- DoHeatmap(cca.merged.obj, features = as.vector(marker_known$Marker), cells = 1:10000)
ggsave("./Seurat_output/cca/marker_Heatmap.png", p6,  width = 42, height = 8, units = "in")

rownames(D1_W1_C0_df)