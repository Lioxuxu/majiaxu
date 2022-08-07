library(Seurat)
library(DT)
library(ggplot2)
library(patchwork)
library(reshape2)
library(dplyr)
library(monocle)
library(reshape2)
clusterd <- readRDS("/Users/lioxuxu/Documents/M-GSGC0246918/clustered.rds")
importCDS(clusterd, import_all = TRUE)
DimPlot(clusterd)
table(clusterd$seurat_clusters)
table(Idents(clusterd))
table(clusterd$orig.ident)
Idents(object = clusterd)
levels(x = clusterd)
head(x = rownames(x = clusterd))
head(x = colnames(x = clusterd))
names(clusterd)
#提取特定细胞群构建新的seurat对象
clusterd29 <- subset(x=clusterd,idents=c("2","4","9","6","8"))
Idents(clusterd29) <- "orig.ident" 
clusterd29 <- subset(x=clusterd29,idents=c("D-1","D-3","D-5","D-7"))
table(clusterd29$orig.ident)
#删除原始对象
rm(clusterd)
gc()
#提取seurat数据转化为cellsetdata
RA_matrix<-as(as.matrix(clusterd29@assays$RNA@data), 'sparseMatrix')
feature_ann<-data.frame(gene_id=rownames(RA_matrix),gene_short_name=rownames(RA_matrix))
rownames(feature_ann)<-rownames(RA_matrix)
RA_fd<-new("AnnotatedDataFrame", data = feature_ann)
sample_ann<-clusterd29@meta.data
rownames(sample_ann)<-colnames(RA_matrix)
RA_pd<-new("AnnotatedDataFrame", data =sample_ann)
RA.cds<-newCellDataSet(RA_matrix,phenoData =RA_pd,featureData =RA_fd,expressionFamily=negbinomial.size())
head(pData(RA.cds))
head(fData(RA.cds))
#删除seurat对象
rm(clusterd29)
gc()
#Estimate size factors and dispersions
RA <- estimateSizeFactors(RA.cds)
RA <- estimateDispersions(RA)
#Ordering based on genes that differ between clusters
RA29_expressed_genes <-  row.names(fData(RA.cds))
clustering_DEG_genes <-
  differentialGeneTest(RA[RA29_expressed_genes,],
                       fullModelFormulaStr = '~seurat_clusters',
                       cores = 1)
RA29_ordering_genes <-
  row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
RA29 <-
  setOrderingFilter(RA,
                    ordering_genes = RA29_ordering_genes)

RA29 <-
  reduceDimension(RA29, method = 'DDRTree')

RA29 <-
  orderCells(RA29)
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$orig.ident)[,"W-1"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
RA29 <-
  orderCells(RA29, root_state = GM_state(RA29))

plot_cell_trajectory(RA29, color_by = "State")
plot_cell_trajectory(RA29, color_by="orig.ident")
plot_cell_trajectory(RA29, color_by="seurat_clusters")
plot_cell_trajectory(RA29, color_by="Pseudotime")
plot_cell_trajectory(RA29, color_by = "State") +
  facet_wrap(~State, nrow = 2)
#Finding Genes that Distinguish Cell Type or State
cds_to_be_test <- RA29[RA29_expressed_genes,]
diff_test_clusters <- differentialGeneTest(cds_to_be_test,
                                      fullModelFormulaStr = "~seurat_clusters")
diff_test_res[,c("gene_short_name", "pval", "qval")]
diff_test_clusters_p0 <-subset(diff_test_res,qval < 0.001)
#展示挑选的基因分布
markergene<- row.names(subset(fData(RA),
                                 gene_short_name %in% c("Arg1", "Mrc1", "Cxcl2")))
marker_subset <- RA[markergene,]
plot_genes_jitter(marker_subset,
                  grouping = "seurat_clusters",
                  color_by = "seurat_clusters",
                  nrow= 1,
                  ncol = NULL,
                  plot_trend = TRUE)

#Finding Genes that Change as a Function of Pseudotime
diff_test_res <- differentialGeneTest(cds_to_be_test,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
diff_test_res[,c("gene_short_name", "pval", "qval")]
#展示挑选基因的拟时序变化
selected_subset <- RA29[markergene,]
plot_genes_in_pseudotime(selected_subset, color_by = "orig.ident")

## 策略2：Selecting genes with high dispersion across cells
disp_table <- dispersionTable(RA)
ordering_genes <- subset(disp_table, 
                         mean_expression >= 0.5 & 
                           dispersion_empirical >= 1 * dispersion_fit)$gene_id

RA <- setOrderingFilter(RA, ordering_genes)
plot_ordering_genes(RA)
## 挑选变异尺度大的基因，如图所示
RA <- reduceDimension(RA, max_components=2)
RA <- orderCells(RA)
#给予拟时序起点#
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$orig.ident)[,"D-1"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
RA_pseudotime<- orderCells(RA, root_state = GM_state(RA))
## 排序好的细胞按照发育顺序可视化
plot_cell_trajectory(RA_pseudotime, color_by="orig.ident")
plot_cell_trajectory(RA_pseudotime, color_by="seurat_clusters")
plot_cell_trajectory(RA_pseudotime, color_by="State")
plot_cell_trajectory(RA_pseudotime, color_by="Pseudotime")
plot_cell_trajectory(RA_pseudotime, color_by = "State") +
  facet_wrap(~State, nrow = 3)
#Clustering Genes by Pseudotemporal Expression Pattern
#select gene names
marker_genes <- row.names(subset(fData(HSMM_myo),
                                 gene_short_name %in% c("MEF2C", "MEF2D", "MYF5",
                                                        "ANPEP", "PDGFRA","MYOG",
                                                        "TPM1",  "TPM2",  "MYH2",
                                                        "MYH3",  "NCAM1", "TNNT1",
                                                        "TNNT2", "TNNC1", "CDK1",
                                                        "CDK2",  "CCNB1", "CCNB2",
                                                        "CCND1", "CCNA1", "ID1")))
#相似趋势基因分组可视化，num_clusters=the cluster of gene,not cell
diff_test_res <- differentialGeneTest(RA[ordering_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.001))
png <- plot_pseudotime_heatmap(RA_c2[sig_gene_names,],
                        num_clusters = 5,
                        cores = 1,return_heatmap=T,
                        show_rownames = T)
ggsave("C/Users/DELL/Documents/monocle",png, width = 210, height = 297, units = "mm") 
#beam分析 分析轨迹分支点
beam_res <- BEAM(RA,branch_point = 3,cores = 1)
beam_res <- beam_res[order(beam_res$qval),]
beam_res <- beam_res[,c("gene_short_name","pval","qval")]
#特殊类型热图，可视化明显依赖于分支基因的变化
#plot_genes_branched_heatmap(RA[row.names(subset(beam_res,qval<1e-5))],
png <- plot_genes_branched_heatmap(RA[row.names(beam_res)[1:50]],
                            branch_point= 3,
                            num_cluster =3,
                            core=1,
                            use_gene_short_name = T,
                            show_rownames = T)
ggsave("/Users/lioxuxu/Desktop/24689D/plan2/branch1.png",png,width = 210, height = 297, units = "mm",dpi = 300)
# 提取热图基因名称画图的方法
p=plot_genes_branched_heatmap(RA[row.names(subset(beam_res,qval<1e-5))],
                            branch_point= 3,
                            num_cluster =3,
                            core=1,
                            return_heatmap=T,
                            show_rownames = T)
#提取热图基因名
p[["ph_res"]][["tree_row"]]
clusters <- cutree(p[["ph_res"]][["tree_row"]], k = 4)
clustering <- data.frame(clusters)
clustering[,2] <- as.character(clustering[,1])
colnames(clustering) <- c("gene","Gene_Clusters")
clustering[,1] <- row.names(clustering)
table(clustering)
write_csv(clustering,"/Users/lioxuxu/Desktop/24689D/249/clustered/D/branch3gene.csv",col_names = T)

