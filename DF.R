library(Seurat)
library(DT)
library(ggplot2)
library(patchwork)
library(reshape2)
library(dplyr)
stringsAsFactors = FALSE
set.seed(12345)
clustered <- readRDS("D:/gk/clustered.rds")
DimPlot(clustered)
table(clustered$seurat_clusters)
table(Idents(clustered))
Idents(object = clustered)
levels(x = clustered)
head(x = rownames(x = clustered))
head(x = colnames(x = clustered))
names(clustered)
# ?????Ƚϱ?ǩ
Idents(clustered) <- clustered$seurat_clusters
clustered$seurat_clusters_day <- paste(Idents(clustered), clustered$orig.ident, sep = "_")
table(clustered$seurat_clusters_day)
Idents(clustered) <- clustered$seurat_clusters_day
table(Idents(clustered)) 
name <- as.vector(unique(clustered$seurat_clusters))
D_7 <- paste0(name, "_D-7")
W_7 <- paste0(name, "_W-7")
#d vs.W
for (i in 1:length(name)){
  cat(W_7[i],"\n")
  if(W_7[i]%in%clustered$seurat_clusters_day==F|W_1[i]%in%clustered$seurat_clusters_day==F) next
  outname <- paste0("D:/gk/output/W7-W1/", name[i], ".csv") #????????·??
  temp <- FindMarkers(clustered, ident.1 = W_7[i], ident.2 = W_1[i], logfc.threshold = 0.25, min.pct = 0.1)
  write.csv(temp, outname, quote=F)
}

#options(scipen = 200) ??ѧ??????ת??
#????????????ͼ??
library(ggpubr)
options(stringsAsFactors = F)
D1_W1_C17 <- read.table("D:/gk/output/D1-W1/17.csv",header=T,sep = ",",row.names = 1)
rownames( D1_W1_C17)
#????????
  nrDEG= D1_W1_C17
  head(nrDEG)
  attach(nrDEG)
  #plot("avg_log2FC",-"log10(p_val)")
  #aBiocManager::install("ggpubr")
  
  df=nrDEG
  df$v= -log10(p_val) #df??????һ??'v',ֵΪ-log10(p_val)
  ggscatter(df, x = "avg_log2FC", y = "v",size=1.0)
  
  
  df$g=ifelse(df$p_val>0.01,'stable', #if ?жϣ???????һ??????P.Value>0.01????Ϊstable????
              ifelse( df$avg_log2FC >1,'up', #???Ͼ?else ???򣺽???��??ʼ?ж???ЩP.Value<0.01?Ļ???????if ?жϣ?????logFC >1.5,??Ϊup???ϵ???????
                      ifelse( df$avg_log2FC < -1,'down','stable') )#???Ͼ?else ???򣺽???��??ʼ?ж???ЩlogFC <1.5 ?Ļ???????if ?жϣ?????logFC <1.5????Ϊdown???µ??????򣬷???Ϊstable????
  )
  table(df$g)
  df$name=rownames(df)
  head(df)
  write.csv(df,"D:/gk/output/dfgene/D-WC17/D1_W1_C17_DF",row.names = T)
  #D5_W5_C4_DF <- read.csv("D:/gk/output/dfgene/D5_W5_C4_DF",header=T,sep = ",")
  ggscatter(df, x = "avg_log2FC", y = "v",size=0.5,color = "g")
  ggscatter(df, 
            x = "avg_log2FC", y = "v", color = "g",size = 0.5,label = "name", repel = T,
            label.select = rownames(df)[df$g != 'stable'] ,
            #label.select = c('Slpi', 'Retnla'), #??ѡһЩ??????ͼ????ʾ??��
            palette = c("#00AFBB","#2E9FDF" , "#FC4E07","#E7B800") )
  ggsave('D:/gk/output/dfgene/D-WC17/D1_W1_C17_1volcano.png')
 #GO????
  #????????????
  library(ggplot2)
  
  library(clusterProfiler)
  library(org.Mm.eg.db)
  #BiocManager#:install("org.Mm.eg.db")
  D1_W1_C0 <- #BiocManager::install("clusterProfiler")
  #BiocManager::install("monocle")
  
  read.table("D:/gk/output/D1-W1/0.csv",header=T,sep = ",",row.names = 1) 
  rownames( D1_W1_C0)
  nrDEG= D1_W1_C0
  head(nrDEG)
  attach(nrDEG)
  
  
  #plot("avg_log2FC",-"log10(p_val)")
  #BiocManager#:install("ggpubr")
  
  df=nrDEG
  df$v= -log10(p_val) #df??????һ??'v',ֵΪ-log10(p_val)
  ggscatter(df, x = "avg_log2FC", y = "v",size=1.0)
  
  
  df$g=ifelse(df$p_val>0.01,'stable', #if ?жϣ???????һ??????P.Value>0.01????Ϊstable????
              ifelse( df$avg_log2FC >1,'up', #???Ͼ?else ???򣺽???��??ʼ?ж???ЩP.Value<0.01?Ļ???????if ?жϣ?????log2FC >1,??Ϊup???ϵ???????
                      ifelse( df$avg_log2FC < -1,'down','stable') )#???Ͼ?else ???򣺽???��??ʼ?ж???Щlog2FC <1 ?Ļ???????if ?жϣ?????log2FC <1????Ϊdown???µ??????򣬷???Ϊstable????
  )
  df$symbol=rownames(df)
  deg <- bitr(unique(df$symbol), fromType = "SYMBOL",
             toType = c( "ENTREZID"),
             OrgDb = org.Mm.eg.db)
  
  #bitr????ΪIDת????
  #bitr(geneID, fromType, toType, OrgDb, drop = TRUE)??
  #geneid ??????ID???? ?? fromtype ?? ????ID?ͣ?toType??????ID?ͣ?orgdb ??ע?????ݿ⣩
  head(df)
  DEG=df#??deg???ݸ?ֵ??DEG????
  head(DEG)
  DEG=merge(DEG,deg,by.y='SYMBOL',by.x='symbol')
  #??????DEG,dfͨ????DEG??'symbol'?У?df??'SYMBOL'??��????һ????ת??ID
  head(DEG)
  save(DEG,file = 'anno_DEG.Rdata')
  gene_up= DEG[DEG$g == 'up','ENTREZID'] #ѡ???ϵ?????ID
  gene_down=DEG[DEG$g == 'down','ENTREZID'] #ѡ???µ?????ID
  gene_diff=c(gene_up,gene_down)#?ó????µ?????ID
  gene_all=as.character(DEG[ ,'ENTREZID'] )#?ó????л???ID
  data(geneList, package="DOSE")#?ó?geneList????
  head(geneList)#?鿴????{
    
   g_list=list(gene_up=gene_up,
               gene_diff=gene_diff,
               gene_down=gene_down)
    
      go_enrich_results <- lapply( g_list , function(gene) {
        lapply( c('BP','MF','CC') , function(ont) {
          cat(paste('Now process ',ont ))
          ego <- enrichGO(gene          = gene,
                          universe      = gene_all,
                          OrgDb         = org.Mm.eg.db,
                          ont           = ont ,
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.99,
                          qvalueCutoff  = 0.99,
                          readable      = TRUE)
          
          print( head(ego) )
          return(ego)
        })
      })
      
      save(go_enrich_results,file = 'H:/gk/rds/D1_W1_C17_go_enrich_results.Rdata')
    load(file = 'H:/gk/rds/D1_W1_C17_go_enrich_results.Rdata')
    n1= c('gene_up','gene_down','gene_diff')
    n2= c('BP','MF','CC') 
    for (i in 1:3){
      for (j in 1:3){
        fn=paste0('dotplot_',n1[i],'_',n2[j],'.png')
        cat(paste0(fn,'\n'))
        png(fn,res=150,width = 1480)
        print( dotplot(go_enrich_results[[i]][[j]] ))
        dev.off()
      }
    }
 #KEGG????????
    ###  ?ֱ???ͼ
    kk.up <- enrichKEGG(gene         = gene_up,
                        organism     = 'mmu',
                        universe     = gene_all,
                        pvalueCutoff = 0.9,
                        qvalueCutoff =0.9)
    head(kk.up)[,1:6]
    #write.csv(as.data.frame(kk.up),"D:/gk/output/dfgene/D-WC7/D1_W1_C7_kk.up.csv",row.names =FALSE)
    barplot(kk.up );ggsave('D:/gk/output/dfgene/D-WC17/D1-W1/D1_W1_C17_kk.up.dotplot.png')
    kk.down <- enrichKEGG(gene         =  gene_down,
                          organism     = 'mmu',
                          universe     = gene_all,
                          pvalueCutoff = 0.9,
                          qvalueCutoff =0.9)
    head(kk.down)[,1:6]
    #write.csv(as.data.frame(kk.down),"D:/gk/output/dfgene/D-WC7/D1_W1_C7_kk.up.csv",row.names =FALSE)
    barplot(kk.down );ggsave('D:/gk/output/dfgene/D-WC17/D1-W1/D1_W1_C17_kk.down.dotplot.png')
    kk.diff <- enrichKEGG(gene         = gene_diff,
                          organism     = 'mmu',
                          universe     = gene_all,
                          pvalueCutoff = 0.9,
                          qvalueCutoff =0.9)
    head(kk.diff)[,1:6]
    #write.csv(as.data.frame(kk.diff),"D:/gk/output/dfgene/D-WC7/D3_W3_C7_kk.up.csv",row.names =FALSE)
    barplot(kk.diff );ggsave('D:/gk/output/dfgene/D-WC17/D1-W1/D1_W1_C17_kk.diff.dotplot.png')
    
    kegg_diff_dt <- as.data.frame(kk.diff)
    kegg_down_dt <- as.data.frame(kk.down)
    kegg_up_dt <- as.data.frame(kk.up)
    down_kegg<-kegg_down_dt[kegg_down_dt$pvalue<0.05,];down_kegg$group=-1
    up_kegg<-kegg_up_dt[kegg_up_dt$pvalue<0.05,];up_kegg$group=1
    source('functions.R')
    g_kegg=kegg_plot(up_kegg,down_kegg)
    print(g_kegg)
    ggsave(g_kegg,filename = 'D:/gk/output/dfgene/D-WC7/D5_W5_C7_
           kegg_up_down.png')
    
    #ͨ·ͼ
    browseKEGG(kk.down,"mmu04929")
 #????Ĭ?ϴ洢·??
  setwd("D:/gk/output/dfgene/")
  ggsave('D1_W1_C4MA.png')
  

#for heatmap





load(file = 'D:/gk/rds/D1-W1-C0_go_enrich_results.Rdata')
table(go_enrich_results$g)
write.table(as.data.frame(go_enrich_results), 'D3-W3-C11_go', row.names = FALSE, quote = FALSE)
D7_W7_C0 <- read.table("D:/gk/output/dfgene/D-WC0/D7_W7_C0_DF",header=T,sep = ",",row.names = 1)
table(list$g)
D1_W1_C17 <- read.table("D:/gk/output/D1-W1/17.csv",header=T,sep = ",",row.names = 1)

table(Idents(clustered))

