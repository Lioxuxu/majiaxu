library(Seurat)
library(DT)
library(ggplot2)
library(patchwork)
library(reshape2)
library(dplyr)
library(reshape2)
library(ggpubr)
stringsAsFactors = FALSE
set.seed(12345)
clusterd <- readRDS("/Users/lioxuxu/Desktop/M-GSGC0246918/clustered.rds")
table(Idents(clusterd))
Idents(clusterd)
levels(clusterd)
head(clusterd@meta.data)
new.cluster.ids <- c("0","1","2","3","2","5","6","7","8","2","10","11","12","13","14","15","16","17")
names(new.cluster.ids) <- levels(clusterd)
clusterd <- RenameIdents(clusterd,new.cluster.ids)

c11 <- FindMarkers(clusterd,ident.1=11,idents.2=2,logfc.threshold = 0.25, min.pct = 0.1)
write.csv(c11,"/Users/lioxuxu/Desktop/c11vs249/c11.csv")
options(stringsAsFactors = F)
c11 <- read.table("/Users/lioxuxu/Desktop/c11vs249/c11.csv",header=T,sep = ",",row.names = 1)
nrDEG= c11
head(nrDEG)
attach(nrDEG)
df=nrDEG
df$log10p_val= -log10(p_val) 
ggscatter(df, x = "avg_log2FC", y = "log10p_val",size=1.0)
df$g=ifelse(df$p_val>0.01,'stable', 
            ifelse( df$avg_log2FC >1,'up', 
                    ifelse( df$avg_log2FC < -1,'down','stable') )
)
table(df$g)
df$name=rownames(df)
head(df)
write.csv(df,"/Users/lioxuxu/Desktop/c11vs249/c11_Df.csv",row.names = T)
#D5_W5_C4_DF <- read.csv("D:/gk/output/dfgene/D5_W5_C4_DF",header=T,sep = ",")
ggscatter(df, x = "avg_log2FC", y = "log10p_val",size=0.5,color = "g")
ggscatter(df, 
          x = "avg_log2FC", y = "log10p_val", color = "g",size = 0.5,label = "name", repel = T,
          label.select = rownames(df)[df$g != 'stable'] ,
          #label.select = c('Slpi', 'Retnla'),
          palette = c("#00AFBB","#2E9FDF" , "#FC4E07","#E7B800") )
ggsave('/Users/lioxuxu/Desktop/c11vs249/c11vs249.png')
#GO????
#????????????
library(ggplot2)

library(clusterProfiler)
library(org.Mm.eg.db)
#BiocManager#:install("org.Mm.eg.db")
df$symbol=rownames(df)
deg <- bitr(unique(df$symbol), fromType = "SYMBOL",
            toType = c( "ENTREZID"),
            OrgDb = org.Mm.eg.db)


head(df)
DEG=df
head(DEG)
DEG=merge(DEG,deg,by.y='SYMBOL',by.x='symbol')
head(DEG)
save(DEG,file = 'anno_DEG.Rdata')
gene_up= DEG[DEG$g == 'up','ENTREZID'] 
gene_down=DEG[DEG$g == 'down','ENTREZID']
gene_diff=c(gene_up,gene_down)
gene_all=as.character(DEG[ ,'ENTREZID'] )
data(geneList, package="DOSE")
head(geneList)

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

save(go_enrich_results,file = '/Users/lioxuxu/Desktop/c11vs249/c11vs249_go_enrich_results.Rdata')
load(file = '/Users/lioxuxu/Desktop/c11vs249/c11vs249_go_enrich_results.Rdata')
n1= c('gene_up','gene_down','gene_diff')
n2= c('BP','MF','CC') 
for (i in 1:3){
  for (j in 1:3){
    fn=paste0('dotplot_',n1[i],'_',n2[j],'.png')
    cat(paste0(fn,'\n'))
    png(fn,res=150,width = 1080)
    print( dotplot(go_enrich_results[[i]][[j]] ))
    dev.off()
  }
}

write.csv(as.data.frame(go_enrich_results[["gene_diff"]][[1]]),"/Users/lioxuxu/Desktop/c11vs249/2/c11vs249_diff_BP.csv",row.names =F) 
kk.up <- enrichKEGG(gene         = gene_up,
                    organism     = 'mmu',
                    universe     = gene_all,
                    pvalueCutoff = 0.9,
                    qvalueCutoff =0.9)
head(kk.up)[,1:6]
#write.csv(as.data.frame(kk.up),"D:/gk/output/dfgene/D-WC7/D1_W1_C7_kk.up.csv",row.names =FALSE)
barplot(kk.up );ggsave('/Users/lioxuxu/Desktop/c11vs249/c11vs249_kk.up.dotplot.png')
kk.down <- enrichKEGG(gene         =  gene_down,
                      organism     = 'mmu',
                      universe     = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff =0.9)
head(kk.down)[,1:6]
#write.csv(as.data.frame(kk.down),"D:/gk/output/dfgene/D-WC7/D1_W1_C7_kk.up.csv",row.names =FALSE)
barplot(kk.down );ggsave('/Users/lioxuxu/Desktop/c11vs249/c11vs249_kk.down.dotplot.png')
kk.diff <- enrichKEGG(gene         = gene_diff,
                      organism     = 'mmu',
                      universe     = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff =0.9)
head(kk.diff)[,1:6]
#write.csv(as.data.frame(kk.diff),"D:/gk/output/dfgene/D-WC7/D3_W3_C7_kk.up.csv",row.names =FALSE)
barplot(kk.diff );ggsave('/Users/lioxuxu/Desktop/c11vs249/c11vs249_kk.diff.dotplot.png')

kegg_diff_dt <- as.data.frame(kk.diff)
kegg_down_dt <- as.data.frame(kk.down)
kegg_up_dt <- as.data.frame(kk.up)
down_kegg<-kegg_down_dt[kegg_down_dt$pvalue<0.05,];down_kegg$group=-1
up_kegg<-kegg_up_dt[kegg_up_dt$pvalue<0.05,];up_kegg$group=1
source('functions.R')
g_kegg=kegg_plot(up_kegg,down_kegg)
print(g_kegg)
ggsave(g_kegg,filename = '/Users/lioxuxu/Desktop/c11vs249/c11vs249_
           kegg_up_down.png')

#ͨ·ͼ
browseKEGG(kk.down,"mmu04929")
#????Ĭ?ϴ洢·??
setwd("D:/gk/output/dfgene/")
ggsave('D1_W1_C4MA.png')
