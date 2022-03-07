# Main R code for artical titled 'High glucose-upregulated PD-L1 expression through RAS signaling-driven downregulation of PTRH1 leads to suppression of T cell cytotoxic function in tumor environment'

#singleR

library(SingleR)
singler <- SingleR(test = scrna, 
                         ref =HPCA=hpca.se,
                         labels = hpca.se$label.main)

#seurat
#T cells
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)

pbmc <- CreateSeuratObject(counts = scrna[,rownames(singler)[singler$label == 'T_cells']],project = "seurat", min.cells = 3, min.features = 200, names.delim = "_",)
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
pdf(file="04.featureViolin.pdf",width=10,height=6) 
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 50 & percent.mt < 10) 
save(pbmc,file = 'SeuratObject.Rdata')
pdf(file="04.featureCor.pdf",width=10,height=6) 
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",,pt.size=1.5)
CombinePlots(plots = list(plot1, plot2))
dev.off()

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)
top10 <- head(x = VariableFeatures(object = pbmc), 10)
pdf(file="04.featureVar.pdf",width=10,height=6)     
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()

pbmc=ScaleData(pbmc)             
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc)) 
pdf(file="05.pcaGene.pdf",width=10,height=8)
VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca",nfeatures = 20)
dev.off()
pdf(file="05.PCA.pdf",width=6.5,height=6)
DimPlot(object = pbmc, reduction = "pca")
dev.off()
pdf(file="05.pcaHeatmap.pdf",width=10,height=8)
DimHeatmap(object = pbmc, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)
dev.off()
pbmc <- JackStraw(object = pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)
pdf(file="05.pcaJackStraw.pdf",width=8,height=6)
JackStrawPlot(object = pbmc, dims = 1:20)
dev.off()

pdf(file="05.ElbowPlot.pdf",width=8,height=6)
ElbowPlot(pbmc)
dev.off()
pcSelect=20
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect) 
pbmc <- FindClusters(object = pbmc, resolution = 0.5) 
pbmc <- RunTSNE(object = pbmc, dims = 1:pcSelect)
pdf(file="06.TSNE.pdf",width=6.5,height=6)
TSNEPlot(object = pbmc, label = TRUE) 
dev.off()
write.table(pbmc$seurat_clusters,file="06.tsneCluster.txt",quote=F,sep="/",col.names=F)
logFCfilter=0.5
adjPvalFilter=0.05
pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
save(pbmc,file = 'pbmc.Rdata')
sig.markers=pbmc.markers[abs(as.numeric(as.vector(pbmc.markers$avg_logFC)))]>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
pdf(file="06.tsneHeatmap.pdf",width=10,height=20)
DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()
dev.off()

pdf(file="06.markerScatter.pdf",width=10,height=9)
FeaturePlot(object = pbmc, features = c('CTLA4','PDCD1','LAG3','TIGIT','HAVCR2','CD96','BTLA'),cols = c("skyblue", "red"))
dev.off()
pdf(file="06.cluster0_markerBubble.pdf",width=8,height=4)
cluster10Marker=c('PTPRC','CD2','CD3D','CD3E','CD3G','CD4','CD8A','CD8B','IL2RA','FOXP3','SELL')
DotPlot(object = pbmc, features = cluster10Marker)
dev.off()
pdf(file="06.cluster_markerBubble.pdf",width=8,height=8)
cluster10Marker=c('FXYD2','ANXA4','KRT18','SLC4A4','DEFB1','KRT8','MMP7','SPP1','TFF2','TFF1','TFF3','LYZ')
DotPlot(object = pbmc, features = cluster10Marker)
dev.off()
pdf(file="06.cluster2_markerBubble.pdf",width=8,height=8)
cluster10Marker=c('RAMP2','PLVAP','CLDN5','RGS5','ADIRF','NDUFA4L2')
DotPlot(object = pbmc, features = cluster10Marker)
dev.off()
pdf(file="06.cluster1_markerBubble.pdf",width=8,height=8)
cluster10Marker=c('COL1A1','COL3A1','COL1A2','C1QC','C1QB','CD1C')
DotPlot(object = pbmc, features = cluster10Marker)
dev.off() 
pdf(file="06.cluster3_markerBubble.pdf",width=8,height=8)
cluster10Marker=c('CD79A','IGJ','MZB1','DERL3','FKBP11','IGLL5','SSR4','VPREB3')
DotPlot(object = pbmc, features = cluster10Marker)
dev.off() 
pdf(file="06.immunecheckpointsBubble.pdf",width=8,height=8)
cluster10Marker=c('CTLA4','PDCD1','LAG3','TIGIT','HAVCR2','CD96','BTLA')
DotPlot(object = pbmc, features = cluster10Marker)
dev.off() 
pdf(file="06.Acinar.pdf",width=8,height=8)
cluster10Marker=c('CLPS','CPA1','CTRB1','PRSS1','CELA3A','CPB1','PNLIP')
DotPlot(object = pbmc, features = cluster10Marker)
dev.off() 
pdf(file="06.gene.pdf",width=8,height=8)
cluster10Marker=c('TCF7','CXCR5','IL7R','NR4A1','NR4A2','NR4A3','TOX','TOX2','TOX3','TOX4')
DotPlot(object = pbmc, features = cluster10Marker)
dev.off() 
pdf(file="06.icpl.pdf",width=8,height=8)
cluster10Marker=c('CD274','PDCD1LG2','CD80','CD86','LGALS3','PVR','TNFRSF14','CEACAM1')
DotPlot(object = pbmc, features = cluster10Marker)
dev.off() 
pdf(file="06.M1M2.pdf",width=8,height=8)
cluster10Marker=c('TRIM2','CD163')
DotPlot(object = pbmc, features = cluster10Marker)
dev.off() 
pdf(file="06.fibro.pdf",width=8,height=8)
cluster10Marker=c('CXCL12','CXCR4','CCL7','TGFB1','HGF','IGF','CTGF','IL6','IL8','IL10','CCL2','CCL5','CXCL9','CXCL10','IL1R','IL1B','LIF','IL11')
DotPlot(object = pbmc, features = cluster10Marker)+theme(axis.text.x=element_text(angle = 90))
dev.off() 

#epithelial cells
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)

pbmc <- CreateSeuratObject(counts = scrna[,rownames(singler)[singler$label == 'Epithelial_cells']],project = "seurat", min.cells = 3, min.features = 200, names.delim = "_",)
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
pdf(file="04.featureViolin.pdf",width=10,height=6) 
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 50 & percent.mt < 10) 
save(pbmc,file = 'SeuratObject.Rdata')
pdf(file="04.featureCor.pdf",width=10,height=6) 
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",,pt.size=1.5)
CombinePlots(plots = list(plot1, plot2))
dev.off()

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)
top10 <- head(x = VariableFeatures(object = pbmc), 10)
pdf(file="04.featureVar.pdf",width=10,height=6)     
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()

pbmc=ScaleData(pbmc)             
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc)) 
pdf(file="05.pcaGene.pdf",width=10,height=8)
VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca",nfeatures = 20)
dev.off()
pdf(file="05.PCA.pdf",width=6.5,height=6)
DimPlot(object = pbmc, reduction = "pca")
dev.off()
pdf(file="05.pcaHeatmap.pdf",width=10,height=8)
DimHeatmap(object = pbmc, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)
dev.off()
pbmc <- JackStraw(object = pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)
pdf(file="05.pcaJackStraw.pdf",width=8,height=6)
JackStrawPlot(object = pbmc, dims = 1:20)
dev.off()

pdf(file="05.ElbowPlot.pdf",width=8,height=6)
ElbowPlot(pbmc)
dev.off()
pcSelect=20
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect) 
pbmc <- FindClusters(object = pbmc, resolution = 0.5) 
pbmc <- RunTSNE(object = pbmc, dims = 1:pcSelect)
pdf(file="06.TSNE.pdf",width=6.5,height=6)
TSNEPlot(object = pbmc, label = TRUE) 
dev.off()
write.table(pbmc$seurat_clusters,file="06.tsneCluster.txt",quote=F,sep="/",col.names=F)
logFCfilter=0.5
adjPvalFilter=0.05
pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
save(pbmc,file = 'pbmc.Rdata')
sig.markers=pbmc.markers[abs(as.numeric(as.vector(pbmc.markers$avg_logFC)))]>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
pdf(file="06.tsneHeatmap.pdf",width=10,height=20)
DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()
dev.off()

pdf(file="06.markerScatter.pdf",width=10,height=9)
FeaturePlot(object = pbmc, features = c('CTLA4','PDCD1','LAG3','TIGIT','HAVCR2','CD96','BTLA'),cols = c("skyblue", "red"))
dev.off()
pdf(file="06.cluster0_markerBubble.pdf",width=8,height=4)
cluster10Marker=c('PTPRC','CD2','CD3D','CD3E','CD3G','CD4','CD8A','CD8B','IL2RA','FOXP3','SELL')
DotPlot(object = pbmc, features = cluster10Marker)
dev.off()
pdf(file="06.cluster_markerBubble.pdf",width=8,height=8)
cluster10Marker=c('FXYD2','ANXA4','KRT18','SLC4A4','DEFB1','KRT8','MMP7','SPP1','TFF2','TFF1','TFF3','LYZ')
DotPlot(object = pbmc, features = cluster10Marker)
dev.off()
pdf(file="06.cluster2_markerBubble.pdf",width=8,height=8)
cluster10Marker=c('RAMP2','PLVAP','CLDN5','RGS5','ADIRF','NDUFA4L2')
DotPlot(object = pbmc, features = cluster10Marker)
dev.off()
pdf(file="06.cluster1_markerBubble.pdf",width=8,height=8)
cluster10Marker=c('COL1A1','COL3A1','COL1A2','C1QC','C1QB','CD1C')
DotPlot(object = pbmc, features = cluster10Marker)
dev.off() 
pdf(file="06.cluster3_markerBubble.pdf",width=8,height=8)
cluster10Marker=c('CD79A','IGJ','MZB1','DERL3','FKBP11','IGLL5','SSR4','VPREB3')
DotPlot(object = pbmc, features = cluster10Marker)
dev.off() 
pdf(file="06.immunecheckpointsBubble.pdf",width=8,height=8)
cluster10Marker=c('CTLA4','PDCD1','LAG3','TIGIT','HAVCR2','CD96','BTLA')
DotPlot(object = pbmc, features = cluster10Marker)
dev.off() 
pdf(file="06.Acinar.pdf",width=8,height=8)
cluster10Marker=c('CLPS','CPA1','CTRB1','PRSS1','CELA3A','CPB1','PNLIP')
DotPlot(object = pbmc, features = cluster10Marker)
dev.off() 
pdf(file="06.gene.pdf",width=8,height=8)
cluster10Marker=c('TCF7','CXCR5','IL7R','NR4A1','NR4A2','NR4A3','TOX','TOX2','TOX3','TOX4')
DotPlot(object = pbmc, features = cluster10Marker)
dev.off() 
pdf(file="06.icpl.pdf",width=8,height=8)
cluster10Marker=c('CD274','PDCD1LG2','CD80','CD86','LGALS3','PVR','TNFRSF14','CEACAM1')
DotPlot(object = pbmc, features = cluster10Marker)
dev.off() 
pdf(file="06.M1M2.pdf",width=8,height=8)
cluster10Marker=c('TRIM2','CD163')
DotPlot(object = pbmc, features = cluster10Marker)
dev.off() 
pdf(file="06.fibro.pdf",width=8,height=8)
cluster10Marker=c('CXCL12','CXCR4','CCL7','TGFB1','HGF','IGF','CTGF','IL6','IL8','IL10','CCL2','CCL5','CXCL9','CXCL10','IL1R','IL1B','LIF','IL11')
DotPlot(object = pbmc, features = cluster10Marker)+theme(axis.text.x=element_text(angle = 90))
dev.off() 


#DEG
library(limma)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

group_list=data$sample
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
exprSet=t(data[,2:15])
rownames(design)=colnames(exprSet)
contrast.matrix<-makeContrasts("HG-NC",
                               levels = design)

deg = function(exprSet,design,contrast.matrix){
  fit <- lmFit(exprSet,design)
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  fit2 <- eBayes(fit2)  
  tempOutput = topTable(fit2, coef=1, n=Inf)
  nrDEG = na.omit(tempOutput) 
  head(nrDEG)
  return(nrDEG)
}

deg = deg(exprSet,design,contrast.matrix)
deg <- deg[rownames(Xcell1),]
#single_cell_seq
Tmean <- apply(data[,rownames(cell)][,which(cell$CD8=='low')], 1, mean)
Nmean <- apply(data[,rownames(cell)][,which(cell$CD8=='high')], 1, mean)
deg <- deg[rownames(data),]
deg <- cbind(deg,Tmean,Nmean)
deg$FC=deg$Tmean/deg$Nmean
write.csv(deg,file = 'CD8_LOW_VS_HIGH_tumorcelldeg.csv')

#KEGG/GO analysis
head(deg)
logFC_t=2
deg$g=ifelse(deg$P.Value>0.05,'stable',
             ifelse( deg$logFC > logFC_t,'UP',
                     ifelse( deg$logFC < -logFC_t,'DOWN','stable')))

deg$g=ifelse(deg$adj.P.Val>0.05,'stable',
             ifelse( deg$FC > logFC_t,'up',
                     ifelse( deg$FC < 1/logFC_t,'down','stable')))


table(deg$g)
head(deg)
deg$symbol=rownames(deg)
df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
head(df)
DEG=deg
head(DEG)

DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
head(DEG)
save(DEG,file = 'anno_DEG.Rdata')


gene_up= DEG[DEG$g == 'UP','ENTREZID'] 
gene_down=DEG[DEG$g == 'DOWN','ENTREZID'] 
gene_diff=c(gene_up,gene_down)
gene_all=as.character(DEG[ ,'ENTREZID'] )
data(geneList, package="DOSE")
head(geneList)
boxplot(geneList)
boxplot(DEG$logFC)

geneList=DEG$logFC
names(geneList)=DEG$ENTREZID
geneList=sort(geneList,decreasing = T)

kk.up <- enrichKEGG(gene         = gene_up,
                    organism     = 'hsa',
                    universe     = gene_all,
                    pvalueCutoff = 0.9,
                    qvalueCutoff =0.9)
head(kk.up)[,1:6]
dotplot(kk.up );ggsave('kk.up.dotplot.png')
kk.down <- enrichKEGG(gene         =  gene_down,
                      organism     = 'hsa',
                      universe     = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff =0.9)
head(kk.down)[,1:6]
dotplot(kk.down );ggsave('kk.down.dotplot.png')
kk.diff <- enrichKEGG(gene         = gene_diff,
                      organism     = 'hsa',
                      pvalueCutoff = 0.05)
head(kk.diff)[,1:6]
dotplot(kk.diff );ggsave('kk.diff.dotplot.png')

kegg_diff_dt <- as.data.frame(kk.diff)
kegg_down_dt <- as.data.frame(kk.down)
kegg_up_dt <- as.data.frame(kk.up)
down_kegg<-kegg_down_dt[kegg_down_dt$pvalue<0.05,];down_kegg$group=-1
up_kegg<-kegg_up_dt[kegg_up_dt$pvalue<0.05,];up_kegg$group=1

kegg_plot <- function(up_kegg,down_kegg){
  dat=rbind(up_kegg,down_kegg)
  colnames(dat)
  dat$pvalue = -log10(dat$pvalue)
  dat$pvalue=dat$pvalue*dat$group 
  
  dat=dat[order(dat$pvalue,decreasing = F),]
  
  g_kegg<- ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
    geom_bar(stat="identity") + 
    scale_fill_gradient(low="blue",high="red",guide = FALSE) + 
    scale_x_discrete(name ="Pathway names") +
    scale_y_continuous(name ="log10P-value") +
    coord_flip() + theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
    ggtitle("Pathway Enrichment") 
}

g_kegg=kegg_plot(up_kegg,down_kegg)
print(g_kegg)

ggsave(g_kegg,filename = 'kegg_up_down.png',width = 6,height = 5,limitsize = F)

kk_gse <- gseKEGG(geneList     = geneList,
                  organism     = 'hsa',
                  nPerm        = 1000,
                  minGSSize    = 120,
                  pvalueCutoff = 0.9,
                  verbose      = FALSE)
head(kk_gse)[,1:6]
kk_gsea=gseaplot(kk_gse, geneSetID = rownames(kk_gse[1,]))
ggsave(kk_gsea,filename = 'kk_gsea.png')

down_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore < 0,];down_kegg$group=-1
up_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore > 0,];up_kegg$group=1

g_kegg=kegg_plot(up_kegg,down_kegg)
print(g_kegg)
ggsave(g_kegg,filename = 'kegg_up_down_gsea.png')

g_list=list(gene_up=gene_up,
            gene_down=gene_down,
            gene_diff=gene_diff)

go_enrich_results <- lapply( g_list , function(gene) {
  lapply( c('BP','MF','CC') , function(ont) {
    cat(paste('Now process ',ont ))
    ego <- enrichGO(gene          = gene,
                    universe      = gene_all,
                    OrgDb         = org.Hs.eg.db,
                    ont           = ont ,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.99,
                    qvalueCutoff  = 0.99,
                    readable      = TRUE)
    
    print( head(ego) )
    return(ego)
  })
})

n1= c('gene_up','gene_down','gene_diff')
n2= c('BP','MF','CC') 
for (i in 1:3){
  for (j in 1:3){
    fn=paste0('dotplot_',n1[i],'_',n2[j],'.png')
    cat(paste0(fn,'/n'))
    png(fn,res=150,width = 1300)
    print( dotplot(go_enrich_results[[i]][[j]] ))
    dev.off()
  }
}

