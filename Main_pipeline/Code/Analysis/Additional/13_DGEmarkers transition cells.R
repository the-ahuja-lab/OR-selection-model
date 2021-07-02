### analysis for sidrah data


setwd("G:/My Drive/projects/WTA/DGE for transition")
library(Seurat)
library(dplyr)
library(Matrix)
library(GSEABase)

data<-read.csv("clustered_an.csv", header= TRUE, sep=',',row.names = as.numeric(TRUE),stringsAsFactors = FALSE)


#data<-t(data) # use  when samples are in rows

pdf("WTA analysis.pdf")

pb <<- CreateSeuratObject(counts = data, min.cells = 3, min.features = 200,project = "scRNAseq")

VlnPlot(pb, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
pb <<- FindVariableFeatures(pb, selection.method = "vst", nfeatures = 2000)
#list_of_variable_features<<-VariableFeatures(pb)
#list_of_variable_features<<-as.data.frame(list_of_variable_features)
top10 <<- head(VariableFeatures(pb), 20)
pb <<- NormalizeData(pb)

pb <<- FindVariableFeatures(pb, selection.method = "vst", nfeatures = 2000)

top10 <<- head(VariableFeatures(pb), 20)
plot1 <- VariableFeaturePlot(pb)
plot2 <- LabelPoints(plot = plot1, points = top10,repel=TRUE,xnudge = 0,ynudge = 0)
plot(plot2) 
all.genes <- rownames(pb)
pb <<- ScaleData(pb, features = all.genes)
pb <<- RunPCA(pb, features = VariableFeatures(object = pb))
############### addition of meta data#############
meta<-read.csv("metadata.csv", header= TRUE, sep=',',row.names = as.numeric(TRUE),stringsAsFactors = FALSE)
meta<-t(meta)
meta<-as.data.frame(meta)
meta$Samples<-rownames(meta)
pb<<-AddMetaData(pb,meta)
Idents(pb)<-pb@meta.data$cluster
##########################

DimPlot(pb, reduction = "pca", pt.size = 4, group.by = "ident")

DimHeatmap(pb, dims = 1:6, cells = 500, balanced = TRUE)
pb <<- JackStraw(pb, num.replicate = 100)
pb <<- ScoreJackStraw(pb, dims = 1:20)
JackStrawPlot(pb, dims = 1:15)
ElbowPlot(pb)
pb <<- FindNeighbors(pb, dims = 1:10)
pb <<- FindClusters(pb, resolution = 0.5,algorithm=2)
Idents(pb)<-pb@meta.data$cluster
head(Idents(pb), 5)

pb <<- RunUMAP(pb, dims = 1:10)
DimPlot(pb, reduction = "umap",pt.size = 3)

dev.off()






############################ for all cell types

pb.markers <<- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
pb.markers$genes<-rownames(pb.markers)
write.table(pb.markers,file="For all cells markers_wilcox logfc 1.csv",
            sep=",",row.names=FALSE,col.names = TRUE,quote = FALSE)


############################# if cluster identity want to check 
cluster_ident<-pb@active.ident

cluster_ident<-as.data.frame(cluster_ident)
write.table(cluster_ident,file="For all cells markers_cluster_info.csv",
            sep=",",row.names=FALSE,col.names = TRUE,quote = FALSE)


pdf("markers_fpr_transition_cells.pdf")
RidgePlot(object = pb, features = c("Ptn","Stmn4","Gnal","Ckb")) 
dev.off()


