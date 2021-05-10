#为了节省时间以及受到设备性能影响，首先先用python将每个类器官样本中的rep1提取出来
#并提取文件第一列作为gene文件，提取文件第一行并拆分成barcode_rep1（含有血管化类器官的barcode）和barcode（含有普通类器官的barcode）
library(data.table)
library(Seurat)
library(scater)
library(Matrix)
#分别导入四个类器官样本（仅含有数据，不含有行名和列名）
huvec_d65=fread("F:/graduate/GSE131094/huvec_d65_rep1.txt",header = F)
huvec_d100=fread("F:/graduate/GSE131094/huvec_d100_rep1.txt",header = F)
h9_d65=fread("F:/graduate/GSE131094/h9_d65_rep1.txt",header = F)
h9_d100=fread("F:/graduate/GSE131094/h9_d100_rep1.txt",header = F)
#导入基因文件
gene=fread("F:/graduate/GSE131094/gene.txt",encoding="UTF-8",header=F)
#将基因文件作为四个样本数据的行名
tmpb=data.frame(gene,huvec_d65)
tmpc=data.frame(gene,huvec_d100)
tmpd=data.frame(gene,h9_d65)
tmpe=data.frame(gene,h9_d100)
rownames(tmpb)=tmpb[,1]
rownames(tmpc)=tmpc[,1]
rownames(tmpd)=tmpd[,1]
rownames(tmpe)=tmpe[,1]
tmpb=tmpb[,-1]
tmpc=tmpc[,-1]
tmpd=tmpd[,-1]
tmpe=tmpe[,-1]
#导入barcode文件
barcode=fread("F:/graduate/GSE131094/barcode_rep1.txt",header = F)
barcode1=fread("F:/graduate/GSE131094/barcode.txt",header = F)
#将barcode作为数据的列名
barcode_b=barcode[,1:2085]
barcode_c=barcode[,2086:6852]
barcode_d=barcode1[,1:1537]
barcode_e=barcode1[,1538:5844]
#至此，构建seurat对象的准备工作已全部完成


#分别构建4个样本的seurat对象
colnames(tmpb)=barcode_b
colnames(tmpc)=barcode_c
colnames(tmpd)=barcode_d
colnames(tmpe)=barcode_e
#构建huvec_d65（培养65天血管化类器官）seurat对象
meta=as.data.frame(colnames(tmpb))
colnames(meta)=c('cell name1')
rownames(meta)=colnames(tmpb)
head(meta)
huvec_d65 <- CreateSeuratObject(counts = tmpb,
                                meta.data = meta,
                                min.cells = 3, min.features = 200,project="huvec_d65")
#构建huvec_d100（培养100天血管化类器官）seurat对象
meta=as.data.frame(colnames(tmpc))
colnames(meta)=c('cell name2')
rownames(meta)=colnames(tmpc)
head(meta)
huvec_d100 <- CreateSeuratObject(counts = tmpc,
                                 meta.data = meta,
                                 min.cells = 3, min.features = 200,project = 'huvec_d100')

#构建h9_d65（培养65天普通类器官）seurat对象
meta=as.data.frame(colnames(tmpd))
colnames(meta)=c('cell name3')
rownames(meta)=colnames(tmpd)
head(meta)
h9_d65 <- CreateSeuratObject(counts = tmpd,
                             meta.data = meta,
                             min.cells = 3, min.features = 200,project="h9_d65")
#构建h9_d100（培养100天普通类器官）seurat对象
meta=as.data.frame(colnames(tmpe))
colnames(meta)=c('cell name4')
rownames(meta)=colnames(tmpe)
head(meta)
h9_d100 <- CreateSeuratObject(counts = tmpe,
                              meta.data = meta,
                              min.cells = 3, min.features = 200,project="h9_d100")

rm(barcode,barcode_d,barcode_e,gene,meta,tmpd,tmpe,barcode1,barcode_b,barcode_c,tmpb,tmpc)

#合并数据
data <- merge(huvec_d65,huvec_d100, add.cell.ids=c("huvec_d65",'huvec_d100'))
data1 <- merge(h9_d65,h9_d100, add.cell.ids=c("h9_d65",'h9_d100'))
data_all <- merge(data,data1, add.cell.ids=c("huvec",'h9'))


#数据的过滤与质量控制
#正则表达式，表示以MT-开头；test.seu[["percent.mt"]]这种写法会在meta.data矩阵加上一列
data_all[["percent.mt"]] <- PercentageFeatureSet(data_all, pattern = "^MT-")  
data_all <- subset(data_all, subset = 
                     nFeature_RNA < 6000 & 
                     percent.mt < 10 & 
                     nFeature_RNA > 200)
#画小提琴图
VlnPlot(data_all, features = c("nFeature_RNA", "percent.mt"), ncol = 3,pt.size = 0)
#画散点图
plot2 <- FeatureScatter(data_all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2
data_all <- FindVariableFeatures(data_all, selection.method = "vst", nfeatures = 2000)
data_all
# 确定10个最高变基因
top10 <- head(VariableFeatures(data_all), 10)
top10
# 筛选可变基因并作图
plot1 <- VariableFeaturePlot(data_all)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 
plot2
#接下来进行数据缩放
all.genes <- rownames(data_all)
data_all <- ScaleData(data_all, features = all.genes)
View(all.genes)
View(data_all[["RNA"]]@scale.data)


#接下来进行线性降维
##对PCA进行可视化操作的方法
data_all <- RunPCA(data_all, features = VariableFeatures(object = data_all))##进行PCA操作
##对PCA进行可视化操作的方法
print(data_all[["pca"]], dims = 1:20, nfeatures = 10)
VizDimLoadings(data_all, dims = 19:20, reduction = "pca",nfeatures=10)##分别查看对PC1和PC2其主要作用的基因
DimPlot(data_all, reduction = "pca")###通过点图展示PC1和PC2

DimHeatmap(data_all, dims = 1:20, cells = 500,balanced = TRUE)##通过热图可以对PC成分进行可视化，允许我们对基因集的异质性进行审阅，从而决定哪个PC可用于后续的下游分析，一般默认可以选择500个细胞，细胞和基因都是按照PCA的分值进行排序的???
VizDimLoadings(data_all, dims = 1:20, reduction = "pca",nfeatures = 30)##分别查看对PC1和PC2其主要作用的基因


#接下来判断选择哪个PC
#way1
ElbowPlot(data_all)
#way2
data_all<-JackStraw(data_all,num.replicate=100)
data_all<-ScoreJackStraw(data_all,dim=1:20)
JackStrawPlot(data_all,dims=1:20)
#选了前10个PC
#细胞的聚类分析 
data_all <- FindNeighbors(data_all, dims = 1:10)
data_all <- FindClusters(data_all, resolution = 1)

#UMAP/tSNE的非线性降维分析 
data_all <- RunUMAP(data_all, dims = 1:10,label = TRUE)
DimPlot(data_all, reduction = "umap")
library(cowplot)
data_all<- RunTSNE(data_all, dims = 1:5)
p1 <- DimPlot(data_tsne, reduction = "tsne", pt.size=0.5)
p1


#发现存在一定的批次效应，于是进行去批次
library(Seurat)
library(tidyverse)
### h9_d65 ----
h9_d65 <- FindVariableFeatures(h9_d65, selection.method = "vst", nfeatures = 2000)
### h9_d100 ----
h9_d100 <- FindVariableFeatures(h9_d100, selection.method = "vst", nfeatures = 2000)
### huvec_d65 ----
huvec_d65 <- FindVariableFeatures(huvec_d65, selection.method = "vst", nfeatures = 2000)
### huvec_d100 ----
huvec_d100 <- FindVariableFeatures(huvec_d100, selection.method = "vst", nfeatures = 2000)
### Integration ----object.list参数是由多个Seurat对象构成的列
data_all.anchors <- FindIntegrationAnchors(object.list = list(h9_d65,h9_d100,huvec_d65,huvec_d100), dims = 1:20)
data_all.integrated <- IntegrateData(anchorset = data_all.anchors, dims = 1:20)
data_all.integrated
#这一步之后就多了一个整合后的assay（原先有一个RNA的assay），整合前后的数据分别存储在这两个assay
DefaultAssay(data_all.integrated) <- "integrated"
#后续仍然是标准流程，基于上面得到的整合data矩阵Run the standard workflow for visualization and clustering
data_all.integrated <- ScaleData(data_all.integrated, features = rownames(data_all.integrated))
data_all.integrated <- RunPCA(data_all.integrated, npcs = 50, verbose = FALSE)
data_all.integrated <- FindNeighbors(data_all.integrated, dims = 1:5)
data_all.integrated <- FindClusters(data_all.integrated, resolution = 0.5)
data_all.integrated <- RunUMAP(data_all.integrated, dims = 1:30)
data_all.integrated <- RunTSNE(data_all.integrated, dims = 1:30)

data_all=data_all.integrated
#再次画图，分cluster
p1 <- DimPlot(data_all, reduction = "umap", pt.size=0.5)
p1
data_tsne<- RunTSNE(data_all, dims = 1:10)
p2<- DimPlot(data_tsne, reduction = "tsne", pt.size=0.5)
p2
#再次画图，图中分样本显示
p3 <- DimPlot(data_tsne, reduction = "tsne", group.by = "orig.ident",label.size = 0, pt.size=0.5)
p3
p4 <- DimPlot(data_all,reduction = "umap", group.by = "orig.ident",   pt.size=0.5, label = F,repel = TRUE)
p4

#注释organoids
celltype_marker=c(
  "SOX2",#RG
  "FAM107A",#oRG
  "MKI67",#cell cycle
  "EOMES",#IPC
  "DCX",#immature neuron
  "NEUROD2",#excitatory neuron
  "GAD1",#interneuron
  "AIF1",#microglia
  "AQP4",#astrocyte
  "OLIG1",#oligodendrocyte
  "RSPO2"#choroid plexus
)
#看表达量进行初步判断
plot1=VlnPlot(data_all,features = celltype_marker,pt.size = 0,ncol = 2)
plot1
0: excitatory neuron
1: immature neuron
2: oRG
3: unknown
4: IPC
5: immature neuron
6: excitatory neuron
7: cell cycle
8: choroid plexus
9: oRG
10: microglia
11: RG
12: excitatory neuron
13: RG
14: RG
15: unknown
16: immature neuron
#列表显示差异基因表达量，验证上述判断
markers <- FindAllMarkers(data_all, logfc.threshold = 0.25, min.pct = 0.1, 
                          only.pos = TRUE, data_all = "wilcox")
rm(h9_d100,h9_d65,huvec_d100,huvec_d65,data1,data)
rm(data_tsne)
library(dplyr)
markers_df = markers %>% group_by(cluster) %>% top_n(n = 500, wt = avg_logFC)
#加入注释信息
excitatory_neuron=c(0,6,12,3,15)
immature_neuron=c(1,5,16)
oRG=c(2,9)
cell_cycle=c(7)
IPC=c(4)
choroid_plexus=c(8)
microglia=c(10)
RG=c(11,13,14)

current.cluster.ids <- c(excitatory_neuron,
                         immature_neuron,
                         oRG,
                         cell_cycle,
                         IPC,
                         choroid_plexus,
                         microglia,
                         RG
)

new.cluster.ids <- c(rep("excitatory_neuron",length(excitatory_neuron)),
                     rep("immature_neuron",length(immature_neuron)),
                     rep("oRG",length(oRG)),
                     rep("cell_cycle",length(cell_cycle)),
                     rep("IPC",length(IPC)),
                     rep("choroid_plexus",length(choroid_plexus)),
                     rep("microglia",length(microglia)),
                     rep("RG",length(RG))
)

data_all@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(data_all@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
#统计一下各种细胞的数目
table(data_all@meta.data$celltype)
head(data_all@meta.data)
#画图
plotCB=as.data.frame(data_all@meta.data)[,"cell.name1"]
data_all
celltype_tsne=DimPlot(data_all, reduction = "tsne", group.by = "celltype", pt.size=0.5)
celltype_tsne

#完成各种细胞分别的计数，在excel中合并sample_type和celltype
write.csv(data_all@meta.data,"F:/graduate/result/data_all_organoids.txt")
data_all@meta.data
data_all_organoids=read.table("F:/graduate/result/data_all_organoids.txt")
#统计一下各种细胞的数目
table(data_all_organoids$sample_type_celltype)
head(data_all@meta.data)
#保存图片成pdf和seurat对象
ggsave(filename = "tsne5_nounknown.pdf", width = 15, height = 12, units = 'cm')
saveRDS(data_all,file = "data_all.seu.rds") #保存test.seu对象，下次可以直接调用，不用重复以上步骤
celltype_umap=DimPlot(data_all, reduction = "umap", group.by = "celltype", pt.size=0.5)
celltype_umap
ggsave(filename = "umap_nounknown.pdf", width = 15, height = 12, units = 'cm')
saveRDS(celltype_tsne,file = "celltype_tsne.rds")
saveRDS(celltype_umap,file = "celltype_umap.rds")

#后续按照相同步骤处理human数据
#随后进行GWEA，获得热图及rank gene list
#使用python将human和organoids的数据分别合并为human_ranked.tsv和organoids_ranked.tsv
human=read.table("F:/graduate/GSE131094/rank_gene/human/human_ranked.tsv",header=T)
organoids=read.table("F:/graduate/GSE131094/rank_gene/organoids/organoids_ranked.tsv",header=T)
#寻找共同基因
keep_genes <- Reduce(intersect, list(human[,1],organoids[,1]))
human_try <- human[match(keep_genes,human[,1] ), ]
organoids_try<- organoids[match(keep_genes, organoids[,1]), ]
write.csv(human_try, file="F:/graduate/GSE131094/rank_gene/human/.txt")
write.csv(organoids_try, file="F:/graduate/GSE131094/rank_gene/organoids/.txt")
data=read.table("F:/graduate/GSE131094/rank_gene/organoids/ranked_gene_list_1_versus_0_1613656798605.tsv",header=T)
data <- data[match(keep_genes,data[,1] ), ]
write.csv(data, file="F:/graduate/GSE131094/rank_gene/organoids_commongene/ranked_gene_list_1_versus_0_1613656798605.txt")
organoid=read.table("F:/graduate/GSE131094/rank_gene/organoids_commongene/ranked_gene_list_1_versus_0_1613656798605.txt",header=T,sep=',')
#合并所有human数据
myfolderpath <- "F:/graduate/GSE131094/rank_gene/human_commongene_tryagain"
setwd(myfolderpath)
allfoldN <- dir(myfolderpath)
foldnewN <- paste0('nF',1:length(allfoldN),".txt")
file.rename(allfoldN,foldnewN)
a=read.table("F:/graduate/GSE131094/rank_gene/nF1_human.txt",header=T,sep=",")
a=a[,-1]
for(i in 2:length(allfoldN)){
  myfilepath <- paste(myfolderpath,foldnewN[i],sep='/')
  tmp=read.table(myfilepath,sep=",",header=T)
  tmp=tmp[,-1:-3]
  a=data.frame(a,tmp)
}
b=read.table("F:/graduate/GSE131094/rank_gene/time_human.txt")
colnames(a)=b[,1]
write.csv(a,"F:/graduate/GSE131094/rank_gene/human_all.txt")

myfolderpath <- "F:/graduate/GSE131094/rank_gene/organoids_commongene"
setwd(myfolderpath)
allfoldN <- dir(myfolderpath)
foldnewN <- paste0('nF',1:length(allfoldN),".txt")
file.rename(allfoldN,foldnewN)
a=read.table("F:/graduate/GSE131094/rank_gene/nF1_organoids.txt",header=T,sep=",")
a=a[,-1]
for(i in 2:length(allfoldN)){
  myfilepath <- paste(myfolderpath,foldnewN[i],sep='/')
  tmp=read.table(myfilepath,sep=",",header=T)
  tmp=tmp[,-1:-3]
  a=data.frame(a,tmp)
}
b=read.table("F:/graduate/GSE131094/rank_gene/time_organoids.txt")
colnames(a)=b[,1]
write.csv(a,"F:/graduate/GSE131094/rank_gene/human_all.txt")
#至此获得human和organoids的分数排布
#在excel中使用Z-score对分数进行标准化

#最后进行MDS
organoids1=read.table("F:/graduate/GSE131094/rank_gene/organoids_all_z_score.txt",header = T,row.names = 1)
human1=read.table("F:/graduate/GSE131094/rank_gene/human_all_z_score.txt",header = T,row.names = 1)
# Before running MDS, we first calculate a distance matrix between all pairs of cells.  Here we use a simple euclidean distance metric on all genes, using scale.data as input
data_all=data.frame(organoids1,human1)
d <- dist(t(data_all))
# Run the MDS procedure, k determines the number of dimensions
mds <- cmdscale(d = d, k =2)
head(mds)
colnames(mds) <- paste0("MDS_", 1:2)
data=data.frame(mds)
write.csv(data,"F:/graduate/GSE131094/rank_gene/mds_all.txt")

data=read.table("F:/graduate/GSE131094/rank_gene/mds_all.txt",header = T,row.names = 1)

ggplot(data,aes(x=data[,1],y=data[,2],color=row.names(data)))+geom_point(size=3.5,shape=data[,3])+xlab("") + ylab("")+
  geom_text(aes(label = row.names(data),color=row.names(data)),nudge_x=0.5,nudge_y=3.5,size=3)+theme(legend.position="none")

distance=read.table("F:/graduate/GSE131094/rank_gene/mds_all_notype.txt",header = T,row.names = 1)
dis<-dist(distance)
dis1=matrix(dis)
write.csv(dis1,"F:/graduate/GSE131094/rank_gene/mds_all_notype_dis.txt")

#制作都从GW8开始的图
data=read.table("F:/graduate/GSE131094/rank_gene/mds_all/mds_all_GW8.txt",header = T,row.names = 1)

ggplot(data,aes(x=data[,1],y=data[,2],color=row.names(data)))+geom_point(size=4,shape=data[,3])+xlab("") + ylab("")+
  geom_text(aes(label = row.names(data),color=row.names(data)),nudge_x=0.5,nudge_y=3.5,size=3)+theme(legend.position="none")

distance=read.table("F:/graduate/GSE131094/rank_gene/mds_all_notype.txt",header = T,row.names = 1)
dis<-dist(distance)
dis1=matrix(dis)
write.csv(dis1,"F:/graduate/GSE131094/rank_gene/mds_all_notype_dis.txt")
