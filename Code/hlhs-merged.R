#hlhs

library(Seurat)
library(ggsci)
library(ggpubr)
library(ggplot2)
library(cowplot)
library(viridis)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(robustSingleCell)

setwd("~/Box/XinZhouFiles/Projects/ZXE7_SingleCell_HLHS/DATA/")
load("./combined/all.Robj")

hlhs.all <- UpdateSeuratObject(all)

DefaultAssay(hlhs.all)
lhv.new <- hlhs.all

lhv.new$oldidents <- Idents(lhv.new)
lhv.new$condition <- lhv.new$orig.ident
lhv.new$condition[lhv.new$orig.ident=="d83_LVp"] <- "Control"
lhv.new$condition[lhv.new$orig.ident=="LVp"] <-  "LHV"

##########################################
lhls.lhv.list <- SplitObject(lhv.new, split.by = "orig.ident")
for (i in 1:length(lhls.lhv.list)) {
  lhls.lhv.list[[i]] <- SCTransform(lhls.lhv.list[[i]], verbose = T , vars.to.regress="nCount_RNA")
}

lhls.lhv.features <- SelectIntegrationFeatures(object.list = lhls.lhv.list, nfeatures = 3000)
lhls.lhv.list <- PrepSCTIntegration(object.list = lhls.lhv.list, anchor.features = lhls.lhv.features, verbose = T)

lhls.lhv.anchors <- FindIntegrationAnchors(object.list = lhls.lhv.list, normalization.method = "SCT", anchor.features = lhls.lhv.features, verbose = T)
lhls.new.merge <- IntegrateData(anchorset = lhls.lhv.anchors, normalization.method = "SCT", verbose = T)

lhls.new.merge <- RunPCA(lhls.new.merge, dim=1:15)
lhls.new.merge <- RunUMAP(lhls.new.merge, dim=1:30)
lhls.new.merge <- FindNeighbors(lhls.new.merge, dims = 1:10)
lhls.new.merge <- FindClusters(lhls.new.merge, resolution = 0.1)
Idents(lhv.new)
pmergeUMAP <- DimPlot(lhls.new.merge, split.by = "condition")
pmergeUMAP

#ggsave2(filename = "../NewAnalysis/F1a.Condition.ulv.pdf", pmergeUMAP, scale=0.77)

pvmerge <- VlnPlot(object = lhls.new.merge, features = c("NPR3","CDH11", "percent.mito","APLN"), ncol=2, combine = T, pt.size=0.05, group.by = "SCT_snn_res.0.1")
pvmerge
ggsave2(filename = "../NewAnalysis/S1B.QC.ulv.pdf", pvmerge, scale=0.77)

lhls.new.merge$celltype2 <- as.character(lhls.new.merge$integrated_snn_res.0.1)
lhls.new.merge$celltype2[lhls.new.merge$celltype2 %in% c(0,2,3)] <- "Endothelium"
lhls.new.merge$celltype2[lhls.new.merge$celltype2 %in% c(1)] <- "Endocardium"
lhls.new.merge$celltype2[lhls.new.merge$celltype2 %in% c(4,5)] <- "Others" 
###########################################


#################################################################################################
Idents(lhls.new.merge) <- "celltype2"
lhv.Endothelium <- subset(lhls.new.merge, idents="Endothelium")

plots4 <- DimPlot(lhv.Endothelium,group.by ="integrated_snn_res.0.1", combine = T) #+ scale_color_aaas()
plots4

ggsave2("../NewAnalysis/before.subset.part.UMAP.pdf", plots4, scale = 0.77)

lhv.Endothelium <- RunPCA(lhv.Endothelium, dim=1:10)
lhv.Endothelium <- RunUMAP(lhv.Endothelium, dim=1:30)
lhv.Endothelium <- FindNeighbors(lhv.Endothelium, dims = 1:15)
lhv.Endothelium <- FindClusters(lhv.Endothelium, resolution = 0.2)

plots5 <- DimPlot(lhv.Endothelium,group.by ="integrated_snn_res.0.2", combine = T) #+ scale_color_aaas()
plots5
ggsave2("../NewAnalysis/after.subset.part.UMAP.pdf", plots5, scale = 0.77)

plots6 <- DimPlot(lhv.Endothelium,group.by ="condition", combine = T) + scale_color_d3()
plots6
ggsave2("../NewAnalysis/after.subset.condition.UMAP.pdf", plots6, scale = 0.77)

lhv.Endothelium.Marker <- FindAllMarkers(lhv.Endothelium)
write.csv(file = "../NewAnalysis/lhv.all.markers.aftermerge.csv", lhv.Endothelium.Marker)

table(lhv.Endothelium$condition)

######
lhv.marker <- read.csv("../NewAnalysis/features/lhv2.markers.csv", header = F)
lhv.marker <- lhv.marker$V1
lhv.marker <- as.character(lhv.marker)

DefaultAssay(lhv.Endothelium) <- "SCT"
for(gene in lhv.marker){
  print(gene)
  p=FeaturePlot(object = lhv.Endothelium, features = c(gene), split.by="condition",cols=c("lightgrey", "#931F21"),pt.size=1.2, sort.cell=T)
  ggsave(p,filename = paste0("../NewAnalysis/features/lhvmerged/",gene,".pdf"),height=6.5,width=12)
}

#run cell cycle analysis 
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

Idents(lhv.Endothelium) <- "integrated_snn_res.0.2"

table(lhv.Endothelium$integrated_snn_res.0.2,lhv.Endothelium$condition)

lhv.Endothelium <-  CellCycleScoring(lhv.Endothelium , s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, assay="RNA")

C0 <- subset(lhv.Endothelium, idents = c(0))
C1 <- subset(lhv.Endothelium, idents = c(1))
C2 <- subset(lhv.Endothelium, idents = c(2))
C3 <- subset(lhv.Endothelium, idents = c(3))
C4 <- subset(lhv.Endothelium, idents = c(4))
C5 <- subset(lhv.Endothelium, idents = c(5))

C0 <-  CellCycleScoring(C0 , s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, assay="RNA")
C1 <-  CellCycleScoring(C1 , s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, assay="RNA")
C2 <-  CellCycleScoring(C2 , s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, assay="RNA")
C3 <-  CellCycleScoring(C3 , s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, assay="RNA")
C4 <-  CellCycleScoring(C4 , s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, assay="RNA")

#Plot cellcycle heatmap
DefaultAssay(lhv.Endothelium) <- "SCT"
lhv.Endothelium <- ScaleData(lhv.Endothelium)
Idents(lhv.Endothelium) <- "old.ident"
cellcyclegene <- c(s.genes,g2m.genes)
p <- DoHeatmap(lhv.Endothelium, features = cellcyclegene, group.by= "old.ident")
p <- p + scale_fill_gradientn(colors=c("white","#ADB6B6FF","grey","#ED0000FF", "#AD002AFF"))
p
#ggsave2("./Heatmap.cellcyclegene.pdf", p, scale=1.2)

table(lhv.Endothelium$celltype3)

lhv.Endothelium$celltype3 <- paste(lhv.Endothelium$old.ident, lhv.Endothelium$condition)

p <- ggplot(C1[[]], aes(x=S.Score, y=G2M.Score)) + geom_point(size=0.2)

#plot <- subset(lhv.Endothelium[[]]) 
          
#cells = sample(x = colnames(seurat), size = 1835)
#seurat <- subset(seurat, cells = sample(x = colnames(seurat), size = 1835))

DimPlot(lhv.Endothelium)

p <- ggplot(subset(lhv.Endothelium[[]],integrated_snn_res.0.2==3), aes(x=S.Score, y=G2M.Score, color=Phase)) + geom_point(size=0.2)
p <- p + facet_wrap(condition~integrated_snn_res.0.2, ncol=1) + ylim(c(-0.25,1.0)) + xlim(-0.3,0.9)
p
p1 <- p + stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="grey", fill="white")
p1 <- p +stat_density_2d(aes(fill = ..density..), geom = "raster", contour = F) +
  scale_fill_gradientn(colors=c("black","white","#FFD31D","#FFD31D","#ED0000FF","#AD002AFF","#AD002AFF","#AD002AFF"))
p1

#scale_fill_gradientn(colors=c("blue","green","green","#FFD31D","#FFD31D","#ED0000FF","#AD002AFF"))
#scale_fill_gradientn(colors=c("#5c2a9d","green","#FFD31D","#FFD31D","#ED0000FF","#AD002AFF"))
#scale_fill_gradientn(colors=c("#5c2a9d","white","#50bda1","#FFD31D","#FFD31D","#ED0000FF","#AD002AFF"))
#scale_fill_gradientn(colors=c("#5c2a9d","white","#50bda1","#FFD31D","#FFD31D","#ED0000FF","#AD002AFF"))
#scale_fill_gradientn(colors=c("black","#50bda1","#FFD31D","#FFD31D","#ED0000FF","#AD002AFF"))
C0 <-  CellCycleScoring(C0 , s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, assay="RNA")
table(C0$condition)
C0.cellcycle <- DimPlot(C0, group.by = "Phase",split.by = "condition")
C0.cellcycle <- C0.cellcycle + ggtitle("Cell Cycle of Cluster 0")
C0.cellcycle 
table(C0$condition, C0$Phase)/c(614,1775,614,1775,614,1775)

Idents(C0) <- "condition"
ID_0 <- c(sample(x = colnames(subset(C0,idents = "LHV")), size = 614),colnames(subset(C0,idents = "Control")))
C0_Even <- subset(C0, cells= ID_0)
table(C0_Even$condition)

table(C0_Even[[]][C0_Even[[]]$condition == "Control",]$S.Score<0.09&C0_Even[[]][C0_Even[[]]$condition == "Control",]$G2M.Score<0.05)
table(C0_Even[[]][C0_Even[[]]$condition == "LHV",]$S.Score<0.09&C0_Even[[]][C0_Even[[]]$condition == "Control",]$G2M.Score<0.05)
418/614
441/614

p <- ggplot(C0_Even[[]], aes(x=S.Score, y=G2M.Score)) + geom_point(size=0.2) 
p <- p + facet_wrap(.~condition, ncol=1) + ylim(c(-0.15,0.4)) + xlim(-0.15,0.4)+ ggtitle("C0_V_A")
p

pp <- ggplot(C0_Even[[]], aes(x=condition, y=G2M.Score)) + geom_violin()
pp <- pp +  stat_compare_means()
pp

pp2 <- ggplot(C0_Even[[]], aes(x=condition, y=S.Score)) + geom_violin()
pp2 <- pp2 +  stat_compare_means()
pp2

pv0 <- pp + pp2 +p
pv0

p0 <- p + stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="grey", fill="white")
p0 <- p +stat_density_2d(aes(fill = ..density..), geom = "raster", contour = F) + ggtitle("C0_V_A") +
  scale_fill_gradientn(colors=c("black","white","white","white","#FFD31D","#FFD31D","#FFD31D","#FFD31D"))
  #scale_fill_gradientn(colors=c("black","white","white","#FFD31D","#FFD31D","#ED0000FF","#ED0000FF","#AD002AFF"))
  #scale_fill_gradientn(colors=c("black","white","#FFD31D","#FFD31D","#ED0000FF","#AD002AFF","#AD002AFF","#AD002AFF"))
#p0 <- p0 + geom_segment(aes(x = 0.09, y = -0.15, xend = 0.09, yend = 0.05), color="blue") +
 # geom_segment(aes(x = -0.14, y = 0.05, xend = 0.09, yend = 0.05), color="blue")
p0
ggsave2(filename = "~/Box/XinZhouFiles/Projects/ZXE7_SingleCell_HLHS/NewAnalysis/cellcycle/C0.dot.pdf",pv0, width = 8, height = 5, units = "in",dpi=300)
ggsave2(filename = "~/Box/XinZhouFiles/Projects/ZXE7_SingleCell_HLHS/NewAnalysis/cellcycle/C0.Seperate_Evennumber_final.pdf",p0, width = 5, height = 8, units = "in",dpi=300)

C1 <-  CellCycleScoring(C1 , s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, assay="RNA")
C1.cellcycle <- DimPlot(C1, group.by = "Phase",split.by = "condition")
C1.cellcycle <- C1.cellcycle + ggtitle("Cell Cycle of Cluster 1")
C1.cellcycle 
table(C1$condition, C1$Phase)/c(243,786,243,786,243,786)

Idents(C1) <- "condition"
ID_1 <- c(sample(x = colnames(subset(C1,idents = "LHV")), size = 243),colnames(subset(C1,idents = "Control")))
C1_Even <- subset(C1, cells= ID_1)
table(C1_Even$condition)

table(C1_Even[[]][C1_Even[[]]$condition == "Control",]$S.Score<0.0&C1_Even[[]][C1_Even[[]]$condition == "Control",]$G2M.Score<0.0)
table(C1_Even[[]][C1_Even[[]]$condition == "LHV",]$S.Score<0.0&C1_Even[[]][C1_Even[[]]$condition == "Control",]$G2M.Score<0.0)
217/243
220/243

p <- ggplot(C1_Even[[]], aes(x=S.Score, y=G2M.Score)) + geom_point(size=0.2)
p <- p + facet_wrap(.~condition, ncol=1) + ylim(c(-0.15,0.5)) + xlim(-0.13,0.3) + ggtitle("C1_Artery")
p

pp <- ggplot(C1_Even[[]], aes(x=condition, y=G2M.Score)) + geom_violin()
pp <- pp +  stat_compare_means()
pp

pp2 <- ggplot(C1_Even[[]], aes(x=condition, y=S.Score)) + geom_violin()
pp2 <- pp2 +  stat_compare_means()
pp2

pv1 <- pp + pp2 + p
pv1

p1 <- p + stat_density_2d(aes(fill = ..density..), geom = "raster", contour = F) + ggtitle("C1_Artery") +
  scale_fill_gradientn(colors=c("black","white","white","white","#FFD31D","#FFD31D","#FFD31D","#FFD31D"))
  #scale_fill_gradientn(colors=c("black","white","white","#FFD31D","#FFD31D","#ED0000FF","#ED0000FF","#AD002AFF"))
  #scale_fill_gradientn(colors=c("black","white","#FFD31D","#ED0000FF","#AD002AFF"))
  #+scale_fill_gradientn(colors=c("black","white","#FFD31D","#FFD31D","#ED0000FF","#AD002AFF","#AD002AFF","#AD002AFF"))
#p1 <- p1 + geom_segment(aes(x = 0.0, y = -0.15, xend = 0.0, yend = 0.0), color="blue") +
  #geom_segment(aes(x = -0.1, y = 0.0, xend = 0.0, yend = 0.0), color="blue")
p1

ggsave2(filename = "~/Box/XinZhouFiles/Projects/ZXE7_SingleCell_HLHS/NewAnalysis/cellcycle/C1.dot.pdf",pv1, width = 8, height = 5, units = "in",dpi=300)
ggsave2(filename = "~/Box/XinZhouFiles/Projects/ZXE7_SingleCell_HLHS/NewAnalysis/cellcycle/C1.Seperate_Evennumber_final.pdf",p1, width = 5, height = 8, units = "in",dpi=300)

C2 <-  CellCycleScoring(C2 , s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, assay="RNA")
C2.cellcycle <- DimPlot(C2, group.by = "Phase",split.by = "condition")
C2.cellcycle <- C2.cellcycle + ggtitle("Cell Cycle of Cluster 2")
C2.cellcycle
table(C2$condition, C2$Phase)/c(404,522,404,522,404,522)

Idents(C2) <- "condition"
ID_2 <- c(sample(x = colnames(subset(C2,idents = "LHV")), size = 404),colnames(subset(C2,idents = "Control")))
C2_Even <- subset(C2, cells= ID_2)
table(C2_Even$condition)

table(C2_Even[[]][C2_Even[[]]$condition == "Control",]$S.Score< (-0.1)&C2_Even[[]][C2_Even[[]]$condition == "Control",]$G2M.Score<(-0.3))
table(C2_Even[[]][C2_Even[[]]$condition == "LHV",]$S.Score < (-0.1)&C2_Even[[]][C2_Even[[]]$condition == "Control",]$G2M.Score<(-0.3))
29/404
43/243

p <- ggplot(C2_Even[[]], aes(x=S.Score, y=G2M.Score)) + geom_point(size=0.2)
p <- p + facet_wrap(.~condition, ncol=1) + ylim(c(-0.8,1.2)) + xlim(-0.6,0.6)
p

pp <- ggplot(C2_Even[[]], aes(x=condition, y=G2M.Score)) + geom_violin()
pp <- pp +  stat_compare_means()
pp

pp2 <- ggplot(C2_Even[[]], aes(x=condition, y=S.Score)) + geom_violin()
pp2 <- pp2 +  stat_compare_means()
pp2

pv2 <- pp + pp2 + p
pv2

p2 <- p + stat_density_2d(aes(fill = ..density..), geom = "raster", contour = F) + ggtitle("C2_Vein") +
  scale_fill_gradientn(colors=c("black","white","white","white","#FFD31D","#FFD31D","#FFD31D","#FFD31D"))
  #scale_fill_gradientn(colors=c("black","white","white","#FFD31D","#FFD31D","#ED0000FF","#ED0000FF","#AD002AFF"))
#+scale_fill_gradientn(colors=c("black","white","#FFD31D","#FFD31D","#ED0000FF","#AD002AFF","#AD002AFF","#AD002AFF"))
#p2 <- p2 + geom_segment(aes(x = -0.1, y = -0.7, xend = -0.1, yend = -0.3), color="blue") +
 #geom_segment(aes(x = -0.4, y = -0.3, xend = -0.1, yend = -0.3), color="blue")
p2
ggsave2(filename = "~/Box/XinZhouFiles/Projects/ZXE7_SingleCell_HLHS/NewAnalysis/cellcycle/C2.dot.pdf",pv2, width = 8, height = 5, units = "in",dpi=300)
ggsave2(filename = "~/Box/XinZhouFiles/Projects/ZXE7_SingleCell_HLHS/NewAnalysis/cellcycle/C2.Seperate_Evennumber_final.pdf",p2, width = 5, height = 8, units = "in",dpi=300)

C3 <-  CellCycleScoring(C3 , s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, assay="RNA")
C3.cellcycle <- DimPlot(C3, group.by = "Phase",split.by = "condition")
C3.cellcycle <- C3.cellcycle + ggtitle("Cell Cycle of Cluster 3")
C3.cellcycle
table(C3$condition, C3$Phase)/c(225,488,225,488,225,488)

Idents(C3) <- "condition"
ID_3 <- c(sample(x = colnames(subset(C3,idents = "LHV")), size = 225),colnames(subset(C3,idents = "Control")))
C3_Even <- subset(C3, cells= ID_3)
table(C3_Even$condition)

table(C3_Even[[]][C3_Even[[]]$condition == "Control",]$S.Score<0.07&C3_Even[[]][C3_Even[[]]$condition == "Control",]$G2M.Score<0.06)
table(C3_Even[[]][C3_Even[[]]$condition == "LHV",]$S.Score<0.07&C3_Even[[]][C3_Even[[]]$condition == "Control",]$G2M.Score<0.06)
163/225
175/225

p <- ggplot(C3_Even[[]], aes(x=S.Score, y=G2M.Score)) + geom_point(size=0.2)
p <- p + facet_wrap(.~condition, ncol=1) + ylim(c(-0.15,0.3)) + xlim(-0.15,0.4)+ ggtitle("C3_V_A")
#p <- p + geom_line(aes(y = ..density..), colour=1, stat = 'density', size = 2, alpha = .6)
p

pp <- ggplot(C3_Even[[]], aes(x=condition, y=G2M.Score)) + geom_violin()
pp <- pp +  stat_compare_means()
pp

pp2 <- ggplot(C3_Even[[]], aes(x=condition, y=S.Score)) + geom_violin()
pp2 <- pp2 +  stat_compare_means()
pp2

pv3 <- pp + pp2 +p
pv3

p3 <- p +stat_density_2d(aes(fill = ..density..), geom = "raster", contour = F) + ggtitle("C3_V_A") +
  scale_fill_gradientn(colors=c("black","white","white","white","#FFD31D","#FFD31D","#FFD31D","#FFD31D"))
 # scale_fill_gradientn(colors=c("black","white","white","#FFD31D","#FFD31D","#ED0000FF","#ED0000FF","#AD002AFF"))
#+scale_fill_gradientn(colors=c("black","white","#FFD31D","#FFD31D","#ED0000FF","#AD002AFF","#AD002AFF","#AD002AFF"))
#p3 <- p3 + geom_segment(aes(x = 0.07, y = -0.11, xend = 0.07, yend = 0.06), color="blue") +
 # geom_segment(aes(x = -0.14, y = 0.06, xend = 0.07, yend = 0.06), color="blue")
p3
p0+p3+p1+p2

ggsave2(filename = "~/Box/XinZhouFiles/Projects/ZXE7_SingleCell_HLHS/NewAnalysis/cellcycle/C3.dot.pdf",pv3, width = 8, height = 5, units = "in",dpi=300)
ggsave2(filename = "~/Box/XinZhouFiles/Projects/ZXE7_SingleCell_HLHS/NewAnalysis/cellcycle/C3.Seperate_Evennumber_Final.pdf",p3, width = 5, height = 8, units = "in",dpi=300)

lhv.Endothelium<- CellCycleScoring(lhv.Endothelium, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, assay="RNA")
Idents(lhv.Endothelium) <- "celltype2"

pcellcycle <- DimPlot(lhv.Endothelium, group.by = "Phase", split.by = "condition")
pcellcycle

p <- ggplot(subset(lhv.Endothelium[[]], celltype2=="Endothelium"), aes(x=S.Score, y=G2M.Score)) + geom_point(size=0.2)
p

########################
#DEG
########################
DefaultAssay(C0)<- "RNA"
DimPlot(C0)
C0 <- NormalizeData(C0,normalization.method = "LogNormalize", scale.factor = 10000)
Idents(C0) <- "condition"
C0.DEG <- FindMarkers(C0, ident.1 = "Control", ident.2 = "LHV", min.pct = 0.2, test.use="wilcox", assay="RNA",logfc.threshold=log(1.2))
C0.DEG$genes <- row.names(C0.DEG)
write.csv(file = "../NewAnalysis/DEG/fetal.CO.DEG.csv",C0.DEG )

p <- ggplot(C0.DEG, aes(x=avg_logFC, y=-log(p_val_adj))) + geom_point() + theme_set(theme_cowplot(12))
p <- p + xlab("Average Log Fold Change") + ylab("Negative Log Adjusted P Value") + ggtitle("C0 Deseq Result on Control/LHV DEGs")
p <- p + geom_label_repel(data=subset(C0.DEG, C0.DEG$genes %in% c("GAPDG","B2M","MICA","FN1","APLN", "MGLL", "NOTCH1","DLL4", "JAK1","ACTB")), aes(label=genes, color="red"))
p <- p + geom_point(data=subset(C0.DEG, C0.DEG$genes %in% c("GAPDG","B2M","MICA","FN1","APLN", "MGLL", "NOTCH1", "DLL4","JAK1","ACTB")), color="red")
p

DefaultAssay(C1) <- "RNA"
DimPlot(C1)
Idents(C1) <- "condition"
C1.DEG <- FindMarkers(C1, ident.1 = "Control", ident.2 = "LHV", min.pct = 0.2, test.use="wilcox", assay="RNA",logfc.threshold=log(1.2))
C1.DEG$genes <- row.names(C1.DEG)
write.csv(file ="../NewAnalysis/DEG/fetal.C1.DEG.csv",C1.DEG)

p1 <- ggplot(C1.DEG, aes(x=avg_logFC, y=-log(p_val_adj))) + geom_point() + theme_set(theme_cowplot(12))
p1 <- p1 + xlab("Average Log Fold Change") + ylab("Negative Log Adjusted P Value") + ggtitle("C1 Deseq Result on Control/LHV DEGs")
p1 <- p1 + geom_label_repel(data=subset(C1.DEG, C1.DEG$genes %in% c("GAPDG","B2M","MICA","FN1","APLN", "MGLL", "NOTCH1","DLL4", "JAK1","ACTB")), aes(label=genes, color="red"))
p1 <- p1 + geom_point(data=subset(C1.DEG, C1.DEG$genes %in% c("GAPDG","B2M","MICA","FN1","APLN", "MGLL", "NOTCH1", "DLL4","JAK1","ACTB")), color="red")
p1

DefaultAssay(C2) <- "RNA"
Idents(C2) <- "condition"
C2.DEG <- FindMarkers(C2, ident.1 = "Control", ident.2 = "LHV", min.pct = 0.2, test.use="wilcox", assay="RNA",logfc.threshold=log(1.2))
C2.DEG$genes <- row.names(C2.DEG)
write.csv(file ="../NewAnalysis/DEG/fetal.C2.DEG.csv",C2.DEG)

DefaultAssay(C3) <- "RNA"
Idents(C3) <- "condition"
C3.DEG <- FindMarkers(C3, ident.1 = "Control", ident.2 = "LHV", min.pct = 0.2, test.use="wilcox", assay="RNA",logfc.threshold=log(1.2))
C3.DEG$genes <- row.names(C3.DEG)
write.csv(file ="../NewAnalysis/DEG/fetal.C3.DEG.csv",C3.DEG)

DefaultAssay(C4) <- "RNA"
Idents(C4) <- "condition"
C4.DEG <- FindMarkers(C4, ident.1 = "Control", ident.2 = "LHV", min.pct = 0.2, test.use="wilcox", assay="RNA",logfc.threshold=log(1.2))
C4.DEG$genes <- row.names(C4.DEG)
write.csv(file ="../NewAnalysis/DEG/fetal.C4.DEG.csv",C4.DEG)

DefaultAssay(C5) <- "RNA"
Idents(C5) <- "condition"
C5.DEG <- FindMarkers(C5, ident.1 = "Control", ident.2 = "LHV", min.pct = 0.2, test.use="wilcox", assay="RNA",logfc.threshold=log(1.2))
C5.DEG$genes <- row.names(C5.DEG)
write.csv(file = "../NewAnalysis/DEG/fetal.C5.DEG.csv",C5.DEG)

#########################################################
uplist <- filter(C1.DEG, avg_logFC > log(1.2) & p_val_adj < 0.05)
downlist <- filter(C1.DEG, avg_logFC < -log(1.2) & p_val_adj < 0.05)

dim(uplist)
dim(downlist)

upgeneList <- uplist[, 2]
names(upgeneList) <- as.character(uplist$gene)
upgeneList <- sort(upgeneList, decreasing = TRUE)
head(upgeneList)

downlistList <- downlist[, 2]
names(downlistList) <- as.character(downlist$genes)
downlistList <- sort(downlistList, decreasing = TRUE)
head(downlistList)

length(upgeneList)
length(downlistList)

gene.df <- bitr(names(upgeneList), fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
ego_up <- enrichGO(gene         = gene.df$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENTREZID',
                   ont           = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable = T)
GO_Up_result <- ego_up@result

write.csv(file = "../NewAnalysis/DEG/fetal.GO.C5.Up.csv",GO_Up_result)

gene.df_2 <- bitr(names(downlistList), fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
ego_down <- enrichGO(gene         = gene.df_2$ENTREZID,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'ENTREZID',
                     ont           = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05, 
                     readable = T)

GO_Down_result <- ego_down@result

#write.csv(file = "../NewAnalysis/DEG/fetal.GO.C1.down.csv",GO_Down_result)

edox <- setReadable(ego_down, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=downlistList)
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=downlistList)
p3 <- cnetplot(edox, foldChange=downlistList, circular = TRUE, colorEdge = TRUE)

subset(edox, gene="9246")

cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))

#########################################################
#total DEG
#########################################################
DefaultAssay(lhv.Endothelium)<- "RNA"
DimPlot(lhv.Endothelium)
lhv.Endothelium <- NormalizeData(lhv.Endothelium,normalization.method = "LogNormalize", scale.factor = 10000)
Idents(lhv.Endothelium) <- "condition"
lhv.Endothelium.DEG <- FindMarkers(lhv.Endothelium, ident.1 = "Control", ident.2 = "LHV", min.pct = 0.2, test.use="wilcox", assay="RNA",logfc.threshold=log(1.2))
lhv.Endothelium.DEG$genes <- row.names(lhv.Endothelium.DEG)
write.csv(file = "../NewAnalysis/DEG/fetal.Endothelium.total.DEG.csv",lhv.Endothelium.DEG )

p <- ggplot(lhv.Endothelium.DEG, aes(x=avg_logFC, y=-log(p_val_adj))) + geom_point() + theme_set(theme_cowplot(12))
p <- p + xlab("Average Log Fold Change") + ylab("Negative Log Adjusted P Value") + ggtitle("C0 Deseq Result on Control/LHV DEGs")
p <- p + geom_label_repel(data=subset(lhv.Endothelium.DEG, lhv.Endothelium.DEG$genes %in% c("MICA","FN1","APLN", "MGLL", "NOTCH1","DLL4", "JAK1","ACTB")), aes(label=genes, color="red"))
p <- p + geom_point(data=subset(lhv.Endothelium.DEG, lhv.Endothelium.DEG$genes %in% c("MICA","FN1","APLN", "MGLL", "NOTCH1", "DLL4","JAK1","ACTB")), color="red")
p

ggsave2(filename = "../NewAnalysis/DEG/fetal.total.EDG.pdf", p, scale = 0.7)

uplist <- filter(lhv.Endothelium.DEG, avg_logFC > log(1.2) & p_val_adj < 0.05)
downlist <- filter(lhv.Endothelium.DEG, avg_logFC < -log(1.2) & p_val_adj < 0.05)

dim(uplist)
dim(downlist)

upgeneList <- uplist[, 2]
names(upgeneList) <- as.character(uplist$gene)
upgeneList <- sort(upgeneList, decreasing = TRUE)
head(upgeneList)

downlistList <- downlist[, 2]
names(downlistList) <- as.character(downlist$genes)
downlistList <- sort(downlistList, decreasing = TRUE)
head(downlistList)

length(upgeneList)
length(downlistList)

gene.df <- bitr(names(upgeneList), fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
ego_up <- enrichGO(gene         = gene.df$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENTREZID',
                   ont           = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable = T)
GO_Up_result <- ego_up@result

write.csv(file = "../NewAnalysis/DEG/fetal.GO.total.Up.csv",GO_Up_result)

gene.df_2 <- bitr(names(downlistList), fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
ego_down <- enrichGO(gene         = gene.df_2$ENTREZID,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'ENTREZID',
                     ont           = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05, 
                     readable = T)
GO_Down_result <- ego_down@result
write.csv(file = "../NewAnalysis/DEG/fetal.GO.total.down.csv",GO_Down_result)

commonGO <- read.csv("../NewAnalysis/DEG/commonGOUP.csv", header = F)
table(GO_Up_result$ID %in% commonGO$V1)
commonUP <- GO_Up_result[GO_Up_result $ID %in% commonGO$V1,]
write.csv(file = "../NewAnalysis/DEG/commonUP.csv",commonUP)

commonGOD <- read.csv("../NewAnalysis/DEG/CommonGOdown.csv", header = F)
table(GO_Down_result$ID %in% commonGOD$V1)
commonDOWN <- GO_Down_result[GO_Down_result $ID %in% commonGOD$V1,]
write.csv(file = "../NewAnalysis/DEG/commonDOWN.csv",commonDOWN)

Idents(lhv.Endothelium)
lhv.Endothelium$integrated_snn_res.0.2 <- factor((lhv.Endothelium$integrated_snn_res.0.2), levels = c(2,0,3,1,4,5))
DefaultAssay(lhv.Endothelium) <- "SCT"
lhv.Endothelium_03 <- subset(lhv.Endothelium, idents=c(0,1,2,3))
Idents(lhv.Endothelium_03)
vl.markers <- c("NR2F2", "GJA5", "AQP1","HEY1","DLL4")
for(gene in vl.markers){
  print(gene)
  p=VlnPlot(object = lhv.Endothelium_03, features = c(gene),pt.size = 0, split.by = "condition",group.by="integrated_snn_res.0.2",split.plot = T) + scale_fill_d3() +
    theme(axis.title.x = element_blank(),axis.text.y  = element_text( vjust=0.5, size=34,face="bold"),
          axis.text.x  = element_text( vjust=0.5, size=34,face="bold"), plot.title = element_text(lineheight=.8, size=40, face="bold.italic"))
  #+VlnPlot(object = all, features.plot = c(gene),group.by="orig.ident")
  #p2=VlnPlot(object = all, features.plot = c(gene),group.by="orig.ident")
  ggsave(p,filename = paste0("../NewAnalysis/Final Violin/lvh0_3.split.",gene,".pdf"),height=6,width=12)
}

#save(list=ls(.GlobalEnv), file = "~/Desktop/lhls.RData")

save(lhv.Endothelium, file = "~/Desktop/lhv.Endothelium.RData")
