# Combine HLHS and normal iPSC-ECs
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggrepel)
library(cowplot)
library(DESeq2)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(robustSingleCell)

setwd("~/Box/XinZhouFiles/Projects/ZXE7_SingleCell_HLHS/DATA/")
load("./iPSC-EC/all.Robj")

#all <- FindClusters(object = all, reduction.type = "pca", dims.use = 1:20, resolution = 0.8,force.recalc=T, print.output = 0, save.SNN = TRUE)
#all <- RunTSNE(object = all, dims.use = 1:20, do.fast = TRUE,check_duplicates=F)

lhls.all <- UpdateSeuratObject(all)

hkgene <- read.table("~/Box/XinZhouFiles/Projects/ZXP1_PAH/Single Cell Related/House_keeping_Gene.txt")
hkgenes <- as.vector(hkgene$V1)
hkgenes.found <- which(toupper(rownames(lhls.all@assays$RNA)) %in% hkgenes)
n.expressed.hkgenes <- Matrix::colSums(lhls.all@assays$RNA[hkgenes.found, ] > 0)
lhls.all <- AddMetaData(object = lhls.all, metadata = n.expressed.hkgenes, col.name = "n.exp.hkgenes")

vplot <- VlnPlot(object = lhls.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "n.exp.hkgenes"), ncol=2, combine = T, pt.size=0.05, group.by = "orig.ident")
vplot
#ggsave2(filename = "./R_work_place/Sup/S1.QC.pdf", vplot, scale=0.77)

vplot2 <- VlnPlot(object = lhls.all, features = c("nGene", "nUMI", "percent.mito", "n.exp.hkgenes"), ncol=2, combine = T, pt.size=0.05, group.by = "orig.ident")
vplot2

vplot3 <- VlnPlot(object = lhls.all, features = c("ACTB", "PPIA", "RPL26", "UBB"), ncol=2, combine = T, pt.size=0.05, group.by = "orig.ident")
vplot3
#ggsave2(filename = "./R_work_place/Sup/S2.HKG.pdf", vplot3, scale=0.77)

############ Combime two dataset
lhls.list <- SplitObject(lhls.all, split.by = "orig.ident")
for (i in 1:length(lhls.list)) {
  lhls.list[[i]] <- SCTransform(lhls.list[[i]], verbose = T,vars.to.regress="nCount_RNA")
}

#for (i in 1:length(lhls.list)) {
#  lhls.list[[i]] <- FindVariableFeatures(lhls.list[[i]],  verbose = T)
#}

#pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
lhls.features <- SelectIntegrationFeatures(object.list = lhls.list, nfeatures = 3000)
lhls.list <- PrepSCTIntegration(object.list = lhls.list, anchor.features = lhls.features, verbose = T)

lhls.anchors <- FindIntegrationAnchors(object.list = lhls.list, normalization.method = "SCT", anchor.features = lhls.features, verbose = T)
lhls.new <- IntegrateData(anchorset = lhls.anchors, normalization.method = "SCT", verbose = T)

lhls.new <- RunPCA(lhls.new, dim=1:15)
lhls.new <- RunUMAP(lhls.new, dim=1:30)

lhls.new$condition <- lhls.new$orig.ident
lhls.new$condition[lhls.new$orig.ident=="BA060"] <- "Control"
lhls.new$condition[lhls.new$orig.ident=="M146"] <-  "HLHS"
table(lhls.new$condition)

#JackStrawPlot
#lhls.new <- JackStraw(lhls.new, num.replicate = 100)
#lhls.new <- ScoreJackStraw(lhls.new, dims = 1:20)
#JackStrawPlot(lhls.new, dims = 1:20)
#ElbowPlot(lhls.new)

lhls.new <- FindNeighbors(lhls.new, dims = 1:10)
lhls.new <- FindClusters(lhls.new, resolution = 0.2)
Idents(lhls.new)


lhls.new <- RenameIdents(lhls.new,  "0"="Endothelium_0", "1"="Endocardium_1", "2"="Endocardium_2", "3"="Endothelium_3")
lhls.new$celltype <- Idents(lhls.new)

#figure1A
plots1 <- DimPlot(lhls.new,group.by ="celltype", combine = T) + scale_color_aaas()
plots1
ggsave2(filename = "./R_work_place/F1.UMAP.pdf", plots1, scale=0.77)
#figure1B
plots2 <- DimPlot(lhls.new, group.by = "condition", combine = T) + scale_color_aaas()
plots2
ggsave2(filename = "./R_work_place/F2.Condition.pdf", plots2, scale=0.77)
#figureS1C
plots3 <- DimPlot(lhls.new, group.by = "Phase", combine = T) + scale_color_d3()
plots3
ggsave2(filename = "./R_work_place/Sup/S3.phase.pdf", plots3, scale=0.77)
#compare with previous plot
plots4 <- DimPlot(lhls.new, group.by = "res.0.8", combine = T)
plots4
ggsave2(filename = "./R_work_place/Sup/S4.withpreviousUMAP.pdf", plots4, scale=0.77)

count <- group_by(lhls.new@meta.data, condition, integrated_snn_res.0.2) %>%
  dplyr::summarise(count = n()) %>% 
  mutate(freq = count / sum(count))
colnames(count)[2] <- "group"

# Barplot
bp<- ggplot(count, aes(x="", y=freq, fill=group)) + theme_set(theme_cowplot(12)) +
  geom_bar(width = 1, stat = "identity") + scale_fill_aaas() + theme(axis.text.y=element_blank())
bp
pie <- bp + coord_polar("y", start=0) + facet_wrap(.~condition)
pie
ggsave2(filename = "./R_work_place/F3.Pie.pdf", pie, scale=0.77)

DefaultAssay(lhls.new)

#look for markers that distinguishi cluster
Idents(lhls.new) <- "integrated_snn_res.0.2"
Idents(lhls.new)
DefaultAssay(lhls.new) <- "integrated"
clustermarker0_3 <- FindMarkers(lhls.new, ident.1 = 0, ident.2=3, verbose = TRUE, test.use = "LR", logfc.threshold = 0.05, min.pct = 0.5)
clustermarker0_3
clustermarker1_2 <- FindMarkers(lhls.new, ident.1 = 1, ident.2=2, verbose = TRUE, test.use = "LR", logfc.threshold = 0.05, min.pct = 0.5)
clustermarker1_2 
clustermarker12_03_prevalence <- FindMarkers(lhls.new, ident.1 = c(1,2), ident.2=c(0,3), verbose = TRUE, test.use = "LR", logfc.threshold = 0.05, min.pct = 0.00, min.diff.pct=0.4)
clustermarker12_03_prevalence

write.csv(clustermarker0_3, file = "./R_work_place/DEG/NewClustermarker0_3.csv")
write.csv(clustermarker1_2, file = "./R_work_place/DEG/Newclustermarker1_2.csv")
write.csv(clustermarker12_03_prevalence, file = "./R_work_place/DEG/Newclustermarker12_03_prevalence.csv")

######################################subset####################################

#########################################################################################################
lhls.new$celltype2 <- as.character(lhls.new$integrated_snn_res.0.2)
lhls.new$celltype2[lhls.new$celltype2 %in% c(0,3)] <- "Endothelium"
lhls.new$celltype2[lhls.new$celltype2 %in% c(1,2)] <- "Endocardium"

Idents(lhls.new) <- "celltype2"
lhls.Endothelium <- subset(lhls.new, idents="Endothelium")

lhls.Endothelium <- RunPCA(lhls.Endothelium, dim=1:15)
lhls.Endothelium <- RunUMAP(lhls.Endothelium, dim=1:30)
lhls.Endothelium <- FindNeighbors(lhls.Endothelium, dims = 1:10)
lhls.Endothelium <- FindClusters(lhls.Endothelium, resolution = 0.2)

table(lhls.Endothelium$condition)

Idents(lhls.Endothelium)
plots5 <- DimPlot(lhls.Endothelium,group.by ="integrated_snn_res.0.2", combine = T) #+ scale_color_aaas()
plots5

ggsave2("../NewAnalysis/UMAP/iPSC.after.subset.part.UMAP.pdf", plots5, scale = 0.77)

plots6 <- DimPlot(lhls.Endothelium,group.by ="condition", combine = T) + scale_color_d3()
plots6

ggsave2("../NewAnalysis/UMAP/iPSC.after.subset.condition.UMAP.pdf", plots6, scale = 0.77)

Idents(lhls.Endothelium)
DefaultAssay(lhls.Endothelium) <- "SCT"
vl.markers <- c("NR2F2", "GJA5", "AQP1", "HEY1", "DLL4")
for(gene in vl.markers){
  print(gene)
  p=VlnPlot(object = lhls.Endothelium, features = c(gene),group.by="integrated_snn_res.0.2",split.by = "condition",pt.size = 0, split.plot = T) + scale_fill_d3() + 
    geom_boxplot(width=.05,outlier.alpha=0, fill="black",colour="gray",coef=0.1) + 
    theme(axis.title.x = element_blank(),axis.text.y  = element_text( vjust=0.5, size=34,face="bold"),
          axis.text.x  = element_text( vjust=0.5, size=34,face="bold"), plot.title = element_text(lineheight=.8, size=40, face="bold.italic"))
  #+VlnPlot(object = all, features.plot = c(gene),group.by="orig.ident")
  #p2=VlnPlot(object = all, features.plot = c(gene),group.by="orig.ident")
  ggsave(p,filename = paste0("../NewAnalysis/Final Violin/ipsc.split.",gene,".pdf"),height=6,width=12)
}


#lhls.Endothelium.Marker <- FindAllMarkers(lhls.Endothelium)
#write.csv(file = "../NewAnalysis/UMAP/iPSC.all.markers.csv", lhls.Endothelium.Marker)
############################################################################################################################################

############################################################################################################################################
#plot differential gene by cluster
lowof0 <- VlnPlot(subset(lhls.new,idents=c(0,3)), features = row.names(subset(clustermarker0_3, avg_logFC<0))[1:9], ncol=3, pt.size=0, cols=c("#3B4992FF","#631879FF"),combine = T)
lowof0 <- lapply(X=lowof0, FUN = function(x) x + scale_fill_d3())
lowof0 <- CombinePlots(lowof0, legend="right")
lowof0
ggsave2(filename = "./R_work_place/Sup/S9.0_3DEG_low.pdf", lowof0, scale=0.77)

highof0 <- VlnPlot(subset(lhls.new,idents=c(0,3)), features = row.names(subset(clustermarker0_3, avg_logFC>0))[1:9], ncol=3,pt.size=0,cols=c("#3B4992FF","#631879FF"),combine = T)
highof0 <- lapply(X=highof0, FUN = function(x) x + scale_fill_d3())
highof0 <- CombinePlots(highof0, legend="right")
highof0
ggsave2(filename = "./R_work_place/Sup/S10.0_3DEG_high.pdf", highof0, scale=0.77)



p=VlnPlot(object = lhls.new, features = "TGFB2", pt.size=1.2, assay = "SCT", group.by = "celltype2", split.by = "condition")
p

DefaultAssay(lhls.new) <-"SCT"
DefaultAssay(lhls.new) <-"RNA"
DefaultAssay(lhls.new) <- "integrated"

#remove sample that are LYVE1 + (Lymphatic EC) or NPR3+ (Endothelium)
p <- FeatureScatter(lhls.new, feature1 = "ROR1",feature2 = "TEAD1", group.by = "integrated_snn_res.0.2") + scale_color_aaas()
p <- p + stat_density2d(aes(fill = ..level..), alpha=0.5, geom="polygon")
#p <- p + facet_wrap(.~colors, ncol=2)
p

lhls.new.allmarkers <- FindAllMarkers(lhls.new, only.pos = TRUE, min.pct = 0.25, test.use = "LR", logfc.threshold = 0.25)

#Find certain markers
Idents(lhls.new) <- "celltype2"
lhls.Endothelium <- subset(lhls.new, idents="Endothelium")
DefaultAssay(lhls.Endothelium) <-"SCT"
marksgrep <- grep(pattern = "^TGFB", x = rownames(x = lhls.new@assays$RNA), value = TRUE)
marksgrep
vplot3 <- VlnPlot(object = lhls.Endothelium, features = cellcycle, ncol=2, combine = T, pt.size=0.05, group.by = "condition", cols =c("#1F77B4FF","#FF7F0EFF") ) 
vplot3

subset(lhls.new, CDH5>0)
lhls.new

#http://www.oncternal.com/Pipeline/Cirmtuzumab
WNT5A_P <- c("ROR1", "ROR2","AKT2", "TEAD1", "RAC1", "RHOA")
Notch_P <- c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4")
Notch_R_P <- c("DLL4", "JAG1", "HEY1", "TNFRSF1B")

#Liange marker
lin <- c("CDH5","PECAM1","CD34","APLN")
cellcycle <- c("CCND1", "MPP1","MPP6", "TP53", "MKI67", "PCNA", "PCNA", "CCNA2", "CCNB1", "CCND2")

ALK1 <- c("ACVRL1", "TMEM100")
P53 <- c("TP53BP2", "PTEN", "MAPK14", "HRAS", "SOD2", "PRKCD", "SMAD2","SOD2")



#Subset Endothelium
DefaultAssay(lhls.Endothelium) <- "SCT"
Idents(lhls.Endothelium) <- "condition"
#DefaultAssay(lhls.Endothelium)
ScaleData(lhls.Endothelium, vars.to.regress = "nUMI")
#Endothelium.markers.SCT.deseq <- FindMarkers(lhls.Endothelium, ident.1 = "Control", ident.2 = "HLHS", min.pct = 0.2, test.use="DESeq2", assay="SCT", logfc.threshold=0.3)
Endothelium.markers.SCT.mast <- FindMarkers(lhls.Endothelium, ident.1 = "Control", ident.2 = "HLHS", min.pct = 0.2, test.use="MAST", assay="SCT",logfc.threshold=0)

DefaultAssay(lhls.Endothelium) <- "RNA"
lhls.Endothelium <- NormalizeData(lhls.Endothelium, normalization.method = "LogNormalize", scale.factor = 10000)
Idents(lhls.Endothelium) <- "condition"
DefaultAssay(lhls.Endothelium)
Idents(lhls.Endothelium)
#ScaleData(lhls.new, vars.to.regress = "nUMI")
Endothelium.markers.RNA.deseq.ori <- FindMarkers(lhls.Endothelium, ident.1 = "Control", ident.2 = "HLHS", min.pct = 0.2, test.use="wilcox", assay="RNA",logfc.threshold=log(1))
#Endothelium.markers.RNA.mast <- FindMarkers(lhls.Endothelium, ident.1 = "Control", ident.2 = "HLHS", min.pct = 0.2, test.use="MAST", assay="RNA", logfc.threshold=0.3)

#EM.SCT.deseq <- Endothelium.markers.SCT.deseq[abs(Endothelium.markers.SCT.deseq$avg_logFC) > 0.3,]
#EM.SCT.mast <- Endothelium.markers.SCT.mast[abs(Endothelium.markers.SCT.mast$avg_logFC) > 0.25,]

#EM.RNA.deseq <- Endothelium.markers.RNA.deseq[abs(Endothelium.markers.RNA.deseq$avg_logFC) > 0.3,]
cluster3.markers <- Endothelium.markers.RNA.deseq.ori
cluster3.markers$genes <- row.names(cluster3.markers)
write.csv(file = "../NewAnalysis/DEG/ipsc.total.EDG.csv",cluster3.markers)

#oldlist <- read.table(file = "./R_work_place/DEG/c1345.M146vsBA060.txt")
#setdiff(row.names(oldlist), row.names(EM.RNA.deseq))

#Check specific genes
cluster3.markers[cluster3.markers$genes == "JAK1",]
cluster3.markers[cluster3.markers$genes == "NOTCH1",]
cluster3.markers[cluster3.markers$genes == "GAPDH",]
cluster3.markers[cluster3.markers$genes == "ACTB",]
cluster3.markers[cluster3.markers$genes == "MGLL",]
cluster3.markers[cluster3.markers$genes == "APLN",]
cluster3.markers[cluster3.markers$genes == "MICA",]
cluster3.markers[cluster3.markers$genes == "TGFB2",]

p <- ggplot(cluster3.markers, aes(x=avg_logFC, y=-log(p_val_adj))) + geom_point() + theme_set(theme_cowplot(12))
p <- p + xlab("Average Log Fold Change") + ylab("Negative Log Adjusted P Value") + ggtitle("Deseq Result on Control/HLHS DEGs")
p <- p + geom_point(data=subset(cluster3.markers, avg_logFC < -log(1.5) | avg_logFC >log(1.5) ), aes(color="red"))
p <- p + geom_label_repel(data=subset(cluster3.markers, cluster3.markers$genes %in% c("MICA","FN1","APLN", "MGLL", "NOTCH1","DLL4", "JAK1","ACTB")), aes(label=genes, color="red"))
p <- p + geom_point(data=subset(cluster3.markers, cluster3.markers$genes %in% c("MICA","FN1","APLN", "MGLL", "NOTCH1", "DLL4","JAK1","ACTB")), color="red")
#p <- p + geom_label_repel(data=subset(cluster3.markers, avg_logFC < -log(1.5) | avg_logFC >0.6), aes(label=genes))
p
ggsave2(filename = "../NewAnalysis/DEG/ipsc.total.EDG.pdf", p, scale = 0.7)


#CD31="PECAM1"; CD144="CDH5"; CD34="CD34"
DefaultAssay(lhls.new) <- "SCT"
p=FeaturePlot(object = lhls.new, features = c("NPR3","APLN"),cols=c("lightgrey","#3C5488B2", "#931F21"),blend = T,
              blend.threshold=c(0.001),pt.size=0.4, sort.cell = T, combine = F)

p[[1]]
p[[3]]
p[[4]]

layout <- c(
  area(t = 1, l = 1, b = 5, r = 5),
  area(t = 1, l = 5, b = 1.5, r = 5.5)
)
pcelltype <- p[[3]] + p[[4]] + plot_layout(design = layout)
#Figure 4
ggsave2(filename = "./R_work_place/F4.NPR3_APLN.pdf",pcelltype, scale = 0.7)

p2=RidgePlot(object = lhls.new, features = "CCND1",pt.size=0, group.by = "integrated_snn_res.0.2")
p2
p / p2

########################################################
Endothelium.markers.RNA.deseq <- Endothelium.markers.RNA.deseq.ori
Endothelium.markers.RNA.deseq$genes <- row.names(Endothelium.markers.RNA.deseq)
Endothelium.markers.RNA.deseq <- filter(Endothelium.markers.RNA.deseq, -log(p_val_adj) > 30)

#determine how to shift X axis
#df <-subset(Endothelium.markers.RNA.deseq, -log(p_val_adj) <150)
#df.lo <- loess(-log(df$p_val_adj)~ df$avg_logFC)
#plot(df.lo$x,df.lo$fitted)
#plot(df.lo$x,df.lo$y)
#min(df.lo$fitted)
#df.lo$x[df.lo$fitted<1.18639,]
#df.lo$fitted[4808]
#df.lo$fitted[5271]

Endothelium.markers.RNA.deseq$avg_logFC2 <-  Endothelium.markers.RNA.deseq$avg_logFC + 0.2783

Endothelium.markers.RNA.deseq <- Endothelium.markers.RNA.deseq[abs(Endothelium.markers.RNA.deseq$avg_logFC2) > log(2, 10),]

Endothelium.markers.RNA.deseq[Endothelium.markers.RNA.deseq$genes == "ACTB",]
Endothelium.markers.RNA.deseq[Endothelium.markers.RNA.deseq$genes == "GAPDH",]

#use DESEQ and RNA assay as result to call DEP
p <- ggplot(subset(Endothelium.markers.RNA.deseq, -log(p_val_adj) <150), aes(x=(avg_logFC2), y=-log(p_val_adj))) + geom_point(size=0.1) + theme_set(theme_cowplot(12))
p <- p + geom_label_repel(data=subset(Endothelium.markers.RNA.deseq, Endothelium.markers.RNA.deseq$genes %in% c("ACTB", "GAPDH")), aes(label=genes),color="#931F21")
p <- p + geom_point(data=subset(Endothelium.markers.RNA.deseq, Endothelium.markers.RNA.deseq$genes %in% c("ACTB","GAPDH")), color="#931F21")
p <- p + geom_smooth()
p
ggsave2(filename = "./R_work_place/Sup/S13.Vocano.plot_ACTB.pdf", p, scale=0.5)

p2 <- ggplot(Endothelium.markers.RNA.deseq, aes(x=(avg_logFC2), y=-log(p_val_adj))) + geom_point(size=0.01) + theme_set(theme_cowplot(12))
p2 <- p2 + geom_label_repel(data=subset(Endothelium.markers.RNA.deseq, Endothelium.markers.RNA.deseq$genes %in% c("TMEM100","SRGN","GREM1","CCND1","GJA4","IFI6","GJA5","CD24","MMP1","FN1","APLN", "MGLL", "NOTCH1", "DLL4")), aes(label=genes),color="#931F21")
p2 <- p2 + geom_point(data=subset(Endothelium.markers.RNA.deseq, Endothelium.markers.RNA.deseq$genes %in% c("TMEM100","SRGN","GREM1","CCND1","GJA4","IFI6","GJA5","CD24","MMP1","FN1","APLN", "MGLL", "NOTCH1","DLL4")), color="#931F21")
p2 <- p2 + xlab("LogE Average fold change") + ylab("- LogE adjusted P value") + ggtitle("Endothelium Gene Expression \n Control/HLHS")
p2
ggsave2(filename = "./R_work_place/Sup/S11.Vocano.plot.pdf", p2, scale=0.77)



Endothelium.markers.RNA.deseq.ori$genes <- row.names(Endothelium.markers.RNA.deseq.ori)

uplist <- filter(Endothelium.markers.RNA.deseq.ori, avg_logFC>log(1.25) & p_val_adj < 0.05)
downlist <- filter(Endothelium.markers.RNA.deseq.ori, avg_logFC< -log(1.25) & p_val_adj < 0.05)

dim(uplist)
dim(downlist)

#Use all gene for analysis (optional)
#alllist <- Endothelium.markers.RNA.deseq[,7]
#names(alllist) <- as.character(Endothelium.markers.RNA.deseq$genes)
#alllist <- sort(alllist, decreasing = TRUE)
#head(alllist)

#gene.df.all <- bitr(names(alllist), fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)

#ego.all <- enrichGO(gene         = gene.df.all$ENSEMBL,
#                 OrgDb         = org.Hs.eg.db,
#                 keyType       = 'ENSEMBL',
#                 ont           = "ALL",
#                 pAdjustMethod = "BH",
#                 pvalueCutoff  = 0.01,
#                 qvalueCutoff  = 0.05,
#                 readable = T)

#edox <- setReadable(ego.all, 'org.Hs.eg.db', keyType = "auto")
#p1 <- cnetplot(edox, foldChange=downlistList)
## categorySize can be scaled by 'pvalue' or 'geneNum'
#p2 <- cnetplot(downlistList, categorySize="pvalue", foldChange=alllist)
#p3 <- cnetplot(downlistList, foldChange=alllist, circular = TRUE, colorEdge = TRUE)
#p1 + p2 + p3

#specifically for upregulated gene


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
#from <- gene.df$SYMBOL
#to <- gene.df$ENTREZID
#map <- setNames(to, from)

#genenomatch <- setdiff(names(upgeneList),gene.df$SYMBOL)
#genenomatch 
#genenomatch <- gsub("MT", "mt", genenomatch)
#genenomatch <- gsub("CO", "Co", genenomatch)
#gene.nomatch  <- bitr(genenomatch, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)

#geneList <- upgeneList
#names(geneList) <- map[names(geneList)]
#genelist1 <- geneList[!duplicated(names(geneList))]

#gene <- unique(names(geneList))

#############################GO pathway analysis############################
#check GO terms in each layer, Notch in level 5
#ggo <- groupGO(gene     = gene, 
#               OrgDb    = org.Hs.eg.db,
#               ont      = "BP",
#               level =5,
#               readable = TRUE)

#table(as.character(ggo$ID) %in% GOdiscover)
#ggo$Description[as.character(ggo$ID) %in% GOnotch]

ego_up <- enrichGO(gene         = gene.df$ENTREZID,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable = T)

GOdiscover <- c("GO:0061621", "GO:0006069", "GO:0007219", "GO:0030198",  "GO:0032963", "GO:0005925",  "GO:0016477",  "GO:0001944", "GO:0048844", "GO:0042127",  "GO:0060485")
GOnotch <- c("GO:0007219")
GOcellcycle <- c("GO:0007049")
GOGly <- c("GO:0070326","GO:0000010", "GO:0060621", "GO:0006069")

GO_Up_result <- ego_up@result


table(GO_Up_result$ID %in% GOcellcycle)

GO_Up_result[GO_Up_result$ID %in% GOGly,]

ECtable <- ego_up[ego_up $ID %in% GOdiscover,]
#write.csv(file = "./R_work_place/Sup/EC_IPS_GOpathway_IPSpvalue.csv",ECtable)
ego_up$Description[1:20]


GO_Up_result_chemo <- filter(GO_Up_result, grepl("chemotaxis",GO_Up_result$Description))
write.csv(file = "../NewAnalysis/DEG/ipsc.GO.result.Up.csv",GO_Up_result)
##cannot use cutted P value
#ego3 <- gseGO(geneList     = subset(geneList,!duplicated(names(geneList))),
#              OrgDb        = org.Hs.eg.db,
#              ont          = "ALL",
#              nPerm        = 1000,
#              minGSSize    = 100,
#              maxGSSize    = 500,
#              pvalueCutoff = 0.05,
#              verbose      = FALSE)
length(downlistList)
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

table(GO_Down_result$ID %in% GOdiscover)

GO_Down_result[GO_Down_result$ID %in% GOdiscover,]

GO_Down_result[1:20,]
#GO_Down_result_chemo <- filter(GO_Down_result, grepl("chemotaxis",GO_Down_result$Description))
write.csv(file = "../NewAnalysis/DEG/ipsc.GO.result.down.csv",GO_Down_result)

dot1 <- dotplot(ego_up, showCategory = 20)
dot1
ggsave2(filename = "./R_work_place/F5.upregulated.dot.pdf", dot1)
#no positive result
dot2 <- dotplot(ego_down,showCategory = 10)
dot2
#ggsave2(filename = "./R_work_place/F6.downregulated.dot.pdf", dot2)

b1 <- barplot(ego_up, showCategory=20,horiz=F) + theme(legend.position="left")
b1
b2 <- barplot(ego_down, showCategory=10) + scale_x_discrete(position = "top") + scale_y_continuous(position="left") + scale_y_reverse() + theme(legend.position="right")
b2
b <- b1 + b2 + plot_annotation(title = "Pathway down-regulated and up-regulated in HLHS")
b

ggsave2(filename = "./R_work_place/F7a.DEPathway.pdf", b1, width =12, height = 7, units = "in", dpi = 300, scale=0.7)
ggsave2(filename = "./R_work_place/F7b.DEPathway.pdf", b2, width =12, height = 7, units = "in", dpi = 300, scale=0.7)
ggsave2(filename = "./R_work_place/F7.DEPathway.pdf", b, width =24, height = 7, units = "in", dpi = 300, scale=0.7)

## convert gene ID to Symbol
edox <- setReadable(ego_up, 'org.Hs.eg.db', keyType = "auto")
p1 <- cnetplot(edox, foldChange=upgeneList,showCategory = 20)
p1
## categorySize can be scaled by 'pvalue' or 'geneNum'
#p2 <- cnetplot(edox, categorySize="pvalue", foldChange=upgeneList)
#p3 <- cnetplot(edox, foldChange=upgeneList, circular = TRUE, colorEdge = TRUE)
ggsave2("./R_work_place/Sup/S14.upregulatedincontrol.pdf", p1,width =12, height = 9, units = "in", dpi = 300, scale=0.77)

edox <- setReadable(ego_down, 'org.Hs.eg.db', keyType = "auto")
p1 <- cnetplot(edox, foldChange=downlistList,showCategory = 15)
p1 <- p1 + coord_flip()
p1
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=downlistList)
#p3 <- cnetplot(edox, foldChange=upgeneList, circular = TRUE, colorEdge = TRUE)
p1 + p2 + p3
ggsave2("./R_work_place/Sup/S15.downregulatedincontrol.pdf", p1,width =6, height = 9, units = "in", dpi = 300)

grep("chemotaxis",GO_Down_result$Description, value = T)

#cell cycle analysis
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
#PMID: 24109597
#G1_S.gene <- c("E2F5", "CCNE1","CCNE2", "CDC25A","CDC45L","CDC6",
#               "CDKN1A","CDKN3","E2F1", "MCM2","MCM6","NPAT",
#               "PCNA","SLBP",  "BRCA1","BRCA2","CCNG2", "CDKN2C", "DHFR","MSH2",
#               "NASP", "RRM1", "RRM2","TYMS")
#G1.S.combine <- unique(c(s.genes, G1_S.gene))
#G2.combine <- unique(c(g2m.genes, G2))

#ss <- read.table("~/Desktop/67.core.cellcycle.gene.txt", header=T,sep = "\t")
#G1_S <- as.character(ss[1:16,2])
#G2 <- as.character(ss[17:67,2])


#lhls.Endothelium<- CellCycleScoring(lhls.Endothelium, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, assay="SCT")

lhls.new<- CellCycleScoring(lhls.new, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, assay="RNA")
Idents(lhls.new) <- "celltype2"
lhls.Endocardium <- subset(lhls.new, idents="Endocardium")
lhls.Endothelium <- subset(lhls.new, idents="Endothelium")

table(lhls.Endothelium$celltype)

#"polygon"
p <- ggplot(subset(lhls.new[[]], celltype2=="Endothelium"), aes(x=S.Score, y=G2M.Score)) + geom_point(size=0.2)
p <- p +  stat_density2d(aes(alpha = ..density..),contour=F, geom="tile")
p <- p + facet_grid(integrated_snn_res.0.2~condition)
p1 <- p + gridExtra::tableGrob(table(lhls.Endothelium$celltype,lhls.Endothelium$Phase))

layout <- "
AABB
AACC
"
p2 <- p + gridExtra::tableGrob(table(hlhs.Endothelium0$Phase, hlhs.Endothelium0$condition)) + gridExtra::tableGrob(table(hlhs.Endothelium3$Phase, hlhs.Endothelium3$condition)) +
   plot_layout(design = layout)
  
p2
ggsave2(filename = "./R_work_place/F8.cellcycle.B.pdf", p1)

coltable <- table(lhls.new$Phase, lhls.new$celltype)
coltable

Idents(lhls.Endothelium) <- "celltype"
hlhs.Endothelium0 <- subset(lhls.Endothelium, idents="Endothelium_0")
hlhs.Endothelium0<- CellCycleScoring(hlhs.Endothelium0, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, assay="RNA")

hlhs.Endothelium3 <- subset(lhls.Endothelium, idents="Endothelium_3")
hlhs.Endothelium3<- CellCycleScoring(hlhs.Endothelium3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, assay="RNA")


lhls.Endothelium

table(hlhs.Endothelium3$Phase, hlhs.Endothelium3$condition)

pcc0 <- ggplot(hlhs.Endothelium0[[]], aes_string(x="S.Score" , y="G2M.Score", color="Phase")) + geom_point(size=0.5) + scale_color_d3() + theme_cowplot()
#pcc0 <- pcc0 + facet_wrap(.~condition, ncol=2) + scale_color_d3()
ppc0t <- pcc0 + gridExtra::tableGrob(table(hlhs.Endothelium0$Phase, hlhs.Endothelium0$condition))
ppc0t 

pcc3 <- ggplot(hlhs.Endothelium3[[]], aes_string(x="S.Score" , y="G2M.Score", color="Phase")) + geom_point(size=0.5) + scale_color_d3()+ theme_cowplot()
#pcc3 <- pcc3 + facet_wrap(.~condition, ncol=2) + scale_color_d3()
ppc3t <- pcc3 + gridExtra::tableGrob(table(hlhs.Endothelium3$Phase, hlhs.Endothelium3$condition))
ppc3t

p <- ppc0t/ppc3t
p
ggsave2(filename = "./R_work_place/F8.cellcycle.pdf",p)

df.endothelium0 <- subset(df, df$celltype=="Endothelium_0")
table(df.endothelium0$Phase,df.endothelium0$condition)

df.endothelium3 <- subset(df, df$celltype=="Endothelium_3")
table(df.endothelium3$Phase,df.endothelium3$condition)

table(df.endothelium0$condition, df.endothelium0$Phase)
table(df.endothelium3$condition, df.endothelium3$Phase)

Idents(lhls.Endothelium) <- "condition"

lhls.Endothelium_C <- subset(lhls.Endothelium, condition=="Control")
lhls.Endothelium_C <- NormalizeData(lhls.Endothelium_C, normalization.method = "LogNormalize",scale.factor = 10000)
lhls.Endothelium_C <- CellCycleScoring(lhls.Endothelium, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, assay="RNA")

lhls.Endothelium_P <- subset(lhls.Endothelium, condition!="Control")
lhls.Endothelium_P <- NormalizeData(lhls.Endothelium_P, normalization.method = "LogNormalize",scale.factor = 10000)
lhls.Endothelium_P <- CellCycleScoring(lhls.Endothelium_P, s.features = G1.S.combine, g2m.features = g2m.genes, set.ident = TRUE, assay="RNA")

Idents(lhls.Endothelium_C) <- "Phase"
Idents(lhls.Endothelium_P) <- "Phase"

p1 <- FeatureScatter(lhls.Endothelium_C, feature2 = "S.Score", feature1 = "G2M.Score") #,group.by="integrated_snn_res.0.2"
p1 <- p1 + xlim(-0.25, 0.75) + ylim(-0.2, 0.4) + coord_flip()
#p1 <- p1 + facet_wrap(.~colors)
p1

p1 <- ggplot(lhls.Endothelium_C[[]], aes(x=S.Score, y=G2M.Score, color=Phase)) + geom_point()
p1

p1 + p1.n
#seurat_clusters
p2 <- FeatureScatter(lhls.Endothelium_P, feature2 = "S.Score", feature1 = "G2M.Score") # group.by="integrated_snn_res.0.2") 
p2 <- p2 + xlim(-0.25, 0.75) + ylim(-0.2, 0.4) + coord_flip()
#p2 <- p2 + facet_wrap(.~colors)
p2

pc <- p1 + p2

table(lhls.Endothelium_C$Phase)/1241
table(lhls.Endothelium_P$Phase)/584

RidgePlot(lhls.Endothelium, features = c("UBE2C", "TOP2A", "MCM6", "MKI67", "TP53", "CDK2", "E2F5", "CDKN1A"), ncol = 2, assay="RNA")

lhls.new <- ScaleData(lhls.new, features=row.names(lhls.new), assay = "RNA")

lhls.new <- RunPCA(lhls.new, features = c(s.genes, g2m.genes), assay = "integrated")
DimPlot(lhls.new, split.by = "condition")

Idents(lhls.new)

p <- FeatureScatter(lhls.new, feature1 = "S.Score", feature2 = "G2M.Score", group.by="condition") 
p+p1

download_LCMV()

LCMV1 <- initialize.project(datasets = "LCMV1", 
                            origins = "CD44+ cells",
                            experiments = "Rep1",
                            data.path = file.path(tempdir(), "LCMV"),
                            work.path = file.path(tempdir(), "LCMV/LCMV_analysis"))

LCMV1 <- read.data(LCMV1, subsample = 500)
LCMV2 <- LCMV1[c("work.path", "baseline.data.path","normalized","genes")]
LCMV3 <- LCMV1[c("normalized","genes")]

names(LCMV1)
LCMVcounts <- LCMV1$counts
LCMVnormal <- LCMV1$normalized
LCMVgene <- LCMV1$genes
hist(LCMVcounts)
hist(LCMVnormal)
#x is count data
hist(x)
#y is normalized count data
hist(y)

callinglist <- list(x,y)

y1 <- as.matrix(lhls.Endothelium_C@assays$RNA@data)
genes_C <- rownames(lhls.Endothelium_C)

y2 <- as.matrix(lhls.Endothelium_P@assays$RNA@data)
genes_P <- rownames(lhls.Endothelium_P)

new.list_C <- list(LCMV1$work.path,LCMV1$baseline.data.path,y1,genes_C)
new.list_P <- list(LCMV1$work.path,LCMV1$baseline.data.path,y2,genes_P)

names(new.list_C) <- c("work.path", "baseline.data.path","normalized","genes")
names(new.list_P) <- c("work.path", "baseline.data.path","normalized","genes")

cell.cycle.score.newlist_C <- cell.cycle.score(new.list_C, knn=5)
cell.cycle.score.newlist_P <- cell.cycle.score(new.list, knn=5)


p <- ggplot(cell.cycle.score.newlist_C, aes(x=S.stage, y=G2M.stage)) + geom_point()
p

rm(p)

identical(lhls.new[[]][16], lhls.new[[]][17])

cbbPalette <- c("#3399CC", "#339999", "#339966", "#993333", "#993366",  "#999933", "#CC6666")


markersUMAP <- c("FN1","NTS","APOE","COLEC11","IGFBP5","C7","NPR3","APCDD1","CFH","ENG","CGNL1","SCD","CCDC80","IGFBP3","TMEM100","LDB2","TMSB4X","SERPINE2","ADGRG6")
markers=c("FN1","CGNL1","ADGRG6","TFPI","SAT1","NPR3","HAPLN1","PLVAP","CDH11","MGLL","APLN","TGFB1","TGFB2","SNAI2","FN1","TAGLN","NOTCH1","DLL4","JAG1","HEY1")
markers_Maturation <- c("TNNI3", "TNNI1", "MYL2", "TNNT2", "TNN-N2B", "TNN-N2BA")

markers_ITGA <- grep(pattern = "^ROR", x = rownames(x = lhls.new@assays$RNA), value = TRUE)
markers_ITGA

DefaultAssay(lhls.new)
p <- VlnPlot(lhls.new, features = markers_ITGA)
p

ECmarker <- read.csv("../NewAnalysis/features/IPSC.Markers.csv", header = F)
ECmarker <- ECmarker$V1
ECmarker <- as.character(ECmarker)
DefaultAssay(lhls.Endothelium) <- "SCT"
for(gene in ECmarker){
  print(gene)
  p=FeaturePlot(object = lhls.Endothelium, features = c(gene), split.by="condition",cols=c("lightgrey", "#931F21"),pt.size=1.2, sort.cell=T)+
    theme(axis.title.x = element_text(vjust=0.5, size=24,face="bold"),axis.text.x  = element_text( vjust=0.5, size=24,face="bold"),
          axis.title.y = element_text(vjust=0.5, size=24,face="bold"),axis.text.y  = element_text( vjust=0.5, size=24,face="bold"),
          legend.text = element_text(size = 20, face = "bold"), plot.title = element_text(lineheight=.8, size=40, face="bold.italic"))
  ggsave(p,filename = paste0("../NewAnalysis/features/IPSC/",gene,".pdf"),height=6.5,width=12)
}


DefaultAssay(lhls.new) <- "SCT"
for(gene in markers){
  print(gene)
  p=FeaturePlot(object = lhls.new, features = c(gene), split.by="condition",cols=c("lightgrey", "#931F21"),pt.size=1.2, sort.cell=T)+
    theme(axis.title.x = element_text( vjust=0.5, size=24,face="bold"),axis.text.x  = element_text( vjust=0.5, size=24,face="bold"),
          axis.title.y = element_text( vjust=0.5, size=24,face="bold"),axis.text.y  = element_text( vjust=0.5, size=24,face="bold"),
          legend.text = element_text( size = 20, face = "bold"), plot.title = element_text(lineheight=.8, size=40, face="bold.italic"))
  ggsave(p,filename = paste0("./R_work_place/features/",gene,".pdf"),height=7.5,width=6)
}


s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

lhls.Endothelium<- CellCycleScoring(lhls.Endothelium, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, assay="RNA")
Idents(lhls.Endothelium) <- "celltype2"

pcellcycle <- DimPlot(lhls.Endothelium, group.by = "Phase", split.by = "condition")
pcellcycle
#lhls.all@meta.data$tmp <- factor(x = lhls.all@active.ident, levels = c(2,0,5,3,6,4,1))

for(gene in markers){
  print(gene)
  p=VlnPlot(object = lhls.new, features = c(gene),group.by="integrated_snn_res.0.2",split.by="condition",pt.size = 0 ) + scale_fill_aaas() +
    geom_boxplot(width=.05,outlier.alpha=0, fill="black",colour="gray",coef=0.1) + 
    theme(axis.title.x = element_blank(),axis.text.y  = element_text( vjust=0.5, size=34,face="bold"),
          axis.text.x  = element_text( vjust=0.5, size=34,face="bold"), plot.title = element_text(lineheight=.8, size=40, face="bold.italic"))
  #+VlnPlot(object = all, features.plot = c(gene),group.by="orig.ident")
  #p2=VlnPlot(object = all, features.plot = c(gene),group.by="orig.ident")
  ggsave(p,filename = paste0("./R_work_place/Final Violin/",gene,".pdf"),height=6,width=12)
}
p
save(lhls.all,file="iPSCEC/all.Robj")

c026=subset(lhls.all, idents = c(0,2,6))
c1345=subset(lhls.all, idents = c(1,3,4,5)) 
Idents(c026) <- c026@meta.data$orig.ident
Idents(c1345) <- c1345@meta.data$orig.ident 

c026.markers <- FindMarkers(object = c026, ident.1="M146", ident.2 = "BA060")
write.table(c026.markers,"DEG/c026.M146vsBA060.txt",quote=F,sep="\t")
c1345.markers <- FindMarkers(object = c1345, ident.1="M146", ident.2 = "BA060")
write.table(c1345.markers,"DEG/c1345.M146vsBA060.txt",quote=F,sep="\t")


genes=read.table("iPSCEC/genes_mutation_expression.txt")
genes=as.vector(genes[,1])
pdf("iPSCEC/mutated.genes.heatmap.pdf",height=12,width=9)
DoHeatmap(object = all, genes.use = genes, slim.col.label = TRUE, remove.key = TRUE,cex.row = 15,group.cex = 20)
DoHeatmap(object = all, genes.use = genes, slim.col.label = TRUE, remove.key = TRUE,cex.row = 15,group.cex = 20, group.by = "orig.ident")
dev.off()
all=SetAllIdent(all,id="orig.ident")
combined.markers <- FindMarkers(object = all, genes.use =genes, ident.1="M146", ident.2 = "BA060", min.pct = 0, logfc.threshold = 0)
write.table(combined.markers,"iPSCEC/mutated.genes.statistics.txt",quote=F,sep="\t")


kkk=sample(which(tmp$group=="others"),25)
tmp[kkk,"group"]="random"
bbb=tmp[which(tmp$group %in% c("Mutation","random")),]
t.test(bbb[which(bbb$group=="Mutation"),"avg_logFC"],bbb[which(bbb$group=="Random"),"avg_logFC"])
pdf("random.pdf",height=7,width=8)
par(mar=c(5,8,5,8))
ggplot(bbb,aes(x=group,y=avg_logFC,colour=group))+geom_boxplot(alpha=0.8)+geom_point(size=3)+
  ylab(bquote('log'['2']*'FC(HLHS/Control)'))+
  theme(axis.text.x = element_text(size=28),axis.text.y = element_text(size=28),axis.title.y = element_text(size=32,face="bold"),legend.position="none")+
  scale_x_discrete(name="")+ ggtitle("iPSC-EC: HLHS v.s. Control")+scale_fill_manual(values=c("red", "blue"))+scale_colour_manual(values=c("red", "blue"))+
  theme(plot.title = element_text(face="bold",size=36,colour="#990000",hjust=1.1))
dev.off()

genes=c("HAPLN1","PLVAP","FN1","LIMCH1","JUNB","PXDN","CD24","PALMD","SPTBN1","TGM2","NREP","LMO2","LDB2","MAP4K4","TM4SF1","CD9")
gene_order=c("FN1","PXDN","SPTBN1","TGM2","LDB2","MAP4K4","LIMCH1","JUNB","HAPLN1","PLVAP","CD24","PALMD","LMO2","NREP","TM4SF1","CD9")
degs=read.table("iPSCEC/DEGs.txt")
degs$avg_logFC=log2(exp(degs$avg_logFC))
degs=degs[which(degs$p_val_adj<0.05),]
degs=degs[which(degs$gene %in% genes),]
degs$absFC=abs(degs$avg_logFC)
degs02=degs[which(degs$cluster %in% c(0,2)),]
degs02$avg_logFC=-degs02$avg_logFC
degs6=degs[which(degs$cluster==6),]
degs=rbind(degs02,degs6)
degs2=degs %>% group_by(gene) %>% top_n(1,absFC)
degs2=degs2[order(degs2$avg_logFC),]
pdf("iPSCEC/iPSC-EC.Endocarduim.FC.barplot.pdf",height=9,width=6)
ggplot(degs2,aes(x=gene,y=avg_logFC))+geom_col(fill="forestgreen",width=0.6)+geom_hline(yintercept = 0)+scale_x_discrete(limits=rev(gene_order))+
  xlab("")+ylab(bquote('log'['2']*'Fold change'))+expand_limits(y=c(-2.5,2))+coord_flip()+
  theme(axis.title.x  = element_text( vjust=0.5, size=24,face="bold"),
        axis.text.x  = element_text( vjust=0.5, size=24,face="bold"),
        axis.text.y  = element_text( vjust=0.5, size=24,face="bold.italic"))
dev.off()

genes=c("CALM2","SRGN","S100A10","HMGB1","CD59","ENO1","CAV1","MT2A","TAGLN","HIST1H4C")
degs=read.table("iPSCEC/DEGs.txt")
degs$avg_logFC=log2(exp(degs$avg_logFC))
degs=degs[which(degs$p_val_adj<0.05),]
degs=degs[which(degs$gene %in% genes),]
degs$absFC=abs(degs$avg_logFC)
degs02=degs[which(degs$cluster %in% c(3,5)),]
degs02$avg_logFC=-degs02$avg_logFC
degs6=degs[which(degs$cluster %in% c(1,4)),]
degs=rbind(degs02,degs6)
degs2=degs %>% group_by(gene) %>% top_n(1,absFC)
degs2=degs2[order(degs2$avg_logFC),]
pdf("iPSCEC/iPSC-EC.Endothelium.FC.barplot.pdf",height=6.5,width=6)
ggplot(degs2,aes(x=gene,y=avg_logFC))+geom_col(fill="orange",width=0.6)+geom_hline(yintercept = 0)+scale_x_discrete(limits=degs1$gene_name)+
  xlab("")+ylab(bquote('log'['2']*'Fold change'))+expand_limits(y=c(-2,2))+coord_flip()+
  theme(axis.title.x  = element_text( vjust=0.5, size=24,face="bold"),
        axis.text.x  = element_text( vjust=0.5, size=24,face="bold"),
        axis.text.y  = element_text( vjust=0.5, size=24,face="bold.italic"))
dev.off()

load("iPSCEC/all.Robj")
genes=c("UBE2S","SMC4","MT2A","CAV1","CD59","HMGB1","SRGN","CALM2")
degs=read.table("iPSCEC/DEGs.txt")
degs$avg_logFC=log2(exp(degs$avg_logFC))
degs=degs[which(degs$p_val_adj<0.05),]
degs=degs[which(degs$gene %in% genes),]
degs$absFC=abs(degs$avg_logFC)
degs02=degs[which(degs$cluster %in% c(3,5)),]
degs02$avg_logFC=-degs02$avg_logFC
degs6=degs[which(degs$cluster %in% c(1,4)),]
degs=rbind(degs02,degs6)
degs2=degs %>% group_by(gene) %>% top_n(1,absFC)
degs2=degs2[order(degs2$avg_logFC),]
pdf("iPSCEC/iPSC-EC.Endothelium.FC.barplot.RV.pdf",height=6.5,width=6)
ggplot(degs2,aes(x=gene,y=avg_logFC))+geom_col(fill="orange",width=0.6)+geom_hline(yintercept = 0)+scale_x_discrete(limits=rev(genes))+
  xlab("")+ylab(bquote('log'['2']*'Fold change'))+expand_limits(y=c(-1.5,2))+coord_flip()+
  theme(axis.title.x  = element_text( vjust=0.5, size=24,face="bold"),
        axis.text.x  = element_text( vjust=0.5, size=24,face="bold"),
        axis.text.y  = element_text( vjust=0.5, size=24,face="bold.italic"))
dev.off()

genes=c("PALMD","HAPLN1","JUNB","CD24","ITM2C","TM4SF1","CD9","HLA-C")
gene_order=c("PALMD","HAPLN1","JUNB","CD24","ITM2C","TM4SF1","CD9","HLA-C")
degs=read.table("iPSCEC/DEGs.txt")
degs$avg_logFC=log2(exp(degs$avg_logFC))
degs=degs[which(degs$p_val_adj<0.05),]
degs=degs[which(degs$gene %in% genes),]
degs$absFC=abs(degs$avg_logFC)
degs02=degs[which(degs$cluster %in% c(0,2)),]
degs02$avg_logFC=-degs02$avg_logFC
degs6=degs[which(degs$cluster==6),]
degs=rbind(degs02,degs6)
degs2=degs %>% group_by(gene) %>% top_n(1,absFC)
degs2=degs2[order(degs2$avg_logFC),]
pdf("iPSCEC/iPSC-EC.Endocarduim.FC.barplot.RV.pdf",height=9,width=6)
ggplot(degs2,aes(x=gene,y=avg_logFC))+geom_col(fill="forestgreen",width=0.6)+geom_hline(yintercept = 0)+scale_x_discrete(limits=rev(gene_order))+
  xlab("")+ylab(bquote('log'['2']*'Fold change'))+expand_limits(y=c(-1.5,2))+coord_flip()+
  theme(axis.title.x  = element_text( vjust=0.5, size=24,face="bold"),
        axis.text.x  = element_text( vjust=0.5, size=24,face="bold"),
        axis.text.y  = element_text( vjust=0.5, size=24,face="bold.italic"))
dev.off()

load("iPSCEC/all.Robj")
newID=rep("c0235",nrow(all@meta.data))
newID[which(all@ident %in% c(1,4,6))]="c146"
all@meta.data$newID=newID
all=SetAllIdent(all,id="newID")
c0235.vs.c146 <- FindMarkers(object = all, ident.1="c0235", ident.2 = "c146")
tmp <- FindMarkers(object = all,genes.use = c("ETS1","NOTCH1","DLL4", "JAG1", "HEY1"), ident.1="c0235", ident.2 = "c146", min.pct=0, logfc.threshold = 0)
c0235.vs.c146=rbind(c0235.vs.c146,tmp)
write.table(c0235.vs.c146,"iPSCEC/c0235.vs.c146.txt",quote=F,sep="\t")
