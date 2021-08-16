
library("Seurat")
library("dplyr")
library("Matrix")
library("ggpubr")
library("DoubletFinder")
library("ggsci")
library("DESeq2")
library("ggrepel")
library("scales")
library("cowplot")
library("ggsci")
library("monocle3")
library("MAST")

options(stringsAsFactors = FALSE)

load("/Users/xzhou7/Box/XinZhouFiles/Projects/ZXE7_SingleCell_HLHS/Manuscript/data/lhv.Endothelium.RData")
lhv.Endothelium 

plots5 <- DimPlot(lhv.Endothelium,group.by ="integrated_snn_res.0.2", split.by = "condition", combine = T) #+ scale_color_aaas()
plots5

table(lhv.Endothelium$condition)

perc <- as.data.frame(table(lhv.Endothelium$integrated_snn_res.0.2,lhv.Endothelium$condition)/c(1835,1835,1835,1835,1835,1835,3742,3742,3742,3742,3742,3742))

p<- ggplot(perc, aes(x=Var2, y=Freq)) + geom_bar(stat="identity")
p <- p + facet_grid(.~Var1)
p

table(lhv.Endothelium$integrated_snn_res.0.2)
table(lhv.Endothelium$condition)

########################################################################
Idents(lhv.Endothelium) <- "condition"
seurat  <- subset(lhv.Endothelium, idents="Control")
Idents(seurat) <- "integrated_snn_res.0.2"
umap.control <- DimPlot(seurat)
umap.control

gene_annotation <- as.data.frame(rownames(seurat@reductions[["pca"]]@feature.loadings), row.names = rownames(seurat@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"

cell_metadata <- as.data.frame(seurat@assays[["integrated"]]@data@Dimnames[[2]], row.names = seurat@assays[["integrated"]]@data@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"

New_matrix <- seurat@assays[["integrated"]]@data
New_matrix <- New_matrix[rownames(seurat@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix

cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)

recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

list_cluster <- seurat@meta.data[["integrated_snn_res.0.2"]]
names(list_cluster) <- seurat@assays[["integrated"]]@data@Dimnames[[2]]

cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"

cds_from_seurat@int_colData@listData$reducedDims@listData[["UMAP"]]<-seurat@reductions[["umap"]]@cell.embeddings
cds_from_seurat@preprocess_aux$gene_loadings <- seurat@reductions[["pca"]]@feature.loadings

cds_from_seurat <- learn_graph(cds_from_seurat,use_partition = F)
cds_from_seurat <- order_cells(cds_from_seurat,reduction_method = "UMAP")# ,root_pr_nodes = "Y_2")

cds_from_seurat_control <- cds_from_seurat

plot_cells(cds_from_seurat_control, 
           color_cells_by = 'cluster',
           label_roots=T,
           label_groups_by_cluster=F,
           label_leaves=F,
           label_branch_points=F,
           graph_label_size=4)

plot_cells(cds_from_seurat_control, color_cells_by = 'pseudotime',trajectory_graph_color="white", 
           label_branch_points=F,label_groups_by_cluster=F,show_trajectory_graph=F)

cluster_control <- plot_cells(cds_from_seurat_control, color_cells_by = "pseudotime",trajectory_graph_color="#AD002AFF")
cluster_control <- cluster_control + ggtitle("Control")
cluster_control

########################################################################
Idents(lhv.Endothelium) <- "condition"
seurat  <- subset(lhv.Endothelium, idents="LHV")
Idents(seurat) <- "integrated_snn_res.0.2"
DimPlot(seurat)

table(lhv.Endothelium$condition)
seurat <- subset(seurat, cells= sample(x = colnames(seurat), size = 1835))
umap.lhv <- DimPlot(seurat)
umap.lhv
seurat

gene_annotation <- as.data.frame(rownames(seurat@reductions[["pca"]]@feature.loadings), row.names = rownames(seurat@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"

cell_metadata <- as.data.frame(seurat@assays[["integrated"]]@data@Dimnames[[2]], row.names = seurat@assays[["integrated"]]@data@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"

New_matrix <- seurat@assays[["integrated"]]@data
New_matrix <- New_matrix[rownames(seurat@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix

cds_from_seurat <- new_cell_data_set(expression_matrix,cell_metadata = cell_metadata, gene_metadata = gene_annotation)

recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

list_cluster <- seurat@meta.data[["integrated_snn_res.0.2"]]
names(list_cluster) <- seurat@assays[["integrated"]]@data@Dimnames[[2]]

cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"

cds_from_seurat@int_colData@listData$reducedDims@listData[["UMAP"]]<-seurat@reductions[["umap"]]@cell.embeddings
cds_from_seurat@preprocess_aux$gene_loadings <- seurat@reductions[["pca"]]@feature.loadings

cds_from_seurat <- learn_graph(cds_from_seurat,use_partition = F)
cds_from_seurat <- order_cells(cds_from_seurat,reduction_method = "UMAP")#,root_pr_nodes = "Y_5")

plot_cells(cds_from_seurat, color_cells_by = 'cluster',
           label_groups_by_cluster=TRUE, label_leaves=FALSE, label_branch_points=TRUE,
           graph_label_size=4,trajectory_graph_color="grey")

plot_cells(cds_from_seurat, color_cells_by = 'pseudotime',trajectory_graph_color="white", 
           label_branch_points=F,label_groups_by_cluster=F,show_trajectory_graph=F)

#cds_from_seurat <- preprocess_cds(cds_from_seurat, num_dim = 50)
#cds_from_seurat <- reduce_dimension(cds_from_seurat,preprocess_method="PCA",reduction_method = "UMAP")
#cds_from_seurat <- cluster_cells(cds_from_seurat, reduction_method ="UMAP",cluster_method = "louvain")

cluster_pah <- plot_cells(cds_from_seurat, color_cells_by = "pseudotime",trajectory_graph_color="#AD002AFF")
cluster_pah <- cluster_pah + ggtitle("LHV")
cluster_pah

p.trajactory <- cluster_control + cluster_pah
p.trajactory


p.umap.samenumber <- umap.control + umap.lhv
p.umap.samenumber

Idents(lhv.Endothelium) <- "condition"
seurat.control  <- subset(lhv.Endothelium, idents="Control")
Idents(seurat.control) <- "integrated_snn_res.0.2"

p4marker <- FindMarkers(seurat.control, ident.1 = c(4))
write.csv(file = "./P4.marker.csv", p4marker)

p.f4 <- FeaturePlot(lhv.Endothelium, features = c("APLN", "NPR3","VEGFA"), min.cutoff = 0, split.by ="condition", order = T)
p.f4
