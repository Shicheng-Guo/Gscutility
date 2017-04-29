# Install related Packages
source("https://bioconductor.org/biocLite.R")
biocLite("ggplot2")
biocLite("VGAM")
biocLite("DDRTree")
biocLite("igraph")
biocLite("HSMMSingleCell")
biocLite("combinat")
biocLite("fastICA")
biocLite("irlba")
biocLite("densityClust")
biocLite("Rtsne")
biocLite("reshape2")
biocLite("dplyr")
biocLite("qlcMatrix")
biocLite("proxy")
biocLite("slam")

library(dplyr)
library(Matrix)
library(BiocParallel)
library(scran)
library(Seurat)
library(monocle)
library(reshape)

source("http://bioconductor.")
#install.packages("https://cran.r-project.org/src/contrib/densityClust_0.2.1.tar.gz", repo=NULL, type="source")
#install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/monocle_2.4.0.tar.gz", repo=NULL, type="source")
#sessionInfo()

load("/media/Home_Raid1/zhl002/NAS1/RNA_seq/hiseq_020617/seurat_analysis/monocle-latest/ipsnt_33k_4.RData")
pdf("monocle_prep.pdf")

#Prepare metadata/phenoData file
pbmc33k.downsampled <- SubsetData(pbmc33k.merged, max.cells.per.ident = 150)

phen=pbmc33k.downsampled@data.info[rownames(pbmc33k.downsampled@data.info) %in% colnames(pbmc33k.downsampled@data),]
phen[,"ident"]<-pbmc33k.downsampled@ident
pheno_data <- new("AnnotatedDataFrame", phen)
head(phen)

namevector <- c("media")
pheno_data@data[,namevector] <- NA
Name<-as.numeric(as.character(unlist(lapply(rownames(pheno_data@data),function(x) unlist(strsplit(x,"-"))[2]))))
pheno_data@data[which(Name==1),12]<-"serum"
pheno_data@data[which(Name==2),12]<-"serum"
pheno_data@data[which(Name==3),12]<-"serum"
pheno_data@data[which(Name==4),12]<-"serum"
pheno_data@data[which(Name==5),12]<-"media2i"
pheno_data@data[which(Name==6),12]<-"media2i"
pheno_data@data[which(Name==7),12]<-"media2i"
pheno_data@data[which(Name==8),12]<-"media2i"
pheno_data@data[,12]<-as.factor(pheno_data@data[,12])

namevector2 <- c("type")
pheno_data@data[,namevector2] <- NA
Name<-as.numeric(as.character(unlist(lapply(rownames(pheno_data@data),function(x) unlist(strsplit(x,"-"))[2]))))
pheno_data@data[which(Name==1),13]<-"iPS"
pheno_data@data[which(Name==2),13]<-"iPS"
pheno_data@data[which(Name==3),13]<-"SCNT"
pheno_data@data[which(Name==4),13]<-"SCNT"
pheno_data@data[which(Name==5),13]<-"iPS"
pheno_data@data[which(Name==6),13]<-"iPS"
pheno_data@data[which(Name==7),13]<-"SCNT"
pheno_data@data[which(Name==8),13]<-"SCNT"
pheno_data@data[,13]<-as.factor(pheno_data@data[,13])

namevector3 <- c("cell")
pheno_data@data[,namevector3] <- NA
Name<-as.numeric(as.character(unlist(lapply(rownames(pheno_data@data),function(x) unlist(strsplit(x,"-"))[2]))))
pheno_data@data[which(Name==1),14]<-"serum_iPS"
pheno_data@data[which(Name==2),14]<-"serum_iPS"
pheno_data@data[which(Name==3),14]<-"serum_SCNT"
pheno_data@data[which(Name==4),14]<-"serum_SCNT"
pheno_data@data[which(Name==5),14]<-"iPS_2i"
pheno_data@data[which(Name==6),14]<-"iPS_2i"
pheno_data@data[which(Name==7),14]<-"SCNT_2i"
pheno_data@data[which(Name==8),14]<-"SCNT_2i"
pheno_data@data[,14]<-as.factor(pheno_data@data[,14])

counts<-pbmc33k.downsampled@data[,(colnames(pbmc33k.downsampled@data) %in% rownames(phen))]
counts<-pbmc33k.downsampled@data[,rownames(phen)]

#QC Metrics and further filtering using scater
sce <- newSCESet(countData=counts, phenoData = pheno_data)
ntips<-convertTo(sce, type="monocle", normalize=FALSE)

#Prepare sample sheet
ntips_sample_sheet <- pData(ntips)

namevector4 <- c("gene_short_name")
fData(ntips)[,namevector4] <- rownames(fData(ntips))

#prepare featureData
ntips_gene_annotation <- fData(ntips)
pd <- pheno_data
fd <- new("AnnotatedDataFrame", data = ntips_gene_annotation)

ntips_matrix <- newCellDataSet(as.matrix(counts),
                               phenoData = pd,
                               featureData = fd)

#Make a new CellDataSet using the RNA counts
ntips <- newCellDataSet(as(as.matrix(ntips_matrix), "sparseMatrix"),
                        phenoData = pd,
                        featureData = fd,
                        lowerDetectionLimit=0.5,
                        expressionFamily=negbinomial())


ntips <- estimateSizeFactors(ntips)
HSMM <- estimateDispersions(ntips)

HSMM <- detectGenes(HSMM, min_expr = 0.01)
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))
print(head(pData(HSMM)))

# Log-transform each value in the expression matrix.
L <- log(exprs(HSMM[expressed_genes,]))
# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily"
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
# Plot the distribution of the standardized gene expression values.
qplot(value, geom="density", data=melted_dens_df) + stat_function(fun = dnorm, size=0.5, color='red') + xlab("Standardized log(read counts)") + ylab("Density")

disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.01)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
plot_ordering_genes(HSMM)
plot_pc_variance_explained(HSMM, return_all = F)
HSMM <- reduceDimension(HSMM, max_components=2, num_dim = 6, reduction_method = 'tSNE', verbose = T)
## Remove noise by PCA ...
## Reduce dimension by tSNE ...
HSMM <- clusterCells(HSMM, num_clusters=2)
## Distance cutoff calculated to 1.243471
plot_cell_clusters(HSMM, 1, 2, color="ident", markers=c("Sox2", "Anxa3"))
plot_cell_clusters(HSMM, 1, 2, color="media")
plot_cell_clusters(HSMM, 1, 2, color="type")

HSMM <- reduceDimension(HSMM, max_components=2, num_dim = 2, reduction_method = 'tSNE',residualModelFormulaStr="~media + num_genes_expressed", verbose = T) 
HSMM <- clusterCells(HSMM, num_clusters=2)
## Distance cutoff calculated to 1.277005
plot_cell_clusters(HSMM, 1, 2, color="ident")
HSMM <- clusterCells(HSMM, num_clusters=2)
## Distance cutoff calculated to 1.277005
plot_cell_clusters(HSMM, 1, 2, color="Cluster") + facet_wrap(~type)

HSMM <- clusterCells(HSMM,num_clusters=2,frequency_thresh=0.1)
plot_cell_clusters(HSMM, 1, 2, color="type", markers = c("Sox2", "Anxa3"))
pie <- ggplot(pData(HSMM), aes(x = factor(1), fill = factor(type))) + geom_bar(width = 1)
pie + coord_polar(theta = "y") + theme(axis.title.x=element_blank(), axis.title.y=element_blank())

ntips <- estimateDispersions(HSMM)

diff_test_res <- differentialGeneTest(ntips[expressed_genes,],fullModelFormulaStr="~media")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
save.image(file="monocle_prep.RData")
dev.off

