##load RData
##extract peaks x cells matrix from snapATAC
load("snapATAC_files/hALL_merged_SnapATAC_counts.RData") #counts
load("x.sp_peakcell_postmacs2_hALL_12122019.RData")
load("peaks.gr_macs2_hALL.RData")
peak.names = paste0(seqnames(peaks.gr), ":", start(peaks.gr), "-", end(peaks.gr))
# length(peak.names) 

counts <- x.sp.mat@pmat
head(rownames(counts))
head(colnames(counts))
rownames(counts) <- x.sp.mat@barcode
colnames(counts) <- peak.names
dim(counts) 
## transpose matrix
counts = t(counts)
colnames(counts) <- gsub("@", "", colnames(counts))
colnames(counts) <- gsub("#", "_", colnames(counts))
colnames(counts) <- gsub("_NOALT", "", colnames(counts))
head(colnames(counts))

filtered.counts <- counts
# keep only peaks from chr 1-22 or chrX
chrom <- sapply(rownames(filtered.counts), function(x) strsplit(x, split = ":")[[1]][[1]])
loc <- sapply(rownames(filtered.counts), function(x) strsplit(x, split = ":")[[1]][[2]])
rownames(filtered.counts) <- paste0(chrom, ":", loc)
chrom.selected <- as.factor(chrom)
levels(chrom.selected) #chr1-22, X only
head(chrom.selected)
chrom.remove <- names(chrom.selected[chrom.selected %in% c("chr17_KI270729v1_random",
                                                           "chr22_KI270733v1_random", "chr5_GL000208v1_random",
                                                           "chr1_KI270709v1_random", "chr14_GL000225v1_random",
                                                           "chr16_KI270728v1_random", "chr17_KI270730v1_random",
                                                           "chr22_KI270735v1_random", "chr9_KI270718v1_random",
                                                           "chr1_KI270713v1_random", "chr14_KI270725v1_random",
                                                           "chr22_KI270736v1_random", "chr4_GL000008v2_random",
                                                           "chr9_KI270719v1_random", "chr1_KI270714v1_random",
                                                           "chr17_GL000205v2_random", "chr22_KI270738v1_random"
                                                           )]) #"chrUn", "chrM", "chrY",
chrom.keep <- setdiff(names(chrom.selected), chrom.remove)
filtered.counts<- counts[chrom.keep,]
dim(filtered.counts) 
counts=filtered.counts
rm(filtered.counts)

## load timepoint info
file.names <- colnames(counts)
timepoint <- sapply(file.names, function(x) strsplit(x, split = "_")[[1]][[5]])

timepoint <- factor(timepoint)
names(timepoint) <- file.names
head(timepoint)
save(timepoint, file = "timepoint.RData")

## initialize cisTopic object from count matrix
library(cisTopic)
cisTopicObject <- createcisTopicObject(counts, min.cells = 25, min.regions = 200, keepCountsMatrix = FALSE)
dim(cisTopicObject@binary.count.matrix) 
# save(cisTopicObject, file = "cisTopic_hALL1500_min25_200.RData")

## run LDA model
topics.range <- c((seq(10,45, by=5)),(seq(50, 60, by=10)))
cisTopicObject <- runModels(cisTopicObject, topic = topics.range, seed = 2018, nCores = 8,
                            burnin = 120, iterations = 200) #burnin=500, iterations = 250
# save(cisTopicObject, file="cisTopicObject_hALL1500snap_sp_runmodels.RData")

## check likelihood stablization
pdf("loglikeli.pdf", width = 6, height = 6)
logLikelihoodByIter(cisTopicObject)
dev.off()

## select for model
pdf("model_sel45.pdf", width = 6, height = 6)
cisTopicObject <- selectModel(cisTopicObject, select=45)
dev.off()

## run UMAP
cisTopicObject <- runUmap(cisTopicObject, target = 'cell',  seed=2018,n_neighbors=10L, min_dist=0.1)
## pull out umap coordinates
umap.coordinates <- cisTopicObject@dr$cell[["Umap"]]
save(umap.coordinates, file = "umap.coords_hALL.RData")

## (Explicitly running umap to get the significance matrix)
topics.emb <- t(modelMatSelection(cisTopicObject, target = "cell", method = 'Z-score'))
save(topics.emb, file = "topics.emb_hALL.RData")

##Cluster by Seurat
library(Seurat)
seurat_obj = CreateSeuratObject(cisTopicObject@binary.count.matrix, project = 'hALLsnap_tp',min.cells = 0)
seurat_obj[["topic"]] <- CreateDimReducObject(embeddings = topics.emb, key = "Topic_", assay = DefaultAssay(seurat_obj))
seurat_obj[["umap"]] <- CreateDimReducObject(embeddings = umap.coordinates, key = "Umap_", assay = DefaultAssay(seurat_obj))

## run clustering
seurat_obj <- FindNeighbors(seurat_obj, reduction = "topic", dims = 1:45, k.param = 30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.8)

## pull out clusters
seurat.clusters <-Idents(seurat_obj); names(seurat.clusters) <- colnames(seurat_obj);
head(Idents(seurat_obj))
table(seurat.clusters);
save(seurat.clusters, file="clusters_snapcisTopic_res08k30_top45_hALL.RData")

## Export clusters as a table
clusters.df <- data.frame(seurat.clusters); colnames(clusters.df) <- NULL;
write.table(clusters.df, file = "hALL_rescueAT2late_k30res08_top40_snapcistop_clusters.txt", sep = "\t", quote = F)

library(swne)
pdf("plots/umap_hALL_k30res08_top45_snapcisTopicseurat_clusters%01d.pdf", onefile = F, width = 6, height = 6)
PlotDims(umap.coordinates, sample.groups = seurat.clusters, x.lab = "umap1", y.lab = "",
        pt.size = 0.3, alpha.plot = 0.3, label.size = 4, show.legend = F, show.axes = T,
        seed = 625) #, font.size = 27
PlotDims(umap.coordinates, sample.groups = seurat.clusters, x.lab = "umap1", y.lab = "",
         pt.size = 0.3, alpha.plot = 0.3, do.label=F, label.size = 4, show.legend = T, show.axes = T,
         seed = 625)
PlotDims(umap.coordinates, sample.groups = seurat.clusters, x.lab = "", y.lab = "",
         pt.size = 0.3, alpha.plot = 0.3, do.label=F, label.size = 4, show.legend = F, show.axes = T,
         seed = 625)
dev.off()

## umap by timepoint
load("timepoint.RData")
vcol=c("D1L"="yellowgreen", "M14L"="darkgreen", "Y3L"="cyan", "Y9L"="darkblue")

pdf("plots/umap_tp_hALL_k30res08_top45_snapcisTopicseurat_clusters%01d.pdf", onefile = F, width = 6, height = 6)
PlotDims(umap.coordinates, sample.groups = timepoint, x.lab = "umap1", y.lab = "umap2",
        pt.size = 0.3, alpha.plot = 0.3, label.size = 8, show.legend = T, show.axes = T,
        seed = 625, font.size = 27, do.label = F, colors.use = vcol)
dev.off()

#----------- Cicero analysis: generate predicted gene activity matrix from chrom data
# library(methods)
# library(cisTopic)
# library(Seurat)
# library(swne)
library(cicero)
# source("chrom.R")

## Run cicero
##Make input_cds
counts=cisTopicObject@binary.count.matrix
dim(counts) #  427057  40427
peak.names <- rownames(counts)
peak.names = gsub(":","_",peak.names)
peak.names = gsub("-","_",peak.names)
rownames(counts)=peak.names

## load timepoint info
file.names <- colnames(counts)
timepoint <- sapply(file.names, function(x) strsplit(x, split = "_")[[1]][[5]])
timepoint <- factor(timepoint)
names(timepoint) <- file.names
head(timepoint)

pData <- data.frame(timepoint, seurat.clusters)
pData <- new("AnnotatedDataFrame", data = pData)
save(pData, file = "pData.RData")

## Get peak metadata
peak.names <- rownames(counts)
peak.names = gsub(":","_",peak.names)
peak.names = gsub("-","_",peak.names)
rownames(counts)=peak.names
chrom <- sapply(peak.names, function(x) strsplit(x, split = "_")[[1]][[1]])
bp1 <- as.integer(sapply(peak.names, function(x) strsplit(x, split = "_")[[1]][[2]]))
bp2 <- as.integer(sapply(peak.names, function(x) strsplit(x, split = "_")[[1]][[3]]))

fData <- data.frame(site_name = peak.names, chromosome = chrom, bp1 = bp1, bp2 = bp2)
fData <- new("AnnotatedDataFrame", data = fData)
save(fData, file = "fData.RData")

input_cds <- newCellDataSet(counts, phenoData = pData, featureData = fData,
                            expressionFamily = VGAM::binomialff())
rm(counts); invisible(gc());

input_cds <- detectGenes(input_cds)
input_cds <- estimateSizeFactors(input_cds)
# save(input_cds, file="input_cds_hALL_snapcistop.RData")

dim(input_cds@assayData$exprs) 
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap.coordinates, k = 50)

print("Done making cicero cds")
# save(cicero_cds, file="cicero_cds_hALL_top45_k30res08.RData")
hg38.chr.lengths <- read.table("hg38/hg38.chr.lengths.txt", header = F, sep = "\t") ## A tab separated text file with chromosome lengths
hg38.chr.lengths[[2]] <- as.numeric(hg38.chr.lengths[[2]])
conns <- run_cicero(cicero_cds, hg38.chr.lengths) # Takes a few minutes to run
print("Done finding connections")
# save(conns, file="conns_hALL_snapcistop.RData")

## Assign peaks to modules
ccan.assigns <- generate_ccans(conns)
# save(ccan.assigns, file="ccan.assigns_hALL.RData")

## Format peak annotation dataframe
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

load("fData.RData")
peak.names <- rownames(fData)
chrom <- sapply(peak.names, ExtractField, field = 1, delim = "_")
loc_start <- sapply(peak.names, ExtractField, field = 2, delim = "_")
loc_end <- sapply(peak.names, ExtractField, field = 3, delim = "_")

peaks.gr <- makeGRangesFromDataFrame(data.frame(seqnames = chrom, start = loc_start, end = loc_end))
peak_anno <- annotatePeak(peaks.gr, tssRegion = c(-5000, 5000),
                          TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                          annoDb = "org.Hs.eg.db")

peak_anno_df <- as.data.frame(peak_anno)
rownames(peak_anno_df) <- paste0(peak_anno_df$seqnames, "_", peak_anno_df$start, "_", peak_anno_df$end)
peak_anno_df[!grepl("Promoter|Exon|Intron", peak_anno_df$annotation), "SYMBOL"] <- NA

peak_anno_df <- peak_anno_df[,c("seqnames", "start", "end", "SYMBOL")]
colnames(peak_anno_df) <- c("chromosome", "start", "end", "gene")

input_cds <- detectGenes(input_cds)
input_cds <- estimateSizeFactors(input_cds)
dim(input_cds@assayData$exprs) #427057  40427

head(rownames(input_cds))
head(rownames(peak_anno_df))

peak_anno_df <- peak_anno_df[intersect(rownames(input_cds), rownames(peak_anno_df)),]
head(peak_anno_df)
sum(!is.na(peak_anno_df)) #1578146 peaks with genes

## Write peak annotations to bed file
# write.table(peak_anno_df, file = "hALL_k30res08_peaks_anno_filtered.bed", sep = "\t",
            # row.names = F, col.names = F)
write.table(peak_anno_df[,1:3], file = "hALL_k30res08_peaks_filtered.bed", sep = "\t",
            row.names = F, col.names = F)

## Annotate sites by gene
input_cds <- input_cds[intersect(rownames(input_cds), rownames(peak_anno_df)),]
input_cds <- annotate_cds_by_site(input_cds, peak_anno_df)
# save(input_cds, file="input_cds_GENES_hALL_k30res08.RData")

## Generate unnormalized gene activity matrix
source("scripts/chrom.R")
library(cicero)

unnorm_ga <- build_gene_activity_matrix(input_cds, conns, coaccess_cutoff = 0.1)
print(dim(unnorm_ga)) 

ga.obj <- CreateSeuratObject(unnorm_ga[,names(seurat.clusters)], project="hALL", assay="RNA")
ga.obj <- NormalizeData(ga.obj, normalization.method = "LogNormalize", scale.factor = median(Matrix::colSums(unnorm_ga)))
ga.obj[["umap"]] <- CreateDimReducObject(embeddings = umap.coordinates, key = "Umap_", assay = DefaultAssay(ga.obj))
ga.obj$clusters <- seurat.clusters
Idents(object = ga.obj) <- "clusters"
table(Idents(ga.obj))
# save(ga.obj, file = "ga.obj_hALL_k30res08.RData")
