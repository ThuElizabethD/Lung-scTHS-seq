# Create snap
library(SnapATAC)
# load("snapATAC_hALL.RData")
# 
# output.file <- "snapATAC_hALL_dec19.RData"

## Step 1 - create and visualize the data within a SNAP object
file.list = "hALL.merged.snap"
sample.list = "hALL"
x.sp = createSnap(file=file.list, sample=sample.list, num.cores=8)
summarySnap(x.sp)

## Step 2 - barcode selection
x.sp.filtered = filterCells(
  obj=x.sp,
  subset.names=c("fragment.num","UMI", "umap.ratio"),
  low.thresholds=c(1500, 0, 0.6),
  high.thresholds=c(Inf,Inf, 1)
)
summarySnap(x.sp.filtered)

## Step 3 - bin size selection
x.sp.filtered = addBmatToSnap(
  obj=x.sp.filtered,
  bin.size=5000,
  num.cores=8
)
# calBmatCor(x.sp.filtered) 
# save(x.sp.filtered, file = "snapATAC_files/x.sp.filtered_spbins_hALL.RData")

## Step 4 - matrix binarization
x.sp.binary = makeBinary(x.sp.filtered, mat="bmat", outlier.filter=1e-3)
summarySnap(x.sp.binary)

## Step 5 - feature selection
black_list <- read.table("hg38/hg38.blacklist.bed")
library(GenomicRanges);

black_list.gr = GRanges(
  black_list[,1], 
  IRanges(black_list[,2], black_list[,3])
);
idy = queryHits(findOverlaps(x.sp.binary@feature, black_list.gr));
if(length(idy) > 0){x.sp.binary = x.sp.binary[,-idy, mat="bmat"]};
x.sp.binary

x.sp=x.sp.binary 
bin.cov = log10(Matrix::colSums(x.sp@bmat)+1);
bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95);
idy = which(bin.cov <= bin.cutoff & bin.cov > 0);
x.sp = x.sp[, idy, mat="bmat"];
x.sp
rm(x.sp.binary, x.sp.filtered) 

## Step 6 - Jaccard index matrix generation
x.sp = runJaccard(
  x.sp,
  tmp.folder=tempdir(),
  mat = "bmat",
  max.var=2000,
  ncell.chunk=1000,
  seed.use=10,
  do.par=FALSE,
  num.cores=1
);

## Step 7 - normalization
x.sp = runNormJaccard(
  obj=x.sp,
  tmp.folder=tempdir(),
  ncell.chunk=1000,
  method="normOVE",
  row.center=TRUE,
  row.scale=TRUE,
  low.threshold=-5,
  high.threshold=5,
  do.par=TRUE,
  num.cores=5,
  seed.use=10
);

# Step 8 - linear dimensionality reduction
x.sp = runDimReduct(
  x.sp,
  pc.num=50,
  input.mat="jmat",
  method="svd",
  center=TRUE,
  scale=FALSE,
  seed.use=10
);

## Step 9 - principle components analysis
plotDimReductElbow(
  obj=x.sp,
  point.size=1,
  point.shape=19,
  point.color="red",
  point.alpha=1,
  pdf.file.name=NULL,
  pdf.height=6,
  pdf.width=6
);

plotDimReductPW(
  obj=x.sp,
  pca.dims=1:50,
  point.size=0.3,
  point.color="grey",
  point.shape=19,
  point.alpha=0.6,
  down.sample=3000,
  pdf.file.name=NULL,
  pdf.height=6,
  pdf.width=6
);

## Step 10 - K-nearest neighbor graph construction
x.sp = runKNN(
  obj=x.sp,
  pca.dims=1:16,
  weight.by.sd=FALSE,
  k=20
)

## Step 11 Clustering
# # NOTE: ledien doesn't work out well when transforming from sparse -> dense matrix
# # Sol. use igraph for clustering instead
# install.packages("leiden")
# library(leiden)
# system("pip install leidenalg python-igraph") # install leiden python package
x.sp = runCluster(
  obj=x.sp,
  tmp.folder=tempdir(),
  louvain.lib="R-igraph",
  #louvain.lib = "leiden",
  seed.use=10,
  resolution=1
  #path.to.snaptools = "/media/Home_Raid1_Voyager/dinh/miniconda3/bin/snaptools"
)

## Step 12 - non-linear dimensionality reduction
#umap
require(umap)
x.sp = runViz(
  obj=x.sp,
  tmp.folder=tempdir(),
  dims=2,
  pca.dims=1:16,
  weight.by.sd=FALSE,
  method="umap",
  fast_tsne_path=NULL,
  Y.init=NULL,
  seed.use=10,
  num.cores=5
)
# #tsne
# x.sp = runViz(
#   obj=x.sp,
#   tmp.folder=tempdir(),
#   dims=2,
#   pca.dims=1:16,
#   weight.by.sd=FALSE,
#   method="Rtsne",
#   fast_tsne_path=NULL,
#   Y.init=NULL,
#   seed.use=10,
#   num.cores=5
# )
# 
## Step 13 - umap and tsne visualization
pdf("plots/umap_hALL_snaptools_clusters_pca16.pdf", width = 6, height = 6)
plotViz(
  obj=x.sp,
  method="umap",
  point.size=0.5,
  point.shape=19,
  point.alpha=0.8,
  point.color="cluster",
  text.add=TRUE,
  text.size=1.2,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=10000,
  pdf.file.name=NULL,
  pdf.width=6,
  pdf.height=6,
  legend.add=FALSE
)
dev.off()
# 
# 
# # tsne
# pdf("plots/tsne_hALL_snaptools_clusters.pdf", width = 6, height = 6)
# plotViz(
#   obj=x.sp,
#   method="tsne",
#   point.size=0.5,
#   point.shape=19,
#   point.alpha=0.8,
#   point.color="cluster",
#   text.add=TRUE,
#   text.size=1.2,
#   text.color="black",
#   text.halo.add=TRUE,
#   text.halo.color="white",
#   text.halo.width=0.2,
#   down.sample=10000,
#   pdf.file.name=NULL,
#   pdf.width=6,
#   pdf.height=6,
#   legend.add=FALSE
# )
# dev.off()
save(x.sp, file="x.sp_preMACS2_pc16_1500_hALL.RData")
# save.image(output.file)

## Step 14 - call peaks for all
library(SnapATAC)
load("x.sp_preMACS2_pc16_1500_hALL.RData")
x.sp@file<-"hALL.merged.snap"

# install snaptools and macs2
system("pip install snaptools")
system("pip install MACS2")

Sys.time()

##pick species for genome size
peaks.gr = runMACSForAll(
  obj=x.sp,
  tmp.folder=tempdir(),
  output.prefix="hALL_clusters",
  path.to.snaptools="/.local/bin/snaptools",
  path.to.macs="/.local/bin/macs2",
  num.cores=16,
  min.cells=75,
  gsize="hs",
  buffer.size=500,
  macs.options="--nomodel --shift 37 --ext 73 --qvalue 5e-2 -B --SPMR --call-summits"
)
save(peaks.gr, file = "peaks.gr_macs2_hALL.RData")
save.image(output.file)
Sys.time()

### Step 15 - create peak by cell matrix
load("peaks.gr_macs2_hALL.RData")
x.sp.mat = createPmat(
  obj=x.sp,
  peaks=peaks.gr,
  ncell.chunk=20,
  do.par=TRUE,
  num.cores=1
)
save(x.sp.mat, file = "x.sp_peakcell_postmacs2_hALL_1core.RData")
Sys.time()

# save.image(output.file)


