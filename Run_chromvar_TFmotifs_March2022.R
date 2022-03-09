#Add in motif matrix
library(Signac)
library(Seurat)
library(chromVAR)
library(chromVARmotifs)
library(chromfunks)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
set.seed(1234)
library(motifmatchr)
library(viridis)
library(stringr)

ths.obj

## Get a list of position frequency matrices from chromvar motifs - more TFs #870
data("human_pwms_v1")
sout <- sapply(strsplit(names(human_pwms_v1), split = "_"), function(s) c(s[3]))
human_pwms_v2 <- human_pwms_v1[match(unique(sout), sout)]

DefaultAssay(ths.obj)<-"THS"
## Scan the DNA sequence of each peak for the presence of each motif
motif.matrix <- CreateMotifMatrix(
  features = StringToGRanges(rownames(ths.obj), sep = c(":", "-")),
  pwm = human_pwms_v2,#pfm,
  genome = 'hg38',
  sep = c(":", "-")
)

## Create a new Mofif object to store the results
motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = human_pwms_v2 #pfm
)
motif

## Add the Motif object to the assay
ths.obj[['THS']] <- AddMotifObject(
  object = ths.obj[['THS']],
  motif.object = motif
)

##Make motif ID to name dictionary
motif.names.dict <- GetMotifData(object = ths.obj, assay = 'THS', slot = 'motif.names')
motif.names.dict=as.data.frame(motif.names.dict)
motif.names.dict=t(motif.names.dict)
colnames(motif.names.dict)<-"motif.name"
motif.id.dict<-rownames(motif.names.dict)
chromvar.dict=cbind(motif.names.dict, motif.id.dict)
chromvar.dict=as.data.frame(chromvar.dict)
# save(chromvar.dict, file = "../chromvar_motif_id2name_HUMAN_dictionary.Rdata")

ths.obj <- RegionStats(
  object = ths.obj,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  sep = c(":", "-")
)

# save(ths.obj, file = "ths.obj_hALL.rda")

##Run Chromvar TF motif activites
BiocParallel::register(BiocParallel::SerialParam())

ths.obj <- RunChromVAR(
  object = ths.obj,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

##TF motif analysis in subclass DEGs
DefaultAssay(ths.obj) <- "chromvar"

tf.markers <- FindAllMarkers(
  object = ths.obj,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_THS'
)

# save(tf.markers, file = "TF_chromotifs_hALL.rda")

motif.id.dict=chromvar.dict$motif.id.dict
head(motif.id.dict)
motif.id.dict=gsub("_", "-", motif.id.dict)
chromvar.dict$motif.id.dict=motif.id.dict
motif.names.dict=as.character(chromvar.dict$motif.name)
head(motif.names.dict)

motif.id=tf.markers$gene
head(motif.id)
tf.markers.motif.names<-str_replace_all(motif.id, setNames(motif.names.dict, motif.id.dict))
head(tf.markers.motif.names)
tf.markers$motif.name=tf.markers.motif.names
# write.csv(tf.markers, file="TF_chromotifs_table_hALL.csv",row.names=TRUE)

##Look at TFs
tf.mark <- tf.markers[tf.markers$p_val < 0.01,]
tf.mark <- tf.mark[tf.mark$avg_logFC > 0.69,]
tf.mark <- distinct(tf.mark, gene, .keep_all = TRUE) 

# tf.mark %>% group_by(cluster) %>% top_n(5, avg_logFC) -> top5

ave.tf <- AverageExpression(ths.obj, assays = "chromvar", features = top2$gene, slot = "data" )
scaled <- t(scale(t(ave.tf$chromvar)))
scaled <- scale(scaled)

##Rename motif.id with motif name using chromvar.dictionary
library(stringr)
motif.id=rownames(scaled)
scaled.motif.names<-str_replace_all(motif.id, setNames(motif.names.dict, motif.id.dict))
rownames(scaled)=scaled.motif.names
head(rownames(scaled))

range(scaled)
scaled[scaled < 0] <- 0
# scaled[scaled > 5] <- 5

ggHeat((scaled), rescaling = "none", clustering = "none", x.lab.size = 12, y.lab.size = 12,
) 

DefaultAssay(ths.obj.epi) <- "chromvar"
FeaturePlot(
  object = ths.obj.hALL,
  features = "ENSG00000136352-LINE2353-NKX21-D", 
  cells= names(Idents(ths.obj.hALL)),
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1, reduction = "umap",
  cols = c("lightgrey", "blue")
)
