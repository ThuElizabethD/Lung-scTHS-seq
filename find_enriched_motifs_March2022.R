##Find motif enrichments for peaks belonging to linked genes
library(Seurat)
library(Signac)
library(dplyr)

##Use cDARs linked to peaks
#Subset cDARs for cluster or gene(s)
# peak_anno_df_ac=subset(peak_anno_df_ac, cluster==1)
peak_anno_df_ac=subset(peak_anno_df_ac, gene=="Myrf" | gene=="Pdgfa")
unique.genes<-peak_anno_df_ac$gene
unique.genes<-unique(unique.genes)
unique.genes<-as.character(unique.genes)
unique.genes
rownames(peak_anno_df_ac) <- paste0(peak_anno_df_ac$chromosome, ":", peak_anno_df_ac$start, "-", peak_anno_df_ac$end) 

# peak_anno_df=subset(peak_anno_df, gene=="Pdgfra")
# peak_anno_df=subset(peak_anno_df, gene=="Acta2")
# unique.genes<-peak_anno_df$gene
# unique.genes<-unique(unique.genes)
# rownames(peak_anno_df) <- paste0(peak_anno_df$chromosome, ":", peak_anno_df$start, "-", peak_anno_df$end) 
# peak_anno_df=distinct(peak_anno_df) #remove duplicate rows

ths.obj.p1
DefaultAssay(ths.obj.p1)<-"THS"

##Take linked genes with more than 1 peak, identify peaks that belong to each cluster and find overrep motifs
gene.tfbs.list <- lapply(unique.genes, function(cl) {
  print(paste("Running for gene:", cl))
  
  cl.sites <- rownames(peak_anno_df_ac) #
  EM <- FindMotifs(object = ths.obj,
                   features = cl.sites)
  EM <- EM[EM$pvalue < 0.05,]
  EM
})

##Make sure regionstats ran and motif added

names(gene.tfbs.list) <- unique.genes
sapply(gene.tfbs.list, nrow) ## Check the number of motifs for each gene

# save(gene.tfbs.list, file = "gene.tfbs.list_pdgfa_P1mLung.RData")

#Add associated cluster for each marker
gene.tfbs <- do.call("rbind", lapply(gene.tfbs.list, as.data.frame))
marker<-unlist(lapply(rownames(gene.tfbs),function(x) unlist(strsplit(x,"[.]"))[1]))
gene.tfbs$gene <- marker

cl.order <- peak_anno_df_ac$cluster
names(cl.order) <- peak_anno_df_ac$gene
cl.order <- factor(cl.order) 
cl.order <- cl.order[gene.tfbs$gene]
gene.tfbs$cluster <- cl.order
# write.csv(gene.tfbs, file = "gene.tfbs_pdgfra_Enrichments_MESm_table.csv")
