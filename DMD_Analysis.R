################################################################################
# Purpose: analyze DMD RNA-seq features                                        #
# Author: Ayush Jain                                          #
################################################################################
# Library Imports --------------------------------------------------------------
library(DESeq2)
library(apeglm)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(ggridges)
library(tidyverse)
library(dplyr)
library(ggfortify)
library(gplots)
library(RColorBrewer)
library(genefilter)
library(org.Mm.eg.db)
library(pathview)
library(gageData)
library(gage)
library(clusterProfiler)
library(enrichplot)
library(cowplot)
# Functions --------------------------------------------------------------------
# Read in sample sheet and append sample name as row labels
#
# Args:
#   input.file: a .txt samplesheet
#
# Returns:
#   table: a data table
ReadTable <- function(input.file) {
  table <- read.table(input.file, header=T, sep="\t")
  rownames(table) <- table$SampleName
  return(table)
}
# Subsets sample sheet to only include specificed tissue
#
# Args:
#   table: sample sheet on samples
#   age: age of animals to be analyzed
#
# Returns:
#   table: a subsetted table
SetTissue <- function(table, tissue) {
  table <- subset(table, Tissue == tissue)
  return(table)
}
# Keep only those genes that have more than 10 summed raw counts across the samples
# Combine two factors of data into one to extract relevant contrasts later
#
# Args:
#   HTseq_object: matrix of normalized gene counts
#
# Returns:
#   HTseq_object: matrix of normalized gene counts, modeled with combined factors
FilterCounts <- function(HTseq_object) {
  HTseq_object1 <- HTseq_object[rowSums(counts(HTseq_object)) > 10, ]
  
  HTseq_object1$group <- factor(paste0(HTseq_object1$Age, HTseq_object1$Genotype))
  design(HTseq_object1) <- ~ group
}
# Creates a heatmap showing sample-to-sample distances
#
# Args:
#   trans_object: transformed object either by VST or rlog
MakeDistMatrix <- function(trans_object, color) {
  sampleDists <- dist(t(assay(trans_object)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(trans_object$Age, trans_object$Genotype, trans_object$Tissue)
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette(rev(brewer.pal(9, color)))(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors, cluster_rows=F, cluster_cols=F)
}
# Creates a PCA plot to examine underlying variance between samples
#
# Args:
#   trans_object: transformed object either by VST or rlog
#   group: vectorof names to use for grouping
MakePCA <- function(trans_object, group) {
  plotPCA(object = trans_object, intgroup = group)
}
# Converts deseq object into a dataframe, labels significant resulst, and orders in Log2Fold change
#
# Args:
#   deseq_res: deseq results object
# Returns:
#   deseq_df: dataframe of ordered and labeled deseq results
CreateDF <- function(deseq_res){
  deseq_df <- as.data.frame(deseq_res)
  deseq_df <- na.omit(deseq_df)
  deseq_df <- rownames_to_column(deseq_df, "GeneID")
  deseq_df$ENTREZ <- mapIds(org.Mm.eg.db,
                            key = deseq_df$GeneID,
                            column = "ENTREZID",
                            keytype = "SYMBOL",
                            multiVals = "first")
  deseq_df <- mutate(deseq_df, sig = ifelse(deseq_df$padj<0.1, "FDR<0.1","Not Sig"))
  deseq_df[which(abs(deseq_df$log2FoldChange)<1.0), "sig"] <- "Not Sig"
  deseq_df <- deseq_df[order(abs(deseq_df$log2FoldChange), decreasing = TRUE),]
  return(deseq_df)
}
# Creates a volcano plot of log2fold change vs significance of expression
#
# Args:
#   result_df: a dataframe from a deseq results object
# Returns:
#   a volcano plot with significant genes labeled in red
MakeVolcano <- function(resdf, title){
  ggplot(resdf, aes(log2FoldChange, -log10(padj))) + 
    geom_point(aes(col=sig)) + scale_color_manual(values=c("red","blue")) +
    ggtitle(title) +
    geom_text_repel(data = head(resdf, 5), aes(label = GeneID))
}
# Creates a heatmap of top differentially expressed genes and clusters samples based on expression
#
# Args:
#   rld_object: rlog transformed count matrix
#   num_genes: number of top differentially expressed genes you want to examine

MakeHeatMap <- function(rld_object, num_genes){
  topVarGenes <- head(order(-rowVars(assay(rld_object))),num_genes)
  mat <- assay(rld_object)[ topVarGenes, ]
  mat <- mat - rowMeans(mat)
  df <- as.data.frame(colData(rld_object)[,c("Genotype","Age","Tissue")])
  rownames(df) <- colnames(mat)
  hclust_rows <- as.dendrogram(hclust(dist(mat)))  # Calculate hclust dendrograms
  hclust_cols <- as.dendrogram(hclust(dist(t(mat))))
  pheatmap(mat, annotation_col=df, fontsize_row = 7, na_col = "white", Rowv = hclust_rows,
           Colv = hclust_cols)
}

#returns a catplot that I have to save somewhere
GOenrich <- function(res, GO){
  all_genes <- as.character(rownames(res))
  signif_res <- res[res$padj < 0.1 & !is.na(res$padj),]
  signif_genes <- as.character(rownames(signif_res))
  lfc <- signif_res$log2FoldChange
  names(lfc) <- as.character(rownames(signif_res))
  ego <- enrichGO(gene = signif_genes, 
                  universe = all_genes,
                  keyType = "SYMBOL",
                  OrgDb = org.Mm.eg.db, 
                  ont = GO,
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH", 
                  qvalueCutoff = 0.01, 
                  readable = FALSE)

  catplot <- cnetplot(ego,
                      categorySize = "pvalue",
                      foldChange = lfc,
                      colorEdge = TRUE)
  #save_plot("catplot.png", catplot, base_height = 10, base_width = 14)
return(ego)
}
#data("go.sets.mm") - common gene dataset that is required to run the function below
data("go.subs.mm")
GOAnalysis <- function(resdf) {
  lfc <- resdf$log2FoldChange
  names(lfc) <- resdf$ENTREZ
  gobpsets <- go.sets.mm[go.subs.mm$BP]
  gage(exprs = lfc, gsets = gobpsets, same.dir = TRUE)
  #gobpres <- gage(exprs = lfc, gsets = gobpsets, same.dir = TRUE)
  #view(gobpres$greater)
 # view(gobpres$less)
}
#data("kegg.sets.mm") - the common gene dataset for kegg need to run this to run the function below. 
KeggAnalysis <- function(resdf) { 
  lfc <- resdf$log2FoldChange 
  names(lfc) <- resdf$ENTREZ
  gage(exprs = lfc, gsets = kegg.sets.mm, same.dir = TRUE)
  #keggres <- gage(exprs = lfc, gsets = kegg.sets.mm, same.dir = TRUE)
  #view(keggres$greater)
 # view(keggres$less)
  # keggrespathways = data.frame(id = rownames(keggres$greater), keggres$greater) %>%
  #   tibble::as_tibble() %>%
  #   filter(row_number() <= 10) %>%
  #   .$id %>%
  #   as.character()
  # keggrespathways
  # 
  # keggresids <- substr(keggrespathways, start = 1, stop = 8)
  # keggresids
  # 
  # tmp <- sapply(keggresids, function(pid) pathview(gene.data = p8medFC, pathway.id = pid, species = "mmu"))
}
# Main -------------------------------------------------------------------------
# read in the sample sheet
sampletable <- ReadTable("/Users/ayushjain/Desktop/Ayush/RNA-Seq Study/gene_counts/sample_sheet.txt")
sampletable_med <- subset(sampletable, Tissue == "med")
sampletable_cerv <- subset(sampletable, Tissue == "cerv")
sampletable_p8_med <- subset(sampletable_med, Age == "p8")
sampletable_p8_cerv <- subset(sampletable_cerv, Age == "p8")
sampletable_6mo_med <- subset(sampletable_med, Age == "6mo")
sampletable_6mo_cerv <- subset(sampletable_cerv, Age == "6mo")
sampletable_12mo_med <- subset(sampletable_med, Age == "12mo")
sampletable_12mo_cerv <- subset(sampletable_cerv, Age == "12mo")
# load in the count data to produce a count matrix
# directory is path to directory where the counts are stored (one per sample)
# design is how we wish to model the data: here we are measuring differences by genotype



ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampletable_12mo_cerv,
                                       directory = "/Users/ayushjain/Desktop/Ayush/RNA-Seq Study/gene_counts",
                                       design = ~ 1)

# filter low counts and combine factors
ddsHTSeq_filtered <- ddsHTSeq[rowSums(counts(ddsHTSeq)) > 10, ]
ddsHTSeq_filtered$Genotype <- relevel(factor(ddsHTSeq_filtered$Genotype), ref = "WT")
ddsHTSeq_filtered$group <- factor(paste0(ddsHTSeq_filtered$Age, ddsHTSeq_filtered$Genotype, ddsHTSeq_filtered$Tissue))
design(ddsHTSeq_filtered) <- ~ group

# transform matrix using variance stabilizing transformation
vsd <- vst(ddsHTSeq_filtered, blind = FALSE)


# create PCA and distance matrix to observe sample similarity - commented out PCA for now to easily do dist_map and heat map
# PCA <- MakePCA(vsd, "group")
# PCA
#dist_matrix <- MakeDistMatrix(vsd, "Blues") #- commented out dist_matrix for now to easily do heatmap
#dist_matrix

#do one with dendrogram suprressed (f,f) and another with the dendrograms present.
#heatmap <- MakeHeatMap(vsd, 40)
#heatmap
# fit statistical model
dds <- DESeq(ddsHTSeq_filtered)

# extract relevant contrasts between groups and make dot plots for each thing
setwd("/Users/ayushjain/Desktop/Ayush/Plots:Analysis/Catplots")
getwd()

#making and printing the catplots
#this uses a temp GOenrich function that retursn the catplot and has echo commented out at the end (doesn't return echo like later in the fucntion)
#if it is commented out there were no values under 0.05 p value of significance
catplot <- GOenrich(p8cerv_res, "MF")
save_plot("catplot_MF_p8cerv.png", catplot, base_height = 10, base_width = 14)

catplot<- GOenrich(p8cerv_res, "BP")
save_plot("catplot_BP_p8cerv.pdf", catplot, base_height = 10, base_width = 14)

catplot<- GOenrich(p8cerv_res, "CC")
save_plot("catplot_CC_p8cerv.pdf", catplot, base_height = 10, base_width = 14)

catplot<- GOenrich(p8med_res, "MF")
save_plot("catplot_MF_p8med.pdf", catplot, base_height = 10, base_width = 14)

catplot<- GOenrich(p8med_res, "BP")
save_plot("catplot_BP_p8med.pdf", catplot, base_height = 10, base_width = 14)

catplot<- GOenrich(p8med_res, "CC")
save_plot("catplot_CC_p8med.pdf", catplot, base_height = 10, base_width = 14)

catplot<- GOenrich(sixcerv_WTmdx_res, "MF")
save_plot("catplot_MF_6mocerv.pdf", catplot, base_height = 10, base_width = 14)

catplot<- GOenrich(sixcerv_WTmdx_res, "BP")
save_plot("catplot_BP_6mocerv.pdf", catplot, base_height = 10, base_width = 14)

catplot<- GOenrich(sixcerv_WTmdx_res, "CC")
save_plot("catplot_CC_6mocerv.pdf", catplot, base_height = 10, base_width = 14)

catplot<- GOenrich(sixmed_WTmdx_res, "MF")
save_plot("catplot_MF_6momed.pdf", catplot, base_height = 10, base_width = 14)

catplot<- GOenrich(sixmed_WTmdx_res, "BP")
save_plot("catplot_BP_6momed.pdf", catplot, base_height = 10, base_width = 14)

catplot<- GOenrich(sixmed_WTmdx_res, "CC")
save_plot("catplot_CC_6momed.pdf", catplot, base_height = 10, base_width = 14)

#catplot<- GOenrich(twelvecerv_WTmdx_res, "MF")
#save_plot("catplot_MF_12mocerv.pdf", catplot, base_height = 10, base_width = 14)

catplot<- GOenrich(twelvecerv_WTmdx_res, "BP")
save_plot("catplot_BP_12mocerv.pdf", catplot, base_height = 10, base_width = 14)

catplot<- GOenrich(twelvecerv_WTmdx_res, "CC")
save_plot("catplot_CC_12mocerv.pdf", catplot, base_height = 10, base_width = 14)

#catplot<- GOenrich(twelvemed_WTmdx_res, "MF")
#save_plot("catplot_MF_12momed.pdf", catplot, base_height = 10, base_width = 14)

catplot<- GOenrich(twelvemed_WTmdx_res, "BP")
save_plot("catplot_BP_12momed.pdf", catplot, base_height = 10, base_width = 14)

catplot<- GOenrich(twelvemed_WTmdx_res, "CC")
save_plot("catplot_CC_12momed.pdf", catplot, base_height = 10, base_width = 14)

#MF ont for both 12mo medulla and cervical spinal cord had no significant results. 
# additionally, BP and CC for 12 mo med also had the same result


#now I uncommentted out the returning of echo to change the GOEnrich function back to what it originally was to make the dotplots

getwd()
setwd("/Users/ayushjain/Desktop/Ayush/Plots:Analysis/Dotplots")
#if function below is commented out, then it didnt give anything significant
GO2 <- GOenrich(p8cerv_res, "MF")
save_plot("dotplot_MF_p8cerv.png", dotplot(GO2, font_size =10) , base_height = 10, base_width = 14)

GO2 <- GOenrich(p8cerv_res, "BP")
save_plot("dotplot_BP_p8cerv.pdf", dotplot(GO2) , base_height = 10, base_width = 14)

GO2 <- GOenrich(p8cerv_res, "CC")
save_plot("dotplot_CC_p8cerv.pdf", dotplot(GO2) , base_height = 10, base_width = 14)

GO2<- GOenrich(p8med_res, "MF")
save_plot("dotplot_MF_p8med.pdf", dotplot(GO2) , base_height = 10, base_width = 14)

GO2<- GOenrich(p8med_res, "BP")
save_plot("dotplot_BP_p8med.pdf", dotplot(GO2) , base_height = 10, base_width = 14)

GO2<- GOenrich(p8med_res, "CC")
save_plot("dotplot_CC_p8med.pdf", dotplot(GO2) , base_height = 10, base_width = 14)

GO2<- GOenrich(sixcerv_WTmdx_res, "MF")
save_plot("dotplot_MF_6mocerv.pdf", dotplot(GO2) , base_height = 10, base_width = 14)

GO2<- GOenrich(sixcerv_WTmdx_res, "BP")
save_plot("dotplot_BP_6mocerv.pdf", dotplot(GO2) , base_height = 10, base_width = 14)

GO2<- GOenrich(sixcerv_WTmdx_res, "CC")
save_plot("dotplot_CC_6mocerv.pdf", dotplot(GO2) , base_height = 10, base_width = 14)

#GO2<- GOenrich(sixmed_WTmdx_res, "MF")
#save_plot("dotplot_MF_6momed.pdf", dotplot(GO2) , base_height = 10, base_width = 14)

GO2<- GOenrich(sixmed_WTmdx_res, "BP")
save_plot("dotplot_BP_6momed.pdf", dotplot(GO2) , base_height = 10, base_width = 14)

GO2<- GOenrich(sixmed_WTmdx_res, "CC")
save_plot("dotplot_CC_6momed.pdf", dotplot(GO2) , base_height = 10, base_width = 14)

#GO2<- GOenrich(twelvecerv_WTmdx_res, "MF")
#save_plot("dotplot_MF_12mocerv.pdf", dotplot(GO2) , base_height = 10, base_width = 14)

GO2<- GOenrich(twelvecerv_WTmdx_res, "BP")
save_plot("dotplot_BP_12mocerv.pdf", dotplot(GO2) , base_height = 10, base_width = 14)

catplot<- GOenrich(twelvecerv_WTmdx_res, "CC")
save_plot("dotplot_CC_12mocerv.pdf", dotplot(GO2) , base_height = 10, base_width = 14)

#GO2<- GOenrich(twelvemed_WTmdx_res, "MF")
#save_plot("dotplot_MF_12momed.pdf", dotplot(GO2) , base_height = 10, base_width = 14)

#GO2<- GOenrich(twelvemed_WTmdx_res, "BP")
#save_plot("dotplot_BP_12momed.pdf", dotplot(GO2) , base_height = 10, base_width = 14)

#GO2<- GOenrich(twelvemed_WTmdx_res, "CC")
#save_plot("dotplot_CC_12momed.pdf", dotplot(GO2) , base_height = 10, base_width = 14)


#GO2 <- GOenrich()
#dotplot(GO2) 
# + facet_grid(ONTOLOGY~., scale = "free")

#now we got to visualize the KEGG
lfc <- p8med_resdf$log2FoldChange 
names(lfc) <- p8med_resdf$ENTREZ

twelvemed_WTmdx_resdf <- CreateDF(twelvemed_WTmdx_res)
twelvemed_WTmdx_resdf <- twelvemed_WTmdx_resdf[!duplicated(twelvemed_WTmdx_resdf$ENTREZ),]
lfc <- p8med_resdf$log2FoldChange 
names(lfc) <- p8med_resdf$ENTREZ
lfc <- sort(lfc, decreasing = TRUE)
#install.packages("ggridges", dependencies = TRUE)
#will need to edit this function for spacing purposes (kinda hard to read when I run it).
#for p8med, ont = MF and CC there was nothing under 0.05 significance
gseaKEGG <- gseKEGG(geneList = lfc, 
                  keyType = "kegg",
                  OrgDb = org.Mm.eg.db, 
                  minGSSize = 10, maxGSSize = 500,
                  pAdjustMethod = "BH",
                  pvalueCutoff = 5,
                  verbose = FALSE
                  )

ridgeplot(gseaKEGG)


p8med_res <- results(dds, contrast = c("group","p8WTmed","p8mdxmed"))
p8cerv_res <- results(dds, contrast = c("group","p8WTcerv","p8mdxcerv"))

sixmed_WTmdx_res <- results(dds, contrast = c("group", "6moWTmed", "6momdxmed"))
#sixmed_WTd52_res <- results(dds, contrast = c("group", "6mWTmed", "6md52med"))
#sixmed_mdxd52_res <- results(dds, contrast = c("group", "6mmdxmed", "6md52med"))

sixcerv_WTmdx_res <- results(dds, contrast = c("group", "6moWTcerv", "6momdxcerv"))
#sixcerv_WTd52_res <- results(dds, contrast = c("group", "6mWTcerv", "6md52cerv"))
#sixcerv_mdxd52_res <- results(dds, contrast = c("group", "6mmdxcerv", "6md52cerv"))

twelvemed_WTmdx_res <- results(dds, contrast = c("group", "12moWTmed", "12momdxmed"))
#twelvemed_WTd52_res <- results(dds, contrast = c("group", "12mWTmed", "12md52med"))
#twelvemed_mdxd52_res <- results(dds, contrast = c("group", "12mmdxmed", "12md52med"))

twelvecerv_WTmdx_res <- results(dds, contrast = c("group", "12moWTcerv", "12momdxcerv"))
#twelvecerv_WTd52_res <- results(dds, contrast = c("group", "12mWTcerv", "12md52cerv"))
#twelvecerv_mdxd52_res <- results(dds, contrast = c("group", "12mmdxcerv", "12md52cerv"))

# convert results into ordered data frames
p8med_resdf <- CreateDF(p8med_res)
p8cerv_resdf <- CreateDF(p8cerv_res)

sixmed_WTmdx_resdf <- CreateDF(sixmed_WTmdx_res)
#sixmed_WTd52_resdf <- CreateDF(sixmed_WTd52_res)
#sixmed_mdxd52_resdf <- CreateDF(sixmed_mdxd52_res)

sixcerv_WTmdx_resdf <- CreateDF(sixcerv_WTmdx_res)
#sixcerv_WTd52_resdf <- CreateDF(sixcerv_WTd52_res)
#sixcerv_mdxd52_resdf <- CreateDF(sixcerv_mdxd52_res)

twelvemed_WTmdx_resdf <- CreateDF(twelvemed_WTmdx_res)
#twelvemed_WTd52_resdf <- CreateDF(twelvemed_WTd52_res)
#twelvemed_mdxd52_resdf <- CreateDF(twelvemed_mdxd52_res)

twelvecerv_WTmdx_resdf <- CreateDF(twelvecerv_WTmdx_res)
#twelvecerv_WTd52_resdf <- CreateDF(twelvecerv_WTd52_res)
#twelvecerv_mdxd52_resdf <- CreateDF(twelvecerv_mdxd52_res)

# make volcano plots
volc1 <- MakeVolcano(p8med_resdf, "p8 mdx vs. p8 WT medulla")
volc2 <- MakeVolcano(p8cerv_resdf,"p8 mdx vs. p8 WT cervical spinal cord" )
volc3 <- MakeVolcano(sixmed_WTmdx_resdf, "6mo mdx vs. 6mo WT medulla")
volc4 <- MakeVolcano(sixcerv_WTmdx_resdf, "6mo mdx vs. 6mo WT cervical spinal cord")
volc5 <- MakeVolcano(twelvemed_WTmdx_resdf, "12mo mdx vs. 12mo WT medulla")
volc6 <- MakeVolcano(twelvecerv_WTmdx_resdf,"12mo mdx vs. 12mo WT cervical spinal cord" )
volc1
volc2
volc3
volc4
volc5
volc6
# Gene ontology and kegg pathway analysis
GOpressp8med<- GOAnalysis(p8med_resdf)
Keggpressp8med <- KeggAnalysis(p8med_resdf)

GOpressp8cerv <- GOAnalysis(p8cerv_resdf)
Keggpressp8cerv <- KeggAnalysis(p8cerv_resdf)

GOpresssixmed <- GOAnalysis(sixmed_WTmdx_resdf)
Keggpresssixmed <- KeggAnalysis(sixmed_WTmdx_resdf)

GOpresssixcerv <- GOAnalysis(sixcerv_WTmdx_resdf)
Keggpresssixcerv <- KeggAnalysis(sixcerv_WTmdx_resdf)

GOpresstwelvemed <- GOAnalysis(twelvemed_WTmdx_resdf)
Keggpresstwelvemed <- KeggAnalysis(twelvemed_WTmdx_resdf)

GOpresstwelvecerv<- GOAnalysis(twelvecerv_WTmdx_resdf)
Keggpresstwelvecerv <- KeggAnalysis(twelvecerv_WTmdx_resdf)

setwd("/Users/ayushjain/Desktop/Ayush/GOrilla_Analysis")
getwd()

#write the GOrilla results to my folder.
write.table(GOpressp8med$greater,"GO_p8med_greater.csv",sep=",")
write.table(Keggpressp8med$greater,"Kegg_p8med_greater.csv",sep=",")

write.table(GOpressp8cerv$greater,"GO_p8cerv_greater.csv",sep=",")
write.table(Keggpressp8cerv$greater,"Kegg_p8cerv_greater.csv",sep=",")

write.table(GOpresssixmed$greater,"GO_6momed_greater.csv",sep=",")
write.table(Keggpresssixmed$greater,"Kegg_6momed_greater.csv",sep=",")

write.table(GOpresssixcerv$greater,"GO_6mocerv_greater.csv",sep=",")
write.table(Keggpresssixcerv$greater,"Kegg_6mocerv_greater.csv",sep=",")

write.table(GOpresstwelvemed$greater,"GO_12momed_greater.csv",sep=",")
write.table(Keggpresstwelvemed$greater,"Kegg_12momed_greater.csv",sep=",")

write.table(GOpresstwelvecerv$greater,"GO_12mocerv_greater.csv",sep=",")
write.table(Keggpresstwelvecerv$greater,"Kegg_12mocerv_greater.csv",sep=",")

write.table(GOpressp8med$less,"GO_p8med_less.csv",sep=",")
write.table(Keggpressp8med$less,"Kegg_p8med_less.csv",sep=",")

write.table(GOpressp8cerv$less,"GO_p8cerv_less.csv",sep=",")
write.table(Keggpressp8cerv$less,"Kegg_p8cerv_less.csv",sep=",")

write.table(GOpresssixmed$less,"GO_6momed_less.csv",sep=",")
write.table(Keggpresssixmed$less,"Kegg_6momed_less.csv",sep=",")

write.table(GOpresssixcerv$less,"GO_6cerv_less.csv",sep=",")
write.table(Keggpresssixcerv$less,"Kegg_6mocerv_less.csv",sep=",")

write.table(GOpresstwelvemed$less,"GO_12momed_less.csv",sep=",")
write.table(Keggpresstwelvemed$less,"Kegg_12momed_less.csv",sep=",")

write.table(GOpresstwelvecerv$less,"GO_12mocerv_less.csv",sep=",")
write.table(Keggpresstwelvecerv$less,"Kegg_12mocerv_less.csv",sep=",")

getwd()
test <- ReadTable("/Users/ayushjain/Desktop/Ayush/RNA-Seq Study/gene_counts/sample_sheet.txt")
write.table(test, "test.csv", sep=',')
