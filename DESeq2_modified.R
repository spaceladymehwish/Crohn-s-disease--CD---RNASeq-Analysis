#RNASeq analysis in R

#Install BiocManager if not already installed
#install.packages("BiocManager")

#Install DESEq2 and other packages

#BiocManager::install("SummarizedExperiment")
#BiocManager::install("DESeq2")
#install.packages("gtable")
#install.packages("stringi")
#install.packages("tidyverse")
#BiocManager::install("EnhancedVolcano")
#BiocManager::install("apeglm")
#install.packages("pheatmap")
#install.packages("PoiClaClu")

#Load required libraries
library(DESeq2)
library(tidyverse)
library(dplyr)
library(EnhancedVolcano)
library(pheatmap)
library(apeglm)
library(RColorBrewer)
library(ggplot2)

setwd("D:/2nd_sem/NGS/FINALS/Project2-RNASeq/Crohn-s-disease--CD---RNASeq-Analysis/Crohn-s-disease--CD---RNASeq-Analysis")

#Read Metadata
metadata <- read.delim("Metadata.txt")

# Import the gene raw counts
counts <- read.delim("Counts.txt")

#Analysis

#Create DESeq2 Object
dds = DESeqDataSetFromMatrix(counts,metadata,~Category)

#vst log transformation
vst <- vst(dds)


#Visualize
head(assay(vst))

#Calculate Sample-to-Sample distances

sampleDists <- dist( t( assay(vst) ) )

#using normalized data

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vst$Category, vst$Sample, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#PCA
#Simple
plotPCA(vst, intgroup = c("Category"))

#Perform DESeq analysis

#Relevel to make "normal" appear first in the data

dds$Category <- relevel(dds$Category, "non inflammatory bowel disease control")

#Run DESeq2
dds <- DESeq(dds)

#Build results table
res <- results(dds)

#Observe important information from results
mcols(res, use.names=TRUE)

#baseMean = average of normalized counts

#log2FoldChange = Effect size estimate: how much gene's expression changed due to cancer
#This is reported on a logarithmic scale to base 2. This means a log2FC of 1.5 represents 
#change in gene expression by ~3x (2^1.5=2.82)

#lfcSE: LFC Standard Error

#stat: names the test used for testing significance

#pvalue: represents the p-value

#padj: corrected p-value or FDR (False Discovery Rate) significant genes that might be false positive
# to increase confidence and reduce false positives

summary(res)

#FDR can be lowered
#res.05 <- results(dds, alpha=.05)
#summary(res.05)


#LFC2 and alpha 0.02
#lfc2res.05 <- results(dds, alpha=.05, lfcThreshold=1)
#summary(lfc2res.05)




#Filter the results to obtain genes with most significant change in expression
#Here, padj < 0.1 which is also the default
#You can lower this value to get only highly significant results

resSig <- subset(res, padj < 0.1)
summary(resSig)

#Observe in ascending order
head(resSig[ order( resSig$log2FoldChange ), ])

#Observe in descending order
head(resSig[ order( -resSig$log2FoldChange ), ]) #- so descending, not - so ascending

#Quick visual of top gene
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene=topGene, intgroup=c("Category"))

#ggplot2
data <- plotCounts(dds, gene=topGene, intgroup=c("Category","Sample"), returnData=TRUE)
ggplot(data, aes(x=Category, y=count, fill=Category)) +
  scale_y_log10() + 
  geom_dotplot(binaxis="y", stackdir="center")


#MA Plot

plotMA(res, ylim=c(-5,5))

plotMA(res, ylim=c(-5,5))
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

#Create publication grade volcano plot with marked genes of interest
EnhancedVolcano(res,
                lab = rownames(counts),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 10e-9,
                FCcutoff = 2,
                xlim = c(-5.7, 5.7),
                ylim = c(0, -log10(10.2e-12)),
                pointSize = 1.3,
                labSize = 2.6,
                title = 'Crohns disease vs non inflammatory bowel disease control',
                subtitle = 'Differential expression analysis',
                caption = 'log2fc cutoff=2; p value cutof=10e-9',
                legendPosition = "right",
                legendLabSize = 14,
                col = c('lightblue', 'orange', 'blue', 'red2'),
                colAlpha = 0.6,
                drawConnectors = TRUE,
                hline = c(10e-9),
                widthConnectors = 0.5)

#Create publication grade volcanoplot with marked genes of interest
EnhancedVolcano(res,
                lab = rownames(counts),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 10e-7,
                FCcutoff = 2.5,
                xlim = c(-5.7, 5.7),
                ylim = c(0, -log10(10.2e-12)),
                pointSize = 1.3,
                labSize = 2.6,
                title = 'The results',
                subtitle = 'Differential expression analysis',
                caption = 'log2fc cutoff=1.333; p value cutof=10e-6',
                legendPosition = "right",
                legendLabSize = 14,
                col = c('lightblue', 'orange', 'blue', 'red2'),
                colAlpha = 0.6,
                drawConnectors = TRUE,
                hline = c(10e-8),
                widthConnectors = 0.5)

#Create the final dataframe consisting of ordered deseq results based on log2fc
resord=as.data.frame(res)
finaltable1=resord[order(resord$padj),]
write.table(finaltable1, file = 'finaltable.tsv', sep = "\t",
            col.names = NA, quote = F)


#BiocManager::install("genefilter")
library("genefilter")


topVarGenes <- head(order(-rowVars(assay(vst))),20)
mat <- assay(vst)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(vst)[,c("Category","Sample")])
pheatmap(mat, annotation_col=df)

#save plot pheatmap as jpg
dev.copy(jpeg, "pheatmap.jpg")
dev.off()

