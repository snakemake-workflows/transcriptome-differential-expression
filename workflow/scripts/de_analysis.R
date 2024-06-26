#!/usr/bin/env Rscript

library(stringr)
suppressMessages(library("DRIMSeq"))
suppressMessages(library("GenomicFeatures"))

cat("Loading counts, conditions and parameters.\n")
#cts <- as.matrix(read.csv("merged/all_counts.tsv", sep="\t", row.names="Reference", stringsAsFactors=FALSE))

cts <- as.matrix(read.csv(snakemake@input$tsv, sep="\t", row.names="Reference", stringsAsFactors=FALSE))

# Set up sample data frame:
coldata <- read.csv(snakemake@input$coldata, row.names="sample", sep="\t", stringsAsFactors=TRUE)
coldata$sample_id <- rownames(coldata)
coldata$condition <- factor(coldata$condition, levels=rev(levels(coldata$condition)))

de_params <- read.csv(snakemake@input$de_params, sep="\t", stringsAsFactors=FALSE)

cat("Loading annotation database.\n")
txdb <- makeTxDbFromGFF(de_params$Annotation[[1]])
# does not work, because our gtf has not enough columns!
#txdb <- makeTxDbFromGFF(snakemake@input$de_params, format = "gff")
txdf <- select(txdb, keys(txdb,"GENEID"), "TXNAME", "GENEID")
tab <- table(txdf$GENEID)
txdf$ntx<- tab[match(txdf$GENEID, names(tab))]

strip_version<-function(x) {
    tmp<-data.frame(strsplit(x,".", fixed=TRUE), stringsAsFactors=FALSE)
    tmp<-as.vector(tmp[1,])
    colnames(tmp) <- c()
    rownames(tmp) <- c()
    return(tmp)
}

rownames(cts) <- strip_version(rownames(cts))


print(head(txdf))
print("----")
print(head(cts))
print("----")
#print(rownames(cts))#
#print("----")
#print(txdf$TXNAME %>% str_subset(cts))
#cts <- cts[rownames(cts) %in% txdf$TXNAME, ] # FIXME: filter for transcripts which are in the annotation. Why they are not all there? 
#print(rownames(cts))#
#print("--x--")
#txdf <- txdf[rownames(cts) %in% txdf$TXNAME]
txdf <- txdf[txdf$TXNAME %in% rownames(cts)]
print(head(txdf))
#txdf <- txdf[txdf$TXNAME %>% str_subset(cts)]
print("xxxxxxx1Gxxxxx")
#print(rownames(cts))#

# Reorder transcript/gene database to match input counts:
txdf <- txdf[match(rownames(cts), txdf$TXNAME), ]
rownames(txdf) <- NULL
print("======")
print(head(txdf))

# Create counts data frame:
print("===dimensions===")
print(length(txdf$GENEID))
print(length(txdf$TXNAME))
counts<-data.frame(gene_id=txdf$GENEID, feature_id=txdf$TXNAME, cts)
print("=======")
print(head(counts$feature_id))
cat("Filtering counts using DRIMSeq.\n")
d <- dmDSdata(counts=counts, samples=coldata)
print("got here")
trs_cts_unfiltered <- counts(d)
#
d <- dmFilter(d, min_samps_gene_expr = de_params$min_samps_gene_expr[[1]], min_samps_feature_expr = de_params$min_samps_feature_expr[[1]],
              min_gene_expr = de_params$min_gene_expr[[1]], min_feature_expr = de_params$min_feature_expr[[1]])

cat("Building model matrix.\n")
design <- model.matrix(~condition, data=DRIMSeq::samples(d))

suppressMessages(library("dplyr"))

print(snakemake@output)

# Sum transcript counts into gene counts:
cat("Sum transcript counts into gene counts.\n")
trs_cts <- counts(d)
write.table(trs_cts, file="merged/all_counts_filtered.tsv",sep="\t")

gene_cts <- trs_cts_unfiltered %>% dplyr::select(c(1, 3:ncol(trs_cts)))  %>% group_by(gene_id) %>% summarise_all(tibble::lst(sum)) %>% data.frame()
rownames(gene_cts) <- gene_cts$gene_id
gene_cts$gene_id <- NULL
write.table(gene_cts, file="merged/all_gene_counts.tsv",sep="\t")

# Differential gene expression using edgeR:
suppressMessages(library("edgeR"))
cat("Running differential gene expression analysis using edgeR.\n")

y <- DGEList(gene_cts)
y <- calcNormFactors(y)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit)
edger_res <- topTags(qlf, n=nrow(y), sort.by="PValue")[[1]]

pdf("de_analysis/results_dge.pdf")
plotMD(qlf)
abline(h=c(-1,1), col="blue")
plotQLDisp(fit)

write.table(as.data.frame(edger_res), file="de_analysis/results_dge.tsv", sep="\t")

# Differential transcript usage using DEXSeq:
suppressMessages(library("DEXSeq"))
cat("Running differential transcript usage analysis using DEXSeq.\n")

sample.data<-DRIMSeq::samples(d)
count.data <- round(as.matrix(counts(d)[,-c(1:2)]))
dxd <- DEXSeqDataSet(countData=count.data, sampleData=sample.data, design=~sample + exon + condition:exon, featureID=trs_cts$feature_id, groupID=trs_cts$gene_id)
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd)
dxd <- testForDEU(dxd, reducedModel=~sample + exon)
dxd <- estimateExonFoldChanges( dxd, fitExpToVar="condition")
dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)

dev.off()
pdf("de_analysis/results_dtu.pdf")
plotMA(dxr, cex=0.8, alpha=0.05) 
plotDispEsts(dxd)

qval <- perGeneQValue(dxr) 
dxr.g<-data.frame(gene=names(qval), qval)
dxr.g <- dxr.g[order(dxr.g$qval),]

dxr_out <- as.data.frame(dxr[,c("featureID", "groupID", "pvalue")])
dxr_out <- dxr_out[order(dxr$pvalue),]

write.table(dxr.g, file="de_analysis/results_dtu_gene.tsv", sep="\t")
write.table(dxr_out, file="de_analysis/results_dtu_transcript.tsv", sep="\t")

# and writing out some of the DEXSeq metrics to accompany EPI2ME Labs tutorial
colnames(dxr)[grep("log2fold", colnames(dxr))] <- "log2fold"
MADTUdata <- data.frame(dxr)[order(dxr$padj),c("exonBaseMean", "log2fold", "pvalue", "padj")]
MADTUdata$exonBaseMean <- log2(MADTUdata$exonBaseMean)
colnames(MADTUdata)[which(colnames(MADTUdata)=="exonBaseMean")] <- "Log2MeanExon"
colnames(MADTUdata)[which(colnames(MADTUdata)=="log2fold")] <- "Log2FC"
write.table(MADTUdata, file="de_analysis/results_dexseq.tsv", sep="\t")

# stageR analysis of DEXSeq results:
cat("stageR analysis\n")
library(stageR)

cat("Running stageR analysis on the differential transcript usage results.\n")
pConfirmation <- matrix(dxr$pvalue, ncol=1)
dimnames(pConfirmation) <- list(dxr$featureID, "transcript")
pScreen <- qval
tx2gene <- as.data.frame(dxr[,c("featureID", "groupID")])

stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation, pScreenAdjusted=TRUE, tx2gene=tx2gene)
# note: the choice of 0.05 here means you can *only* threshold at 5% OFDR later
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05)
suppressWarnings({dex.padj <- getAdjustedPValues(stageRObj, order=FALSE, onlySignificantGenes=FALSE)})

write.table(dex.padj, file="de_analysis/results_dtu_stageR.tsv", sep="\t")

