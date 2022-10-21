#! /usr/bin/env Rscript

library(DESeq2)
library(genefilter)
library(ggplot2)

source("./my_PCA.R")

exLim <- function(values, ratio) {
  lim.min <- min(values)
  lim.max <- max(values)
  lim.mid <- (lim.max-lim.min)/2 + lim.min
  lim.min.new <- lim.mid - (lim.max-lim.min)*(1+ratio)/2
  lim.max.new <- lim.mid + (lim.max-lim.min)*(1+ratio)/2
  lims <- list(min=lim.min.new, max=lim.max.new)
  lims
}

FC.cutoff <- 2.0

tag <- "MingTao_mRNA"
res.path <- "./"
setwd(res.path)
reads.cnt.tbl <- read.table("D:\\Stanford\\Hongchao_Input.txt",
                            stringsAsFactors=FALSE,
                            header=TRUE, sep="\t")

phenotypes <- c("WT","WT"."ALDH2","ALDH2")

colData <- as.data.frame(cbind( phenotypes=phenotypes))

rownames(reads.cnt.tbl) <- reads.cnt.tbl[ ,1]
reads.cnt.tbl <- reads.cnt.tbl[ , -1]

rownames(colData) <- names(reads.cnt.tbl)
print(colData)

cds <- DESeqDataSetFromMatrix(countData = reads.cnt.tbl,
                              colData = colData,
                              design = ~ 1 + phenotypes)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
vsd = varianceStabilizingTransformation(cds)
vsd.exp <- assay(vsd)
write.table(vsd.exp, file="vsd_exp.txt", sep="\t", quote=FALSE)

cds <- DESeq(cds)
print(resultsNames(cds))


print(colData(cds))
cds.2 <- DESeq(cds, test="LRT", reduced=~1)
print(resultsNames(cds.2))
res.LRT.2 <- results(cds.2)
write.table(res.LRT.2, file="res_WTvsALDH2.txt", sep="\t", quote=FALSE)
