library("tximport")
library("readr")

library(DESeq2)



dataDir <- dataDir <- "/nfs/cds-peta/exports/biol_micro_cds_gr_sunagawa/scratch/ansintsova/fly/12_21_results/salmon"
sampleData <- read.csv("/nfs/cds-peta/exports/biol_micro_cds_gr_sunagawa/scratch/ansintsova/fly/12_21_results/bdgp6_samples.csv")
files <- file.path(dataDir, paste0(sampleData$sample, "_quant"), "quant.sf")
names(files) <- sampleData$sample
tx2gene_file <- "/nfs/cds-peta/exports/biol_micro_cds_gr_sunagawa/scratch/Projects_NCCR/ref/BDGP6/Drosophila_melanogaster.BDGP6.32.104.tx2gene.csv"
tx2gene <- read_csv(tx2gene_file)
txi.tx <- tximport(files, type = "salmon", txOut = TRUE)
txi.gene <- summarizeToGene(txi.tx, tx2gene)

sampleTable <-

dds <- DESeq(dds)

createCountTable <- function(files, dataType='salmon', tx=FALSE){
  if (dataType == 'salmon'){
    txi.tx <- tximport(files, type = "salmon", txOut = TRUE)
    txi.gene <- summarizeToGene(txi.tx, tx2gene)
    if (tx==TRUE) exprData <- txi.tx else txi.gene
  }
  else{
    csv <- lapply(fnames, read.csv)
    exprData <- do.call(rbind, csv)
  }
}


createDds <- function(exprData, sampleData, dataType='salmon', condition1,
                      condition2=NULL){
  require(DESeq2)
  if (condition2 != NULL){
    sampleData[condition2] <- as.factor(sampleData[[condition2]])
    formula <- paste("~", condition1, condition2)
  }
  else {
    formula <- paste("~", condition1)
  }
  sampleData[condition1] <- as.factor(sampleData[[condition1]])
  edf <- exprData[, rownames(sampleData)]
  if (dataType == 'salmon'){
    dds <- DESeqDataSetFromTximport(countData = edf,
                                    colData = sampleData,
                                    desingn = as.formula(formula))
  }
  else {
    dds <- DESeqDataSetFromMatrix(
      countData = edf, colData = sampleData,
      design = as.formula(formula)
    )
  }
  return (dds)
}


runDeseq2 <- function(dds, sampleGroup, refGroup, condition){
  require(DESeq2)
  dds[[condition]] <- relevel(dds[[condition]], ref = refGroup)
  dds <- DESeq(dds)
  res <- results(dds, contrast=c(condition, sampleGroup, refGroup))
}


writeDeseqResults <- function(results, outFile){
  
}