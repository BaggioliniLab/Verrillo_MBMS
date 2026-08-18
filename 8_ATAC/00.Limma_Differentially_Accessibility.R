# Description: Limma analysis to calculate Differentially accessible regions (DARs)
# Date: August 18th, 2026
# Author: Luciano Cascione

rm(list = ls())
setwd("~/OneDrive - USI/ABA/Antonietta_Verrillo/ATAC_Seq/Data/")

library(DESeq2)
library(edgeR)
library(ggplot2)
library(ggrepel)

#load("/data/ATACSeq_Verrilio/results/bwa/merged_library/macs2/broad_peak/consensus/deseq2/consensus_peaks.mLb.clN.dds.RData")
load("/home/cascione/cascione@Bertoni/Cascione/ATACSeq_Verrilio/results/bwa/merged_library/macs2/broad_peak/consensus/deseq2/consensus_peaks.mLb.clN.dds.RData")
df <- read.table("/home/cascione/cascione@Bertoni/Cascione/ATACSeq_Verrilio/results/bwa/merged_library/macs2/broad_peak/consensus/consensus_peaks.mLb.clN.boolean.txt", header = TRUE)

class(dds)
head((dds@assays@data$counts))

counts_mtx <- (dds@assays@data$counts)
head(counts_mtx)
tmp <- read.delim("/home/cascione/cascione@Bertoni/Cascione/ATACSeq_Verrilio/results/DA_Regions_OneMonth_vs_PreExposure.txt", sep = "\t", header = T)

table(rownames(counts_mtx) %in% tmp$DAR_ID)
counts_mtx <- counts_mtx[rownames(counts_mtx) %in% tmp$DAR_ID,]
dim(counts_mtx)

sample_annotation <- data.frame(colnames(counts_mtx)); colnames(sample_annotation)[1] <- "Label"
sample_annotation$Group <- c("ONE_MONTH", "TWO_WEEKS", "ONE_MONTH", "ONE_MONTH", "TWO_WEEKS", "PRE_EXPOSURE", "TWO_WEEKS", "PRE_EXPOSURE", "PRE_EXPOSURE")

dar <- DGEList(counts = (counts_mtx), samples = sample_annotation)
dar <- calcNormFactors(dar)

group <- factor(dar$samples$Group)
designATAC <- model.matrix(~ 0 + group)  # No intercept for full contrast control
colnames(designATAC) <- levels(group)
#colnames(designATAC) <- c("ONE_MONTH", "TWO_WEEKS", "PRE_EXPOSURE")

dim(dar$counts)
v <- voomLmFit(dar, designATAC, plot = TRUE, sample.weights = T)
contrast.matrix <- makeContrasts(
  C1M = ONE_MONTH - PRE_EXPOSURE,
  C2W = TWO_WEEKS - PRE_EXPOSURE,
  PRE_EXPOSURE = PRE_EXPOSURE, ONE_MONTH = ONE_MONTH, TWO_WEEKS = TWO_WEEKS, 
  levels = designATAC)

fitATAC <- contrasts.fit(v, contrast.matrix)
fitATAC <- eBayes(fitATAC)

# Summary of significant changes
summary(decideTests(fitATAC, adjust.method = "BH", p.value = 0.001))

# Get results
resultsTwoWeeks <- topTable(fitATAC, coef = "C2W", number = Inf, adjust = "BH", sort.by = "none")
resultsTwoWeeks$DAR_ID <- rownames(resultsTwoWeeks)
head(resultsTwoWeeks)

resultsOneMonth <- topTable(fitATAC, coef = "C1M", number = Inf, adjust = "BH", sort.by = "none")
resultsOneMonth$DAR_ID <- rownames(resultsOneMonth)


resultsOneMonth$OneMonth <- topTable(fitATAC, coef = "ONE_MONTH", number = Inf, adjust = "BH", sort.by = "none")$logFC
resultsOneMonth$PRE_EXPOSURE <- topTable(fitATAC, coef = "PRE_EXPOSURE", number = Inf, adjust = "BH", sort.by = "none")$logFC
resultsOneMonth$TWO_WEEKS <- topTable(fitATAC, coef = "TWO_WEEKS", number = Inf, adjust = "BH", sort.by = "none")$logFC

resultsTwoWeeks$TWO_WEEKS <- topTable(fitATAC, coef = "TWO_WEEKS", number = Inf, adjust = "BH", sort.by = "none")$logFC
resultsTwoWeeks$PRE_EXPOSURE <- topTable(fitATAC, coef = "PRE_EXPOSURE", number = Inf, adjust = "BH", sort.by = "none")$logFC
resultsTwoWeeks$OneMonth <- topTable(fitATAC, coef = "ONE_MONTH", number = Inf, adjust = "BH", sort.by = "none")$logFC

dar_annotation <- read.delim("/home/cascione/cascione@Bertoni/Cascione/ATACSeq_Verrilio/results/bwa/merged_library/macs2/broad_peak/consensus/consensus_peaks.mLb.clN.annotatePeaks.txt",
                             heade=T)
colnames(dar_annotation)[1] <- "DAR_ID"
head(dar_annotation)

resultsOneMonth <- left_join(resultsOneMonth, dar_annotation, by="DAR_ID")
resultsTwoWeeks <- left_join(resultsTwoWeeks, dar_annotation, by="DAR_ID")

head(resultsOneMonth[resultsOneMonth$DAR_ID == "Interval_264746",])

write.table(resultsOneMonth, file = "/home/cascione/ATACSEQ_VERRILLO_DA_Regions_OneMonth_vs_PreExposure.tsv", quote = F, sep = "\t", row.names = F)
write.table(resultsTwoWeeks, file = "/home/cascione/ATACSEQ_VERRILLO_DA_Regions_TwoWeeks_vs_PreExposure.tsv", quote = F, sep = "\t", row.names = F)

saveRDS(object = resultsOneMonth, file = "/home/cascione/ATACSEQ_VERRILLO_DA_Regions_OneMonth_vs_PreExposure.rds")
saveRDS(object = resultsTwoWeeks, file = "/home/cascione/ATACSEQ_VERRILLO_DA_Regions_TwoWeeks_vs_PreExposure.rds")
