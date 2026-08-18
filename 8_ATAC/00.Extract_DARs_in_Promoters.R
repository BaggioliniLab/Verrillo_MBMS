# Description: Script to extract and save the regions close to TSS (genes' promoter)
# Date: August 18th, 2026
# Author: Luciano Cascione

rm(list = ls())
setwd("~/OneDrive - USI/ABA/Antonietta_Verrillo/ATAC_Seq/Data/")

library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(dplyr)
library(msigdbr)
library(GSVA)
library(fgsea)
  
DARs <- as.data.frame(readRDS("./ATACSEQ_VERRILLO_DARs_TwoWeeks_vs_PreExposure.rds"))
TwoWeeks_Promoter_DARs <- DARs[grepl(DARs$Annotation, pattern="promoter"),]; dim(TwoWeeks_Promoter_DARs)

DARs <- as.data.frame(readRDS("./ATACSEQ_VERRILLO_DARs_OneMonth_vs_PreExposure.rds"))
OneMonth_Promoter_DARs <- DARs[grepl(DARs$Annotation, pattern="promoter"),]; dim(OneMonth_Promoter_DARs)

TwoWeeks_Promoter_DARs$RankStats <- -log10(TwoWeeks_Promoter_DARs$P.Value)*sign(TwoWeeks_Promoter_DARs$logFC)
OneMonth_Promoter_DARs$RankStats <- -log10(OneMonth_Promoter_DARs$P.Value)*sign(OneMonth_Promoter_DARs$logFC)

TwoWeeks_Promoter_DARs <- TwoWeeks_Promoter_DARs %>% filter(!is.na(Gene.Name), Gene.Name != "") %>% distinct(Gene.Name, .keep_all = TRUE)
write.table(TwoWeeks_Promoter_DARs, file = "./DARs_in_Promoters_TwoWeeks_vs_PreExposure.tsv", quote = F, sep = "\t")
saveRDS(object = TwoWeeks_Promoter_DARs, file = "./DARs_in_Promoters_TwoWeeks_vs_PreExposure.rds")
  
OneMonth_Promoter_DARs <- OneMonth_Promoter_DARs %>% filter(!is.na(Gene.Name), Gene.Name != "") %>% distinct(Gene.Name, .keep_all = TRUE)
write.table(OneMonth_Promoter_DARs, file = "./DARs_in_Promoters_OneMonth_vs_PreExposure.tsv", quote = F, sep = "\t")
saveRDS(object = OneMonth_Promoter_DARs, file = "./DARs_in_Promoters_OneMonth_vs_PreExposure.rds")