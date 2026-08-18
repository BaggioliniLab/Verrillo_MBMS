# Date: Aug 2026
# Author: Luciano Cascione
# Description: Code to perform the Gene-Set Enrichment Analysis

rm(list = ls())
setwd("~/OneDrive - USI/ABA/Antonietta_Verrillo/ATAC_Seq/")

if(TRUE){
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(dplyr)
  library(msigdbr)
  library(GSVA)
  library(fgsea)
  
  TwoWeeks_Promoter_DARs <- readRDS("./Data/DARs_in_Promoters_TwoWeeks_vs_PreExposure.rds")
  OneMonth_Promoter_DARs <- readRDS("./Data/DARs_in_Promoters_OneMonth_vs_PreExposure.rds")
  
  gs_list <- msigdbr(species = "Homo sapiens", category = "H")
  gs_list <- split(gs_list$gene_symbol, gs_list$gs_name)

  tmp <- TwoWeeks_Promoter_DARs %>% filter(!is.na(Gene.Name), Gene.Name != "") %>% distinct(Gene.Name, .keep_all = TRUE)
  gene_rank <- tmp$RankStats
  names(gene_rank) <- tmp$Gene.Name
  gene_rank <- sort(gene_rank, decreasing = TRUE)
  fgseaRes <- fgsea(pathways=gs_list, gene_rank)
  write.table(as.data.frame(fgseaRes[,1:7]), file = "./GSEA/GSEA_ATACSEQ_TwoWeeks_vs_PreExposure.Hallmarks.txt", quote = F, sep = "\t")
  saveRDS(object = fgseaRes, file = "./GSEA/GSEA_ATACSEQ_TwoWeeks_vs_PreExposure.rds")
  
  tmp <- OneMonth_Promoter_DARs %>% filter(!is.na(Gene.Name), Gene.Name != "") %>% distinct(Gene.Name, .keep_all = TRUE)
  gene_rank <- tmp$RankStats
  names(gene_rank) <- tmp$Gene.Name
  gene_rank <- sort(gene_rank, decreasing = TRUE)
  fgseaRes <- fgsea(pathways=gs_list, gene_rank)
  write.table(as.data.frame(fgseaRes[,1:7]), file = "./GSEA/GSEA_ATACSEQ_OneMonth_vs_PreExposure.Hallmarks.txt", quote = F, sep = "\t")
  saveRDS(object = fgseaRes, file = "./GSEA/GSEA_ATACSEQ_OneMonth_vs_PreExposure.rds")
  

  gs_list <- msigdbr(species = "Homo sapiens", category = "C5")
  gs_list <- dplyr::select(gs_list %>% dplyr::filter(grepl(gs_subcollection, pattern = "GO:BP")), gs_name, gene_symbol)
  gs_list <- split(gs_list$gene_symbol, gs_list$gs_name)

  tmp <- TwoWeeks_Promoter_DARs %>% filter(!is.na(Gene.Name), Gene.Name != "") %>% distinct(Gene.Name, .keep_all = TRUE)
  gene_rank <- tmp$RankStats
  names(gene_rank) <- tmp$Gene.Name
  gene_rank <- sort(gene_rank, decreasing = TRUE)
  fgseaRes <- fgsea(pathways=gs_list, gene_rank)
  write.table(as.data.frame(fgseaRes[,1:7]), file = "./GSEA/GSEA_ATACSEQ_TwoWeeks_vs_PreExposure.GOBP.txt", quote = F, sep = "\t")
  saveRDS(object = fgseaRes, file = "./GSEA/GSEA_ATACSEQ_TwoWeeks_vs_PreExposure.GOBP.rds")
  
  tmp <- OneMonth_Promoter_DARs %>% filter(!is.na(Gene.Name), Gene.Name != "") %>% distinct(Gene.Name, .keep_all = TRUE)
  gene_rank <- tmp$RankStats
  names(gene_rank) <- tmp$Gene.Name
  gene_rank <- sort(gene_rank, decreasing = TRUE)
  fgseaRes <- fgsea(pathways=gs_list, gene_rank)
  write.table(as.data.frame(fgseaRes[,1:7]), file = "./GSEA/GSEA_ATACSEQ_OneMonth_vs_PreExposure.GOBP.txt", quote = F, sep = "\t")
  saveRDS(object = fgseaRes, file = "./GSEA/GSEA_ATACSEQ_OneMonth_vs_PreExposure.GOBP.rds")

  gs_list <- msigdbr(species = "Homo sapiens", category = "C2")
  gs_list <- dplyr::select(gs_list %>% dplyr::filter(grepl(gs_subcollection, pattern = "CP:KEGG_LEGACY")), gs_name, gene_symbol)
  gs_list <- split(gs_list$gene_symbol, gs_list$gs_name)

  tmp <- TwoWeeks_Promoter_DARs %>% filter(!is.na(Gene.Name), Gene.Name != "") %>% distinct(Gene.Name, .keep_all = TRUE)
  gene_rank <- tmp$RankStats
  names(gene_rank) <- tmp$Gene.Name
  gene_rank <- sort(gene_rank, decreasing = TRUE)
  fgseaRes <- fgsea(pathways=gs_list, gene_rank)
  write.table(as.data.frame(fgseaRes[,1:7]), file = "./GSEA/GSEA_ATACSEQ_TwoWeeks_vs_PreExposure.KEGG.txt", quote = F, sep = "\t")
  saveRDS(object = fgseaRes, file = "./GSEA/GSEA_ATACSEQ_TwoWeeks_vs_PreExposure.KEGG.rds")
  
  tmp <- OneMonth_Promoter_DARs %>% filter(!is.na(Gene.Name), Gene.Name != "") %>% distinct(Gene.Name, .keep_all = TRUE)
  gene_rank <- tmp$RankStats
  names(gene_rank) <- tmp$Gene.Name
  gene_rank <- sort(gene_rank, decreasing = TRUE)
  fgseaRes <- fgsea(pathways=gs_list, gene_rank)
  write.table(as.data.frame(fgseaRes[,1:7]), file = "./GSEA/GSEA_ATACSEQ_OneMonth_vs_PreExposure.KEGG.txt", quote = F, sep = "\t")
  saveRDS(object = fgseaRes, file = "./GSEA/GSEA_ATACSEQ_OneMonth_vs_PreExposure.KEGG.rds")
}
