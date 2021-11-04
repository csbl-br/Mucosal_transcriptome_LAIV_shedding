library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(reshape)
library(ggrepel)
library(RColorBrewer)
library(hrbrthemes)
library(viridis)
library(circlize)
library(extrafont)
library(ggpubr)
library(data.table)
library(clusterProfiler)
library(fgsea)
library(tidyr)

read.gmt <- function(fname){
  res <- list(genes=list(), 
              desc=list())
  gmt <- file(fname)
  gmt.lines <- readLines(gmt)
  close(gmt)
  gmt.list <- lapply(gmt.lines, 
                     function(x) unlist(strsplit(x, split="\t")))
  gmt.names <- sapply(gmt.list, '[', 1)
  gmt.desc <- lapply(gmt.list, '[', 2)
  gmt.genes <- lapply(gmt.list,
                      function(x){x[3:length(x)]})
  names(gmt.desc) <- names(gmt.genes) <- gmt.names
  return(gmt.genes)
}
ssGSEA <- function(gmtfile=gmtfile,fileranks=fileranks,Ptype=Ptype,pval_cutoff=pval_cutoff){
  
  gmt <- read.gmt(gmtfile)
  ranks_df  <- read.delim(fileranks)
  names <- colnames(ranks_df)
  ranks_df$geneID <- row.names(ranks_df)
  ranks_df <- ranks_df %>% select(geneID, all_of(names)) 
  ranks_df <- ranks_df[complete.cases(ranks_df), ]
  ranks_ch <- colnames(ranks_df)[-1]
  ranks_ch <- setNames(ranks_ch, ranks_ch)
  
  #run fastGSEA
  tmp_ranks <- lapply(ranks_ch, function(rankname){
    tmpdf <- ranks_df[,c('geneID', rankname)]
    tmpdf <- tmpdf[complete.cases(tmpdf),]
    tmpranks <- tmpdf[[rankname]]
    names(tmpranks) <- tmpdf$geneID
    fgseaRes <- fgsea(pathways = gmt, stats = tmpranks, 
                      minSize = 5, maxSize = 2000, nperm = 1000)
    fgseaRes <- as.data.frame(fgseaRes)
    fgseaRes <- fgseaRes[,c('pathway', 'pval', 'padj', 'NES','size','leadingEdge')]
    fgseaRes
  })
  
  
  # Remove ranks without enrichment
  tmp_ranks <- Filter(function(x) nrow(x) > 1, tmp_ranks)
  
  # Write output NES
  rm(df)
  for(name_rank in names(tmp_ranks)){
    rank_out <- tmp_ranks[[name_rank]] 
    xxx   <- rank_out[,"NES"]
    names(xxx) <- rank_out[,"pathway"]
    df <- cbind(df,xxx)
    colnames(df)[ncol(df)] <- name_rank
  }
  df <- df[,-1]
  rank_nameout <- paste0('NES_', fileranks,sep="")
  write.table(df, file = rank_nameout, sep = '\t', quote = F, row.names = T,col.names=NA)
  
  # Write output AdjP
  rm(df)
  for(name_rank in names(tmp_ranks)){
    rank_out <- tmp_ranks[[name_rank]] 
    xxx   <- rank_out[,"padj"]
    names(xxx) <- rank_out[,"pathway"]
    df <- cbind(df,xxx)
    colnames(df)[ncol(df)] <- name_rank
  }
  df <- df[,-1]
  rank_nameout <- paste0('padj_', fileranks,sep="")
  write.table(df, file = rank_nameout, sep = '\t', quote = F, row.names = T,col.names=NA)
  
  # Write output pval
  rm(df)
  for(name_rank in names(tmp_ranks)){
    rank_out <- tmp_ranks[[name_rank]] 
    xxx   <- rank_out[,"pval"]
    names(xxx) <- rank_out[,"pathway"]
    df <- cbind(df,xxx)
    colnames(df)[ncol(df)] <- name_rank
  }
  df <- df[,-1]
  rank_nameout <- paste0('pval_', fileranks,sep="")
  write.table(df, file = rank_nameout, sep = '\t', quote = F, row.names = T,col.names=NA)
  
  # Write output Leading Edge genes
  rm(df)
  for(name_rank in names(tmp_ranks)){
    rank_out <- tmp_ranks[[name_rank]] 
    xxx   <- rank_out[,"leadingEdge"]
    names(xxx) <- rank_out[,"pathway"]
    df <- cbind(df,xxx)
    colnames(df)[ncol(df)] <- name_rank
  }
  df <- df[,-1]
  rank_nameout <- paste0('LE_', fileranks,sep="")
  write.table(df, file = rank_nameout, sep = '\t', quote = F, row.names = T,col.names=NA)
  
  
  nes <- data.table::fread(paste0('NES_', fileranks,sep=""))
  colnames(nes)[1] <- "pathway"
  pval <- data.table::fread(paste0(Ptype,'_', fileranks,sep=""))
  colnames(pval)[1] <- "pathway"
  
  nes_melt <- nes %>%
    gather(sample, nes, -pathway)
  
  pval_melt <- pval %>%
    gather(sample, pval, -pathway)
  
  result <- full_join(nes_melt, pval_melt, by=c("pathway", "sample")) %>%
    filter(pval <= pval_cutoff) %>%
    select(pathway, sample, nes) %>%
    spread(sample, nes)
  
  rank_nameout <- paste0('NES_',Ptype,pval_cutoff,"_", fileranks,sep="")
  write.table(result, file = rank_nameout, sep = '\t', quote = F, row.names = T,col.names=NA)
}
# ------------------- Function Block  ------------------- #

# The next steeps perfom the GSEA analysis using Spearman correlation as a rank and cell markers as gene set
RhoB <- read.delim("data/GSEA/B_Rho_RlogViralLoad.tsv")
RhoB$id <- row.names(RhoB)
RhoH1 <- read.delim("data/GSEA/H1_Rho_RlogViralLoad.tsv")
RhoH1$id <- row.names(RhoH1)
RhoH3 <- read.delim("data/GSEA/H3_Rho_RlogViralLoad.tsv")
RhoH3$id <- row.names(RhoH3)

RhoRank <- full_join(RhoH3,RhoH1, by = "id")
RhoRank <- full_join(RhoRank,RhoB, by = "id")
row.names(RhoRank) <-  RhoRank$id
RhoRank <- RhoRank %>% select(h1_v2_logeidml.Spearman_R, h3_v2_logeidml.Spearman_R, b_v2_logeidml.Spearman_R,
                              h1_v7_logeidml.Spearman_R, h3_v7_logeidml.Spearman_R, b_v7_logeidml.Spearman_R)
RhoRank$mean <- rowMeans(RhoRank)

RhoRank <- RhoRank %>% arrange(mean) %>% select(h1_v2_logeidml.Spearman_R, h3_v2_logeidml.Spearman_R, b_v2_logeidml.Spearman_R,
                                                 h1_v7_logeidml.Spearman_R, h3_v7_logeidml.Spearman_R, b_v7_logeidml.Spearman_R)


write.table(RhoRank, "data/deconvolution_InfluenzaMarkers/RhoRank_RlogViralload.tsv", sep = "\t", quote = F )


setwd("data/deconvolution_InfluenzaMarkers/")

gmtfile <- "../raw/Influenza_markers.gmt"
fileranks <- "RhoRank_RlogViralload.tsv"
Ptype <- "padj"
pval_cutoff <- 0.25

#Run ssGSEA
ssGSEA(gmtfile = gmtfile, 
       fileranks = fileranks,
       Ptype = Ptype,
       pval_cutoff = pval_cutoff)

rm(list = ls())

NES <- read.delim("NES_RhoRank_RlogViralload.tsv")
Padj <- read.delim("padj_RhoRank_RlogViralload.tsv")

library(reshape2)
library(dplyr)
library(ggplot2)

NES_melted <- melt(NES, id.vars = "X")
NES_melted$ID <- paste0(NES_melted$X, "_", NES_melted$variable)

Padj_melted <- melt(Padj, id.vars = "X")
Padj_melted$ID <- paste0(Padj_melted$X, "_", Padj_melted$variable)
colnames(Padj_melted)[3] <- "padj"

combined <- full_join(NES_melted,Padj_melted[,c(3,4)], by="ID")

min(-log10(combined$padj))
max(-log10(combined$padj))

# Raw plot of figure 3B
pdf("dotplot_GSEA_Celltype_RhoRank_RlogViralload.pdf")
ggplot(data = combined, aes(x = variable, y = X, color = value, size = -log10(padj))) +
  geom_point() +
  theme_bw() +
  theme(text = element_text(size = 8)) + 
  scale_color_gradient2(midpoint=0, low="blue3", mid="white", high="red3",breaks=c(-2,0,2)) +
  scale_size_continuous(range = c(2,12), breaks = c(0.5,1.0,2,2.5))
dev.off()








