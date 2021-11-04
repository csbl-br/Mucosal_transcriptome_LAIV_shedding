library(ComplexHeatmap)
library(corrplot)
library(circlize)
library(fgsea)
library(data.table)
library(clusterProfiler)
library(ggplot2)
library(fastDummies)
library(tidyr)
library(dplyr)
options(stringsAsFactors = F)

# ------------------- Function Block  ------------------- #
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
                      minSize = 15, maxSize = 2000, nperm = 1000)
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
# ------------------- End of function block  ------------------- #

### This steep is important to perform Spearman correlation between the NES (interferon signaling pathway) 
### and cell type scores for each children.


# ----------- Prepare input for GSEA ------------------- #

# Preparing rlog table to be used as a rank
rlogrank <- read.delim("data/processed/Baseline_nasal-2018_selected_children_rlogNorm.tsv")
kidslist <- colnames(rlogrank)
rlogrank$means <- rowMeans(rlogrank)
rlogrank <- rlogrank %>% arrange(means) %>% select(all_of(kidslist))
write.table(rlogrank, "data/GSEA/RlogRanks_baseline.tsv", sep ="\t", quote = F)

setwd("data/GSEA/")

fileranks <- "RlogRanks_baseline.tsv"
gmtfile <- "ReactomePathwaysLevel3.gmt"
Ptype <- "padj"
pval_cutoff <- 0.25

#Run ssGSEA
ssGSEA(gmtfile = gmtfile, 
       fileranks = fileranks,
       Ptype = Ptype,
       pval_cutoff = pval_cutoff)

