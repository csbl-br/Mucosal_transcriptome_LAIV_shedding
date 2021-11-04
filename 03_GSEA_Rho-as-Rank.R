library(ComplexHeatmap)
library(corrplot)
library(circlize)
library(fgsea)
library(data.table)
#library(clusterProfiler)
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


# ----------- Prepare input for GSEA ------------------- #
system("mkdir data/GSEA")


# Preparing Rho table to be used as a rank
CortableB <- read.delim("data/correlation/spearman_table_rlogViralLoadB.tsv")
CortableB$mean <-  apply(X =CortableB[,1:2], 1, mean )
CortableB <- CortableB %>% arrange(mean)
CortableB <- CortableB %>% select(b_v2_logeidml.Spearman_R, b_v7_logeidml.Spearman_R)
write.table(CortableB, "data/GSEA/B_Rho_RlogViralLoad.tsv", sep = "\t", quote = F)

CortableH1 <- read.delim("data/correlation/spearman_table_rlogViralLoadH1.tsv")
CortableH1$mean <-  apply(X =CortableH1[,1:2], 1, mean )
CortableH1 <- CortableH1 %>% arrange(mean)
CortableH1 <- CortableH1 %>% select(h1_v2_logeidml.Spearman_R, h1_v7_logeidml.Spearman_R)
write.table(CortableH1, "data/GSEA/H1_Rho_RlogViralLoad.tsv",sep = "\t", quote = F)

CortableH3 <- read.delim("data/correlation/spearman_table_rlogViralLoadH3.tsv")
CortableH3$mean <-  apply(X =CortableH3[,1:2], 1, mean )
CortableH3 <- CortableH3 %>% arrange(mean)
CortableH3 <- CortableH3 %>% select(h3_v2_logeidml.Spearman_R, h3_v7_logeidml.Spearman_R)
write.table(CortableH3, "data/GSEA/H3_Rho_RlogViralLoad.tsv",sep = "\t", quote = F)
setwd("data/GSEA/")
system("ln -s ../raw/ReactomePathwaysLevel3.gmt")

# ----------- Reading rank data - Reapeat the the GSEA function for each rank ------------------- #
gmtfile <- "ReactomePathwaysLevel3.gmt"
Ptype <- "padj"
pval_cutoff <- 0.25

toremove <- setdiff(ls(), c("ssGSEA", "read.gmt", "gmtfile", "Ptype", "pval_cutoff"))
rm(list = toremove)

fileranks <- "H3_Rho_RlogViralLoad.tsv"
#Run ssGSEA
ssGSEA(gmtfile = gmtfile, 
       fileranks = fileranks,
       Ptype = Ptype,
       pval_cutoff = pval_cutoff)

toremove <- setdiff(ls(), c("ssGSEA", "read.gmt", "gmtfile", "Ptype", "pval_cutoff"))
rm(list = toremove)

fileranks <- "B_Rho_RlogViralLoad.tsv"
ssGSEA(gmtfile = gmtfile, 
       fileranks = fileranks,
       Ptype = Ptype,
       pval_cutoff = pval_cutoff)


toremove <- setdiff(ls(), c("ssGSEA", "read.gmt", "gmtfile", "Ptype", "pval_cutoff"))
rm(list = toremove)

fileranks <- "H1_Rho_RlogViralLoad.tsv"
ssGSEA(gmtfile = gmtfile, 
       fileranks = fileranks,
       Ptype = Ptype,
       pval_cutoff = pval_cutoff)


# output files were inspected manually and pathways selected to produced file is "Pathways_to_display.tsv"
# The file Pathways_to_display.tsv is included in "data/raw"

setwd("../../")
rm(list = ls())

# ----------- create Corplots using GSEA curated output ------------------- #
selPathways <- read.delim("data/raw/Pathways_to_display.tsv")
selPathways <- selPathways$Reactome.pathway
#selPathways <- c("Nucleotide-binding domain, leucine rich repeat containing receptor NLR signaling pathways", selPathways)


# GSEA p values and NES has changes since last version - using 2020 tables to recreate exactly values from 2020 
B_NESpval <- read.delim("data/GSEA/padj_B_Rho_RlogViralLoad.tsv")
B_NES <- read.delim("data/GSEA/NES_B_Rho_RlogViralLoad.tsv")

H1_NESpval <- read.delim("data/GSEA/padj_H1_Rho_RlogViralLoad.tsv")
H1_NES <- read.delim("data/GSEA/NES_H1_Rho_RlogViralLoad.tsv")

H3_NESpval <- read.delim("data/GSEA/padj_H3_Rho_RlogViralLoad.tsv")
H3_NES <- read.delim("data/GSEA/NES_H3_Rho_RlogViralLoad.tsv")

NES <- full_join(H1_NES, H3_NES, by="X" )
NES <- full_join(NES, B_NES, by="X" )
row.names(NES) <-  NES$X

Padj <- full_join(H1_NESpval, H3_NESpval, by="X" )
Padj <- full_join(Padj, B_NESpval, by="X" )
row.names(Padj) <-  Padj$X

all(row.names(Padj) %in% row.names(NES))
all(row.names(Padj)  == row.names(NES) )

library(reshape2)
NES_melted <- melt(NES, id.vars = "X")
NES_melted$ID <- paste0(NES_melted$X, "_", NES_melted$variable)

Padj_melted <- melt(Padj, id.vars = "X")
Padj_melted$ID <- paste0(Padj_melted$X, "_", Padj_melted$variable)
colnames(Padj_melted)[3] <- "padj"

combined <- full_join(NES_melted,Padj_melted[,c(3,4)], by="ID")
combined <-  combined %>% filter(X %in% selPathways )

colnames(combined)

# Checking the max and min P adjusted values to set the size scale 
min(-log10(combined$padj))
max(-log10(combined$padj))

head(combined)

# Plot figure 2A
pdf("data/GSEA/dotplot_GSEA_RhoRank_rlogViralLoad_selectedPathways.pdf")
ggplot(data = combined, aes(x = variable, y = X, color = value, size = -log10(padj))) +
  geom_point() +
  theme_bw() +
  theme(text = element_text(size = 8)) + 
  scale_color_gradient2(midpoint=0, low="blue3", mid="white", high="red3",breaks=c(-2,0,2)) +
  scale_size_continuous(range = c(1, 10), breaks = c(0.2,1,1.8))
dev.off()

# Leading edge genes for each pathway and strain are save as "LE_STRAIN_Rho_RlogViralLoad.tsv"
# LE tables has has unwanted line breaks in few cases before "(" "," and "Double quotes" which need to be fixed before proceed.
system("./Processing_LE-tables.sh")

# ----------- get LE_genes from selected pathways ------------------- #

LEtable_h1 <- read.delim("data/GSEA/LE_H1_Rho_RlogViralLoad.tsv")
LEtable_h1 <- LEtable_h1 %>% filter(X %in% selPathways ) %>% select(X, h1_v2_logeidml.Spearman_R, h1_v7_logeidml.Spearman_R)
write.table(LEtable_h1, "data/GSEA/LE_H1_Rho_RlogViralLoad_selectedPathways.tsv", sep = "\t", quote = F, row.names = F)

LEtable_h3 <- read.delim("data/GSEA/LE_H3_Rho_RlogViralLoad.tsv")
LEtable_h3 <- LEtable_h3 %>% filter(X %in% selPathways )  %>% select(X, h3_v2_logeidml.Spearman_R, h3_v7_logeidml.Spearman_R)
write.table(LEtable_h3, "data/GSEA/LE_H3_Rho_RlogViralLoad_selectedPathways.tsv", sep = "\t", quote = F, row.names = F)

LEtable_b <- read.delim("data/GSEA/LE_B_Rho_RlogViralLoad.tsv")
LEtable_b <- LEtable_b %>% filter(X %in% selPathways )  %>% select(X, b_v2_logeidml.Spearman_R, b_v7_logeidml.Spearman_R)
write.table(LEtable_b, "data/GSEA/LE_B_Rho_RlogViralLoad_selectedPathways.tsv",  sep = "\t", quote = F, row.names = F)

LEtable <- full_join(LEtable_h1, LEtable_h3, by="X")
LEtable <- full_join(LEtable, LEtable_b, by="X")
write.table(LEtable, "data/GSEA/LE_AllStrains_Rho_RlogViralLoad_selectedPathways.tsv",  sep = "\t", quote = F, row.names = F)

