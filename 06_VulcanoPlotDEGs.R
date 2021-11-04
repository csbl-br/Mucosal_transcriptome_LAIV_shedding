library(EnhancedVolcano)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(reshape)
library(ggrepel)
library(RColorBrewer)

DEGs <- read.delim("data/DEGs/DEG_viral_baseline.tsv")
DEGs <- DEGs %>% filter(padj != "", !is.na(padj))

threshold_DEGs <-DEGs$padj < 0.001 
DEGs$threshold <- threshold_DEGs

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

# Preparing the file to import
system('sed -i "s/c(//g" data/GSEA/LE_AllStrains_Rho_RlogViralLoad_selectedPathways.tsv')
system('sed -i "s/)//g" data/GSEA/LE_AllStrains_Rho_RlogViralLoad_selectedPathways.tsv')
system("awk -F'\t' 'BEGIN{OFS=FS} NR == 1 {for (col = 2; col <= NF; col++) {colname[col]=$col}} NR > 1 {for (col2 = 2; col2 <= NF; col2++) {print $1,colname[col2],$col2}}' data/GSEA/LE_AllStrains_Rho_RlogViralLoad_selectedPathways.tsv > data/GSEA/LE_AllStrains_Rho_RlogViralLoad_selectedPathways_melted.tsv")
system('sed -i "s/\t/|/" data/GSEA/LE_AllStrains_Rho_RlogViralLoad_selectedPathways_melted.tsv')
system('sed -i "s/, /\t/g" data/GSEA/LE_AllStrains_Rho_RlogViralLoad_selectedPathways_melted.tsv')
system('sed -i "s/ /_/g" data/GSEA/LE_AllStrains_Rho_RlogViralLoad_selectedPathways_melted.tsv')

DEGfup <- DEGs %>% filter(log2FoldChange > 0.5, padj <= 0.001 )
DEGfup <- DEGfup$geneID
DEGfdown <- DEGs %>% filter(log2FoldChange < -0.5, padj <= 0.001 )
DEGfdown <- DEGfdown$geneID


LE <- read.gmt("data/GSEA/LE_AllStrains_Rho_RlogViralLoad_selectedPathways_melted.tsv")

# Create a non redundant list of LE genes for IFN and DDX58 pathways
LEgenes <- unique(c(c(LE$`Interferon_Signaling|h1_v2_logeidml.Spearman_R`,
       LE$`Interferon_Signaling|h3_v2_logeidml.Spearman_R`,
       LE$`Interferon_Signaling|b_v2_logeidml.Spearman_R`),
       c(LE$`Interferon_Signaling|h1_v7_logeidml.Spearman_R`,
         LE$`Interferon_Signaling|h3_v7_logeidml.Spearman_R`,
         LE$`Interferon_Signaling|b_v7_logeidml.Spearman_R`),
       c(LE$`DDX58/IFIH1-mediated_induction_of_interferon-alpha/beta|b_v2_logeidml.Spearman_R`,
         LE$`DDX58/IFIH1-mediated_induction_of_interferon-alpha/beta|h1_v2_logeidml.Spearman_R`,
         LE$`DDX58/IFIH1-mediated_induction_of_interferon-alpha/beta|h3_v2_logeidml.Spearman_R`),
       c(LE$`DDX58/IFIH1-mediated_induction_of_interferon-alpha/beta|b_v7_logeidml.Spearman_R`,
         LE$`DDX58/IFIH1-mediated_induction_of_interferon-alpha/beta|h1_v7_logeidml.Spearman_R`,
         LE$`DDX58/IFIH1-mediated_induction_of_interferon-alpha/beta|h3_v7_logeidml.Spearman_R`)))

IFN_DDX58_overlapDEGs  <- intersect(DEGfup,LEgenes ) # Genes highlighted in figure 2 B

keyvals <- ifelse(
  rownames(DEGs) %in% DEGfup ,'red', 'gray30')
keyvals <- as.data.frame(keyvals)
keyvals$geneName <- row.names(DEGs)

keyvals$keyvals <- ifelse(
  keyvals$geneName %in% DEGfdown,'cornflowerblue', keyvals$keyvals)

keyvals$keyvals <- ifelse(
  keyvals$geneName %in% IFN_DDX58_overlapDEGs,'gold', keyvals$keyvals)

unique(keyvals$keyvals)

keyvals <- keyvals$keyvals
names(keyvals)[keyvals == 'red'] <- 'UPreg'
names(keyvals)[keyvals == 'cornflowerblue'] <- 'Downreg'
names(keyvals)[keyvals == 'gold'] <- 'LE'
names(keyvals)[keyvals == 'gray30'] <- 'Notsig'

p1  <- EnhancedVolcano(DEGs,
                       lab = rownames(DEGs),
                       x = 'log2FoldChange',
                       y = 'padj',
                       pCutoff = 0.001, 
                       FCcutoff = 0.5,
                       colCustom = keyvals,
                       max.overlaps = 28,
                       axisLabSize = 20,
                       colAlpha = 0.8,
                       drawConnectors = TRUE,
                       widthConnectors = 0.5,
                       typeConnectors = "closed",
                       colConnectors = 'black',
                       xlim = c(-4, 4),
                       ylim = c(0, 6))

pdf(file = "data/DEGs/VulcanPlot.pdf", 12,12)
plot(p1)
dev.off()

system("mkdir data/network")

# Figure 2C: the following tables are used as input for protein-protein interaction network (https://www.networkanalyst.ca/)
networkInput <- DEGs %>% filter(geneID  %in% IFN_DDX58_overlapDEGs) %>% select(geneID, log2FoldChange)
write.table(networkInput, "data/network/IFN_DDX58_DEGs_overlap.tsv", quote = F, sep = "\t", row.names = F)


