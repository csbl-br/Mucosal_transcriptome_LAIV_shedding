options(stringsAsFactors = F)
library(dplyr)
library(corrplot)
library(DESeq2)
library(BiocParallel)
library(fastDummies)

#--------------- Read and prepare raw data  ------------- #
master_data <- read.delim("data/raw/masterdatafile.tsv")
count_2018 <- read.delim("data/raw/GSE159884_RawCounts.txt")
cildren_list <- read.delim("data/raw/children.list", header = F)
cildren_list <- cildren_list$V1 # List of selected children: Children that are naive for at least one flu strain

# Create a table with row mean for collapsing isoforms
mean.table <- count_2018
mean.table$Mean <- rowMeans(count_2018[,4:length(count_2018)])
table2collapse <- select(mean.table, c("ID", "nameofG", "Mean"))

# Collapsing using the isoform with max row mean
nodup <- table2collapse  %>%  group_by(nameofG) %>% 
  arrange(desc(Mean)) %>% 
  filter(!duplicated(nameofG))
nodupIDs <- as.vector(nodup$ID) 
collapDF  <- subset(count_2018 , ID %in% nodupIDs)
collapDF$ID <- NULL
baseline_2018 <- collapDF
baseline_2018$geneLength <- NULL
rm(collapDF,table2collapse,mean.table, nodupIDs, nodup)

# Match children list with count table
all(cildren_list %in% colnames(baseline_2018)) # Expecting TRUE
all(cildren_list  == colnames(baseline_2018)[2:ncol(baseline_2018)]) # Expecting TRUE

system("mkdir data/processed")
write.table(baseline_2018, "data/processed/Baseline_nasal-2018_selected_children.tsv", quote = F, sep = "\t", row.names = F)

# ------------- Prepare coldata and cts for DEseq2 ------------- #

#Selecting children with information about other virus, 2018 and has nasal RNAseq data
phenodata <- master_data %>% 
  select(subid1, group, year, v0_resp_virus_positive_correct) %>% 
  filter(subid1  %in% cildren_list) %>%
  filter(v0_resp_virus_positive_correct != "NA") 
coldata <- phenodata %>% select(subid1, v0_resp_virus_positive_correct) 
colnames(coldata) <-  c("id", "viral")

# From 82 selected children, 79 has viral status in the baseline 
write.table(coldata, "data/processed/Coldata_OtherVirus_baseline_status.tsv", row.names = F, quote = F, sep = "\t")

row.names(coldata) <- coldata$id
cts <- baseline_2018 %>% select(nameofG, all_of(coldata$id))
row.names(cts) <- cts$nameofG
cts$nameofG <- NULL
coldata$id <- factor(coldata$id, levels = unique(coldata$id))
coldata$viral <- factor(coldata$viral, levels=c(0, 1))

# Checking cts samples and coldata samples
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))

# -------------DEG using DEseq2 ------------- #
outprefix <- "data/DEGs/DEG"
register(MulticoreParam(4))
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ viral)
keep <- rowSums(counts(dds)) >= 30 # Filter low counts

dds <- dds[keep,]
dds <- DESeq(dds, parallel = TRUE)
res <- results(dds, parallel = TRUE)

rld <- rlog(dds, blind=FALSE)
vsd <- vst(dds, blind=FALSE)

colData(dds)
plotPCA(vsd, intgroup=c("viral"))
mat <- assay(vsd)

# ------------- Savig normalized tables ------------- #

# VSD will be used to deconvolution as it takes in to accound the design formula 
write.table(file = "data/processed/Baseline_nasal-2018_selected_children_VSDNorm.tsv", assay(vsd), 
            sep = "\t", quote = F)

# rlog will be used to correlations against viral load as it takes it do not use the design formula 
write.table(file = "data/processed/Baseline_nasal-2018_selected_children_rlogNorm.tsv", assay(rld), 
            sep = "\t", quote = F)

# ------------- Get DEGs using different cutoffs ------------- #
resultsNames(dds)
# [1] "Intercept"    "viral_1_vs_0"  ----- > viral positive versus viral negative (control)

resOrdered <- res[order(res$log2FoldChange),] 

# Preparing the complete table and tables with three different thresholds 
resSig <- subset(resOrdered, padj < 0.1)
resSig_001 <- subset(resOrdered, padj < 0.01)
resSig_0001 <- subset(resOrdered, padj < 0.001)

outfile_padj <- paste0(outprefix, "_viral_baseline_padj_01.tsv")
outfile_padj001 <- paste0(outprefix, "_viral_baseline_padj_001.tsv")
outfile_padj0001 <- paste0(outprefix, "_viral_baseline_padj_0001.tsv")
outfile <- paste0(outprefix, "_viral_baseline.tsv")

DEGs_padj0.1 <- as.data.frame(resSig)
DEGs_padj0.01 <- as.data.frame(resSig_001)
DEGs_padj0.001 <- as.data.frame(resSig_0001)
DEGs <- as.data.frame(resOrdered)

DEGs$geneID <- rownames(DEGs)
DEGs_padj0.1$geneID <- rownames(DEGs_padj0.1)
DEGs_padj0.01$geneID <- rownames(DEGs_padj0.01)
DEGs_padj0.001$geneID <- rownames(DEGs_padj0.001)


DEGs_padj0.1 <- select(DEGs_padj0.1, c("geneID","log2FoldChange", 
                                       "pvalue", "padj", "baseMean", "lfcSE", "stat" ))
DEGs_padj0.01 <- select(DEGs_padj0.01, c("geneID","log2FoldChange", 
                                         "pvalue", "padj", "baseMean", "lfcSE", "stat" ))
DEGs_padj0.001 <- select(DEGs_padj0.001, c("geneID","log2FoldChange", 
                                           "pvalue", "padj", "baseMean", "lfcSE", "stat" ))

DEGs <- select(DEGs, c("geneID","log2FoldChange", "pvalue", "padj", 
                       "baseMean", "lfcSE", "stat" ))

system("mkdir data/DEGs")


write.table(DEGs_padj0.1,file = outfile_padj, sep = '\t', 
            row.names = TRUE, col.names = TRUE)
write.table(DEGs_padj0.01,file = outfile_padj001, sep = '\t', 
            row.names = TRUE, col.names = TRUE)
write.table(DEGs_padj0.001,file = outfile_padj0001, sep = '\t', 
            row.names = TRUE, col.names = TRUE)
write.table(DEGs,file = outfile, 
            sep = '\t', row.names = TRUE, col.names = TRUE)



