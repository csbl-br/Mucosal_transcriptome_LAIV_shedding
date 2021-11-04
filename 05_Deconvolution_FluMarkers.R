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
library(fmsb)

font_import(pattern="[A/a]rial", prompt=FALSE)

vsd <- read.delim("data/processed/Baseline_nasal-2018_selected_children_VSDNorm.tsv")
vsd$geneID <- row.names(vsd)
flumarkers <- read.delim("data/raw/influenza_markers.tsv")
otherVirusOrder <- read.delim("data/processed/Coldata_OtherVirus_baseline_status.tsv")
v0_response <-  otherVirusOrder
IFNgenes <- c("R-HSA-913531","AAAS","ABCE1","ADAR","ARIH1","B2M","BST2","CAMK2A","CAMK2B","CAMK2D","CAMK2G","CD44","CIITA","DDX58",
              "EGR1","EIF2AK2","EIF4A1","EIF4A2","EIF4A3","EIF4E","EIF4E2","EIF4E3","EIF4G1","EIF4G2","EIF4G3","FCGR1A","FCGR1B",
              "FLNA","FLNB","GBP1","GBP2","GBP3","GBP4","GBP5","GBP6","GBP7","HERC5","HLA-A","HLA-B","HLA-C","HLA-DPA1","HLA-DPB1",
              "HLA-DQA1","HLA-DQA2","HLA-DQB1","HLA-DQB2","HLA-DRA","HLA-DRB1","HLA-DRB3","HLA-DRB4","HLA-DRB5","HLA-E","HLA-F","HLA-G",
              "HLA-H","ICAM1","IFI27","IFI30","IFI35","IFI6","IFIT1","IFIT2","IFIT3","IFITM1","IFITM2","IFITM3","IFNA1","IFNA10","IFNA14",
              "IFNA16","IFNA17","IFNA2","IFNA21","IFNA4","IFNA5","IFNA6","IFNA7","IFNA8","IFNAR1","IFNAR2","IFNB1","IFNG","IFNGR1","IFNGR2",
              "IP6K2","IRF1","IRF2","IRF3","IRF4","IRF5","IRF6","IRF7","IRF8","IRF9","ISG15","ISG20","JAK1","JAK2","KPNA1","KPNA2","KPNA3",
              "KPNA4","KPNA5","KPNA7","KPNB1","MAPK3","MID1","MT2A","MX1","MX2","NCAM1","NDC1","NEDD4","NS","NUP107","NUP133","NUP153","NUP155",
              "NUP160","NUP188","NUP205","NUP210","NUP214","NUP35","NUP37","NUP43","NUP50","NUP54","NUP58","NUP62","NUP85","NUP88","NUP93","NUP98",
              "NUPL2","OAS1","OAS2","OAS3","OASL","PDE12","PIAS1","PIN1","PLCG1","PML","POM121","POM121C","PPM1B","PRKCD","PSMB8","PTAFR","PTPN1",
              "PTPN11","PTPN2","PTPN6","RAE1","RANBP2","RNASEL","RPS27A","RSAD2","SAMHD1","SEC13","SEH1L","SOCS1","SOCS3","SP100","STAT1","STAT2",
              "SUMO1","TPR","TRIM10","TRIM14","TRIM17","TRIM2","TRIM21","TRIM22","TRIM25","TRIM26","TRIM29","TRIM3","TRIM31","TRIM34","TRIM35",
              "TRIM38","TRIM45","TRIM46","TRIM48","TRIM5","TRIM6","TRIM62","TRIM68","TRIM8","TYK2","UBA52","UBA7","UBB","UBC","UBE2E1","UBE2L6",
              "UBE2N","USP18","USP41","VCAM1","VP3","XAF1","gag")

markerslist <- list()
markerslist_noIF <- list()
cellIFNgenes <- list()
for(celltype in 1:length(colnames(flumarkers))) {
  markerslist[[celltype]] <- c(intersect(flumarkers[,celltype], row.names(vsd)))
  cellIFNgenes [[celltype]] <-intersect(IFNgenes, markerslist[[celltype]])
  names(cellIFNgenes)[celltype] <- colnames(flumarkers)[celltype]
  markerslist_noIF[[celltype]] <- setdiff(markerslist[[celltype]], IFNgenes)
  names(markerslist)[celltype] <- colnames(flumarkers)[celltype]
  names(markerslist_noIF)[celltype] <- colnames(flumarkers)[celltype]
}

objectstokeep <- ls()
v0_response <-  otherVirusOrder
row.names(v0_response) <- v0_response$id
v0_response$id <- NULL
otherVirusOrder <- otherVirusOrder$id
signaltable <- v0_response
signaltablenotlog <- v0_response
system("mkdir data/deconvolution_InfluenzaMarkers")

# Creating the tables with cell type signal (log10 and raw), but removing IFN genes
for(celltype in 1:length(markerslist_noIF)) {
  genes <- unlist(markerslist_noIF[celltype])
  expressed_markers <- vsd %>% filter(geneID %in% all_of(genes))
  expressed_markers$geneID <- NULL
  expressed_markers <- expressed_markers %>% select( all_of(otherVirusOrder))
  CellType_Signalnotlog <- colSums(expressed_markers)
  signaltablenotlog <- cbind(signaltablenotlog, CellType_Signalnotlog)
  names(signaltablenotlog)[(celltype + 1)] <- names(markerslist_noIF)[celltype]
  CellType_Signal <- log(colSums(expressed_markers),10)
  signaltable <- cbind(signaltable, CellType_Signal)
  names(signaltable)[(celltype + 1)] <- names(markerslist_noIF)[celltype]
  v0_response <- cbind(v0_response, CellType_Signal)
  expressed_markers <- as.matrix(expressed_markers)
  expressed_markers <- t(expressed_markers)
}

write.table(signaltable, "data/deconvolution_InfluenzaMarkers/Celltype-Log10Signal_withoutIFNgenes.tsv", quote = F, sep = "\t")
write.table(signaltablenotlog, "data/deconvolution_InfluenzaMarkers/Celltype-RawSignal_withoutIFNgenes.tsv", quote = F, sep = "\t")

# Creating the tables with cell type signal (log10 and raw) using all markes (including IFN genes)
rm(list= setdiff(ls(), objectstokeep))
otherVirusOrder <- read.delim("data/processed/Coldata_OtherVirus_baseline_status.tsv")
v0_response <-  otherVirusOrder
row.names(v0_response) <- v0_response$id
v0_response$id <- NULL
otherVirusOrder <- otherVirusOrder$id
signaltable <- v0_response
signaltablenotlog <- v0_response
for(celltype in 1:length(markerslist)) {
  genes <- unlist(markerslist[celltype])
  expressed_markers <- vsd %>% filter(geneID %in% all_of(genes))
  expressed_markers$geneID <- NULL
  expressed_markers <- expressed_markers %>% select( all_of(otherVirusOrder))
  CellType_Signalnotlog <- colSums(expressed_markers)
  signaltablenotlog <- cbind(signaltablenotlog, CellType_Signalnotlog)
  names(signaltablenotlog)[(celltype + 1)] <- names(markerslist)[celltype]
  CellType_Signal <- log(colSums(expressed_markers),10)
  signaltable <- cbind(signaltable, CellType_Signal)
  names(signaltable)[(celltype + 1)] <- names(markerslist)[celltype]
  v0_response <- cbind(v0_response, CellType_Signal)
  expressed_markers <- as.matrix(expressed_markers)
  expressed_markers <- t(expressed_markers)
}

write.table(signaltable, "data/deconvolution_InfluenzaMarkers/Celltype-Log10Signal_allMarkers.tsv", quote = F, sep = "\t")
write.table(signaltablenotlog, "data/deconvolution_InfluenzaMarkers/Celltype-RawSignal_allMarkers.tsv", quote = F, sep = "\t")
rm(list = ls())

# Ploting Barplots and BoxPlots with INF Genes
signaltable <- read.delim("data/deconvolution_InfluenzaMarkers/Celltype-Log10Signal_allMarkers.tsv")
signaltable <- cbind(row.names(signaltable), signaltable)
colnames(signaltable)[1] <- "IDs" 
signaltable_melted <- reshape2::melt(signaltable, id.vars=c("IDs","viral"), measure.vars=3:ncol(signaltable))
signaltable_melted$variable2 <- paste(signaltable_melted$viral, signaltable_melted$variable, sep = "_")
signaltable_melted <- signaltable_melted %>% arrange(variable, viral)

signaltablenotlog <- read.delim("data/deconvolution_InfluenzaMarkers/Celltype-RawSignal_allMarkers.tsv")
signaltablenotlog <- cbind(row.names(signaltablenotlog), signaltablenotlog)
colnames(signaltablenotlog)[1] <- "IDs" 
signaltablenotlog_melted <- reshape2::melt(signaltablenotlog, id.vars=c("IDs","viral"), measure.vars=3:ncol(signaltablenotlog))
signaltablenotlog_melted$variable2 <- paste(signaltablenotlog_melted$viral, signaltablenotlog_melted$variable, sep = "_")
signaltablenotlog_melted <- signaltablenotlog_melted %>% arrange(variable, viral)

fctable <- signaltablenotlog_melted %>%
  group_by(variable2) %>%
  summarise(median = median(value), iqr = IQR(value), mad = mad(value) )

g0 <- grep("0_", fctable$variable2, value = T )
g1 <- grep("1_", fctable$variable2, value = T )
fctable$viral <- ifelse(fctable$variable2 %in% g0, "negative", "positive")
fctable$celltype <- sapply(strsplit(fctable$variable2,"_"), `[`, 2)
celllist <- unique(fctable$celltype)

fctable_1 <- fctable %>% filter(viral == "positive") %>% select(median, celltype, viral)
fctable_1 <- as.data.frame(fctable_1)
row.names(fctable_1) <- fctable_1$celltype
fctable_1 <- fctable_1[celllist,]
fctable_0 <- fctable %>% filter(viral == "negative")  %>% select(median, celltype, viral)
fctable_0 <- as.data.frame(fctable_0)
row.names(fctable_0) <- fctable_0$celltype
fctable_0 <- fctable_0[celllist,]

# Checking if rows are the same and in the same order
all(rownames(fctable_1) %in% rownames(fctable_0))
all(rownames(fctable_1) == rownames(fctable_0))

celltype_FC <- log(fctable_1$median, 2) -  log(fctable_0$median,2)
names(celltype_FC) <- celllist
celltype_FC <- as.data.frame(celltype_FC)
celltype_FC$cell <- row.names(celltype_FC)
colnames(celltype_FC) <- c("value", "cell")

celltype_FC$cell <- factor(celltype_FC$cell, levels = c("NEU","B","CD4T", "EOS", "NK",
                                                        "cDC","MAC", "gdT", "pDC", "CD8T",
                                                        "proliferating","MAST","CEP", "GOB", "MHCIIEPI",
                                                        "BasalEPI", "Squamous"))

write.table(celltype_FC, "data/deconvolution_InfluenzaMarkers/celltype_Log2Fold_change_withIFNgenes.tsv", sep = "\t", quote = F, row.names = F)


# Raw plot --> Figure 3A
pdf("data/deconvolution_InfluenzaMarkers/Barplot_celltype_Log2Fold_change_withIFNgenes.pdf")
g1 <- ggplot(data = celltype_FC,
             aes(x = cell, y = value))+
  geom_bar(stat = "identity") +
  ylim(c(-0.2, 0.2))
g1
dev.off()

unique(signaltable_melted$variable)

my_comparisons <- list( c("0_B_Markers", "1_B_Markers"),
                        c("0_BasalEPI_Markers", "1_BasalEPI_Markers"),
                        c("0_CD4T_Markers", "1_CD4T_Markers"),
                        c("0_CD8T_Markers", "1_CD8T_Markers"),
                        c("0_cDC_Markers", "1_cDC_Markers"),
                        c("0_CEP_Markers", "1_CEP_Markers"),
                        c("0_EOS_Markers", "1_EOS_Markers"),
                        c("0_gdT_Markers", "1_gdT_Markers"),
                        c("0_GOB_Markers", "1_GOB_Markers"),
                        c("0_MAC_Markers", "1_MAC_Markers"),
                        c("0_MAST_Markers", "1_MAST_Markers"),
                        c("0_MHCIIEPI_Markers", "1_MHCIIEPI_Markers"),
                        c("0_NEU_Markers", "1_NEU_Markers"),
                        c("0_NK_Markers", "1_NK_Markers"),
                        c("0_pDC_Markers", "1_pDC_Markers"),
                        c("0_proliferating_CD8T_Markers", "1_proliferating_CD8T_Markers"),
                        c("0_Squamous_Markers", "1_Squamous_Markers"))

rm(tableWilcox) # making sure there no previous data in the data frame
for(celltype in 1:length(my_comparisons)) {
  print(my_comparisons[[celltype]])
  grouptable <- signaltable_melted %>% filter(variable2 %in% my_comparisons[[celltype]]) %>%
    select(variable2, value)
  res <- wilcox.test(value ~ variable2, data = grouptable,
                     exact = FALSE)
  if(exists("tableWilcox")) {
    values <- c(my_comparisons[[celltype]][1], my_comparisons[[celltype]][2],res$p.value)
    tableWilcox <- rbind(tableWilcox, c(my_comparisons[[celltype]][1], my_comparisons[[celltype]][2],res$p.value)) }
  if(!exists("tableWilcox")) {
    tableWilcox <- data.frame( 
      group1 = c(my_comparisons[[celltype]][1]),
      group2 = c(my_comparisons[[celltype]][2]),
      pvalue = c(res$p.value), stringsAsFactors = FALSE)
  }
}

write.table(tableWilcox, "data/deconvolution_InfluenzaMarkers/Celltype_wilcoxTest_withIFNgenes.tsv", quote = F, sep = "\t", row.names = F)

signaltable_melted$viral <- ifelse(signaltable_melted$viral == 0, "negative", "positive")

color <- c("grey50", "grey100")

# Raw grouped boxplot - supplementary figure
box1 <- ggplot(signaltable_melted, aes(x=variable2, y=value, fill=viral, colo)) +
  geom_violin(width=1, color="black", alpha=0.5) +
  geom_boxplot(width=0.1, color="black", alpha=0.5) +
  scale_fill_manual(values=c("grey50", "grey10")) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  facet_wrap(~ variable, scale = "free" ) +
  ggtitle("Influenza - Cell type markers | pertubation score per sample ") +
  xlab("")

box1

pdf(file = "data/deconvolution_InfluenzaMarkers/Boxplot_celltypes_WithIFNgenes.pdf", 10, 10)
box1
dev.off()

rm(list = ls())

# Ploting Barplots and BoxPlots without INF Genes
signaltable <- read.delim("data/deconvolution_InfluenzaMarkers/Celltype-Log10Signal_withoutIFNgenes.tsv")
signaltable <- cbind(row.names(signaltable), signaltable)
colnames(signaltable)[1] <- "IDs" 
signaltable_melted <- reshape2::melt(signaltable, id.vars=c("IDs","viral"), measure.vars=3:ncol(signaltable))
signaltable_melted$variable2 <- paste(signaltable_melted$viral, signaltable_melted$variable, sep = "_")
signaltable_melted <- signaltable_melted %>% arrange(variable, viral)

signaltablenotlog <- read.delim("data/deconvolution_InfluenzaMarkers/Celltype-RawSignal_withoutIFNgenes.tsv")
signaltablenotlog <- cbind(row.names(signaltablenotlog), signaltablenotlog)
colnames(signaltablenotlog)[1] <- "IDs" 
signaltablenotlog_melted <- reshape2::melt(signaltablenotlog, id.vars=c("IDs","viral"), measure.vars=3:ncol(signaltablenotlog))
signaltablenotlog_melted$variable2 <- paste(signaltablenotlog_melted$viral, signaltablenotlog_melted$variable, sep = "_")
signaltablenotlog_melted <- signaltablenotlog_melted %>% arrange(variable, viral)

fctable <- signaltablenotlog_melted %>%
  group_by(variable2) %>%
  summarise(median = median(value), iqr = IQR(value), mad = mad(value) )

g0 <- grep("0_", fctable$variable2, value = T )
g1 <- grep("1_", fctable$variable2, value = T )
fctable$viral <- ifelse(fctable$variable2 %in% g0, "negative", "positive")
fctable$celltype <- sapply(strsplit(fctable$variable2,"_"), `[`, 2)
celllist <- unique(fctable$celltype)

fctable_1 <- fctable %>% filter(viral == "positive") %>% select(median, celltype, viral)
fctable_1 <- as.data.frame(fctable_1)
row.names(fctable_1) <- fctable_1$celltype
fctable_1 <- fctable_1[celllist,]
fctable_0 <- fctable %>% filter(viral == "negative")  %>% select(median, celltype, viral)
fctable_0 <- as.data.frame(fctable_0)
row.names(fctable_0) <- fctable_0$celltype
fctable_0 <- fctable_0[celllist,]

# Checking if rows are the same and in the same order
all(rownames(fctable_1) %in% rownames(fctable_0))
all(rownames(fctable_1) == rownames(fctable_0))

celltype_FC <- log(fctable_1$median, 2) -  log(fctable_0$median,2) # ratio of log values
names(celltype_FC) <- celllist
celltype_FC <- as.data.frame(celltype_FC)
celltype_FC$cell <- row.names(celltype_FC)
colnames(celltype_FC) <- c("value", "cell")

celltype_FC$cell <- factor(celltype_FC$cell, levels = c("NEU","B","CD4T", "EOS", "NK",
                                                        "cDC","MAC", "gdT", "pDC", "CD8T",
                                                        "proliferating","MAST","CEP", "GOB", "MHCIIEPI",
                                                        "BasalEPI", "Squamous"))

write.table(celltype_FC, "data/deconvolution_InfluenzaMarkers/celltype_Log2Fold_change_withoutIFNgenes.tsv", sep = "\t", quote = F, row.names = F)


pdf("data/deconvolution_InfluenzaMarkers/Barplot_celltype_Log2Fold_change_withoutIFNgenes.pdf")
g1 <- ggplot(data = celltype_FC,
             aes(x = cell, y = value))+
  geom_bar(stat = "identity") +
  ylim(c(-0.2, 0.2))
g1
dev.off()

my_comparisons <- list( c("0_B_Markers", "1_B_Markers"),
                        c("0_BasalEPI_Markers", "1_BasalEPI_Markers"),
                        c("0_CD4T_Markers", "1_CD4T_Markers"),
                        c("0_CD8T_Markers", "1_CD8T_Markers"),
                        c("0_cDC_Markers", "1_cDC_Markers"),
                        c("0_CEP_Markers", "1_CEP_Markers"),
                        c("0_EOS_Markers", "1_EOS_Markers"),
                        c("0_gdT_Markers", "1_gdT_Markers"),
                        c("0_GOB_Markers", "1_GOB_Markers"),
                        c("0_MAC_Markers", "1_MAC_Markers"),
                        c("0_MAST_Markers", "1_MAST_Markers"),
                        c("0_MHCIIEPI_Markers", "1_MHCIIEPI_Markers"),
                        c("0_NEU_Markers", "1_NEU_Markers"),
                        c("0_NK_Markers", "1_NK_Markers"),
                        c("0_pDC_Markers", "1_pDC_Markers"),
                        c("0_proliferating_CD8T_Markers", "1_proliferating_CD8T_Markers"),
                        c("0_Squamous_Markers", "1_Squamous_Markers"))

rm(tableWilcox)
for(celltype in 1:length(my_comparisons)) {
  print(my_comparisons[[celltype]])
  grouptable <- signaltable_melted %>% filter(variable2 %in% my_comparisons[[celltype]]) %>%
    select(variable2, value)
  res <- wilcox.test(value ~ variable2, data = grouptable,
                     exact = FALSE)
  if(exists("tableWilcox")) {
    values <- c(my_comparisons[[celltype]][1], my_comparisons[[celltype]][2],res$p.value)
    tableWilcox <- rbind(tableWilcox, c(my_comparisons[[celltype]][1], my_comparisons[[celltype]][2],res$p.value)) }
  if(!exists("tableWilcox")) {
    tableWilcox <- data.frame( 
      group1 = c(my_comparisons[[celltype]][1]),
      group2 = c(my_comparisons[[celltype]][2]),
      pvalue = c(res$p.value), stringsAsFactors = FALSE)
  }
}

write.table(tableWilcox, "data/deconvolution_InfluenzaMarkers/Celltype_wilcoxTest_withoutIFNgenes.tsv", quote = F, sep = "\t", row.names = F)


signaltable_melted$viral <- ifelse(signaltable_melted$viral == 0, "negative", "positive")

color <- c("grey50", "grey100")

# grouped boxplot
box1 <- ggplot(signaltable_melted, aes(x=variable2, y=value, fill=viral, colo)) +
  geom_violin(width=1, color="black", alpha=0.5) +
  geom_boxplot(width=0.1, color="black", alpha=0.5) +
  scale_fill_manual(values=c("grey50", "grey10")) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  facet_wrap(~ variable, scale = "free" ) +
  ggtitle("Influenza (Jose - Nature) - Cell type markers | pertubation score per sample ") +
  xlab("")

box1

pdf(file = "data/deconvolution_InfluenzaMarkers/Boxplot_celltypes_WithoutIFNgenes.pdf", 10, 10)
box1
dev.off()
