library(dplyr)
library(ggplot2)
library(stringr)
options(stringsAsFactors = F)

# ------------- Reading and preparing tables ------------- #
cellsig_complete <- read.delim("data/deconvolution_InfluenzaMarkers/Celltype-Log10Signal_withoutIFNgenes.tsv") # From steep 8
NES_IFNsign <- read.delim("data/GSEA/NES_RlogRanks_baseline.tsv") # From steep 11
cellsig <- cellsig_complete
cellsig$viral <- NULL
cellsig <- as.data.frame(t(cellsig))
NES_IFNsign <- NES_IFNsign %>% filter(X =="Interferon Signaling")
row.names(NES_IFNsign) <- NES_IFNsign$X
NES_IFNsign$X <- NULL

colddata <- read.delim("data/processed/Coldata_OtherVirus_baseline_status.tsv")  # From steep 1
cellsig <- cellsig %>% select(all_of(colnames(NES_IFNsign)))
colddata <- colddata %>% filter(id %in% colnames(cellsig) )
colddata$viralstatus <- ifelse(colddata$viral == 0, "Negative", "Positive")
viral <- c(colddata$viralstatus)
names(viral) <- colddata$id 

# ------------- check if samples are the same and in the same order ------------- #
all(colnames(NES_IFNsign) %in% colnames(cellsig))  # expecting TRUE
all(colnames(NES_IFNsign) == colnames(cellsig)) # expecting TRUE
all(colnames(NES_IFNsign) %in% names(viral)) # expecting TRUE
all(colnames(NES_IFNsign) == names(viral)) # expecting TRUE


cells <- row.names(cellsig)

pdf("data/correlation/Interferon_Signaling_XYplot_CelltypeScore.pdf")  # Saving in PDF

for(x in 1:nrow(cellsig)) {
  pvalS_vec <- vector() # Create an empty vector for receiving data 
  pvalP_vec <- vector() # Create an empty vector for receiving data 
  rvalS_vec <- vector() # Create an empty vector for receiving data 
  rvalP_vec <- vector() # Create an empty vector for receiving data 
  
  for (y in 1:nrow(NES_IFNsign)) {
    
    # ------------- Correlation block  ------------- #
    
    print(paste(row.names(cellsig)[x],row.names(NES_IFNsign)[y], sep = "|"))
    pearsonTest  <- cor.test(as.numeric(cellsig[x,]), as.numeric(NES_IFNsign[y,]), method = "pearson",  use = "pairwise.complete.obs")
    spearmanTest <- cor.test(as.numeric(cellsig[x,]), as.numeric(NES_IFNsign[y,]), method = "spearman",  use = "pairwise.complete.obs")
    pPvalue  <- pearsonTest$p.value
    pRvalue  <- pearsonTest$estimate
    sPvalue  <- spearmanTest$p.value
    sRvalue  <- spearmanTest$estimate
    output <- c(pRvalue, pPvalue, sRvalue, sPvalue)
    rlog <- row.names(cellsig)[x]
    viralload <- row.names(NES_IFNsign)[y]
    dot <- cbind(as.numeric(cellsig[x,]),  as.numeric(NES_IFNsign[y,]))
    dot <-  as.data.frame(dot)
    colnames(dot) <- c("rlog", "viralload")
    row.names(dot) <- names(cellsig[x,])
    
    # ------------- Scatter Plot ussing ggplot  ------------- #
    
    if(rlog %in% cells) {
      p1 <- ggplot(dot, aes(rlog, viralload,  color=viral )) +
        geom_point() +
        geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
        xlab(rlog) +
        ylab(viralload) +
        theme_classic() +
        labs(title = paste(rlog,viralload, sep = "_VS_" ) ) + 
        scale_x_continuous(breaks = seq(round(min(dot$rlog), digits = 0), round(max(dot$rlog), digits = 0), by=1)) +
        scale_y_continuous(breaks = seq(round(min(dot$viralload - 2.5), digits = 0), round(max(dot$viralload+0.5), digits = 0), by=1)) +
        expand_limits(x=c(min(dot$rlog), max(dot$rlog)), y=c(min(dot$viralload), max(dot$viralload))) +
        annotate(geom="text", x=min(dot$rlog),y=max(dot$viralload),hjust=-1,
                 label= paste(paste0("p=", sPvalue),paste0("R²=", sRvalue),sep = "\n"), color="red")
      
      filename <- paste0(str_replace_all(paste(rlog,viralload, sep = "_VS_"), " ", "_"), ".pdf")
      print(p1) 
    }
    # ------------- Saving Plot file  ------------- #
    
    # ------------- Preparing corralation and Pvalue tables for saving  ------------- #
    
    name <- paste(rlog,viralload, sep = "|" )
    names(output)[1] <- paste(viralload, "Pearson_R", sep = "|")
    names(output)[2] <- paste(viralload, "Pearson_p", sep = "|")
    names(output)[3] <- paste(viralload, "Spearman_R", sep = "|")
    names(output)[4] <- paste(viralload, "Spearman_p", sep = "|")
    rvalP_vec <- c(rvalP_vec,output[1])
    pvalP_vec <- c(pvalP_vec,output[2])
    rvalS_vec <- c(rvalS_vec,output[3])
    pvalS_vec <- c(pvalS_vec,output[4])
    
  }  
  rvalP_vec <- as.data.frame(rvalP_vec)
  colnames(rvalP_vec) <- rlog
  pvalP_vec <- as.data.frame(pvalP_vec)
  colnames(pvalP_vec) <- rlog
  rvalS_vec <- as.data.frame(rvalS_vec)
  colnames(rvalS_vec) <- rlog
  pvalS_vec <- as.data.frame(pvalS_vec)
  colnames(pvalS_vec) <- rlog
  
  if(exists("RcortableS")) {
    RcortableS <- cbind(RcortableS, rvalS_vec )
    PcortableS <- cbind(PcortableS, pvalS_vec )
    RcortableP <- cbind(RcortableP, rvalP_vec )
    PcortableP <- cbind(PcortableP, pvalP_vec )
  }
  else {
    RcortableS <- rvalS_vec
    PcortableS <- pvalS_vec
    RcortableP <- rvalP_vec
    PcortableP <- pvalP_vec
  }
}    

dev.off() 

PcortableS <- t(PcortableS)
RcortableS <- t(RcortableS)
PcortableP <- t(PcortableP)
RcortableP <- t(RcortableP)

spearman_table <- cbind(RcortableS, PcortableS)
write.table(spearman_table, "data/correlation/spearman_table_IFNsignaling-CelltypeScore.tsv", quote = F, row.names = T, sep = "\t")
