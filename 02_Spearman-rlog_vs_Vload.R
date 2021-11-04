library(dplyr)
library(ggplot2)
library(stringr)
options(stringsAsFactors = F)

# ------------- Reading and preparing tables ------------- #
masterfile <- read.delim("data/raw/masterdatafile.tsv")
shedtable <-  masterfile %>% select(subid1,
                                     h1_v2_logeidml,
                                     h1_v7_logeidml,
                                     h1_v0_seropositive,
                                     h3_v2_logeidml,
                                     h3_v7_logeidml,
                                     h3_v0_seropositive,
                                     b_v2_logeidml,
                                     b_v7_logeidml,
                                     b_v0_seropositive)
generlog_complete <- read.delim("data/processed/Baseline_nasal-2018_selected_children_rlogNorm.tsv") # From steep 1
write.table(shedtable, "data/processed/ViralLoadD0ata.tsv", sep = "\t", quote = F, row.names = F)

shedH3 <- shedtable %>% filter(h3_v0_seropositive == 0) %>% select(subid1, 
                                                                     h3_v2_logeidml,
                                                                     h3_v7_logeidml,
                                                                     h3_v0_seropositive)
# Correlation for H3 strain
generlog <- generlog_complete
row.names(shedH3) <- shedH3$subid1
shedH3 <- shedH3 %>% select(h3_v2_logeidml,h3_v7_logeidml)
generlog <- generlog %>% select(all_of(row.names(shedH3)))
shedH3 <- as.data.frame(t(shedH3))
colddata <- read.delim("data/processed/Coldata_OtherVirus_baseline_status.tsv")  # From steep 1
generlog <- generlog %>% select(all_of(colnames(shedH3)))
colddata <- colddata %>% filter(id %in% colnames(generlog) )
colddata$viralstatus <- ifelse(colddata$viral == 0, "Negative", "Positive")
viral <- c(colddata$viralstatus)
names(viral) <- colddata$id 

# ------------- check if samples are the same and in the same order ------------- #
all(colnames(shedH3) %in% colnames(generlog))  # expecting TRUE
all(colnames(shedH3) == colnames(generlog)) # expecting TRUE
all(colnames(shedH3) %in% names(viral)) # expecting TRUE
all(colnames(shedH3) == names(viral)) # expecting TRUE


# ------------- Use for selecting few genes to plot or bypass for all genes ------------- #

# Interferon signaling pathway - Reactome Gene set
genes <- c("R-HSA-913531","AAAS","ABCE1","ADAR","ARIH1","B2M","BST2","CAMK2A","CAMK2B","CAMK2D","CAMK2G","CD44","CIITA","DDX58",
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


# ------------- Spearman and Person correlation and do a scatterplots  ------------- #
system("mkdir data/correlation")

pdf("data/correlation/IFNsamples_XYplot_rlogViralLoadH3.pdf")  # Saving in PDF


for(x in 1:nrow(generlog)) {
  pvalS_vec <- vector() # Create an empty vector for receiving data 
  pvalP_vec <- vector() # Create an empty vector for receiving data 
  rvalS_vec <- vector() # Create an empty vector for receiving data 
  rvalP_vec <- vector() # Create an empty vector for receiving data 
  
  for (y in 1:nrow(shedH3)) {
    
    # ------------- Correlation block  ------------- #
    
    print(paste(row.names(generlog)[x],row.names(shedH3)[y], sep = "|"))
    pearsonTest  <- cor.test(as.numeric(generlog[x,]), as.numeric(shedH3[y,]), method = "pearson",  use = "pairwise.complete.obs")
    spearmanTest <- cor.test(as.numeric(generlog[x,]), as.numeric(shedH3[y,]), method = "spearman",  use = "pairwise.complete.obs")
    pPvalue  <- pearsonTest$p.value
    pRvalue  <- pearsonTest$estimate
    sPvalue  <- spearmanTest$p.value
    sRvalue  <- spearmanTest$estimate
    output <- c(pRvalue, pPvalue, sRvalue, sPvalue)
    rlog <- row.names(generlog)[x]
    viralload <- row.names(shedH3)[y]
    dot <- cbind(as.numeric(generlog[x,]),  as.numeric(shedH3[y,]))
    dot <-  as.data.frame(dot)
    colnames(dot) <- c("rlog", "viralload")
    row.names(dot) <- names(generlog[x,])
    
    # ------------- Scatter Plot ussing ggplot  ------------- #
    
    if(rlog %in% genes) {
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
write.table(spearman_table, "data/correlation/spearman_table_rlogViralLoadH3.tsv", quote = F, row.names = T, sep = "\t")

#####################################################################################################################
#####################################################################################################################
#####################################################################################################################


toremove <- setdiff(ls(), c("generlog_complete", "masterfile", "genes","shedtable"))
rm(list = toremove)

# Correlation for H1 strain
shedH1 <- shedtable %>% filter(h1_v0_seropositive == 0) %>% select(subid1, 
                                                                     h1_v2_logeidml,
                                                                     h1_v7_logeidml,
                                                                     h1_v0_seropositive)

colddata <- read.delim("data/processed/Coldata_OtherVirus_baseline_status.tsv")  # From steep 1
generlog <- generlog_complete
row.names(shedH1) <- shedH1$subid1
shedH1 <- shedH1 %>% filter(subid1 %in% colnames(generlog)) # two samples were romoved because dont have Viral status at base line
shedH1 <- shedH1 %>% select(h1_v2_logeidml,h1_v7_logeidml)

generlog <- generlog %>% select(all_of(row.names(shedH1)))

shedH1 <- as.data.frame(t(shedH1))
generlog <- generlog %>% select(all_of(colnames(shedH1)))
colddata <- colddata %>% filter(id %in% colnames(generlog) )
colddata$viralstatus <- ifelse(colddata$viral == 0, "Negative", "Positive")
viral <- c(colddata$viralstatus)
names(viral) <- colddata$id 

# ------------- check if samples are the same and in the same order ------------- #
all(colnames(shedH1) %in% colnames(generlog))  # expecting TRUE
all(colnames(shedH1) == colnames(generlog)) # expecting TRUE
all(colnames(shedH1) %in% names(viral)) # expecting TRUE
all(colnames(shedH1) == names(viral)) # expecting TRUE

# ------------- Spearman and Person correlation and do a scatterplots  ------------- #

pdf("data/correlation/IFNsamples_XYplot_rlogViralLoadH1.pdf")  # Saving in PDF


for(x in 1:nrow(generlog)) {
  pvalS_vec <- vector() # Create an empty vector for receiving data 
  pvalP_vec <- vector() # Create an empty vector for receiving data 
  rvalS_vec <- vector() # Create an empty vector for receiving data 
  rvalP_vec <- vector() # Create an empty vector for receiving data 
  
  for (y in 1:nrow(shedH1)) {
    
    # ------------- Correlation block  ------------- #
    
    print(paste(row.names(generlog)[x],row.names(shedH1)[y], sep = "|"))
    pearsonTest  <- cor.test(as.numeric(generlog[x,]), as.numeric(shedH1[y,]), method = "pearson",  use = "pairwise.complete.obs")
    spearmanTest <- cor.test(as.numeric(generlog[x,]), as.numeric(shedH1[y,]), method = "spearman",  use = "pairwise.complete.obs")
    pPvalue  <- pearsonTest$p.value
    pRvalue  <- pearsonTest$estimate
    sPvalue  <- spearmanTest$p.value
    sRvalue  <- spearmanTest$estimate
    output <- c(pRvalue, pPvalue, sRvalue, sPvalue)
    rlog <- row.names(generlog)[x]
    viralload <- row.names(shedH1)[y]
    dot <- cbind(as.numeric(generlog[x,]),  as.numeric(shedH1[y,]))
    dot <-  as.data.frame(dot)
    colnames(dot) <- c("rlog", "viralload")
    row.names(dot) <- names(generlog[x,])
    
    # ------------- Scatter Plot ussing ggplot  ------------- #
    
    if(rlog %in% genes) {
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
write.table(spearman_table, "data/correlation/spearman_table_rlogViralLoadH1.tsv", quote = F, row.names = T, sep = "\t")

#####################################################################################################################
#####################################################################################################################
#####################################################################################################################

toremove <- setdiff(ls(), c("generlog_complete", "masterfile", "genes","shedtable"))
rm(list = toremove)

# Correlation for B strain

shedB <- shedtable %>% filter(b_v0_seropositive == 0) %>% select(subid1, 
                                                                     b_v2_logeidml,
                                                                     b_v7_logeidml,
                                                                     b_v0_seropositive)

generlog <- generlog_complete
row.names(shedB) <- shedB$subid1
shedB <- shedB %>% filter(subid1 %in% colnames(generlog))
shedB <- shedB %>% select(b_v2_logeidml,b_v7_logeidml)
generlog <- generlog %>% select(all_of(row.names(shedB)))
shedB <- as.data.frame(t(shedB))
colddata <- read.delim("data/processed/Coldata_OtherVirus_baseline_status.tsv")  # From steep 1
generlog <- generlog %>% select(all_of(colnames(shedB)))
colddata <- colddata %>% filter(id %in% colnames(generlog) )
colddata$viralstatus <- ifelse(colddata$viral == 0, "Negative", "Positive")
viral <- c(colddata$viralstatus)
names(viral) <- colddata$id 

# ------------- check if samples are the same and in the same order ------------- #
all(colnames(shedB) %in% colnames(generlog))  # expecting TRUE
all(colnames(shedB) == colnames(generlog)) # expecting TRUE
all(colnames(shedB) %in% names(viral)) # expecting TRUE
all(colnames(shedB) == names(viral)) # expecting TRUE

# ------------- Spearman and Person correlation and do a scatterplots  ------------- #

pdf("data/correlation/IFNsamples_XYplot_rlogViralLoadB.pdf")  # Saving in PDF


for(x in 1:nrow(generlog)) {
  pvalS_vec <- vector() # Create an empty vector for receiving data 
  pvalP_vec <- vector() # Create an empty vector for receiving data 
  rvalS_vec <- vector() # Create an empty vector for receiving data 
  rvalP_vec <- vector() # Create an empty vector for receiving data 
  
  for (y in 1:nrow(shedB)) {
    
    # ------------- Correlation block  ------------- #
    
    print(paste(row.names(generlog)[x],row.names(shedB)[y], sep = "|"))
    pearsonTest  <- cor.test(as.numeric(generlog[x,]), as.numeric(shedB[y,]), method = "pearson",  use = "pairwise.complete.obs")
    spearmanTest <- cor.test(as.numeric(generlog[x,]), as.numeric(shedB[y,]), method = "spearman",  use = "pairwise.complete.obs")
    pPvalue  <- pearsonTest$p.value
    pRvalue  <- pearsonTest$estimate
    sPvalue  <- spearmanTest$p.value
    sRvalue  <- spearmanTest$estimate
    output <- c(pRvalue, pPvalue, sRvalue, sPvalue)
    rlog <- row.names(generlog)[x]
    viralload <- row.names(shedB)[y]
    dot <- cbind(as.numeric(generlog[x,]),  as.numeric(shedB[y,]))
    dot <-  as.data.frame(dot)
    colnames(dot) <- c("rlog", "viralload")
    row.names(dot) <- names(generlog[x,])
    
    # ------------- Scatter Plot ussing ggplot  ------------- #
    
    if(rlog %in% genes) {
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
write.table(spearman_table, "data/correlation/spearman_table_rlogViralLoadB.tsv", quote = F, row.names = T, sep = "\t")




