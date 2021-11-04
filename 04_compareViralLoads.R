library(ComplexHeatmap)
library(dplyr)
library(corrplot)
library(circlize)
library(fastDummies)
library(reshape2)
library(matrixStats)
options(stringsAsFactors = F)

# ----------- Figure 1B viral shedding heatmap (all 82 children) ------------------- #
childrens <- read.delim("data/processed/Baseline_nasal-2018_selected_children.tsv")
childrens$nameofG <- NULL
childrens <- colnames(childrens) # Get children codes 82
otherVirusOrder <- read.delim("data/processed/Coldata_OtherVirus_baseline_status.tsv")
phenodata <- read.delim("data/raw/masterdatafile.tsv")
viraload <-  phenodata %>% select(subid1, b_v2_logeidml,b_v7_logeidml, h1_v2_logeidml,h1_v7_logeidml, h3_v2_logeidml,h3_v7_logeidml) %>%
  filter(subid1 %in% childrens)
phenodata <- phenodata %>% select(subid1,b_v0_seropositive,h3_v0_seropositive, h1_v0_seropositive, v0_baseline_virus_correct ) %>%
  filter(subid1 %in% childrens)

# Include the three samples without baseline data for other viruses
otherVirusOrder <- rbind(otherVirusOrder, matrix(c("A134F", "A212B", "B220C", NA, NA,NA ), nrow = 3, ncol = 2, 
                                                 dimnames = list(c("A134F", "A212B", "B220C"), c("id","viral")) ))

othertax <- phenodata %>% select(subid1, v0_baseline_virus_correct) %>% filter(subid1 %in% otherVirusOrder$id)
othertax$v0_baseline_virus_correct <- ifelse(othertax$v0_baseline_virus_correct == "", "Negative", othertax$v0_baseline_virus_correct)

unique(othertax$v0_baseline_virus_correct)
phenodata$v0_baseline_virus_correct <- NULL

colnames(otherVirusOrder)[1] <- "subid1"
otherVirusOrder <- full_join(otherVirusOrder, othertax, by="subid1")
otherVirusOrder$v0_baseline_virus_correct <- ifelse(is.na(otherVirusOrder$viral), NA, otherVirusOrder$v0_baseline_virus_correct )
otherVirusOrder$v0_baseline_virus_correct <- ifelse(otherVirusOrder$v0_baseline_virus_correct == "Rhinovirus + adenovirus", "Rhinovirus & Adenovirus", otherVirusOrder$v0_baseline_virus_correct )


viraload <- full_join(phenodata, viraload, by = "subid1")
viraload <- full_join(otherVirusOrder, viraload, by = "subid1")
viraload <- viraload %>% filter(subid1  %in% otherVirusOrder$subid1)

viraload_b <- viraload %>% select(subid1, viral, v0_baseline_virus_correct, b_v0_seropositive, b_v2_logeidml, b_v7_logeidml )
viraload_b$mean <- rowMeans(as.matrix(viraload_b[,5:6]))
viraload_b <- viraload_b %>%  arrange(desc(mean))%>% arrange(b_v0_seropositive)
viraload_b$mean <- NULL
data_b <- viraload_b %>% select(subid1, b_v2_logeidml,b_v7_logeidml)
row.names(data_b) <- data_b$subid1
data_b$subid1 <- NULL
serob <- viraload_b %>% select(subid1, b_v0_seropositive, viral,v0_baseline_virus_correct)
row.names(serob) <- serob$subid1
serob$subid1 <- NULL

viraload_h1 <- viraload %>% select(subid1, viral,v0_baseline_virus_correct, h1_v0_seropositive, h1_v2_logeidml, h1_v7_logeidml)
viraload_h1$mean <- rowMeans(as.matrix(viraload_h1[,5:6]))
viraload_h1 <- viraload_h1 %>%  arrange(desc(mean))%>% arrange(h1_v0_seropositive)
viraload_h1$mean <- NULL
data_h1 <- viraload_h1 %>% select(subid1, h1_v2_logeidml,h1_v7_logeidml)
row.names(data_h1) <- data_h1$subid1
data_h1$subid1 <- NULL
seroh1 <- viraload_h1 %>% select(subid1, h1_v0_seropositive, viral, v0_baseline_virus_correct)
row.names(seroh1) <- seroh1$subid1
seroh1$subid1 <- NULL

viraload_h3 <- viraload %>% select(subid1, viral,v0_baseline_virus_correct, h3_v0_seropositive, h3_v2_logeidml, h3_v7_logeidml)
viraload_h3$mean <- rowMeans(as.matrix(viraload_h3[,5:6]))
viraload_h3 <- viraload_h3 %>%  arrange(desc(mean))%>% arrange(h3_v0_seropositive)
viraload_h3$mean <- NULL
data_h3 <- viraload_h3 %>% select(subid1, h3_v2_logeidml,h3_v7_logeidml)
row.names(data_h3) <- data_h3$subid1
data_h3$subid1 <- NULL
seroh3 <- viraload_h3 %>% select(subid1, h3_v0_seropositive, viral, v0_baseline_virus_correct)
row.names(seroh3) <- seroh3$subid1
seroh3$subid1 <- NULL

system("mkdir data/viralLoad")

taxa <- unique(viraload$v0_baseline_virus_correct)

col_fun <- colorRamp2(c(-1,7), c("white", "red"))

botHMannob <- HeatmapAnnotation(df = serob, which = "row",
                                col = list(viral = c("1" = "black", "0" = "cornsilk2", "NA" = "gray53"),
                                           b_v0_seropositive = c("1" = "black", "0" = "cornsilk2"),
                                           v0_baseline_virus_correct = c("Negative" =  "black", "Rhinovirus" = "#6495edff", 
                                                                         "Para 1" = "#de854cff", "Coronavirus" = "#5cb196ff",
                                                                         "Adenovirus" = "#f08080ff","Para 1 & Coronavirus" = "#9491c2ff",
                                                                         "Para 1 & Rhinovirus" = "#e9bb4dff", "Rhinovirus & Adenovirus" = "#89b75eff",
                                                                         "Rhinovirus & Coronavirus" = "#ea63a3ff", "NA" = "gray53")))
botHMannoh1 <- HeatmapAnnotation(df = seroh1, which = "row",
                                 col = list(viral = c("1" = "black", "0" = "cornsilk2", "NA" = "gray53"),
                                            h1_v0_seropositive = c("1" = "black", "0" = "cornsilk2"),
                                            v0_baseline_virus_correct = c("Negative" =  "black", "Rhinovirus" = "#6495edff", 
                                                                          "Para 1" = "#de854cff", "Coronavirus" = "#5cb196ff",
                                                                          "Adenovirus" = "#f08080ff","Para 1 & Coronavirus" = "#9491c2ff",
                                                                          "Para 1 & Rhinovirus" = "#e9bb4dff", "Rhinovirus & Adenovirus" = "#89b75eff",
                                                                          "Rhinovirus & Coronavirus" = "#ea63a3ff", "NA" = "gray53")))
botHMannoh3 <- HeatmapAnnotation(df = seroh3, which = "row",
                                 col = list(viral = c("1" = "black", "0" = "cornsilk2", "NA" = "gray53"),
                                            h3_v0_seropositive = c("1" = "black", "0" = "cornsilk2"),
                                            v0_baseline_virus_correct = c("Negative" =  "black", "Rhinovirus" = "#6495edff", 
                                                                          "Para 1" = "#de854cff", "Coronavirus" = "#5cb196ff",
                                                                          "Adenovirus" = "#f08080ff","Para 1 & Coronavirus" = "#9491c2ff",
                                                                          "Para 1 & Rhinovirus" = "#e9bb4dff", "Rhinovirus & Adenovirus" = "#89b75eff",
                                                                          "Rhinovirus & Coronavirus" = "#ea63a3ff", "NA" = "gray53")))

b <- Heatmap(as.matrix(data_b),
             name = "Heatmap\nName",
             cluster_rows = FALSE,
             cluster_columns = FALSE,
             col = col_fun,
             left_annotation = botHMannob, 
             heatmap_height = unit(25, "cm"),
             heatmap_width = unit(15, "cm"),
             row_names_gp = gpar(fontsize = 10),
             column_names_gp = gpar(fontsize = 10))

h1 <- Heatmap(as.matrix(data_h1),
              name = "Heatmap\nName",
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              col = col_fun,
              left_annotation = botHMannoh1, 
              heatmap_height = unit(25, "cm"),
              heatmap_width = unit(15, "cm"),
              row_names_gp = gpar(fontsize = 10),
              column_names_gp = gpar(fontsize = 10))

h3 <- Heatmap(as.matrix(data_h3),
              name = "Heatmap\nName",
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              col = col_fun,
              left_annotation = botHMannoh3, 
              heatmap_height = unit(25, "cm"),
              heatmap_width = unit(15, "cm"),
              row_names_gp = gpar(fontsize = 10),
              column_names_gp = gpar(fontsize = 10))

pdf(file = "data/viralLoad/viralshedding_heatmap_sep.pdf", 15, 15)
draw(b, heatmap_legend_side = "left", annotation_legend_side = "bottom")
draw(h1, heatmap_legend_side = "left", annotation_legend_side = "bottom")
draw(h3, heatmap_legend_side = "left", annotation_legend_side = "bottom")
dev.off()
