########################################################
#phyloseq data, create curated ASV table, Figure 1, 2, SRS curves

#Import libraries
library(phyloseq)
library(vegan)
library(ranacapa)
library(ggpubr)
library(dplyr)
library(SRS)
library(ggplot2)
library(forcats)
library(stats)
library(ade4)
library(ggord)
library(adespatial)
library(sf)
library(readxl)
library(imputeTS)
library(betapart)
library(randomForest)
library(caret)
library(rfPermute)
library(dttr2)
library(agricolae)
library("ranger")
library("ggrepel")
library("PMCMRplus")
library("dunn.test")
library("steprf")
library(VennDiagram)
library(ggtext)
#to ensure reproducibility
set.seed(123)
################################################################
######################################################################

#Add ASV table
asv <- read.csv("~/Euk_18S_data/ASVs.csv", header=1, row.names=1)

#Add taxa table
tax_a_meso <- read.csv("~/Euk_18S_data/taxa_a_meso_filtered_1205.csv", header=1, row.names=1)
tax_a_macro <- read.csv("~/Euk_18S_data/taxa_a_macro_filtered_1205.csv", header=1, row.names=1)
tax_a_micro <- read.csv("~/Euk_18S_data/taxa_a_micro_filtered_1205.csv", header=1, row.names=1)
tax_a_arthro <- tax_a_meso %>% filter(phylum == "Arthropoda")
tax_a_arthro <- tax_a_arthro[,-9]
tax_a_anneli <- tax_a_meso %>% filter(phylum == "Annelida")
tax_a_anneli <- tax_a_anneli[,-9]
tax_a_anneli1 <- tax_a_macro %>% filter(phylum == "Annelida")
tax_a_anneli1 <- tax_a_anneli1[,-9]
tax_a_anneli <- rbind(tax_a_anneli,tax_a_anneli1)
tax_a_nema <- tax_a_micro %>% filter(phylum == "Nematozoa")
tax_a_nema <- tax_a_nema[,-9]
tax_a_tardi <- tax_a_micro %>% filter(phylum == "Tardigrada")
tax_a_tardi <-tax_a_tardi[,-9]
tax_a_roti <- tax_a_micro %>% filter(phylum == "Rotifera")
tax_a_roti <-tax_a_roti[,-9]
tax_a_mollu <- tax_a_macro %>% filter(phylum == "Mollusca")
tax_a_mollu <-tax_a_mollu[,-9]
tax_a_platy <- tax_a_macro %>% filter(phylum == "Platyhelminthes")
tax_a_platy <-tax_a_platy[,-9]
tax_a_platy <- tax_a_macro %>% filter(phylum == "Gastrotricha")
tax_a_platy <-tax_a_platy[,-9]
tax_f <- read.csv("~/Euk_18S_data/taxa_f_filtered_1205.csv", header=1, row.names=1)
tax_p <- read.csv("~/Euk_18S_data/taxa_p_filtered_1205.csv", header=1, row.names=1)

#load sample_data table
sample_data <- read.csv(file="~/Euk_18S_data/Sample_data_18052023.csv",header = 1, row.names=1)
####################################
#Merge taxonomy table + ASVs
tax_asv_f <- merge(tax_f, asv, by=0, all=F)
tax_asv_p <- merge(tax_p, asv, by=0, all=F)
tax_asv_a_arthro <- merge(tax_a_arthro, asv, by=0, all=F)
tax_asv_a_anneli <- merge(tax_a_anneli, asv, by=0, all=F)
tax_asv_a_nema <- merge(tax_a_nema, asv, by=0, all=F)
tax_asv_a_tardi <- merge(tax_a_tardi, asv, by=0, all=F)
tax_asv_a_roti <- merge(tax_a_roti, asv, by=0, all=F)
tax_asv_a_mollu <- merge(tax_a_mollu, asv, by=0, all=F)
tax_asv_a_platy <- merge(tax_a_platy, asv, by=0, all=F)


#cut it to the counts only
asv_p <- tax_asv_p[,10:1059]
asv_f <- tax_asv_f[,10:1059]
asv_a_arthro <- tax_asv_a_arthro[,10:1059]
asv_a_anneli <- tax_asv_a_anneli[,10:1059]
asv_a_nema <- tax_asv_a_nema[,10:1059]
asv_a_tardi <- tax_asv_a_tardi[,10:1059]
asv_a_roti <- tax_asv_a_roti[,10:1059]
asv_a_mollu <- tax_asv_a_mollu[,10:1059]
asv_a_platy <- tax_asv_a_platy[,10:1059]

asv_p1 <-apply(asv_p, 2, FUN=function(x){as.numeric(x)})
asv_p <- as.data.frame(asv_p1)
rownames(asv_p) <- tax_asv_p[,1]

asv_f1 <-apply(asv_f, 2, FUN=function(x){as.numeric(x)})
asv_f <- as.data.frame(asv_f1)
rownames(asv_f) <- tax_asv_f[,1]

asv_a_arthro1 <-apply(asv_a_arthro, 2, FUN=function(x){as.numeric(x)})
asv_a_arthro <- as.data.frame(asv_a_arthro1)
rownames(asv_a_arthro) <- tax_asv_a_arthro[,1]

asv_a_anneli1 <-apply(asv_a_anneli, 2, FUN=function(x){as.numeric(x)})
asv_a_anneli <- as.data.frame(asv_a_anneli1)
rownames(asv_a_anneli) <- tax_asv_a_anneli[,1]

asv_a_nema1 <-apply(asv_a_nema, 2, FUN=function(x){as.numeric(x)})
asv_a_nema <- as.data.frame(asv_a_nema1)
rownames(asv_a_nema) <- tax_asv_a_nema[,1]

asv_a_tardi1 <-apply(asv_a_tardi, 2, FUN=function(x){as.numeric(x)})
asv_a_tardi <- as.data.frame(asv_a_tardi1)
rownames(asv_a_tardi) <- tax_asv_a_tardi[,1]

asv_a_roti1 <-apply(asv_a_roti, 2, FUN=function(x){as.numeric(x)})
asv_a_roti <- as.data.frame(asv_a_roti1)
rownames(asv_a_roti) <- tax_asv_a_roti[,1]

asv_a_platy1 <-apply(asv_a_platy, 2, FUN=function(x){as.numeric(x)})
asv_a_platy <- as.data.frame(asv_a_platy1)
rownames(asv_a_platy) <- tax_asv_a_platy[,1]

asv_a_mollu1 <-apply(asv_a_mollu, 2, FUN=function(x){as.numeric(x)})
asv_a_mollu <- as.data.frame(asv_a_mollu1)
rownames(asv_a_mollu) <- tax_asv_a_mollu[,1]

##############################
#CUT ASVs TO 885 SAMpLES (Curated sample data)
asv_f_t <- t(asv_f)
ss_f <- merge(sample_data, asv_f_t, by = 0, all = F)

asv_p_t <- t(asv_p)
ss_p <- merge(sample_data, asv_p_t, by = 0, all = F)

asv_a_arthro_t <- t(asv_a_arthro)
ss_a_arthro <- merge(sample_data, asv_a_arthro_t, by = 0, all = F)

asv_a_anneli_t <- t(asv_a_anneli)
ss_a_anneli <- merge(sample_data, asv_a_anneli_t, by = 0, all = F)

asv_a_tardi_t <- t(asv_a_tardi)
ss_a_tardi <- merge(sample_data, asv_a_tardi_t, by = 0, all = F)

asv_a_roti_t <- t(asv_a_roti)
ss_a_roti <- merge(sample_data, asv_a_roti_t, by = 0, all = F)

asv_a_nema_t <- t(asv_a_nema)
ss_a_nema <- merge(sample_data, asv_a_nema_t, by = 0, all = F)

asv_a_mollu_t <- t(asv_a_mollu)
ss_a_mollu <- merge(sample_data, asv_a_mollu_t, by = 0, all = F)

asv_a_platy_t <- t(asv_a_platy)
ss_a_platy <- merge(sample_data, asv_a_platy_t, by = 0, all = F)

#cut back to read_counts only: 
asv_f <- ss_f[,c(81:26138)]
asv_f1 <-apply(asv_f, 2, FUN=function(x){as.numeric(x)})
asv_f <- as.data.frame(asv_f1)

asv_p <- ss_p[,81:32138]
asv_p1 <-apply(asv_p, 2, FUN=function(x){as.numeric(x)})
asv_p <- as.data.frame(asv_p1)

asv_a_arthro <- ss_a_arthro[,81:2104]
asv_a_arthro1 <-apply(asv_a_arthro, 2, FUN=function(x){as.numeric(x)})
asv_a_arthro <- as.data.frame(asv_a_arthro1)

asv_a_anneli <- ss_a_anneli[,81:129]
asv_a_anneli1 <-apply(asv_a_anneli, 2, FUN=function(x){as.numeric(x)})
asv_a_anneli <- as.data.frame(asv_a_anneli1)

asv_a_roti <- ss_a_roti[,81:732]
asv_a_roti1 <-apply(asv_a_roti, 2, FUN=function(x){as.numeric(x)})
asv_a_roti <- as.data.frame(asv_a_roti1)

asv_a_tardi <- ss_a_tardi[,81:292]
asv_a_tardi1 <-apply(asv_a_tardi, 2, FUN=function(x){as.numeric(x)})
asv_a_tardi <- as.data.frame(asv_a_tardi1)

asv_a_nema <- ss_a_nema[,81:3290]
asv_a_nema1 <-apply(asv_a_nema, 2, FUN=function(x){as.numeric(x)})
asv_a_nema <- as.data.frame(asv_a_nema1)

asv_a_mollu <- ss_a_mollu[,81:85]
asv_a_mollu1 <-apply(asv_a_mollu, 2, FUN=function(x){as.numeric(x)})
asv_a_mollu <- as.data.frame(asv_a_mollu1)

rownames(asv_f) <- ss_f[,1]
rownames(asv_p) <- ss_p[,1]
rownames(asv_a_arthro) <- ss_a_arthro[,1]
rownames(asv_a_anneli) <- ss_a_anneli[,1]
rownames(asv_a_nema) <- ss_a_nema[,1]
rownames(asv_a_roti) <- ss_a_roti[,1]
rownames(asv_a_tardi) <- ss_a_tardi[,1]
rownames(asv_a_mollu) <- ss_a_mollu[,1]

#change back to ASVs as rows
asv_p1 <- as.data.frame(t(asv_p))
asv_f1 <- as.data.frame(t(asv_f))
asv_a1_arthro <- as.data.frame(t(asv_a_arthro))
asv_a1_anneli <- as.data.frame(t(asv_a_anneli))
asv_a1_nema <- as.data.frame(t(asv_a_nema))
asv_a1_tardi <- as.data.frame(t(asv_a_tardi))
asv_a1_roti <- as.data.frame(t(asv_a_roti))
asv_a1_mollu <- as.data.frame(t(asv_a_mollu))


#################
#Filtering / Scaling with ranked subsamples (SRS)
#################

#filter ASVs with less than 50 count
ab1 <- rowSums(asv_a1_arthro)
asv_a2_arthro <- cbind(asv_a1_arthro, ab1)
asv_a2_arthro <- as.data.frame(asv_a2_arthro)
asv_a2_arthro <- asv_a2_arthro %>% dplyr::filter(asv_a2_arthro[,788] > 1)
asv_a3_arthro <- asv_a2_arthro[,-788]
asv_a3_arthro <- as.data.frame((asv_a3_arthro))

ab1 <- rowSums(asv_a1_anneli)
asv_a2_anneli <- cbind(asv_a1_anneli, ab1)
asv_a2_anneli <- as.data.frame(asv_a2_anneli)
asv_a2_anneli <- asv_a2_anneli %>% dplyr::filter(asv_a2_anneli[,788] > 1)
asv_a3_anneli <- asv_a2_anneli[,-788]
asv_a3_anneli <- as.data.frame((asv_a3_anneli))

ab3 <- rowSums(asv_a1_nema)
asv_a2_nema <- cbind(asv_a1_nema, ab3)
asv_a2_nema <- as.data.frame(asv_a2_nema) 
asv_a2_nema <- asv_a2_nema %>% dplyr::filter(asv_a2_nema[,788] > 1)
asv_a3_nema <- asv_a2_nema[,-788]
asv_a3_nema <- as.data.frame((asv_a3_nema))

ab3 <- rowSums(asv_a1_roti)
asv_a2_roti <- cbind(asv_a1_roti, ab3)
asv_a2_roti <- as.data.frame(asv_a2_roti) 
asv_a2_roti <- asv_a2_roti %>% dplyr::filter(asv_a2_roti[,788] > 1)
asv_a3_roti <- asv_a2_roti[,-788]
asv_a3_roti <- as.data.frame((asv_a3_roti))

ab3 <- rowSums(asv_a1_tardi)
asv_a2_tardi <- cbind(asv_a1_tardi, ab3)
asv_a2_tardi <- as.data.frame(asv_a2_tardi) 
asv_a2_tardi <- asv_a2_tardi %>% dplyr::filter(asv_a2_tardi[,788] > 1)
asv_a3_tardi <- asv_a2_tardi[,-788]
asv_a3_tardi <- as.data.frame((asv_a3_tardi))

ab4 <- rowSums(asv_f1)
asv_f2 <- cbind(asv_f1, ab4)
asv_f2 <- as.data.frame(asv_f2) 
asv_f2 <- asv_f2 %>% dplyr::filter(asv_f2[,788] > 1)
asv_f3 <- asv_f2[,-788]
asv_f3 <- as.data.frame((asv_f3))

ab5 <- rowSums(asv_p1)
asv_p2 <- cbind(asv_p1, ab5)
asv_p2 <- as.data.frame(asv_p2) 
asv_p2 <- asv_p2 %>% dplyr::filter(asv_p2[,788] > 1)
asv_p3 <- asv_p2[,-788]
asv_p3 <- as.data.frame((asv_p3))

ab5 <- rowSums(asv_a1_mollu)
asv_mollu2 <- cbind(asv_a1_mollu, ab5)
asv_mollu2 <- as.data.frame(asv_mollu2) 
asv_mollu2 <- asv_mollu2 %>% dplyr::filter(asv_mollu2[,788] > 1)
asv_mollu3 <- asv_mollu2[,-788]
asv_mollu3 <- as.data.frame((asv_mollu3))

#preparation for normalization via SRS.shiny.app
Cmin_f <-2000 
Cmin_p <- 2000 
Cmin_a_arthro <- 1600 
Cmin_a_anneli <- 330 
Cmin_a_nema <- 2500 
Cmin_a_tardi <- 80
Cmin_a_roti <- 540 

#Normalization with scaling with ranked sub-sampling   
asv_p_SRS <- SRS(asv_p3, Cmin_p,set_seed = TRUE, seed = 1)
asv_f_SRS <- SRS(asv_f3, Cmin_f,set_seed = TRUE, seed = 1)
asv_arthro_SRS <- SRS(asv_a3_arthro, Cmin_a_arthro,set_seed = TRUE, seed = 1)
asv_anneli_SRS <- SRS(asv_a3_anneli, Cmin_a_anneli,set_seed = TRUE, seed = 1)
asv_tardi_SRS <- SRS(asv_a3_tardi, Cmin_a_tardi,set_seed = TRUE, seed = 1)
asv_roti_SRS <- SRS(asv_a3_roti, Cmin_a_roti,set_seed = TRUE, seed = 1)
asv_nema_SRS <- SRS(asv_a3_nema, Cmin_a_nema,set_seed = TRUE, seed = 1)

##SRS erases rownames, to get them back:
rownames(asv_f_SRS) <- rownames(asv_f3)
rownames(asv_p_SRS) <- rownames(asv_p3)
rownames(asv_arthro_SRS) <- rownames(asv_a3_arthro)
rownames(asv_anneli_SRS) <- rownames(asv_a3_anneli)
rownames(asv_tardi_SRS) <- rownames(asv_a3_tardi)
rownames(asv_roti_SRS) <- rownames(asv_a3_roti)
rownames(asv_nema_SRS) <- rownames(asv_a3_nema)

##################
#to phyloseq file
##################
asv_f_SRS_m <- (as.matrix(asv_f_SRS))
asv_p_SRS_m <- (as.matrix(asv_p_SRS))
asv_arthro_SRS_m <- (as.matrix(asv_arthro_SRS))
asv_anneli_SRS_m <- (as.matrix(asv_anneli_SRS))
asv_roti_SRS_m <- (as.matrix(asv_roti_SRS))
asv_tardi_SRS_m <- (as.matrix(asv_tardi_SRS))
asv_nema_SRS_m <- (as.matrix(asv_nema_SRS))

tax_f_m <- as.matrix(tax_f)
tax_p_m <- as.matrix(tax_p)
tax_a_arthro_m <- as.matrix(tax_a_arthro)
tax_a_anneli_m <- as.matrix(tax_a_anneli)
tax_a_tardi_m <- as.matrix(tax_a_tardi)
tax_a_roti_m <- as.matrix(tax_a_roti)
tax_a_nema_m <- as.matrix(tax_a_nema)



TAX_p <- tax_table(tax_p_m)
TAX_f <- tax_table(tax_f_m)
TAX_a_arthro <- tax_table(tax_a_arthro_m)
TAX_a_anneli <- tax_table(tax_a_anneli_m)
TAX_a_tardi <- tax_table(tax_a_tardi_m)
TAX_a_roti <- tax_table(tax_a_roti_m)
TAX_a_nema <- tax_table(tax_a_nema_m)

ASV_p <- otu_table(asv_p_SRS_m, taxa_are_rows = TRUE) 
ASV_f <- otu_table(asv_f_SRS_m, taxa_are_rows = TRUE) 
ASV_a_arthro <- otu_table(asv_arthro_SRS_m, taxa_are_rows = TRUE)
ASV_a_anneli <- otu_table(asv_anneli_SRS_m, taxa_are_rows = TRUE)
ASV_a_tardi <- otu_table(asv_tardi_SRS_m, taxa_are_rows = TRUE)
ASV_a_roti <- otu_table(asv_roti_SRS_m, taxa_are_rows = TRUE)
ASV_a_nema <- otu_table(asv_nema_SRS_m, taxa_are_rows = TRUE)

#add diversity to sample_data data
asv_f_SRS_m1 <- t(as.matrix(asv_f_SRS))
asv_p_SRS_m1 <- t(as.matrix(asv_p_SRS))
asv_arthro_SRS_m1 <- t(as.matrix(asv_arthro_SRS))
asv_anneli_SRS_m1 <- t(as.matrix(asv_anneli_SRS))
asv_tardi_SRS_m1 <- t(as.matrix(asv_tardi_SRS))
asv_roti_SRS_m1 <- t(as.matrix(asv_roti_SRS))
asv_nema_SRS_m1 <- t(as.matrix(asv_nema_SRS))

a_arthro_vegan_shannon <- as.data.frame(vegan::diversity(asv_arthro_SRS_m1, index ="shannon"))
a_arthro_vegan_observed <- as.data.frame(specnumber(asv_arthro_SRS_m1))

a_anneli_vegan_shannon <- as.data.frame(vegan::diversity(asv_anneli_SRS_m1, index ="shannon"))
a_anneli_vegan_observed <- as.data.frame(specnumber(asv_anneli_SRS_m1))

a_tardi_vegan_shannon <- as.data.frame(vegan::diversity(asv_tardi_SRS_m1, index ="shannon"))
a_tardi_vegan_observed <- as.data.frame(specnumber(asv_tardi_SRS_m1))

a_roti_vegan_shannon <- as.data.frame(vegan::diversity(asv_roti_SRS_m1, index ="shannon"))
a_roti_vegan_observed <- as.data.frame(specnumber(asv_roti_SRS_m1))

a_nema_vegan_shannon <- as.data.frame(vegan::diversity(asv_nema_SRS_m1, index ="shannon"))
a_nema_vegan_observed <- as.data.frame(specnumber(asv_nema_SRS_m1))

f_vegan_shannon <- as.data.frame(vegan::diversity(asv_f_SRS_m1, index ="shannon"))
f_vegan_observed <- as.data.frame(specnumber(asv_f_SRS_m1))

p_vegan_shannon <- as.data.frame(vegan::diversity(asv_p_SRS_m1, index ="shannon"))
p_vegan_observed <- as.data.frame(specnumber(asv_p_SRS_m1))

#sample_data data for different groups
p_alpha <- cbind(p_vegan_observed,p_vegan_shannon)
p_sample_data <- merge(p_alpha,sample_data,by=0)
rownames(p_sample_data) <- p_sample_data$BARCODE_ID
p_sample_data <- p_sample_data[,-1]
colnames(p_sample_data)[1] <- "Richness"
colnames(p_sample_data)[2] <- "Shannon"

f_alpha <- cbind(f_vegan_observed,f_vegan_shannon)
f_sample_data <- merge(f_alpha,sample_data,by=0)
rownames(f_sample_data) <- f_sample_data$BARCODE_ID
f_sample_data <- f_sample_data[,-1]
colnames(f_sample_data)[1] <- "Richness"
colnames(f_sample_data)[2] <- "Shannon"

a_arthro_alpha <- cbind(a_arthro_vegan_observed,a_arthro_vegan_shannon)
a_arthro_sample_data <- merge(a_arthro_alpha,sample_data,by=0)
rownames(a_arthro_sample_data) <- a_arthro_sample_data$BARCODE_ID
a_arthro_sample_data <- a_arthro_sample_data[,-1]
colnames(a_arthro_sample_data)[1] <- "Richness"
colnames(a_arthro_sample_data)[2] <- "Shannon"

a_anneli_alpha <- cbind(a_anneli_vegan_observed,a_anneli_vegan_shannon)
a_anneli_sample_data <- merge(a_anneli_alpha,sample_data,by=0)
rownames(a_anneli_sample_data) <- a_anneli_sample_data$BARCODE_ID
a_anneli_sample_data <- a_anneli_sample_data[,-1]
colnames(a_anneli_sample_data)[1] <- "Richness"
colnames(a_anneli_sample_data)[2] <- "Shannon"

a_tardi_alpha <- cbind(a_tardi_vegan_observed,a_tardi_vegan_shannon)
a_tardi_sample_data <- merge(a_tardi_alpha,sample_data,by=0)
rownames(a_tardi_sample_data) <- a_tardi_sample_data$BARCODE_ID
a_tardi_sample_data <- a_tardi_sample_data[,-1]
colnames(a_tardi_sample_data)[1] <- "Richness"
colnames(a_tardi_sample_data)[2] <- "Shannon"

a_roti_alpha <- cbind(a_roti_vegan_observed,a_roti_vegan_shannon)
a_roti_sample_data <- merge(a_roti_alpha,sample_data,by=0)
rownames(a_roti_sample_data) <- a_roti_sample_data$BARCODE_ID
a_roti_sample_data <- a_roti_sample_data[,-1]
colnames(a_roti_sample_data)[1] <- "Richness"
colnames(a_roti_sample_data)[2] <- "Shannon"

a_nema_alpha <- cbind(a_nema_vegan_observed,a_nema_vegan_shannon)
a_nema_sample_data <- merge(a_nema_alpha,sample_data,by=0)
rownames(a_nema_sample_data) <- a_nema_sample_data$BARCODE_ID
a_nema_sample_data <- a_nema_sample_data[,-1]
colnames(a_nema_sample_data)[1] <- "Richness"
colnames(a_nema_sample_data)[2] <- "Shannon"

SAMPLE_DATA_f = sample_data(f_sample_data)
SAMPLE_DATA_p = sample_data(p_sample_data)
SAMPLE_DATA_a_nema = sample_data(a_nema_sample_data)
SAMPLE_DATA_a_arthro = sample_data(a_arthro_sample_data)
SAMPLE_DATA_a_anneli = sample_data(a_anneli_sample_data)
SAMPLE_DATA_a_tardi = sample_data(a_tardi_sample_data)
SAMPLE_DATA_a_roti = sample_data(a_roti_sample_data)


#create phyloseq data
lc_f <- phyloseq(ASV_f, TAX_f, SAMPLE_DATA_f) #create phyloseq file fungi
lc_p <- phyloseq(ASV_p, TAX_p, SAMPLE_DATA_p) #create phyloseq file protists
lc_a_arthro <- phyloseq(ASV_a_arthro, TAX_a_arthro, SAMPLE_DATA_a_arthro) #create phyloseq file arthrofauna
lc_a_anneli <- phyloseq(ASV_a_anneli, TAX_a_anneli, SAMPLE_DATA_a_anneli) #create phyloseq file arthrofauna
lc_a_tardi <- phyloseq(ASV_a_tardi, TAX_a_tardi, SAMPLE_DATA_a_tardi) #create phyloseq file arthrofauna
lc_a_roti <- phyloseq(ASV_a_roti, TAX_a_roti, SAMPLE_DATA_a_roti) #create phyloseq file arthrofauna
lc_a_nema <- phyloseq(ASV_a_nema, TAX_a_nema, SAMPLE_DATA_a_nema) #create phyloseq file nemafauna 

#create relative abundance phyloseq data
lc_f_rel_abund = phyloseq::transform_sample_counts(lc_f, function(x){x / sum(x)})
lc_p_rel_abund = phyloseq::transform_sample_counts(lc_p, function(x){x / sum(x)})
lc_a_arthro_rel_abund = phyloseq::transform_sample_counts(lc_a_arthro, function(x){x / sum(x)})
lc_a_anneli_rel_abund = phyloseq::transform_sample_counts(lc_a_anneli, function(x){x / sum(x)})
lc_a_roti_rel_abund = phyloseq::transform_sample_counts(lc_a_roti, function(x){x / sum(x)})
lc_a_tardi_rel_abund = phyloseq::transform_sample_counts(lc_a_tardi, function(x){x / sum(x)})
lc_a_nema_rel_abund = phyloseq::transform_sample_counts(lc_a_nema, function(x){x / sum(x)})


