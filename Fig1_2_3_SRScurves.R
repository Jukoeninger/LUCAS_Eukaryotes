##################################################
#PLOTS

#Figure 1
#Correlation Season + Land Cover
df1 <- sample_data[,c("Sample_season","lat","lon","LC1_2018","BARCODE_ID","LC1_2018_num")]
df1_long <- df1 %>% tidyr::gather(group, value, -lon, -lat, -LC1_2018_num,-Sample_season)

df1_long$LC1_2018_num[df1_long$LC1_2018_num =="1"] <- "CL1"
df1_long$LC1_2018_num[df1_long$LC1_2018_num =="2"] <- "CL2"
df1_long$LC1_2018_num[df1_long$LC1_2018_num =="3"] <- "GL1"
df1_long$LC1_2018_num[df1_long$LC1_2018_num =="4"] <- "GL2"
df1_long$LC1_2018_num[df1_long$LC1_2018_num =="5"] <- "WL1"
df1_long$LC1_2018_num[df1_long$LC1_2018_num =="6"] <- "WL2"

df1_long$LC1_2018_num <- factor(df1_long$LC1_2018_num, levels=c("CL1","CL2","GL1","GL2","WL1","WL2"))

library("rnaturalearth")
library("rnaturalearthdata")
world <- ne_countries(scale = 50, returnclass = "sf")
EU <- c("Austria","Belgium","Bulgaria","Croatia","Cyprus",
        "Czech Rep.","Denmark","Estonia","Finland","France",
        "Germany","Greece","Hungary","Ireland","Italy","Latvia",
        "Lithuania","Luxembourg","Malta","Netherlands","poland",
        "portugal","Romania","Slovakia","Slovenia","Spain",
        "Sweden","United Kingdom")
EU_map <- 
  world %>% 
  filter(name %in% EU)
svg("1305_Fig1_A.svg",width=12)# SVG graphics device
EU_ecosystem_sites <- ggplot(data = EU_map) +
  geom_sf() +
  geom_point(data = df1_long, aes(x = lon, y = lat, color=LC1_2018_num), shape=20,size =3) +
  coord_sf(xlim = c(-15, 35), ylim = c(35, 70))+theme_bw() +
  ylab("Latitude")+
  xlab("Longitude")+
  scale_colour_manual(values=c("#FCE205","orange","#2AAA8A","#98FB98","#FE036A","#7d1969"))+
  guides(colour = guide_legend("Ecosystem type",override.aes = list(size=7)))+
  theme_bw() +
  theme(text = element_text(size=12, family = "sans"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.grid  = element_blank(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18, hjust = 0.5,vjust = 3),
        plot.margin = unit(c(1, 1, 1, 1), units =, "cm"),
        legend.title = element_text(size=18,face = "bold"),
        legend.text = element_text(size = 18),
        legend.key.size = unit(1,"cm"),
        legend.key.width= unit(0, "cm"),
        legend.background = element_rect(color = "black",
                                         fill = "transparent",
                                         size = 0.5, linetype = "blank"))
dev.off() # Close the graphics device


df1_long$Sample_season <- factor(df1_long$Sample_season, levels=c("Spring","Summer","Autumn"))
#svg("1305_Fig1_B.svg",width=12)# SVG graphics device
EU_season_sites <- ggplot(data = EU_map) +
  geom_sf() +
  geom_point(data = df1_long, aes(x = lon, y = lat, color=Sample_season,), shape=20,size =3) +
  
  guides(colour = guide_legend("Season",override.aes = list(size=7)))+
  coord_sf(xlim = c(-15, 35), ylim = c(35, 70))+theme_bw() +
  ylab("Latitude")+
  xlab("Longitude")+
  scale_colour_manual(values=c("#808000","#FFAE42","brown")) +theme_bw() +
  theme(text = element_text(size=12, family = "sans"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.grid  = element_blank(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18, hjust = 0.5,vjust = 3),
        plot.margin = unit(c(1, 1, 1, 1), units =, "cm"),
        legend.title = element_text(size=18,face = "bold"),
        legend.text = element_text(size = 18),
        legend.key.size = unit(1,"cm"),
        legend.key.width= unit(0, "cm"),
        legend.background = element_rect(color = "black",
                                         fill = "transparent",
                                         size = 0.5, linetype = "blank"))
#dev.off() # Close the graphics device


#Site distribution season, land cover
df <- sample_data[,c("Sample_season","LC1_2018","BARCODE_ID","LC1_2018_num")]
df_long <- df %>% tidyr::gather(group, value, -BARCODE_ID, -LC1_2018_num,-Sample_season)
df_long$Sample_season <- factor(df_long$Sample_season, levels=c("Spring","Summer","Autumn"))
g <- ggplot(df_long,aes(LC1_2018_num)) 
#svg("1305_Fig1_C.svg",width=12)# SVG graphics device
ecosystem_season_barplot <- g + geom_bar(aes(fill=Sample_season)) + 
  xlab("\nEcosystem type")+
  ylab("\nNumber of samples")+
  scale_x_discrete(name="Ecosystem type", limits = c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  scale_fill_manual(name="Sampling season",values=c("#708238","orange","#E47250","#5A4A6F","#EBB261", "#9D5A6C","#4393C3", "#2166AC", "#053061","#E69F00", "#56B4E9")) +
  theme(text = element_text(size=12, family = "sans"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.grid  = element_blank(),
        axis.text = element_text(size = 18),
        plot.title = element_text(size = 12, vjust = 1, hjust = 0),
        axis.title = element_text(size = 18, hjust = 0.5,vjust = 0.5),
        plot.margin = unit(c(1, 1, 1, 1), units =, "cm"),
        legend.text = element_text(size = 18),
        legend.key.size = unit(1,"cm"),
        legend.key.width= unit(1, "cm"),
        legend.title = element_text(size=18,face = "bold"),
        legend.background = element_rect(color = "black",
                                         fill = "transparent",
                                         size = 0.5, linetype = "blank"))
#dev.off() # Close the graphics device

svg("1806_Fig1ABC.svg",width=14,height=15)
#plot
ggarrange(ecosystem_season_barplot,ggarrange(EU_ecosystem_sites,EU_season_sites, legend="bottom",common.legend =F,ncol=2,widths=c(0.9,0.7),labels = c("(b)","(c)")),nrow=2,legend="right",labels="(a)")
dev.off() # Close the graphics device


#Taxa pie charts
tax <- read.csv("~/Euk_18S_data/tax_uncurated.csv", header=1,row.names=1)
rownames(tax) <- 1:nrow(tax)
rownames(tax) <- sub("^","zotu_",rownames(tax))
tax$domain[tax$domain =="ukn"] <- "Unknown" #change ukn to unknown
table_tax_domain <- table(tax$domain)
table_tax_domain <- as.data.frame(table_tax_domain)
table_tax_domain <- table_tax_domain %>%
  arrange(desc(Var1)) %>%
  mutate(percentage=paste0(round(Freq/sum(Freq)*100,1),"%"))
legend_title<- "Domain"
pos <- position_jitter(width = 0.5, seed=1)
#svg("1305_Fig1_A.svg",width=12)# SVG graphics device
euk_pie <- ggplot(table_tax_domain, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white")+
  coord_polar("y")+
  geom_label_repel(size=7.0,aes(label = percentage, y = Freq),show.legend = F,position=position_stack(vjust = 0.5),label.size=0,fontface = 'bold',color = 'white')  +
  scale_fill_manual(legend_title,values = c("#2166AC","#900C3F","#CD5C5C","#999999"))+
  ggtitle("18S read results")+ theme_void() + theme(text = element_text(size=12, family = "sans"),
                                                    panel.grid.major = element_blank(), 
                                                    panel.grid.minor = element_blank(),
                                                    panel.background = element_blank(),
                                                    panel.grid  = element_blank(),
                                                    axis.text = element_blank(), 
                                                    axis.title = element_blank(),
                                                    plot.margin = unit(c(1, 1, 1, 1), units =, "cm"),
                                                    plot.title = element_blank(),
                                                    legend.text = element_text(size = 18), 
                                                    legend.key.size = unit(1,"cm"),
                                                    legend.key.width= unit(1, "cm"),
                                                    legend.title = element_text(size=20),
                                                    legend.background = element_rect(color = "black", 
                                                                                     fill = "transparent", 
                                                                                     size = 0.5, linetype = "blank"))

#dev.off() # Close the graphics device


tax_domain_ukn <- tax %>% filter(domain == "Unknown")
#Filter out eukaryotes
tax_euk_king <- tax %>% filter(domain == "Eukaryota")
#Add kingdom raw to taxa
tax <- mutate(tax_euk_king, kingdom = if_else(tax_euk_king$phylum == "Ascomycota", "Fungi", if_else(tax_euk_king$phylum =="Zoopagomycota","Fungi", if_else(tax_euk_king$phylum == "Basidiomycota", "Fungi", if_else(tax_euk_king$phylum == "Mucoromycota", "Fungi", if_else(tax_euk_king$phylum == "Cryptomycota", "Fungi", if_else(tax_euk_king$phylum == "Chytridiomycota", "Fungi", if_else(tax_euk_king$phylum == "Aphelidea", "Fungi", if_else(tax_euk_king$phylum == "Blastocladiomycota", "Fungi", if_else(tax_euk_king$phylum == "Dikarya_Incertae_Sedis", "Fungi", if_else(tax_euk_king$phylum == "LKM15", "Fungi", if_else(tax_euk_king$phylum == "Neocallimastigomycota", "Fungi", if_else(tax_euk_king$phylum == "undef_Glomeromycotina_sp.", "Fungi", if_else(tax_euk_king$phylum == "undef_Sanchytrium_tribonematis", "Fungi", if_else(tax_euk_king$phylum == "undef_uncultured_fungus", "Fungi", if_else(tax_euk_king$phylum == "undef_uncultured_Tremellaceae", "Fungi", if_else(tax_euk_king$phylum == "Rotifera", "Animalia", if_else(tax_euk_king$phylum == "Platytricha", "Animalia", if_else(tax_euk_king$phylum == "Ctenophora", "Animalia", if_else(tax_euk_king$phylum == "Arthropoda", "Animalia", if_else(tax_euk_king$phylum == "Nematozoa", "Animalia", if_else(tax_euk_king$phylum == "Placozoa", "Animalia", if_else(tax_euk_king$phylum == "Tardigrada", "Animalia", if_else(tax_euk_king$phylum == "Porifera", "Animalia", if_else(tax_euk_king$phylum =="Platyhelminthes", "Animalia", if_else(tax_euk_king$phylum == "Vertebrata", "Animalia", if_else(tax_euk_king$phylum == "Micrognathozoa", "Animalia", if_else(tax_euk_king$phylum == "Mollusca", "Animalia", if_else(tax_euk_king$phylum == "Platytricha", "Animalia", if_else(tax_euk_king$phylum == "undef_Daphnia_magna", "Animalia", if_else(tax_euk_king$phylum == "Annelida", "Animalia", if_else(tax_euk_king$phylum == "Rotifera", "Animalia", "Protista"))))))))))))))))))))))))))))))))

table_tax_kingdom <- table(tax$kingdom)
table_tax_kingdom <- as.data.frame(table_tax_kingdom)
table_tax_kingdom <- table_tax_kingdom %>%
  arrange(desc(Var1)) %>%
  mutate(percentage=paste0(round(Freq/sum(Freq)*100,1),"%"))
legend_title <- "Eukaryota"
#svg("1305_Fig1_B.svg",width=12)# SVG graphics device
kingdom_pie <- ggplot(table_tax_kingdom, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white")+
  coord_polar("y")+
  geom_label_repel(size=7.0,aes(label = percentage, y = Freq),show.legend = F,position=position_stack(vjust = 0.5),label.size=0,fontface = 'bold',color = 'white')  +
  scale_fill_manual(legend_title,values = c("#9D5A6C","#708238","#4393C3","#2166AC"))+
  ggtitle("18S Eukaryotes")+ theme_void() + theme(text = element_text(size=12, family = "sans"),
                                                  panel.grid.major = element_blank(), 
                                                  panel.grid.minor = element_blank(),
                                                  panel.background = element_blank(),
                                                  panel.grid  = element_blank(),
                                                  axis.text = element_blank(), 
                                                  axis.title = element_blank(),
                                                  plot.title = element_blank(),
                                                  plot.margin = unit(c(1, 1, 1, 1), units =, "cm"),
                                                  legend.text = element_text(size = 18), 
                                                  legend.key.size = unit(1,"cm"),
                                                  legend.key.width= unit(1, "cm"),
                                                  legend.title = element_text(size=20),
                                                  legend.background = element_rect(color = "black", 
                                                                                   fill = "transparent", 
                                                                                   size = 0.5, linetype = "blank"))
#dev.off() # Close the graphics device

tax_a <- rbind(tax_a_macro,tax_a_meso,tax_a_micro)
table_tax_fauna <- table(tax_a$fauna)
table_tax_fauna <- as.data.frame(table_tax_fauna)
table_tax_fauna <- table_tax_fauna %>%
  arrange(desc(Var1)) %>%
  mutate(percentage=paste0(round(Freq/sum(Freq)*100,1),"%"))
legend_title<- "Kingdom Animalia"
table_tax_fauna$Var1 <- factor(table_tax_fauna$Var1, levels=c("Microfauna", "Mesofauna", "Macrofauna","Unknown"))

#svg("1305_Fig1_C.svg")# SVG graphics device
fauna_pie <- ggplot(table_tax_fauna, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white")+
  coord_polar("y")+
  geom_label_repel(size=7.0,aes(label = percentage, y = Freq),show.legend = F,position=position_stack(vjust = 0.5),label.size=0,fontface = 'bold',color = 'white')  +
  scale_fill_manual(legend_title,values = c( "#EBB261", "#E47250","#5A4A6F","#FE036A"))+
  ggtitle("Kingdom Animalia")+ theme_void() + theme(text = element_text(size=12, family = "sans"),
                                                    panel.grid.major = element_blank(), 
                                                    panel.grid.minor = element_blank(),
                                                    panel.background = element_blank(),
                                                    panel.grid  = element_blank(),
                                                    axis.text = element_blank(), 
                                                    axis.title = element_blank(),
                                                    plot.title = element_blank(),
                                                    plot.margin = unit(c(1, 1, 1, 1), units =, "cm"),
                                                    legend.text = element_text(size = 18), 
                                                    legend.key.size = unit(1,"cm"),
                                                    legend.key.width= unit(1, "cm"),
                                                    legend.title = element_text(size=20),
                                                    legend.background = element_rect(color = "black", 
                                                                                     fill = "transparent", 
                                                                                     size = 0.5, linetype = "blank"))

tax_fauna <- rbind(tax_a_roti,tax_a_tardi,tax_a_nema,tax_a_arthro,tax_a_anneli)
table_tax_fauna <- table(tax_fauna$phylum)
table_tax_fauna <- as.data.frame(table_tax_fauna)
table_tax_fauna <- table_tax_fauna %>%
  arrange(desc(Var1)) %>%
  mutate(percentage=paste0(round(Freq/sum(Freq)*100,1),"%"))
legend_title<- "Animalia"
table_tax_fauna$Var1 <- factor(table_tax_fauna$Var1, levels=c("Rotifera", "Tardigrada", "Nematozoa","Arthropoda","Annelida"))

svg("1305_Fig1_C.svg")# SVG graphics device
fauna_pie <- ggplot(table_tax_fauna, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white")+
  coord_polar("y")+
  geom_label_repel(size=7.0,aes(label = percentage, y = Freq),show.legend = F,position=position_stack(vjust = 0.5),label.size=0,fontface = 'bold',color = 'white')  +
  scale_fill_manual(legend_title,values = c( "#EBB261", "#E47250","#5A4A6F","#cb99c9","#B22222"))+
  ggtitle("Animalia")+ theme_void() + theme(text = element_text(size=12, family = "sans"),
                                            panel.grid.major = element_blank(), 
                                            panel.grid.minor = element_blank(),
                                            panel.background = element_blank(),
                                            panel.grid  = element_blank(),
                                            axis.text = element_blank(), 
                                            axis.title = element_blank(),
                                            plot.margin = unit(c(1, 1, 1, 1), units =, "cm"),
                                            plot.title = element_blank(),
                                            legend.text = element_text(size = 18), 
                                            legend.key.size = unit(1,"cm"),
                                            legend.key.width= unit(1, "cm"),
                                            legend.title = element_text(size=20),
                                            legend.background = element_rect(color = "black", 
                                                                             fill = "transparent", 
                                                                             size = 0.5, linetype = "blank"))

dev.off() # Close the graphics device

svg("2305_Fig1ABC.svg",width=16,height=14)
#plot
ggarrange(euk_pie,kingdom_pie,fauna_pie,nrow=1,labels = c("(a)","(b)","(c)"),legend="right",ncol=3,vjust=42)
dev.off() # Close the graphics device



#SRS curves, Figure S12
#Fungi
svg("1305_FigS12_SRS_f.svg",width=8,height=4)# SVG graphics device
SRScurve(asv_f3, metric = "richness",  step = 1000, sample =0,xlim=c(0,5000),rarefy.repeats =10,max.sample.size = 6000, rarefy.comparison = FALSE, rarefy.comparison.legend = FALSE, xlab = "sample size",ylab = "richness", label = F,col = c("#5A4A6F", "#E47250",  "#EBB261"))
dev.off() # Close the graphics device
#protists
svg("1305_FigS12_SRS_p.svg",width=8,height=4)# SVG graphics device
SRScurve(asv_p3, metric = "richness",  step = 1000, sample = 0,xlim=c(0,5000),max.sample.size = 6000, rarefy.comparison = FALSE,rarefy.comparison.legend = FALSE, xlab = "sample size", ylab = "richness", label = FALSE,col = c("#5A4A6F", "#E47250",  "#EBB261"))
dev.off() # Close the graphics device
#rotifers
svg("1305_FigS12_SRS_roti.svg",width=8,height=4)# SVG graphics device
SRScurve(asv_a3_roti, metric = "richness",  step = 1000, sample = 0, xlim=c(0,5000),max.sample.size = 6000, rarefy.comparison = FALSE, rarefy.repeats = 10, rarefy.comparison.legend = FALSE, xlab = "sample size", ylab = "richness", label = FALSE,col = c("#5A4A6F", "#E47250",  "#EBB261"))
dev.off() # Close the graphics device
#tardigrades
svg("1305_FigS12_SRS_tardi.svg",width=8,height=4)# SVG graphics device
SRScurve(asv_a3_tardi, metric = "richness", step = 1000, sample = 0,xlim=c(0,5000), max.sample.size = 6000, rarefy.comparison = FALSE, rarefy.repeats = 10, rarefy.comparison.legend = FALSE, xlab = "sample size", ylab = "richness", label = FALSE,col = c("#5A4A6F", "#E47250",  "#EBB261"))
dev.off() # Close the graphics device
#nematodes
svg("1305_FigS12_SRS_nema.svg",width=8,height=4)# SVG graphics device
SRScurve(asv_a3_nema, metric = "richness", step = 1000, sample = 0,xlim=c(0,5000), max.sample.size = 6000, rarefy.comparison = FALSE, rarefy.repeats = 10, rarefy.comparison.legend = FALSE, xlab = "sample size", ylab = "richness", label = FALSE,col = c("#5A4A6F", "#E47250",  "#EBB261"))
dev.off() # Close the graphics device
#arthropods
svg("1305_FigS12_SRS_arthro.svg",width=8,height=4)# SVG graphics device
SRScurve(asv_a3_arthro, metric = "richness", step = 1000, sample = 0,xlim=c(0,5000), max.sample.size = 6000, rarefy.comparison = FALSE, rarefy.repeats = 10, rarefy.comparison.legend = FALSE, xlab = "sample size", ylab = "richness", label = FALSE,col = c("#5A4A6F", "#E47250",  "#EBB261"))
dev.off() # Close the graphics device
#annelids
svg("1305_FigS12_SRS_anneli.svg",width=8,height=4)# SVG graphics device
SRScurve(asv_a3_anneli, metric = "richness", step = 1000, sample = 0,xlim=c(0,5000), max.sample.size = 6000, rarefy.comparison = FALSE, rarefy.repeats = 10, rarefy.comparison.legend = FALSE, xlab = "sample size", ylab = "richness", label = FALSE,col = c("#5A4A6F", "#E47250",  "#EBB261"))
dev.off() # Close the graphics device





#############################################
#Overlapping taxa, Figure 3
#subsampling to ecosystem type
#anneli
lc_a_anneli_cl <- subset_samples(lc_a_anneli, LC_simpl_2018 == "Cropland")
lc_a_anneli_wl <- subset_samples(lc_a_anneli, LC_simpl_2018 == "Woodland")
lc_a_anneli_gl <- subset_samples(lc_a_anneli, LC_simpl_2018 == "Grassland")

sampledf_anneli_cl <- data.frame(sample_data(lc_a_anneli_cl))
sampledf_anneli_wl <- data.frame(sample_data(lc_a_anneli_wl))
sampledf_anneli_gl <- data.frame(sample_data(lc_a_anneli_gl))

asv_anneli_cl <- data.frame(as(otu_table(lc_a_anneli_cl),"matrix"))
asv_anneli_wl <- data.frame(as(otu_table(lc_a_anneli_wl),"matrix"))
asv_anneli_gl <- data.frame(as(otu_table(lc_a_anneli_gl),"matrix"))

#arthro
lc_a_arthro_cl <- subset_samples(lc_a_arthro, LC_simpl_2018 == "Cropland")
lc_a_arthro_wl <- subset_samples(lc_a_arthro, LC_simpl_2018 == "Woodland")
lc_a_arthro_gl <- subset_samples(lc_a_arthro, LC_simpl_2018 == "Grassland")

sampledf_arthro_cl <- data.frame(sample_data(lc_a_arthro_cl))
sampledf_arthro_wl <- data.frame(sample_data(lc_a_arthro_wl))
sampledf_arthro_gl <- data.frame(sample_data(lc_a_arthro_gl))

asv_arthro_cl <- data.frame(as(otu_table(lc_a_arthro_cl),"matrix"))
asv_arthro_wl <- data.frame(as(otu_table(lc_a_arthro_wl),"matrix"))
asv_arthro_gl <- data.frame(as(otu_table(lc_a_arthro_gl),"matrix"))


#tardigrades
lc_a_tardi_cl <- subset_samples(lc_a_tardi, LC_simpl_2018 == "Cropland")
lc_a_tardi_wl <- subset_samples(lc_a_tardi, LC_simpl_2018 == "Woodland")
lc_a_tardi_gl <- subset_samples(lc_a_tardi, LC_simpl_2018 == "Grassland")

sampledf_tardi_cl <- data.frame(sample_data(lc_a_tardi_cl))
sampledf_tardi_wl <- data.frame(sample_data(lc_a_tardi_wl))
sampledf_tardi_gl <- data.frame(sample_data(lc_a_tardi_gl))

asv_tardi_cl <- data.frame(as(otu_table(lc_a_tardi_cl),"matrix"))
asv_tardi_wl <- data.frame(as(otu_table(lc_a_tardi_wl),"matrix"))
asv_tardi_gl <- data.frame(as(otu_table(lc_a_tardi_gl),"matrix"))

#rotifers
lc_a_roti_cl <- subset_samples(lc_a_roti, LC_simpl_2018 == "Cropland")
lc_a_roti_wl <- subset_samples(lc_a_roti, LC_simpl_2018 == "Woodland")
lc_a_roti_gl <- subset_samples(lc_a_roti, LC_simpl_2018 == "Grassland")

sampledf_roti_cl <- data.frame(sample_data(lc_a_roti_cl))
sampledf_roti_wl <- data.frame(sample_data(lc_a_roti_wl))
sampledf_roti_gl <- data.frame(sample_data(lc_a_roti_gl))

asv_roti_cl <- data.frame(as(otu_table(lc_a_roti_cl),"matrix"))
asv_roti_wl <- data.frame(as(otu_table(lc_a_roti_wl),"matrix"))
asv_roti_gl <- data.frame(as(otu_table(lc_a_roti_gl),"matrix"))

#nematodes
lc_a_nema_cl <- subset_samples(lc_a_nema, LC_simpl_2018 == "Cropland")
lc_a_nema_wl <- subset_samples(lc_a_nema, LC_simpl_2018 == "Woodland")
lc_a_nema_gl <- subset_samples(lc_a_nema, LC_simpl_2018 == "Grassland")

sampledf_nema_cl <- data.frame(sample_data(lc_a_nema_cl))
sampledf_nema_wl <- data.frame(sample_data(lc_a_nema_wl))
sampledf_nema_gl <- data.frame(sample_data(lc_a_nema_gl))

asv_nema_cl <- data.frame(as(otu_table(lc_a_nema_cl),"matrix"))
asv_nema_wl <- data.frame(as(otu_table(lc_a_nema_wl),"matrix"))
asv_nema_gl <- data.frame(as(otu_table(lc_a_nema_gl),"matrix"))

#Protists
lc_p_cl <- subset_samples(lc_p, LC_simpl_2018 == "Cropland")
lc_p_wl <- subset_samples(lc_p, LC_simpl_2018 == "Woodland")
lc_p_gl <- subset_samples(lc_p, LC_simpl_2018 == "Grassland")

sampledf_p_cl <- data.frame(sample_data(lc_p_cl))
sampledf_p_wl <- data.frame(sample_data(lc_p_wl))
sampledf_p_gl <- data.frame(sample_data(lc_p_gl))

asv_p_cl <- data.frame(as(otu_table(lc_p_cl),"matrix"))
asv_p_wl <- data.frame(as(otu_table(lc_p_wl),"matrix"))
asv_p_gl <- data.frame(as(otu_table(lc_p_gl),"matrix"))

#Fungi
lc_f_cl <- subset_samples(lc_f, LC_simpl_2018 == "Cropland")
lc_f_wl <- subset_samples(lc_f, LC_simpl_2018 == "Woodland")
lc_f_gl <- subset_samples(lc_f, LC_simpl_2018 == "Grassland")

sampledf_f_cl <- data.frame(sample_data(lc_f_cl))
sampledf_f_wl <- data.frame(sample_data(lc_f_wl))
sampledf_f_gl <- data.frame(sample_data(lc_f_gl))

asv_f_cl <- data.frame(as(otu_table(lc_f_cl),"matrix"))
asv_f_wl <- data.frame(as(otu_table(lc_f_wl),"matrix"))
asv_f_gl <- data.frame(as(otu_table(lc_f_gl),"matrix"))

############################
#Filter for ASVs not present/only present in one site
ab1_anneli_cl <- rowSums(asv_anneli_cl)
asv_a2_anneli_cl <- cbind(asv_anneli_cl, ab1_anneli_cl)
asv_a2_anneli_cl <- as.data.frame(asv_a2_anneli_cl)
asv_a2_anneli_cl <- asv_a2_anneli_cl %>% dplyr::filter(asv_a2_anneli_cl[,305] > 1)
asv_a3_anneli_cl <- asv_a2_anneli_cl[,-305]
asv_a3_anneli_cl <- as.data.frame((asv_a3_anneli_cl))

ab1_anneli_wl <- rowSums(asv_anneli_wl)
asv_a2_anneli_wl <- cbind(asv_anneli_cl, ab1_anneli_wl)
asv_a2_anneli_wl <- as.data.frame(asv_a2_anneli_wl)
asv_a2_anneli_wl <- asv_a2_anneli_wl %>% dplyr::filter(asv_a2_anneli_wl[,305] > 1)
asv_a3_anneli_wl <- asv_a2_anneli_wl[,-305]
asv_a3_anneli_wl <- as.data.frame((asv_a3_anneli_wl))

ab1_anneli_gl <- rowSums(asv_anneli_gl)
asv_a2_anneli_gl <- cbind(asv_anneli_cl, ab1_anneli_gl)
asv_a2_anneli_gl <- as.data.frame(asv_a2_anneli_gl)
asv_a2_anneli_gl <- asv_a2_anneli_gl %>% dplyr::filter(asv_a2_anneli_gl[,305] > 1)
asv_a3_anneli_gl <- asv_a2_anneli_gl[,-305]
asv_a3_anneli_gl <- as.data.frame((asv_a3_anneli_gl))

#Filter for ASVs not present/only present in one site
ab1_arthro_cl <- rowSums(asv_arthro_cl)
asv_a2_arthro_cl <- cbind(asv_arthro_cl, ab1_arthro_cl)
asv_a2_arthro_cl <- as.data.frame(asv_a2_arthro_cl)
asv_a2_arthro_cl <- asv_a2_arthro_cl %>% dplyr::filter(asv_a2_arthro_cl[,277] > 1)
asv_a3_arthro_cl <- asv_a2_arthro_cl[,-277]
asv_a3_arthro_cl <- as.data.frame((asv_a3_arthro_cl))

ab1_arthro_wl <- rowSums(asv_arthro_wl)
asv_a2_arthro_wl <- cbind(asv_arthro_cl, ab1_arthro_wl)
asv_a2_arthro_wl <- as.data.frame(asv_a2_arthro_wl)
asv_a2_arthro_wl <- asv_a2_arthro_wl %>% dplyr::filter(asv_a2_arthro_wl[,277] > 1)
asv_a3_arthro_wl <- asv_a2_arthro_wl[,-277]
asv_a3_arthro_wl <- as.data.frame((asv_a3_arthro_wl))

ab1_arthro_gl <- rowSums(asv_arthro_gl)
asv_a2_arthro_gl <- cbind(asv_arthro_cl, ab1_arthro_gl)
asv_a2_arthro_gl <- as.data.frame(asv_a2_arthro_gl)
asv_a2_arthro_gl <- asv_a2_arthro_gl %>% dplyr::filter(asv_a2_arthro_gl[,277] > 1)
asv_a3_arthro_gl <- asv_a2_arthro_gl[,-277]
asv_a3_arthro_gl <- as.data.frame((asv_a3_arthro_gl))

#Filter for ASVs not present/only present in one site
ab1_roti_cl <- rowSums(asv_roti_cl)
asv_a2_roti_cl <- cbind(asv_roti_cl, ab1_roti_cl)
asv_a2_roti_cl <- as.data.frame(asv_a2_roti_cl)
asv_a2_roti_cl <- asv_a2_roti_cl %>% dplyr::filter(asv_a2_roti_cl[,168] > 1)
asv_a3_roti_cl <- asv_a2_roti_cl[,-168]
asv_a3_roti_cl <- as.data.frame((asv_a3_roti_cl))

ab1_roti_wl <- rowSums(asv_roti_wl)
asv_a2_roti_wl <- cbind(asv_roti_cl, ab1_roti_wl)
asv_a2_roti_wl <- as.data.frame(asv_a2_roti_wl)
asv_a2_roti_wl <- asv_a2_roti_wl %>% dplyr::filter(asv_a2_roti_wl[,168] > 1)
asv_a3_roti_wl <- asv_a2_roti_wl[,-168]
asv_a3_roti_wl <- as.data.frame((asv_a3_roti_wl))

ab1_roti_gl <- rowSums(asv_roti_gl)
asv_a2_roti_gl <- cbind(asv_roti_cl, ab1_roti_gl)
asv_a2_roti_gl <- as.data.frame(asv_a2_roti_gl)
asv_a2_roti_gl <- asv_a2_roti_gl %>% dplyr::filter(asv_a2_roti_gl[,168] > 1)
asv_a3_roti_gl <- asv_a2_roti_gl[,-168]
asv_a3_roti_gl <- as.data.frame((asv_a3_roti_gl))


#tardi
ab1_tardi_cl <- rowSums(asv_tardi_cl)
asv_a2_tardi_cl <- cbind(asv_tardi_cl, ab1_tardi_cl)
asv_a2_tardi_cl <- as.data.frame(asv_a2_tardi_cl)
asv_a2_tardi_cl <- asv_a2_tardi_cl %>% dplyr::filter(asv_a2_tardi_cl[,282] > 1)
asv_a3_tardi_cl <- asv_a2_tardi_cl[,-282]
asv_a3_tardi_cl <- as.data.frame((asv_a3_tardi_cl))

ab1_tardi_wl <- rowSums(asv_tardi_wl)
asv_a2_tardi_wl <- cbind(asv_tardi_cl, ab1_tardi_wl)
asv_a2_tardi_wl <- as.data.frame(asv_a2_tardi_wl)
asv_a2_tardi_wl <- asv_a2_tardi_wl %>% dplyr::filter(asv_a2_tardi_wl[,282] > 1)
asv_a3_tardi_wl <- asv_a2_tardi_wl[,-282]
asv_a3_tardi_wl <- as.data.frame((asv_a3_tardi_wl))

ab1_tardi_gl <- rowSums(asv_tardi_gl)
asv_a2_tardi_gl <- cbind(asv_tardi_cl, ab1_tardi_gl)
asv_a2_tardi_gl <- as.data.frame(asv_a2_tardi_gl)
asv_a2_tardi_gl <- asv_a2_tardi_gl %>% dplyr::filter(asv_a2_tardi_gl[,282] > 1)
asv_a3_tardi_gl <- asv_a2_tardi_gl[,-282]
asv_a3_tardi_gl <- as.data.frame((asv_a3_tardi_gl))


#nema
ab1_nema_cl <- rowSums(asv_nema_cl)
asv_a2_nema_cl <- cbind(asv_nema_cl, ab1_nema_cl)
asv_a2_nema_cl <- as.data.frame(asv_a2_nema_cl)
asv_a2_nema_cl <- asv_a2_nema_cl %>% dplyr::filter(asv_a2_nema_cl[,282] > 1)
asv_a3_nema_cl <- asv_a2_nema_cl[,-282]
asv_a3_nema_cl <- as.data.frame((asv_a3_nema_cl))

ab1_nema_wl <- rowSums(asv_nema_wl)
asv_a2_nema_wl <- cbind(asv_nema_cl, ab1_nema_wl)
asv_a2_nema_wl <- as.data.frame(asv_a2_nema_wl)
asv_a2_nema_wl <- asv_a2_nema_wl %>% dplyr::filter(asv_a2_nema_wl[,282] > 1)
asv_a3_nema_wl <- asv_a2_nema_wl[,-282]
asv_a3_nema_wl <- as.data.frame((asv_a3_nema_wl))

ab1_nema_gl <- rowSums(asv_nema_gl)
asv_a2_nema_gl <- cbind(asv_nema_cl, ab1_nema_gl)
asv_a2_nema_gl <- as.data.frame(asv_a2_nema_gl)
asv_a2_nema_gl <- asv_a2_nema_gl %>% dplyr::filter(asv_a2_nema_gl[,282] > 1)
asv_a3_nema_gl <- asv_a2_nema_gl[,-282]
asv_a3_nema_gl <- as.data.frame((asv_a3_nema_gl))


#Protists
ab1_p_cl <- rowSums(asv_p_cl)
asv_a2_p_cl <- cbind(asv_p_cl, ab1_p_cl)
asv_a2_p_cl <- as.data.frame(asv_a2_p_cl)
asv_a2_p_cl <- asv_a2_p_cl %>% dplyr::filter(asv_a2_p_cl[,349] > 1)
asv_a3_p_cl <- asv_a2_p_cl[,-349]
asv_a3_p_cl <- as.data.frame((asv_a3_p_cl))

ab1_p_wl <- rowSums(asv_p_wl)
asv_a2_p_wl <- cbind(asv_p_cl, ab1_p_wl)
asv_a2_p_wl <- as.data.frame(asv_a2_p_wl)
asv_a2_p_wl <- asv_a2_p_wl %>% dplyr::filter(asv_a2_p_wl[,349] > 1)
asv_a3_p_wl <- asv_a2_p_wl[,-349]
asv_a3_p_wl <- as.data.frame((asv_a3_p_wl))

ab1_p_gl <- rowSums(asv_p_gl)
asv_a2_p_gl <- cbind(asv_p_cl, ab1_p_gl)
asv_a2_p_gl <- as.data.frame(asv_a2_p_gl)
asv_a2_p_gl <- asv_a2_p_gl %>% dplyr::filter(asv_a2_p_gl[,349] > 1)
asv_a3_p_gl <- asv_a2_p_gl[,-349]
asv_a3_p_gl <- as.data.frame((asv_a3_p_gl))


#f
ab1_f_cl <- rowSums(asv_f_cl)
asv_a2_f_cl <- cbind(asv_f_cl, ab1_f_cl)
asv_a2_f_cl <- as.data.frame(asv_a2_f_cl)
asv_a2_f_cl <- asv_a2_f_cl %>% dplyr::filter(asv_a2_f_cl[,350] > 1)
asv_a3_f_cl <- asv_a2_f_cl[,-350]
asv_a3_f_cl <- as.data.frame((asv_a3_f_cl))

ab1_f_wl <- rowSums(asv_f_wl)
asv_a2_f_wl <- cbind(asv_f_cl, ab1_f_wl)
asv_a2_f_wl <- as.data.frame(asv_a2_f_wl)
asv_a2_f_wl <- asv_a2_f_wl %>% dplyr::filter(asv_a2_f_wl[,350] > 1)
asv_a3_f_wl <- asv_a2_f_wl[,-350]
asv_a3_f_wl <- as.data.frame((asv_a3_f_wl))

ab1_f_gl <- rowSums(asv_f_gl)
asv_a2_f_gl <- cbind(asv_f_cl, ab1_f_gl)
asv_a2_f_gl <- as.data.frame(asv_a2_f_gl)
asv_a2_f_gl <- asv_a2_f_gl %>% dplyr::filter(asv_a2_f_gl[,350] > 1)
asv_a3_f_gl <- asv_a2_f_gl[,-350]
asv_a3_f_gl <- as.data.frame((asv_a3_f_gl))


##################################
#ovrelapping taxa
#Fungi
GL_CL_f <- merge(asv_a3_f_gl,asv_a3_f_cl,by=0,all=F) #5641
nrow(GL_CL_f)
CL_WL_f <- merge(asv_a3_f_wl,asv_a3_f_cl,by=0,all=F) #5774
nrow(CL_WL_f)
GL_WL_f <- merge(asv_a3_f_gl,asv_a3_f_wl,by=0,all=F) #4865
nrow(GL_WL_f)
CL_WL_GL_f <- merge(CL_WL_f,GL_WL_f,by=0,all=F) #4865
nrow(CL_WL_GL_f)
GL_CL_WL_f <- merge(GL_CL_f,GL_WL_f,by=0,all=F) #4865
nrow(GL_CL_WL_f)
all_f <- merge(GL_CL_WL_f,CL_WL_GL_f,by=0,all=F) #4865
nrow(all_f)

nrow(asv_a3_f_cl)#CL:13308
nrow(asv_a3_f_wl)#WL: 10906
nrow(asv_a3_f_gl)#GL:8708
total <- (nrow(asv_a3_f_gl)+nrow(asv_a3_f_cl)+nrow(asv_a3_f_wl))-(nrow(GL_CL_f)+nrow(CL_WL_f))
overlapping <- 100/total*nrow(all_f) #22.62%

svg("Fig3_f.svg",width=10,height=6)
triple_f <- draw.triple.venn(area1=13308, area2=10906,area3=8708,n12=5774,n23=4865,n13=5641,n123=4865,category=c("CL","WL","GL"),fill=c("#ffff00","red","green"),cat.dist=0.1,cat.cex=2,cat.fontface=2,cat.fontfamily=2,fontface=2,fontfamily=2,cex=2,alpha=0.5,scaled=T)
dev.off() # Close the graphics device

#Protist
GL_CL_p <- merge(asv_a3_p_gl,asv_a3_p_cl,by=0,all=F) 
nrow(GL_CL_p) #8397
CL_WL_p <- merge(asv_a3_p_wl,asv_a3_p_cl,by=0,all=F) 
nrow(CL_WL_p)#8682
GL_WL_p <- merge(asv_a3_p_gl,asv_a3_p_wl,by=0,all=F)
nrow(GL_WL_p) #7105
GL_WL_CL_p <- merge(GL_WL_p,GL_CL_p,by=0,all=F)
nrow(GL_WL_CL_p)#7105
CL_WL_GL_p <- merge(CL_WL_p,GL_CL_p,by=0,all=F) 
nrow(CL_WL_GL_p) #8397
all_p <- merge(GL_WL_CL_p,CL_WL_GL_p,by=0,all=F) 
nrow(all_p)#7105
nrow(asv_a3_p_cl)#CL:18188
nrow(asv_a3_p_wl)#WL: 15391
nrow(asv_a3_p_gl)#GL:13571
total <- (nrow(asv_a3_p_gl)+nrow(asv_a3_p_cl)+nrow(asv_a3_p_wl))-(nrow(GL_CL_p)+nrow(CL_WL_p))
overlapping <- 100/total*nrow(all_p)
#overlapping:24.07% 

svg("Fig3_p.svg",width=10,height=6)
triple_p <- draw.triple.venn(area1=18188, area2=15391,area3=13571,n12=8682,n23=7105,n13=8397,n123=7105,category=c("CL","WL","GL"),fill=c("#ffff00","red","green"),cat.dist=0.1,cat.cex=2,cat.fontface=2,cat.fontfamily=2,fontface=2,fontfamily=2,cex=2,alpha=0.5,scaled=T)
dev.off() # Close the graphics device




######################
#tardigrades
GL_CL_tardi <- merge(asv_a3_tardi_gl,asv_a3_tardi_cl,by=0,all=F) 
nrow(GL_CL_tardi) #98
CL_WL_tardi <- merge(asv_a3_tardi_wl,asv_a3_tardi_cl,by=0,all=F) 
nrow(CL_WL_tardi) #97
GL_WL_tardi <- merge(asv_a3_tardi_gl,asv_a3_tardi_wl,by=0,all=F) 
nrow(GL_WL_tardi) #84
GL_WL_CL_tardi <- merge(GL_WL_tardi,GL_CL_tardi,by=0,all=F)
nrow(GL_WL_CL_tardi) #84
CL_WL_GL_tardi <- merge(CL_WL_tardi,GL_CL_tardi,by=0,all=F) 
nrow(CL_WL_GL_tardi) #97
all_tardi <- merge(GL_WL_CL_tardi,CL_WL_GL_tardi,by=0,all=F) 
nrow(all_tardi)#84
nrow(asv_a3_tardi_cl)#CL:143
nrow(asv_a3_tardi_wl)#WL: 137
nrow(asv_a3_tardi_gl)#GL:119
total <- (nrow(asv_a3_tardi_gl)+nrow(asv_a3_tardi_cl)+nrow(asv_a3_tardi_wl))-(nrow(GL_CL_tardi)+nrow(CL_WL_tardi))
overlapping <- (100/total)*nrow(all_tardi)
#overlapping:41.18% 


svg("Fig3_tardi.svg",width=10,height=6)
triple_tardi <- draw.triple.venn(area1=143, area2=137,area3=119,n12=97,n23=84,n13=98,n123=84,category=c("CL","WL","GL"),fill=c("#ffff00","red","green"),cat.dist=0.1,cat.cex=2,cat.fontface=2,cat.fontfamily=2,fontface=2,fontfamily=2,cex=2,alpha=0.5,scaled=T)
dev.off() # Close the graphics device


######################
#annelids
GL_CL_anneli <- merge(asv_a3_anneli_gl,asv_a3_anneli_cl,by=0,all=F)
nrow(GL_CL_anneli)#31
CL_WL_anneli <- merge(asv_a3_anneli_wl,asv_a3_anneli_cl,by=0,all=F) 
nrow(CL_WL_anneli)#31
GL_WL_anneli <- merge(asv_a3_anneli_gl,asv_a3_anneli_wl,by=0,all=F) 
nrow(GL_WL_anneli)#30
GL_WL_CL_anneli <- merge(GL_WL_anneli,GL_CL_anneli,by=0,all=F) 
nrow(GL_WL_CL_anneli)#30
CL_WL_GL_anneli <- merge(CL_WL_anneli,GL_CL_anneli,by=0,all=F) 
nrow(CL_WL_GL_anneli)#31
all_anneli <- merge(GL_WL_CL_anneli,CL_WL_GL_anneli,by=0,all=F) 
nrow(all_anneli)#30
nrow(asv_a3_anneli_cl)#CL:39
nrow(asv_a3_anneli_wl)#WL: 37
nrow(asv_a3_anneli_gl)#GL:33
total <- (nrow(asv_a3_anneli_gl)+nrow(asv_a3_anneli_cl)+nrow(asv_a3_anneli_wl))-(nrow(GL_CL_anneli)+nrow(CL_WL_anneli))
overlapping <- (100/total)*nrow(all_anneli)
#overlapping: 63.83%


svg("Fig3_anneli.svg",width=10,height=6)
triple_anneli <- draw.triple.venn(area1=39, area2=37,area3=33,n12=31,n23=30,n13=31,n123=30,category=c("WL","GL","CL"),fill=c("red","green","#ffff00"),cat.cex=2,cat.fontface=2,cat.fontfamily=2,cat.dist=0.1,fontface=2,fontfamily=2,cex=2,scaled=T)
dev.off() # Close the graphics device



######################
#nematodes
GL_CL_nema <- merge(asv_a3_nema_gl,asv_a3_nema_cl,by=0,all=F)
nrow(GL_CL_nema) #1056
CL_WL_nema <- merge(asv_a3_nema_wl,asv_a3_nema_cl,by=0,all=F) 
nrow(CL_WL_nema) #1086
GL_WL_nema <- merge(asv_a3_nema_gl,asv_a3_nema_wl,by=0,all=F) 
nrow(GL_WL_nema) #898
CL_WL_GL_nema <- merge(CL_WL_nema,GL_CL_nema,by=0,all=F) 
nrow(CL_WL_GL_nema)#1056
GL_WL_CL_nema <- merge(GL_WL_nema,GL_CL_nema,by=0,all=F) 
nrow(CL_WL_GL_nema) #1056
all_nema <- merge(GL_WL_CL_nema,CL_WL_GL_nema,by=0,all=F) 
nrow(all_nema) #898
nrow(asv_a3_nema_cl)#CL:2025
nrow(asv_a3_nema_wl)#WL: 1655
nrow(asv_a3_nema_gl)#GL:1382
total <- (nrow(asv_a3_nema_gl)+nrow(asv_a3_nema_cl)+nrow(asv_a3_nema_wl))-(nrow(GL_CL_nema)+nrow(CL_WL_nema))
overlapping <- (100/total)*nrow(all_nema)
#overlapping:30.75%

svg("Fig3_nema.svg",width=10,height=6)
triple_nema <- draw.triple.venn(area1=2025, area2=1655,area3=1382,n12=1086,n23=898,n13=1056,n123=898,category=c("CL","WL","GL"),fill=c("#ffff00","red","green"),cat.dist=0.1,cat.cex=2,cat.fontface=2,cat.fontfamily=2,fontface=2,fontfamily=2,cex=2,alpha=0.5,scaled=T)
dev.off() # Close the graphics device



####################
#arthropods
GL_CL_arthro <- merge(asv_a3_arthro_gl,asv_a3_arthro_cl,by=0,all=F) 
nrow(GL_CL_arthro)#767
CL_WL_arthro <- merge(asv_a3_arthro_wl,asv_a3_arthro_cl,by=0,all=F) 
nrow(CL_WL_arthro)#784
GL_WL_arthro <- merge(asv_a3_arthro_gl,asv_a3_arthro_wl,by=0,all=F) 
nrow(GL_WL_arthro)#665
CL_WL_GL_arthro <- merge(CL_WL_arthro,GL_CL_arthro,by=0,all=F) 
nrow(CL_WL_GL_arthro)#767
GL_WL_CL_arthro <- merge(GL_WL_arthro,GL_CL_arthro,by=0,all=F) 
nrow(CL_WL_GL_arthro) #767
all_arthro <- merge(GL_WL_CL_arthro,CL_WL_GL_arthro,by=0,all=F) 
nrow(all_arthro) #665
nrow(asv_a3_arthro_cl)#CL:1341
nrow(asv_a3_arthro_wl)#WL:1130
nrow(asv_a3_arthro_gl)#GL:980
total <- (nrow(asv_a3_arthro_gl)+nrow(asv_a3_arthro_cl)+nrow(asv_a3_arthro_wl))-(nrow(GL_CL_arthro)+nrow(CL_WL_arthro))
overlapping <- (100/total)*nrow(all_arthro)
#overlapping:35%

svg("Fig3_arthro.svg",width=10,height=6)
triple_arthro <- draw.triple.venn(area1=1341, area2=1130,area3=980,n12=784,n23=665,n13=767,n123=665,category=c("CL","WL","GL"),fill=c("#ffff00","red","green"),cat.dist=0.1,cat.cex=2,cat.fontface=2,cat.fontfamily=2,fontface=2,fontfamily=2,cex=2,alpha=0.5,scaled=T)
# dev.off() # Close the graphics device


######################
#rotifers
GL_CL_roti <- merge(asv_a3_roti_gl,asv_a3_roti_cl,by=0,all=F)
nrow(GL_CL_roti) #156
CL_WL_roti <- merge(asv_a3_roti_wl,asv_a3_roti_cl,by=0,all=F) 
nrow(CL_WL_roti) #159
GL_WL_roti <- merge(asv_a3_roti_gl,asv_a3_roti_wl,by=0,all=F)
nrow(GL_WL_roti) #112
CL_WL_GL_roti <- merge(CL_WL_roti,GL_CL_roti,by=0,all=F) 
nrow(GL_WL_roti) #112
GL_WL_CL_roti <- merge(GL_WL_roti,GL_CL_roti,by=0,all=F) 
nrow(GL_WL_CL_roti) #112
all_roti <- merge(GL_WL_CL_roti,CL_WL_GL_roti,by=0,all=F) 
nrow(all_roti) #112
nrow(asv_a3_roti_cl)#CL:371
nrow(asv_a3_roti_wl)#WL:272
nrow(asv_a3_roti_gl)#GL:216
total <- (nrow(asv_a3_roti_gl)+nrow(asv_a3_roti_cl)+nrow(asv_a3_roti_wl))-(nrow(GL_CL_roti)+nrow(CL_WL_roti))
overlapping <- (100/total)*nrow(all_roti)
#overlapping:30.59%

svg("Fig3_roti.svg",width=10,height=6)
triple_roti <- draw.triple.venn(area1=371, area2=272,area3=216,n12=159,n23=112,n13=156,n123=112,category=c("CL","WL","GL"),fill=c("#ffff00","red","green"),cat.dist=0.1,cat.cex=2,cat.fontface=2,cat.fontfamily=2,fontface=2,fontfamily=2,cex=2,alpha=0.5,scaled=T)
dev.off() # Close the graphics device



