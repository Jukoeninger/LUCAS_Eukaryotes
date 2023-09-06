#############################################
#Alpha diversity, Figure 3

#Set graphical themes
theme_jk2 <- function(){
  theme_bw() +
    theme(text = element_text(size=12, family = "sans"),
          axis.text = element_text(size = 8,angle = 90), 
          axis.title = element_text(size = 12, vjust = -3),
          plot.margin = unit(c(1, 1, 1, 1), units =, "cm"),
          plot.title = element_text(size = 12, vjust = 1, hjust = 0),
          legend.text = element_text(size = 8), 
          legend.key.size = unit(1,"cm"),
          legend.key.width= unit(0, "cm"),
          legend.title = element_text(size=10),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 0.5, linetype = "blank"))
}

theme_jk3 <- function(){
  theme_bw() +
    theme(text = element_text(size=8, family = "sans"),
          axis.text.x =  element_blank(),          
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), units =, "cm"),
          axis.title = element_blank(),
          plot.title = element_text(size = 8, vjust = 1, hjust = 0),
          legend.text = element_text(size = 8),legend.title=element_blank())
}

#with horizontal axis.text/for alpha diversity
theme_jk4 <- function(){
  theme_bw() +
    theme(text = element_text(size=14, family = "sans"),
          axis.text = element_text(size = 16), 
          axis.title.y = element_text(size = 16), 
          axis.text.y = element_text(size = 16),
          axis.title.x = element_blank(),
          plot.margin = unit(c(1, 1, 1, 1), units =, "cm"),
          plot.title = element_text(size = 18, vjust = 1, hjust = 0),
          legend.text = element_text(size = 18), 
          legend.key.size = unit(0.5,"cm"),
          legend.key.width= unit(0.5, "cm"),
          legend.title = element_text(size=18),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 0.5, linetype = "blank"))
}

#Alpha diversity 
#ONE TOTAL DATAFRAME WITH SHANNON INDEX/RICHNESS
#extract shannon index, richness, ecosystem type
tardi_div_obs <- as.data.frame(a_tardi_sample_data[,c("Richness","BARCODE_ID","LC1_2018")])
nema_div_obs <- as.data.frame(a_nema_sample_data[,c("Richness","BARCODE_ID","LC1_2018")])
anneli_div_obs <- as.data.frame(a_anneli_sample_data[,c("Richness","BARCODE_ID","LC1_2018")])
arthro_div_obs <- as.data.frame(a_arthro_sample_data[,c("Richness","BARCODE_ID","LC1_2018")])
p_div_obs <- as.data.frame(p_sample_data[,c("Richness","BARCODE_ID","LC1_2018")])
f_div_obs <- as.data.frame(f_sample_data[,c("Richness","BARCODE_ID","LC1_2018")])
roti_div_obs <- as.data.frame(a_roti_sample_data[,c("Richness","BARCODE_ID","LC1_2018")])

#add column with eukaryotic group
tardi_div_obs$group <-"Tardigrades"
anneli_div_obs$group <-"Annelids"
arthro_div_obs$group <-"Arthropods"
nema_div_obs$group <-"Nematodes"
p_div_obs$group <-"Protists"
f_div_obs$group <-"Fungi"
roti_div_obs$group <-"Rotifers"

tardi_div_H <- as.data.frame(a_tardi_sample_data[,c("Shannon","BARCODE_ID","LC1_2018")])
arthro_div_H <- as.data.frame(a_arthro_sample_data[,c("Shannon","BARCODE_ID","LC1_2018")])
anneli_div_H <- as.data.frame(a_anneli_sample_data[,c("Shannon","BARCODE_ID","LC1_2018")])
nema_div_H <- as.data.frame(a_nema_sample_data[,c("Shannon","BARCODE_ID","LC1_2018")])
p_div_H <- as.data.frame(p_sample_data[,c("Shannon","BARCODE_ID","LC1_2018")])
f_div_H <- as.data.frame(f_sample_data[,c("Shannon","BARCODE_ID","LC1_2018")])
roti_div_H <- as.data.frame(a_roti_sample_data[,c("Shannon","BARCODE_ID","LC1_2018")])

#add column with eukaryotic group
tardi_div_H$group <-"Tardigrades"
anneli_div_H$group <-"Annelids"
arthro_div_H$group <-"Arthropods"
nema_div_H$group <-"Nematodes"
p_div_H$group <-"Protists"
f_div_H$group <-"Fungi"
roti_div_H$group <-"Rotifers"


#single plots with Kruskal-Wallis letters
#Richness
kk_a_arthro_groups_LC1 <- kruskal(a_arthro_sample_data$Richness,a_arthro_sample_data$LC1_2018,group=TRUE,p.adj="bonferroni")
kk_a_arthro_groups_LC1
kk_a_anneli_groups_LC1 <- kruskal(a_anneli_sample_data$Richness,a_anneli_sample_data$LC1_2018,group=TRUE,p.adj="bonferroni")
kk_a_anneli_groups_LC1
kk_a_nema_groups_LC1 <- kruskal(a_nema_sample_data$Richness,a_nema_sample_data$LC1_2018,group=TRUE,p.adj="bonferroni")
kk_a_nema_groups_LC1
kk_a_roti_groups_LC1 <- kruskal(a_roti_sample_data$Richness,a_roti_sample_data$LC1_2018,group=TRUE,p.adj="bonferroni")
kk_a_roti_groups_LC1
kk_a_tardi_groups_LC1 <- kruskal(a_tardi_sample_data$Richness,a_tardi_sample_data$LC1_2018,group=TRUE,p.adj="bonferroni")
kk_a_tardi_groups_LC1
kk_f_groups_LC1 <- kruskal(f_sample_data$Richness,f_sample_data$LC1_2018,group=TRUE,p.adj="bonferroni")
kk_f_groups_LC1
kk_p_groups_LC1 <- kruskal(p_sample_data$Richness,p_sample_data$LC1_2018,group=TRUE,p.adj="bonferroni")
kk_p_groups_LC1

obs_div <- rbind(tardi_div_obs,roti_div_obs,nema_div_obs,roti_div_obs,arthro_div_obs,anneli_div_obs,p_div_obs,f_div_obs)

#rename ecosystem type
obs_div$LC1_2018[obs_div$LC1_2018 =="CL_annual"] <- "CL1"
obs_div$LC1_2018[obs_div$LC1_2018 =="CL_permanent"] <- "CL2"
obs_div$LC1_2018[obs_div$LC1_2018 =="GL_managed"] <- "GL1"
obs_div$LC1_2018[obs_div$LC1_2018 =="GL_unmanaged"] <- "GL2"
obs_div$LC1_2018[obs_div$LC1_2018 =="WL_broadleaved"] <-"WL1" 
obs_div$LC1_2018[obs_div$LC1_2018 =="WL_coniferous"] <- "WL2"
#set order of variables
obs_div$LC1_2018 <- factor(obs_div$LC1_2018, levels=c("CL1", "CL2", "GL1", "GL2","WL1", "WL2"))
obs_div$group <- factor(obs_div$group, levels=c("Fungi","Protists","Rotifers","Tardigrades","Nematodes","Arthropods","Annelids"))
obs_div_f_p <- obs_div %>% filter(group==c("Fungi","Protists"))
la_f_p2 <- c("a","ab","bc","bcd","cd","d","a","ab","b","bc","bc","c")
obs_value_max1 = obs_div_f_p %>% group_by(LC1_2018) %>% summarize(max_value = max(Richness))
la_animals2 <- c("a","ab","ab","ab","ab","b","a","ab","ab","ab","b","b","a","ab","ab","b","c","c","a","ab","bc","bcd","cd","d","a","ab","b","b","b","b")
obs_div_animals <- obs_div %>% filter(group==c("Rotifers","Tardigrades","Nematodes","Arthropods","Annelids"))
obs_value_max2 = obs_div_animals %>% group_by(LC1_2018) %>% summarize(max_value = max(Richness))

#PLOT
svg("1312_Fig3A.svg")# SVG graphics device
obs_div_f_p_plot2 <-  ggplot(obs_div_f_p, aes(x = LC1_2018, y = Richness)) +
  geom_boxplot(aes(fill =LC1_2018),width = 0.7, outlier.shape = NA, alpha = 0.7,fatten = 6) +
  geom_richtext(data = obs_value_max1, aes(x=LC1_2018, y = 490, label = la_f_p2), vjust=-0.3,fill = "#D3D3D3", label.color = NA,size=4.5)+
  facet_grid(~group)+
  stat_compare_means(aes(label = paste0("  p = ", p.format,after_stat(p.signif)), vjust=-0.1,hjust=0.01), label.y=580)+
  scale_colour_manual(name='Ecosystem type',values = c("#FCE205","orange","#2AAA8A","#98FB98","#FE036A","#7d1969")) +
  scale_fill_manual(name='Ecosystem type',values = c("#FCE205","orange","#2AAA8A","#98FB98","#FE036A","#7d1969")) +
  labs(y="ASV richness") +
  ylim(0,595)+
  theme_jk4()+theme(axis.text.x=element_blank(),
                    axis.ticks.x=element_blank(),
                                        axis.text = element_text(size=12),
                                        axis.title = element_text(size=14),
                                        strip.text.x = element_text(size = 16))
dev.off() # Close the graphics device

#PLOT
svg("1312_Fig3B.svg")# SVG graphics device
obs_div_animalia_plot2 <- ggplot(obs_div_animals, aes(x = LC1_2018, y = Richness)) +
  geom_boxplot(aes(fill=LC1_2018),width = 0.7, outlier.shape = NA, alpha = 0.7,fatten = 6) +
  geom_richtext(data = obs_value_max2, aes(x=LC1_2018, y = 75, label = la_animals2),vjust=-0.3,fill = "#D3D3D3", label.color = NA,size=3.5)+
  facet_grid(~group)+
  scale_colour_manual(name='Ecosystem type',values = c("#FCE205","orange","#2AAA8A","#98FB98","#FE036A","#7d1969")) +
  scale_fill_manual(name='Ecosystem type',values = c("#FCE205","orange","#2AAA8A","#98FB98","#FE036A","#7d1969")) +
  stat_compare_means(aes(label = paste0("   p = ", p.format,after_stat(p.signif)), hjust=0.2),label.y=88)+
  labs(color='Ecosystem type',y="ASV richness") +
  ylim(0,90)+
  theme_jk4()+theme(axis.text.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    axis.text = element_text(size=14),
                    axis.title = element_text(size=16),
                    strip.text.x = element_text(size = 12))
dev.off() # Close the graphics device


#Shannon index
kk_a_arthro_groups_LC1_H <- kruskal(a_arthro_sample_data$Shannon,a_arthro_sample_data$LC1_2018,group=TRUE,p.adj="bonferroni")
kk_a_arthro_groups_LC1_H
kk_a_anneli_groups_LC1_H <- kruskal(a_anneli_sample_data$Shannon,a_anneli_sample_data$LC1_2018,group=TRUE,p.adj="bonferroni")
kk_a_anneli_groups_LC1_H
kk_a_nema_groups_LC1_H <- kruskal(a_nema_sample_data$Shannon,a_nema_sample_data$LC1_2018,group=TRUE,p.adj="bonferroni")
kk_a_nema_groups_LC1_H
kk_a_tardi_groups_LC1_H <- kruskal(a_tardi_sample_data$Shannon,a_tardi_sample_data$LC1_2018,group=TRUE,p.adj="bonferroni")
kk_a_tardi_groups_LC1_H
kk_a_roti_groups_LC1_H <- kruskal(a_roti_sample_data$Shannon,a_roti_sample_data$LC1_2018,group=TRUE,p.adj="bonferroni")
kk_a_roti_groups_LC1_H
kk_f_groups_LC1_H <- kruskal(f_sample_data$Shannon,f_sample_data$LC1_2018,group=TRUE,p.adj="bonferroni")
kk_f_groups_LC1_H
kk_p_groups_LC1_H <- kruskal(p_sample_data$Shannon,p_sample_data$LC1_2018,group=TRUE,p.adj="bonferroni")
kk_p_groups_LC1_H

H_div <- rbind(tardi_div_H,roti_div_H,nema_div_H,arthro_div_H,anneli_div_H,p_div_H,f_div_H)

#rename ecosystem type
H_div$LC1_2018[H_div$LC1_2018 =="CL_annual"] <- "CL1"
H_div$LC1_2018[H_div$LC1_2018 =="CL_permanent"] <- "CL2"
H_div$LC1_2018[H_div$LC1_2018 =="GL_managed"] <- "GL1"
H_div$LC1_2018[H_div$LC1_2018 =="GL_unmanaged"] <- "GL2"
H_div$LC1_2018[H_div$LC1_2018 =="WL_broadleaved"] <-"WL1" 
H_div$LC1_2018[H_div$LC1_2018 =="WL_coniferous"] <- "WL2"

#rearrange order
H_div$LC1_2018 <- factor(H_div$LC1_2018, levels=c("CL1", "CL2", "GL1", "GL2","WL1", "WL2"))
H_div$group <- factor(H_div$group, levels=c("Fungi","Protists","Rotifers","Tardigrades","Nematodes","Arthropods","Annelids"))
H_div_f_p <- H_div %>% filter(group==c("Fungi","Protists"))
H_div_animals <- H_div %>% filter(group==c("Rotifers","Tardigrades","Nematodes","Arthropods","Annelids"))
la_f_p2_H <- c("a","a","ab","ab","b","b","a","a","ab","abc","bc","c")
la_animals2_H <- c("a","a","a","a","a","a","a","ab","b","b","b","b","a","ab","ab","b","b","b","a","ab","bc","bc","bc","c","a","ab","b","b","b","b")
H_value_max1 = H_div_f_p %>% group_by(LC1_2018) %>% summarize(max_value = max(Shannon))
H_value_max2 = H_div_animals %>% group_by(LC1_2018) %>% summarize(max_value = max(Shannon))

#PLOT
svg("1312_Fig3C.svg")# SVG graphics device
H_div_f_p_plot <- ggplot(H_div_f_p, aes(x = LC1_2018, y = Shannon)) +
  geom_boxplot(aes(fill =LC1_2018),width = 0.7, outlier.shape = NA, alpha = 0.7,fatten = 6) +
  geom_richtext(data = H_value_max1, aes(x=LC1_2018, y = 5.5, label = la_f_p2_H),vjust=-0.3,fill = "#D3D3D3", label.color = NA,size=4.5)+
  facet_grid(~group)+
  scale_fill_manual(name='Ecosystem type',values = c("#FCE205","orange","#2AAA8A","#98FB98","#FE036A","#7d1969")) +
  scale_colour_manual(name='Ecosystem type',values = c("#FCE205","orange","#2AAA8A","#98FB98","#FE036A","#7d1969")) +   
  ylim(0, 7) +
  stat_compare_means(aes(label = paste0("  p = ", p.format,after_stat(p.signif)) ,vjust=-0.1,hjust=0.01), label.y=6.7)+
  labs(y="Shannon index") +
  theme_jk4()+theme(axis.text.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    axis.text = element_text(size=15),
                    axis.title = element_text(size=17),
                    strip.text.x = element_text(size = 16))
dev.off() # Close the graphics device

#PLOT
svg("1312_Fig3D.svg")# SVG graphics device
H_div_animalia_plot <- ggplot(H_div_animals, aes(x = LC1_2018, y = Shannon)) +
  geom_boxplot(aes(fill =LC1_2018),width = 0.7, outlier.shape = NA, alpha = 0.7,fatten = 6) +
  geom_richtext(data = H_value_max2, aes(x=LC1_2018, y = 3.5, label = la_animals2_H), vjust=-0.3,fill = "#D3D3D3", label.color = NA,size=3.5)+
  facet_grid(~group)+
  scale_fill_manual(name='Ecosystem type',values = c("#FCE205","orange","#2AAA8A","#98FB98","#FE036A","#7d1969")) +
  scale_colour_manual(name='Ecosystem type',values = c("#FCE205","orange","#2AAA8A","#98FB98","#FE036A","#7d1969")) +
  ylim(0,4.2)+
  stat_compare_means(aes(label = paste0("   p = ", p.format,after_stat(p.signif))),label.y=4.1)+
  labs(y="Shannon index") +
  theme_jk4()+theme(axis.text.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    axis.text = element_text(size=15),
                    axis.title = element_text(size=17),
                    strip.text.x = element_text(size = 12))
dev.off() # Close the graphics device


#combine plots, Figure 4
svg("1505_Fig3.svg",width=17.5,height=9)# SVG graphics device
ggarrange(obs_div_f_p_plot2,obs_div_animalia_plot2,H_div_f_p_plot,H_div_animalia_plot,ncol = 2,nrow=2,
          labels = c( "(a)", "(b)","(c)","(d)"),common.legend = TRUE,legend="bottom")
dev.off() # Close the graphics device










###########################################
###########################################
###########################################
#TESTING

#Test for normal distribution
shapiro.test(p_sample_data$Richness) #normal distribution
shapiro.test(p_sample_data$Shannon) #no normal distribution

shapiro.test(f_sample_data$Richness) #normal distribution
shapiro.test(f_sample_data$Shannon) #no normal distribution

shapiro.test(a_tardi_sample_data$Richness) #no normal distribution
shapiro.test(a_tardi_sample_data$Shannon) #no normal distribution

shapiro.test(a_roti_sample_data$Richness) #no normal distribution
shapiro.test(a_roti_sample_data$Shannon) #no no rmal distribution

shapiro.test(a_nema_sample_data$Richness) #no normal distribution
shapiro.test(a_nema_sample_data$Shannon) #no normal distribution

shapiro.test(a_anneli_sample_data$Richness) #no normal distribution
shapiro.test(a_anneli_sample_data$Shannon) #no normal distribution

shapiro.test(a_arthro_sample_data$Richness) #no normal distribution
shapiro.test(a_arthro_sample_data$Shannon) #no normal distribution


#Kruskal-Wallis test +post-hoc Dunn test

#Protists
#Richness
kruskal.test(Richness ~ LC1_2018, data= p_sample_data)
kk_p_groups_LC1 <- kruskal(p_sample_data$Richness,p_sample_data$LC1_2018,group=TRUE,p.adj="bonferroni")
kk_p_groups_LC1
dunn_Richness_p_LC1 <- kwAllPairsDunnTest(Richness ~ as.factor(LC1_2018), data= p_sample_data, p.adjust.method = "BH")
dunn_Richness_p_LC1
table_dunn_p_richness_LC1 <- dunn.test(p_sample_data$Richness, g= p_sample_data$LC1_2018, method = "bh",list=TRUE)
dunn_richness_p_LC1_df <- as.data.frame(table_dunn_p_richness_LC1)
dunn_richness_p_LC1_df$Variable <- "LC1"


kruskal.test(Richness ~ LU2, data= p_sample_data)
kk_p_groups_LU2 <- kruskal(p_sample_data$Richness,p_sample_data$LU2,group=TRUE,p.adj="bonferroni")
kk_p_groups_LU2
dunn_Richness_p_LU2 <- kwAllPairsDunnTest(Richness ~ as.factor(LU2), data= p_sample_data, p.adjust.method = "BH")
dunn_Richness_p_LU2
table_dunn_p_richness_LU2 <- dunn.test(p_sample_data$Richness, g= p_sample_data$LU2, method = "bh",list=TRUE)
dunn_richness_p_LU2_df <- as.data.frame(table_dunn_p_richness_LU2)
dunn_richness_p_LU2_df$Variable <- "LU2"

kruskal.test(Richness ~ Sample_season, data= p_sample_data)
kk_p_groups_Season <- kruskal(p_sample_data$Richness,p_sample_data$Sample_season,group=TRUE,p.adj="bonferroni")
kk_p_groups_Season
dunn_Richness_p_season <- kwAllPairsDunnTest(Richness ~ as.factor(Sample_season), data= p_sample_data, p.adjust.method = "BH")
dunn_Richness_p_season
table_dunn_p_richness_season <- dunn.test(p_sample_data$Richness, g= p_sample_data$Sample_season, method = "bh",list=TRUE)
dunn_richness_p_season_df <- as.data.frame(table_dunn_p_richness_season)
dunn_richness_p_season_df$Variable <- "Season"

kruskal.test(Richness ~ pH_grouped, data= p_sample_data)
kk_p_groups_pH <- kruskal(p_sample_data$Richness,p_sample_data$pH_grouped,group=TRUE,p.adj="bonferroni")
kk_p_groups_pH
dunn_Richness_p_pH <- kwAllPairsDunnTest(Richness ~ as.factor(pH_grouped), data= p_sample_data, p.adjust.method = "BH")
dunn_Richness_p_pH
table_dunn_p_richness_pH <- dunn.test(p_sample_data$Richness, g= p_sample_data$pH_grouped, method = "bh",list=TRUE)
dunn_richness_p_pH_df <- as.data.frame(table_dunn_p_richness_pH)
dunn_richness_p_pH_df$Variable <- "pH"

kruskal.test(Richness ~ depth_grouped, data= p_sample_data)
kk_p_groups_depth <- kruskal(p_sample_data$Richness,p_sample_data$depth_grouped,group=TRUE,p.adj="bonferroni")
kk_p_groups_depth
dunn_Richness_p_depth <- kwAllPairsDunnTest(Richness ~ as.factor(depth_grouped), data= p_sample_data, p.adjust.method = "BH")
dunn_Richness_p_depth
table_dunn_p_richness_depth <- dunn.test(p_sample_data$Richness, g= p_sample_data$depth_grouped, method = "bh",list=TRUE)
dunn_richness_p_depth_df <- as.data.frame(table_dunn_p_richness_depth)
dunn_richness_p_depth_df$Variable <- "Depth"

kruskal.test(Richness ~ Erosion_grouped, data= p_sample_data)
kk_p_groups_erosion <- kruskal(p_sample_data$Richness,p_sample_data$Erosion_grouped,group=TRUE,p.adj="bonferroni")
kk_p_groups_erosion
dunn_Richness_p_erosion <- kwAllPairsDunnTest(Richness ~ as.factor(Erosion_grouped), data= p_sample_data, p.adjust.method = "BH")
dunn_Richness_p_erosion
table_dunn_p_richness_erosion <- dunn.test(p_sample_data$Richness, g= p_sample_data$Erosion_grouped, method = "bh",list=TRUE)
dunn_richness_p_erosion_df <- as.data.frame(table_dunn_p_richness_erosion)
dunn_richness_p_erosion_df$Variable <- "Erosion"

kruskal.test(Richness ~  cluster_temp_range, data= p_sample_data)
kk_p_groups_cluster_temp_range <- kruskal(p_sample_data$Richness,p_sample_data$cluster_temp_range,group=TRUE,p.adj="bonferroni")
kk_p_groups_cluster_temp_range
dunn_Richness_p_cluster_temp_range <- kwAllPairsDunnTest(Richness ~ as.factor(cluster_temp_range), data= p_sample_data, f.adjust.method = "BH")
dunn_Richness_p_cluster_temp_range
table_dunn_p_richness_cluster_temp_range <- dunn.test(p_sample_data$Richness, g= p_sample_data$cluster_temp_range, method = "bh",list=TRUE)
dunn_p_richness_cluster_temp_range_df <- as.data.frame(table_dunn_p_richness_cluster_temp_range)
dunn_p_richness_cluster_temp_range_df 
dunn_p_richness_cluster_temp_range_df$Variable <- "Temp_range"

kruskal.test(Richness ~  cluster_temp, data= p_sample_data)
kk_p_groups_cluster_temp <- kruskal(p_sample_data$Richness,p_sample_data$cluster_temp,group=TRUE,p.adj="bonferroni")
kk_p_groups_cluster_temp
dunn_Richness_p_cluster_temp <- kwAllPairsDunnTest(Richness ~ as.factor(cluster_temp), data= p_sample_data, f.adjust.method = "BH")
dunn_Richness_p_cluster_temp
table_dunn_p_richness_cluster_temp <- dunn.test(p_sample_data$Richness, g= p_sample_data$cluster_temp, method = "bh",list=TRUE)
dunn_p_richness_cluster_temp_df <- as.data.frame(table_dunn_p_richness_cluster_temp)
dunn_p_richness_cluster_temp_df$Variable <- "Temp"

kruskal.test(Richness ~  cluster_prec, data= p_sample_data)
#not significant

#Shannon
kruskal.test(Shannon ~ LC1_2018, data= p_sample_data)
kk_H_p_groups_LC1 <- kruskal(p_sample_data$Shannon,p_sample_data$LC1_2018,group=TRUE,p.adj="bonferroni")
kk_H_p_groups_LC1
dunn_shannon_p_LC1 <- kwAllPairsDunnTest(Shannon ~ as.factor(LC1_2018), data= p_sample_data, p.adjust.method = "BH")
dunn_shannon_p_LC1
table_dunn_p_shannon_LC1 <- dunn.test(p_sample_data$Shannon, g= p_sample_data$LC1_2018, method = "bh",list=TRUE)
dunn_shannon_p_LC1_df <- as.data.frame(table_dunn_p_shannon_LC1)
dunn_shannon_p_LC1_df$Variable <- "LC1"

kruskal.test(Shannon ~ LU2, data= p_sample_data)
kk_H_p_groups_LU2 <- kruskal(p_sample_data$Shannon,p_sample_data$LU2,group=TRUE,p.adj="bonferroni")
kk_H_p_groups_LU2
dunn_Shannon_p_LU2 <- kwAllPairsDunnTest(Shannon ~ as.factor(LU2), data= p_sample_data, p.adjust.method = "BH")
dunn_Shannon_p_LU2
table_dunn_p_Shannon_LU2 <- dunn.test(p_sample_data$Shannon, g= p_sample_data$LU2, method = "bh",list=TRUE)
dunn_Shannon_p_LU2_df <- as.data.frame(table_dunn_p_Shannon_LU2)
dunn_Shannon_p_LU2_df$Variable <- "LU2"

kruskal.test(Shannon ~ Sample_season, data= p_sample_data)
kk_H_p_groups_ss <- kruskal(p_sample_data$Shannon,p_sample_data$Sample_season,group=TRUE,p.adj="bonferroni")
kk_H_p_groups_ss
dunn_Shannon_p_season <- kwAllPairsDunnTest(Shannon ~ as.factor(Sample_season), data= p_sample_data, p.adjust.method = "BH")
dunn_Shannon_p_season
table_dunn_p_shannon_season <- dunn.test(p_sample_data$Shannon, g= p_sample_data$Sample_season, method = "bh",list=TRUE)
dunn_shannon_p_season_df <- as.data.frame(table_dunn_p_shannon_season)
dunn_shannon_p_season_df$Variable <- "Season"

kruskal.test(Shannon ~ pH_grouped, data= p_sample_data)
kk_H_p_groups_pH <- kruskal(p_sample_data$Shannon,p_sample_data$pH_grouped,group=TRUE,p.adj="bonferroni")
kk_H_p_groups_pH 
dunn_Shannon_p_pH <- kwAllPairsDunnTest(Shannon ~ as.factor(pH_grouped), data= p_sample_data, p.adjust.method = "BH")
dunn_Shannon_p_pH
table_dunn_p_shannon_pH <- dunn.test(p_sample_data$Shannon, g= p_sample_data$pH_grouped, method = "bh",list=TRUE)
dunn_shannon_p_pH_df <- as.data.frame(table_dunn_p_shannon_pH)
dunn_shannon_p_pH_df$Variable <- "pH"

kruskal.test(Shannon ~ depth_grouped, data= p_sample_data)
kk_H_p_groups_depth <- kruskal(p_sample_data$Shannon,p_sample_data$depth_grouped,group=TRUE,p.adj="bonferroni")
kk_H_p_groups_depth
dunn_Shannon_p_depth <- kwAllPairsDunnTest(Shannon ~ as.factor(depth_grouped), data= p_sample_data, p.adjust.method = "BH")
dunn_Shannon_p_depth
table_dunn_p_shannon_depth <- dunn.test(p_sample_data$Shannon, g= p_sample_data$depth_grouped, method = "bh",list=TRUE)
dunn_shannon_p_depth_df <- as.data.frame(table_dunn_p_shannon_depth)
dunn_shannon_p_depth_df$Variable <- "Depth"

kruskal.test(Shannon ~ Erosion_grouped, data= p_sample_data)
kk_H_p_groups_erosion <- kruskal(p_sample_data$Shannon,p_sample_data$Erosion_grouped,group=TRUE,p.adj="bonferroni")
kk_H_p_groups_erosion
dunn_Shannon_p_erosion <- kwAllPairsDunnTest(Shannon ~ as.factor(Erosion_grouped), data= p_sample_data, p.adjust.method = "BH")
dunn_Shannon_p_erosion
table_dunn_p_shannon_erosion <- dunn.test(p_sample_data$Shannon, g= p_sample_data$Erosion_grouped, method = "bh",list=TRUE)
dunn_shannon_p_erosion_df <- as.data.frame(table_dunn_p_shannon_erosion)
dunn_shannon_p_erosion_df$Variable <- "Erosion"

kruskal.test(Shannon ~  cluster_temp, data= p_sample_data)
kk_p_groups_cluster_temp <- kruskal(p_sample_data$Shannon,p_sample_data$cluster_temp,group=TRUE,p.adj="bonferroni")
kk_p_groups_cluster_temp
dunn_Shannon_p_cluster_temp <- kwAllPairsDunnTest(Shannon ~ as.factor(cluster_temp), data= p_sample_data, f.adjust.method = "BH")
dunn_Shannon_p_cluster_temp
table_dunn_p_Shannon_cluster_temp <- dunn.test(p_sample_data$Shannon, g= p_sample_data$cluster_temp, method = "bh",list=TRUE)
dunn_p_Shannon_cluster_temp_df <- as.data.frame(table_dunn_p_Shannon_cluster_temp)
dunn_p_Shannon_cluster_temp_df$Variable <- "Temp"

kruskal.test(Shannon ~  cluster_temp_range, data= p_sample_data)
kk_p_groups_cluster_temp_range <- kruskal(p_sample_data$Shannon,p_sample_data$cluster_temp_range,group=TRUE,p.adj="bonferroni")
kk_p_groups_cluster_temp_range
dunn_Shannon_p_cluster_temp_range <- kwAllPairsDunnTest(Shannon ~ as.factor(cluster_temp_range), data= p_sample_data, f.adjust.method = "BH")
dunn_Shannon_p_cluster_temp_range
table_dunn_p_Shannon_cluster_temp_range <- dunn.test(p_sample_data$Shannon, g= p_sample_data$cluster_temp_range, method = "bh",list=TRUE)
dunn_p_Shannon_cluster_temp_range_df <- as.data.frame(table_dunn_p_Shannon_cluster_temp_range)
dunn_p_Shannon_cluster_temp_range_df$Variable <- "Temp_range"

kruskal.test(Shannon ~  cluster_prec, data= p_sample_data) 
#not sign

#Fungi
#Richness
kruskal.test(Richness ~ LC1_2018, data= f_sample_data)
kk_f_groups_LC1 <- kruskal(f_sample_data$Shannon,f_sample_data$LC1_2018,group=TRUE,p.adj="bonferroni")
kk_f_groups_LC1
dunn_Richness_f_LC1 <- kwAllPairsDunnTest(Richness ~ as.factor(LC1_2018), data= f_sample_data, f.adjust.method = "BH")
dunn_Richness_f_LC1
table_dunn_f_richness_LC1 <- dunn.test(f_sample_data$Richness, g= f_sample_data$LC1_2018, method = "bh",list=TRUE)
dunn_f_richness_LC1_df <- as.data.frame(table_dunn_f_richness_LC1)
dunn_f_richness_LC1_df$Variable <- "LC1"

kruskal.test(Richness ~ LU2, data= f_sample_data)
kk_f_groups_LU2 <- kruskal(f_sample_data$Richness,f_sample_data$LU2,group=TRUE,p.adj="bonferroni")
dunn_Richness_f_LU2 <- kwAllPairsDunnTest(Richness ~ as.factor(LU2), data= f_sample_data, p.adjust.method = "BH")
dunn_Richness_f_LU2
table_dunn_f_richness_LU2 <- dunn.test(f_sample_data$Richness, g= f_sample_data$LU2, method = "bh",list=TRUE)
dunn_richness_f_LU2_df <- as.data.frame(table_dunn_f_richness_LU2)
dunn_richness_f_LU2_df$Variable <- "LU2"

kruskal.test(Richness ~ Sample_season, data= f_sample_data)
#not significant

kruskal.test(Richness ~ pH_grouped, data= f_sample_data)
kk_f_groups_pH <- kruskal(f_sample_data$Richness,f_sample_data$pH_grouped,group=TRUE,p.adj="bonferroni")
kk_f_groups_pH
dunn_Richness_f_pH <- kwAllPairsDunnTest(Richness ~ as.factor(pH_grouped), data= f_sample_data, f.adjust.method = "BH")
dunn_Richness_f_pH
table_dunn_f_richness_pH <- dunn.test(f_sample_data$Richness, g= f_sample_data$pH_grouped, method = "bh",list=TRUE)
dunn_f_richness_pH_df <- as.data.frame(table_dunn_f_richness_pH)
dunn_f_richness_pH_df$Variable <- "pH"

kruskal.test(Richness ~ depth_grouped, data= f_sample_data)
kk_f_groups_depth <- kruskal(f_sample_data$Richness,f_sample_data$depth_grouped,group=TRUE,p.adj="bonferroni")
kk_f_groups_depth
dunn_Richness_f_depth <- kwAllPairsDunnTest(Richness ~ as.factor(depth_grouped), data= f_sample_data, f.adjust.method = "BH")
dunn_Richness_f_depth
table_dunn_f_richness_depth <- dunn.test(f_sample_data$Richness, g= f_sample_data$depth_grouped, method = "bh",list=TRUE)
dunn_f_richness_depth_df <- as.data.frame(table_dunn_f_richness_depth)
dunn_f_richness_depth_df$Variable <- "Depth"

kruskal.test(Richness ~ Erosion_grouped, data= f_sample_data)
kk_f_groups_Erosion <- kruskal(f_sample_data$Richness,f_sample_data$Erosion_grouped,group=TRUE,p.adj="bonferroni")
kk_f_groups_Erosion
dunn_Richness_f_erosion <- kwAllPairsDunnTest(Richness ~ as.factor(Erosion_grouped), data= f_sample_data, f.adjust.method = "BH")
dunn_Richness_f_erosion
table_dunn_f_richness_erosion <- dunn.test(f_sample_data$Richness, g= f_sample_data$Erosion_grouped, method = "bh",list=TRUE)
dunn_f_richness_erosion_df <- as.data.frame(table_dunn_f_richness_erosion)
dunn_f_richness_erosion_df$Variable <- "Erosion"

kruskal.test(Richness ~  cluster_temp_range, data= f_sample_data)
kk_f_groups_cluster_temp_range <- kruskal(f_sample_data$Richness,f_sample_data$cluster_temp_range,group=TRUE,p.adj="bonferroni")
kk_f_groups_cluster_temp_range
dunn_Richness_f_cluster_temp_range <- kwAllPairsDunnTest(Richness ~ as.factor(cluster_temp_range), data= f_sample_data, f.adjust.method = "BH")
dunn_Richness_f_cluster_temp_range
table_dunn_f_richness_cluster_temp_range <- dunn.test(f_sample_data$Richness, g= f_sample_data$cluster_temp_range, method = "bh",list=TRUE)
dunn_f_richness_cluster_temp_range_df <- as.data.frame(table_dunn_f_richness_cluster_temp_range)
dunn_f_richness_cluster_temp_range_df$Variable <- "Temp_range"

kruskal.test(Richness ~  cluster_temp, data= f_sample_data)
kk_f_groups_cluster_temp <- kruskal(f_sample_data$Richness,f_sample_data$cluster_temp,group=TRUE,p.adj="bonferroni")
kk_f_groups_cluster_temp
dunn_Richness_f_cluster_temp <- kwAllPairsDunnTest(Richness ~ as.factor(cluster_temp), data= f_sample_data, f.adjust.method = "BH")
dunn_Richness_f_cluster_temp
table_dunn_f_richness_cluster_temp <- dunn.test(f_sample_data$Richness, g= f_sample_data$cluster_temp, method = "bh",list=TRUE)
dunn_f_richness_cluster_temp_df <- as.data.frame(table_dunn_f_richness_cluster_temp)
dunn_f_richness_cluster_temp_df$Variable <- "Temp"

kruskal.test(Richness ~  cluster_prec, data= f_sample_data)
kk_f_groups_cluster_prec <- kruskal(f_sample_data$Richness,f_sample_data$cluster_prec,group=TRUE,p.adj="bonferroni")
kk_f_groups_cluster_prec
dunn_Richness_f_cluster_prec <- kwAllPairsDunnTest(Richness ~ as.factor(cluster_prec), data= f_sample_data, f.adjust.method = "BH")
dunn_Richness_f_cluster_prec
table_dunn_f_richness_cluster_prec <- dunn.test(f_sample_data$Richness, g= f_sample_data$cluster_prec, method = "bh",list=TRUE)
dunn_f_richness_cluster_prec_df <- as.data.frame(table_dunn_f_richness_cluster_prec)
dunn_f_richness_cluster_prec_df$Variable <- "Prec"


#Shannon
kruskal.test(Shannon ~ LC1_2018, data= f_sample_data)
kk_H_f_groups_LC1 <- kruskal(f_sample_data$Shannon,f_sample_data$LC1_2018,group=TRUE,p.adj="bonferroni")
kk_H_f_groups_LC1 
dunn_Shannon_f_LC1 <- kwAllPairsDunnTest(Shannon ~ as.factor(LC1_2018), data= f_sample_data, f.adjust.method = "BH")
dunn_Shannon_f_LC1
table_dunn_f_shannon_LC1 <- dunn.test(f_sample_data$Shannon, g= f_sample_data$LC1_2018, method = "bh",list=TRUE)
dunn_f_shannon_LC1_df <- as.data.frame(table_dunn_f_shannon_LC1)
dunn_f_shannon_LC1_df$Variable <- "LC1"

kruskal.test(Shannon ~ LU2, data= f_sample_data)
kk_H_f_groups_LU2 <- kruskal(f_sample_data$Shannon,f_sample_data$LU2,group=TRUE,p.adj="bonferroni")
kk_H_f_groups_LU2
dunn_Shannon_f_LU2 <- kwAllPairsDunnTest(Shannon ~ as.factor(LU2), data= f_sample_data, p.adjust.method = "BH")
dunn_Shannon_f_LU2
table_dunn_f_Shannon_LU2 <- dunn.test(f_sample_data$Shannon, g= f_sample_data$LU2, method = "bh",list=TRUE)
dunn_Shannon_f_LU2_df <- as.data.frame(table_dunn_f_Shannon_LU2)
dunn_Shannon_f_LU2_df$Variable <- "LU2"

kruskal.test(Shannon ~ Sample_season, data= f_sample_data)
kk_H_f_groups_ss <- kruskal(f_sample_data$Shannon,f_sample_data$Sample_season,group=TRUE,p.adj="bonferroni")
kk_H_f_groups_ss
dunn_Shannon_f_season <- kwAllPairsDunnTest(Shannon ~ as.factor(Sample_season), data= f_sample_data, f.adjust.method = "BH")
dunn_Shannon_f_season
table_dunn_f_shannon_season <- dunn.test(f_sample_data$Shannon, g= f_sample_data$Sample_season, method = "bh",list=TRUE)
dunn_f_shannon_season_df <- as.data.frame(table_dunn_f_shannon_season)
dunn_f_shannon_season_df$Variable <- "Season"

kruskal.test(Shannon ~ pH_grouped, data= f_sample_data)
kk_H_f_groups_pH <- kruskal(f_sample_data$Shannon,f_sample_data$pH_grouped,group=TRUE,p.adj="bonferroni")
kk_H_f_groups_pH 
dunn_Shannon_f_pH <- kwAllPairsDunnTest(Shannon ~ as.factor(pH_grouped), data= f_sample_data, f.adjust.method = "BH")
dunn_Shannon_f_pH
table_dunn_f_shannon_pH <- dunn.test(f_sample_data$Shannon, g= f_sample_data$pH_grouped, method = "bh",list=TRUE)
dunn_f_shannon_pH_df <- as.data.frame(table_dunn_f_shannon_pH)
dunn_f_shannon_pH_df$Variable <- "pH"

kruskal.test(Shannon ~ depth_grouped, data= f_sample_data)
kk_H_f_groups_depth <- kruskal(f_sample_data$Shannon,f_sample_data$depth_grouped,group=TRUE,p.adj="bonferroni")
kk_H_f_groups_depth
dunn_Shannon_f_depth <- kwAllPairsDunnTest(Shannon ~ as.factor(depth_grouped), data= f_sample_data, f.adjust.method = "BH")
dunn_Shannon_f_depth
table_dunn_f_shannon_depth <- dunn.test(f_sample_data$Shannon, g= f_sample_data$depth_grouped, method = "bh",list=TRUE)
dunn_f_shannon_depth_df <- as.data.frame(table_dunn_f_shannon_depth)
dunn_f_shannon_depth_df$Variable <- "Depth"

kruskal.test(Shannon ~ Erosion_grouped, data= f_sample_data)
kk_H_f_groups_Erosion <- kruskal(f_sample_data$Shannon,f_sample_data$Erosion_grouped,group=TRUE,p.adj="bonferroni")
kk_H_f_groups_Erosion
dunn_Shannon_f_erosion <- kwAllPairsDunnTest(Shannon ~ as.factor(Erosion_grouped), data= f_sample_data, f.adjust.method = "BH")
dunn_Shannon_f_erosion
table_dunn_f_shannon_erosion <- dunn.test(f_sample_data$Shannon, g= f_sample_data$Erosion_grouped, method = "bh",list=TRUE)
dunn_f_shannon_erosion_df <- as.data.frame(table_dunn_f_shannon_erosion)
dunn_f_shannon_erosion_df$Variable <- "Erosion"

kruskal.test(Shannon ~  cluster_temp, data= f_sample_data)
kk_f_groups_cluster_temp <- kruskal(f_sample_data$Shannon,f_sample_data$cluster_temp,group=TRUE,p.adj="bonferroni")
kk_f_groups_cluster_temp
dunn_Shannon_f_cluster_temp <- kwAllPairsDunnTest(Shannon ~ as.factor(cluster_temp), data= f_sample_data, f.adjust.method = "BH")
dunn_Shannon_f_cluster_temp
table_dunn_f_Shannon_cluster_temp <- dunn.test(f_sample_data$Shannon, g= f_sample_data$cluster_temp, method = "bh",list=TRUE)
dunn_f_Shannon_cluster_temp_df <- as.data.frame(table_dunn_f_Shannon_cluster_temp)
dunn_f_Shannon_cluster_temp_df$Variable <- "Temp"

kruskal.test(Shannon ~  cluster_temp_range, data= f_sample_data)
kk_f_groups_cluster_temp_range <- kruskal(f_sample_data$Shannon,f_sample_data$cluster_temp_range,group=TRUE,p.adj="bonferroni")
kk_f_groups_cluster_temp_range
dunn_Shannon_f_cluster_temp_range <- kwAllPairsDunnTest(Shannon ~ as.factor(cluster_temp_range), data= f_sample_data, f.adjust.method = "BH")
dunn_Shannon_f_cluster_temp_range
table_dunn_f_Shannon_cluster_temp_range <- dunn.test(f_sample_data$Shannon, g= f_sample_data$cluster_temp_range, method = "bh",list=TRUE)
dunn_f_Shannon_cluster_temp_range_df <- as.data.frame(table_dunn_f_Shannon_cluster_temp_range)
dunn_f_Shannon_cluster_temp_range_df$Variable <- "Temp_range"

kruskal.test(Shannon ~  cluster_prec, data= f_sample_data)
kk_f_groups_cluster_prec <- kruskal(f_sample_data$Shannon,f_sample_data$cluster_prec,group=TRUE,p.adj="bonferroni")
kk_f_groups_cluster_prec
dunn_Shannon_f_cluster_prec <- kwAllPairsDunnTest(Shannon ~ as.factor(cluster_prec), data= f_sample_data, f.adjust.method = "BH")
dunn_Shannon_f_cluster_prec
table_dunn_f_Shannon_cluster_prec <- dunn.test(f_sample_data$Shannon, g= f_sample_data$cluster_prec, method = "bh",list=TRUE)
dunn_f_Shannon_cluster_prec_df <- as.data.frame(table_dunn_f_Shannon_cluster_prec)
dunn_f_Shannon_cluster_prec_df$Variable <- "Prec"



#rotifers
#Richness
kruskal.test(Richness ~ LC1_2018, data= a_roti_sample_data)
kk_a_roti_groups_LC1 <- kruskal(a_roti_sample_data$Richness,a_roti_sample_data$LC1_2018,group=TRUE,p.adj="bonferroni")
kk_a_roti_groups_LC1
dunn_richness_roti_LC1 <- kwAllPairsDunnTest(Richness ~ as.factor(LC1_2018), data= a_roti_sample_data, f.adjust.method = "BH")
dunn_richness_roti_LC1
table_dunn_roti_richness_LC1 <- dunn.test(a_roti_sample_data$Richness, g= a_roti_sample_data$LC1_2018, method = "bh",list=TRUE)
dunn_roti_richness_LC1_df <- as.data.frame(table_dunn_roti_richness_LC1)
dunn_roti_richness_LC1_df$Variable <- "LC1"

kruskal.test(Richness ~ LU2, data= a_roti_sample_data)
kk_a_roti_groups_LU2 <- kruskal(a_roti_sample_data$Richness,a_roti_sample_data$LU2,group=TRUE,p.adj="bonferroni")
kk_a_roti_groups_LU2
dunn_Richness_a_roti_LU2 <- kwAllPairsDunnTest(Richness ~ as.factor(LU2), data= a_roti_sample_data, p.adjust.method = "BH")
dunn_Richness_a_roti_LU2
table_dunn_a_roti_richness_LU2 <- dunn.test(a_roti_sample_data$Richness, g= a_roti_sample_data$LU2, method = "bh",list=TRUE)
dunn_richness_a_roti_LU2_df <- as.data.frame(table_dunn_a_roti_richness_LU2)
dunn_richness_a_roti_LU2_df$Variable <- "LU2"

kruskal.test(Richness ~ Sample_season, data= a_roti_sample_data) #not significant

kruskal.test(Richness ~ pH_grouped, data= a_roti_sample_data)
#not significant 

kruskal.test(Richness ~ depth_grouped, data= a_roti_sample_data)
#not significant 

kruskal.test(Richness ~ Erosion_grouped, data= a_roti_sample_data)
#not significant 

kruskal.test(Richness ~  cluster_temp_range, data= a_roti_sample_data)
#not significant 

kruskal.test(Richness ~  cluster_temp, data= a_roti_sample_data)
kk_roti_groups_cluster_temp <- kruskal(a_roti_sample_data$Richness,a_roti_sample_data$cluster_temp,group=TRUE,p.adj="bonferroni")
kk_roti_groups_cluster_temp
dunn_Richness_roti_cluster_temp <- kwAllPairsDunnTest(Richness ~ as.factor(cluster_temp), data= a_roti_sample_data, f.adjust.method = "BH")
dunn_Richness_roti_cluster_temp
table_dunn_roti_richness_cluster_temp <- dunn.test(a_roti_sample_data$Richness, g= a_roti_sample_data$cluster_temp, method = "bh",list=TRUE)
dunn_roti_richness_cluster_temp_df <- as.data.frame(table_dunn_roti_richness_cluster_temp)
dunn_roti_richness_cluster_temp_df$Variable <- "Temp"

kruskal.test(Richness ~  cluster_temp_range, data= a_roti_sample_data)
#not significant 



#Shannon
kruskal.test(Shannon ~ LC1_2018, data= a_roti_sample_data)
#not significant 

kruskal.test(Shannon ~ LU2, data= a_roti_sample_data)
#not significant 

kruskal.test(Shannon ~ Sample_season, data= a_roti_sample_data)
#not significant 

kruskal.test(Shannon ~ pH_grouped, data= a_roti_sample_data)
#not significant 

kruskal.test(Shannon ~ depth_grouped, data= a_roti_sample_data) 
#not significant

kruskal.test(Shannon ~ Erosion_grouped, data= a_roti_sample_data)
kk_H_a_roti_groups_Erosion <- kruskal(a_roti_sample_data$Shannon,a_roti_sample_data$Erosion_grouped,group=TRUE,p.adj="bonferroni")
kk_H_a_roti_groups_Erosion
dunn_Shannon_roti_erosion <- kwAllPairsDunnTest(Shannon ~ as.factor(Erosion_grouped), data= a_roti_sample_data, f.adjust.method = "BH")
dunn_Shannon_roti_erosion
table_dunn_roti_shannon_erosion <- dunn.test(a_roti_sample_data$Shannon, g= a_roti_sample_data$Erosion_grouped, method = "bh",list=TRUE)
dunn_roti_shannon_erosion_df <- as.data.frame(table_dunn_roti_shannon_erosion)
dunn_roti_shannon_erosion_df$Variable <- "Erosion"

kruskal.test(Shannon ~  cluster_temp, data= a_roti_sample_data) 
#not sign

kruskal.test(Shannon ~  cluster_prec, data= a_roti_sample_data) 
#not sign

kruskal.test(Shannon ~  cluster_temp_range, data= a_roti_sample_data)
#not significant 


#nematodes
#Richness
kruskal.test(Richness ~ LC1_2018, data= a_nema_sample_data)
kk_H_a_nema_groups_LC1 <- kruskal(a_nema_sample_data$Richness,a_nema_sample_data$LC1_2018,group=TRUE,p.adj="bonferroni")
kk_H_a_nema_groups_LC1
dunn_Richness_nema_LC1 <- kwAllPairsDunnTest(Richness ~ as.factor(LC1_2018), data= a_nema_sample_data, f.adjust.method = "BH")
dunn_Richness_nema_LC1
table_dunn_nema_Richness_LC1 <- dunn.test(a_nema_sample_data$Richness, g= a_nema_sample_data$LC1_2018, method = "bh",list=TRUE)
dunn_nema_Richness_LC1_df <- as.data.frame(table_dunn_nema_Richness_LC1)
dunn_nema_Richness_LC1_df$Variable <- "LC1"

kruskal.test(Richness ~ LU2, data= a_nema_sample_data)
kk_H_a_nema_groups_LU2 <- kruskal(a_nema_sample_data$Richness,a_nema_sample_data$LU2,group=TRUE,p.adj="bonferroni")
kk_H_a_nema_groups_LU2
dunn_Richness_nema_LU2 <- kwAllPairsDunnTest(Richness ~ as.factor(LU2), data= a_nema_sample_data, f.adjust.method = "BH")
dunn_Richness_nema_LU2
table_dunn_nema_Richness_LU2 <- dunn.test(a_nema_sample_data$Richness, g= a_nema_sample_data$LU2, method = "bh",list=TRUE)
dunn_nema_Richness_LU2_df <- as.data.frame(table_dunn_nema_Richness_LU2)
dunn_nema_Richness_LU2_df$Variable <- "LU2"

kruskal.test(Richness ~ Sample_season, data= a_nema_sample_data) 
#not significant

kruskal.test(Richness ~ pH_grouped, data= a_nema_sample_data) 
kk_H_a_nema_groups_pH <- kruskal(a_nema_sample_data$Shannon,a_nema_sample_data$pH_grouped,group=TRUE,p.adj="bonferroni")
kk_H_a_nema_groups_pH
dunn_Richness_nema_pH <- kwAllPairsDunnTest(Richness ~ as.factor(pH_grouped), data= a_nema_sample_data, f.adjust.method = "BH")
dunn_Richness_nema_pH
table_dunn_nema_Richness_pH <- dunn.test(a_nema_sample_data$Richness, g= a_nema_sample_data$pH_grouped, method = "bh",list=TRUE)
dunn_nema_Richness_pH_df <- as.data.frame(table_dunn_nema_Richness_pH)
dunn_nema_Richness_pH_df$Variable <- "pH"

kruskal.test(Richness ~ depth_grouped, data= a_nema_sample_data) 
#not significant

kruskal.test(Richness ~ Erosion_grouped, data= a_nema_sample_data) 
#not significant

kruskal.test(Richness ~  cluster_prec, data= a_nema_sample_data) 
#not significant

kruskal.test(Richness ~  cluster_temp_range, data= a_nema_sample_data) 
kk_H_a_nema_groups_cluster_temp_range <- kruskal(a_nema_sample_data$Richness,a_nema_sample_data$cluster_temp_range,group=TRUE,p.adj="bonferroni")
kk_H_a_nema_groups_cluster_temp_range
dunn_Richness_nema_cluster_temp_range <- kwAllPairsDunnTest(Richness ~ as.factor(cluster_temp_range), data= a_nema_sample_data, f.adjust.method = "BH")
dunn_Richness_nema_cluster_temp_range
table_dunn_nema_Richness_cluster_temp_range <- dunn.test(a_nema_sample_data$Richness, g= a_nema_sample_data$cluster_temp_range, method = "bh",list=TRUE)
dunn_nema_Richness_cluster_temp_range_df <- as.data.frame(table_dunn_nema_Richness_cluster_temp_range)
dunn_nema_Richness_cluster_temp_range_df$Variable <- "Temp_range"

kruskal.test(Richness ~  cluster_temp, data= a_nema_sample_data) 
kk_H_a_nema_groups_cluster_temp <- kruskal(a_nema_sample_data$Richness,a_nema_sample_data$cluster_temp,group=TRUE,p.adj="bonferroni")
kk_H_a_nema_groups_cluster_temp
dunn_Richness_nema_cluster_temp <- kwAllPairsDunnTest(Richness ~ as.factor(cluster_temp), data= a_nema_sample_data, f.adjust.method = "BH")
dunn_Richness_nema_cluster_temp
table_dunn_nema_Richness_cluster_temp <- dunn.test(a_nema_sample_data$Richness, g= a_nema_sample_data$cluster_temp, method = "bh",list=TRUE)
dunn_nema_Richness_cluster_temp_df <- as.data.frame(table_dunn_nema_Richness_cluster_temp)
dunn_nema_Richness_cluster_temp_df$Variable <- "Temp"

#Shannon
kruskal.test(Shannon ~ LC1_2018, data= a_nema_sample_data)
kk_H_a_nema_groups_LC1 <- kruskal(a_nema_sample_data$Shannon,a_nema_sample_data$LC1_2018,group=TRUE,p.adj="bonferroni")
kk_H_a_nema_groups_LC1
dunn_Shannon_nema_LC1 <- kwAllPairsDunnTest(Shannon ~ as.factor(LC1_2018), data= a_nema_sample_data, f.adjust.method = "BH")
dunn_Shannon_nema_LC1
table_dunn_nema_shannon_LC1 <- dunn.test(a_nema_sample_data$Shannon, g= a_nema_sample_data$LC1_2018, method = "bh",list=TRUE)
dunn_nema_shannon_LC1_df <- as.data.frame(table_dunn_nema_shannon_LC1)
dunn_nema_shannon_LC1_df$Variable <- "LC1"

kruskal.test(Shannon ~ LU2, data= a_nema_sample_data)
kk_H_a_nema_groups_LU2 <- kruskal(a_nema_sample_data$Shannon,a_nema_sample_data$LU2,group=TRUE,p.adj="bonferroni")
kk_H_a_nema_groups_LU2 
dunn_Shannon_a_nema_LU2 <- kwAllPairsDunnTest(Shannon ~ as.factor(LU2), data= a_nema_sample_data, p.adjust.method = "BH")
dunn_Shannon_a_nema_LU2
table_dunn_a_nema_Shannon_LU2 <- dunn.test(a_nema_sample_data$Shannon, g= a_nema_sample_data$LU2, method = "bh",list=TRUE)
dunn_Shannon_a_nema_LU2_df <- as.data.frame(table_dunn_a_nema_Shannon_LU2)
dunn_Shannon_a_nema_LU2_df$Variable <- "LU2"

kruskal.test(Shannon ~ Sample_season, data= a_nema_sample_data) 
#not significant

kruskal.test(Shannon ~ Erosion_grouped, data= a_nema_sample_data) 
#not significant

kruskal.test(Shannon ~ pH_grouped, data= a_nema_sample_data) 
kk_H_a_nema_groups_pH <- kruskal(a_nema_sample_data$Shannon,a_nema_sample_data$pH_grouped,group=TRUE,p.adj="bonferroni")
kk_H_a_nema_groups_pH
dunn_Shannon_nema_pH <- kwAllPairsDunnTest(Shannon ~ as.factor(pH_grouped), data= a_nema_sample_data, f.adjust.method = "BH")
dunn_Shannon_nema_pH
table_dunn_nema_shannon_pH <- dunn.test(a_nema_sample_data$Shannon, g= a_nema_sample_data$pH_grouped, method = "bh",list=TRUE)
dunn_nema_Shannon_pH_df <- as.data.frame(table_dunn_nema_shannon_pH)
dunn_nema_Shannon_pH_df$Variable <- "pH"

kruskal.test(Shannon ~ depth_grouped, data= a_nema_sample_data) 
#not significant

kruskal.test(Shannon ~  cluster_prec, data= a_nema_sample_data) 
#not sign

kruskal.test(Shannon ~  cluster_temp_range, data= a_nema_sample_data) 
kk_nema_groups_cluster_temp_range <- kruskal(a_nema_sample_data$Shannon,a_nema_sample_data$cluster_temp_range,group=TRUE,p.adj="bonferroni")
kk_nema_groups_cluster_temp_range
dunn_Shannon_nema_cluster_temp_range <- kwAllPairsDunnTest(Shannon ~ as.factor(cluster_temp_range), data= a_nema_sample_data, f.adjust.method = "BH")
dunn_Shannon_nema_cluster_temp_range
table_dunn_nema_Shannon_cluster_temp_range <- dunn.test(a_nema_sample_data$Shannon, g= a_nema_sample_data$cluster_temp_range, method = "bh",list=TRUE)
dunn_nema_Shannon_cluster_temp_range_df <- as.data.frame(table_dunn_nema_Shannon_cluster_temp_range)
dunn_nema_Shannon_cluster_temp_range_df$Variable <- "Temp_range"

kruskal.test(Shannon ~  cluster_temp, data= a_nema_sample_data) 
#not sign


#tardigrades
#Richness
kruskal.test(Richness ~ LC1_2018, data= a_tardi_sample_data)
kk_a_tardi_groups_LC1 <- kruskal(a_tardi_sample_data$Richness,a_tardi_sample_data$LC1_2018,group=TRUE,p.adj="bonferroni")
kk_a_tardi_groups_LC1
dunn_richness_tardi_LC1 <- kwAllPairsDunnTest(Richness ~ as.factor(LC1_2018), data= a_tardi_sample_data, f.adjust.method = "BH")
dunn_richness_tardi_LC1
table_dunn_tardi_richness_LC1 <- dunn.test(a_tardi_sample_data$Richness, g= a_tardi_sample_data$LC1_2018, method = "bh",list=TRUE)
dunn_tardi_richness_LC1_df <- as.data.frame(table_dunn_tardi_richness_LC1)
dunn_tardi_richness_LC1_df$Variable <- "LC1"

kruskal.test(Richness ~ LU2, data= a_tardi_sample_data)
kk_a_tardi_groups_LU2 <- kruskal(a_tardi_sample_data$Richness,a_tardi_sample_data$LU2,group=TRUE,p.adj="bonferroni")
dunn_Richness_a_tardi_LU2 <- kwAllPairsDunnTest(Richness ~ as.factor(LU2), data= a_tardi_sample_data, p.adjust.method = "BH")
dunn_Richness_a_tardi_LU2
table_dunn_a_tardi_richness_LU2 <- dunn.test(a_tardi_sample_data$Richness, g= a_tardi_sample_data$LU2, method = "bh",list=TRUE)
dunn_richness_a_tardi_LU2_df <- as.data.frame(table_dunn_a_tardi_richness_LU2)
dunn_richness_a_tardi_LU2_df$Variable <- "LU2"

kruskal.test(Richness ~ Sample_season, data= a_tardi_sample_data) 
#not significant

kruskal.test(Richness ~ pH_grouped, data= a_tardi_sample_data)
kk_a_tardi_groups_pH <- kruskal(a_tardi_sample_data$Richness,a_tardi_sample_data$pH_grouped,group=TRUE,p.adj="bonferroni")
kk_a_tardi_groups_pH 
dunn_richness_tardi_pH <- kwAllPairsDunnTest(Richness ~ as.factor(pH_grouped), data= a_tardi_sample_data, f.adjust.method = "BH")
dunn_richness_tardi_pH
table_dunn_tardi_richness_pH <- dunn.test(a_tardi_sample_data$Richness, g= a_tardi_sample_data$pH_grouped, method = "bh",list=TRUE)
dunn_tardi_richness_pH_df <- as.data.frame(table_dunn_tardi_richness_pH)
dunn_tardi_richness_pH_df$Variable <- "pH"

kruskal.test(Richness ~ depth_grouped, data= a_tardi_sample_data)
kk_a_tardi_groups_depth <- kruskal(a_tardi_sample_data$Richness,a_tardi_sample_data$depth_grouped,group=TRUE,p.adj="bonferroni")
kk_a_tardi_groups_depth
dunn_richness_tardi_depth <- kwAllPairsDunnTest(Richness ~ as.factor(depth_grouped), data= a_tardi_sample_data, f.adjust.method = "BH")
dunn_richness_tardi_depth
table_dunn_tardi_richness_depth <- dunn.test(a_tardi_sample_data$Richness, g= a_tardi_sample_data$depth_grouped, method = "bh",list=TRUE)
dunn_tardi_richness_depth_df <- as.data.frame(table_dunn_tardi_richness_depth)
dunn_tardi_richness_depth_df$Variable <- "Depth"

kruskal.test(Richness ~ Erosion_grouped, data= a_tardi_sample_data)
#not significant

kruskal.test(Richness ~  cluster_temp_range, data= a_tardi_sample_data)
kk_tardi_groups_cluster_temp_range <- kruskal(a_tardi_sample_data$Richness,a_tardi_sample_data$cluster_temp_range,group=TRUE,p.adj="bonferroni")
kk_tardi_groups_cluster_temp_range
dunn_Richness_tardi_cluster_temp_range <- kwAllPairsDunnTest(Richness ~ as.factor(cluster_temp_range), data= a_tardi_sample_data, f.adjust.method = "BH")
dunn_Richness_tardi_cluster_temp_range
table_dunn_tardi_richness_cluster_temp_range <- dunn.test(a_tardi_sample_data$Richness, g= a_tardi_sample_data$cluster_temp_range, method = "bh",list=TRUE)
dunn_tardi_richness_cluster_temp_range_df <- as.data.frame(table_dunn_tardi_richness_cluster_temp_range)
dunn_tardi_richness_cluster_temp_range_df$Variable <- "Temp_range"

kruskal.test(Richness ~  cluster_prec, data= a_tardi_sample_data)
#not sign

kruskal.test(Richness ~  cluster_temp, data= a_tardi_sample_data)
#not sign

#Shannon
kruskal.test(Shannon ~ LC1_2018, data= a_tardi_sample_data)
kk_H_a_tardi_groups_LC1 <- kruskal(a_tardi_sample_data$Shannon,a_tardi_sample_data$LC1_2018,group=TRUE,p.adj="bonferroni")
kk_H_a_tardi_groups_LC1
dunn_Shannon_tardi_LC1 <- kwAllPairsDunnTest(Shannon ~ as.factor(LC1_2018), data= a_tardi_sample_data, f.adjust.method = "BH")
dunn_Shannon_tardi_LC1
table_dunn_tardi_shannon_LC1 <- dunn.test(a_tardi_sample_data$Shannon, g= a_tardi_sample_data$LC1_2018, method = "bh",list=TRUE)
dunn_tardi_shannon_LC1_df <- as.data.frame(table_dunn_tardi_shannon_LC1)
dunn_tardi_shannon_LC1_df$Variable <- "LC1"

kruskal.test(Shannon ~ LU2, data= a_tardi_sample_data)
kk_H_a_tardi_groups_LU2 <- kruskal(a_tardi_sample_data$Shannon,a_tardi_sample_data$LU2,group=TRUE,p.adj="bonferroni")
kk_H_a_tardi_groups_LU2 
dunn_Shannon_a_tardi_LU2 <- kwAllPairsDunnTest(Shannon ~ as.factor(LU2), data= a_tardi_sample_data, p.adjust.method = "BH")
dunn_Shannon_a_tardi_LU2
table_dunn_a_tardi_Shannon_LU2 <- dunn.test(a_tardi_sample_data$Shannon, g= a_tardi_sample_data$LU2, method = "bh",list=TRUE)
dunn_Shannon_a_tardi_LU2_df <- as.data.frame(table_dunn_a_tardi_Shannon_LU2)
dunn_Shannon_a_tardi_LU2_df$Variable <- "LU2"

kruskal.test(Shannon ~ Sample_season, data= a_tardi_sample_data) 
#not significant

kruskal.test(Shannon ~ Erosion_grouped, data= a_tardi_sample_data)
#not significant

kruskal.test(Shannon ~ pH_grouped, data= a_tardi_sample_data)
kk_H_a_tardi_groups_pH <- kruskal(a_tardi_sample_data$Shannon,a_tardi_sample_data$pH_grouped,group=TRUE,p.adj="bonferroni")
kk_H_a_tardi_groups_pH
dunn_Shannon_tardi_pH <- kwAllPairsDunnTest(Shannon ~ as.factor(pH_grouped), data= a_tardi_sample_data, f.adjust.method = "BH")
dunn_Shannon_tardi_pH
table_dunn_tardi_shannon_pH <- dunn.test(a_tardi_sample_data$Shannon, g= a_tardi_sample_data$pH_grouped, method = "bh",list=TRUE)
dunn_tardi_shannon_pH_df <- as.data.frame(table_dunn_tardi_shannon_pH)
dunn_tardi_shannon_pH_df$Variable <- "pH"

kruskal.test(Shannon ~ depth_grouped, data= a_tardi_sample_data)
kk_H_a_tardi_groups_depth <- kruskal(a_tardi_sample_data$Shannon,a_tardi_sample_data$depth_grouped,group=TRUE,p.adj="bonferroni")
kk_H_a_tardi_groups_depth
dunn_Shannon_tardi_depth <- kwAllPairsDunnTest(Shannon ~ as.factor(depth_grouped), data= a_tardi_sample_data, f.adjust.method = "BH")
dunn_Shannon_tardi_depth
table_dunn_tardi_shannon_depth <- dunn.test(a_tardi_sample_data$Shannon, g= a_tardi_sample_data$depth_grouped, method = "bh",list=TRUE)
dunn_tardi_shannon_depth_df <- as.data.frame(table_dunn_tardi_shannon_depth)
dunn_tardi_shannon_depth_df$Variable <- "Depth"


kruskal.test(Shannon ~  cluster_prec, data= a_tardi_sample_data)
#not significant

kruskal.test(Shannon ~  cluster_temp_range, data= a_tardi_sample_data)
kk_tardi_groups_cluster_temp_range <- kruskal(a_tardi_sample_data$Shannon,a_tardi_sample_data$cluster_temp_range,group=TRUE,p.adj="bonferroni")
kk_tardi_groups_cluster_temp_range
dunn_Shannon_tardi_cluster_temp_range <- kwAllPairsDunnTest(Shannon ~ as.factor(cluster_temp_range), data= a_tardi_sample_data, f.adjust.method = "BH")
dunn_Shannon_tardi_cluster_temp_range
table_dunn_tardi_Shannon_cluster_temp_range <- dunn.test(a_tardi_sample_data$Shannon, g= a_tardi_sample_data$cluster_temp_range, method = "bh",list=TRUE)
dunn_tardi_Shannon_cluster_temp_range_df <- as.data.frame(table_dunn_tardi_Shannon_cluster_temp_range)
dunn_tardi_Shannon_cluster_temp_range_df$Variable <- "Temp_range"

kruskal.test(Shannon ~  cluster_temp, data= a_tardi_sample_data)
#not significant


#arthropods
#Richness
kruskal.test(Richness ~ LC1_2018, data= a_arthro_sample_data)
kk_a_arthro_groups_LC1 <- kruskal(a_arthro_sample_data$Richness,a_arthro_sample_data$LC1_2018,group=TRUE,p.adj="bonferroni")
kk_a_arthro_groups_LC1
dunn_richness_arthro_LC1 <- kwAllPairsDunnTest(Richness ~ as.factor(LC1_2018), data= a_arthro_sample_data, f.adjust.method = "BH")
dunn_richness_arthro_LC1
table_dunn_arthro_richness_LC1 <- dunn.test(a_arthro_sample_data$Richness, g= a_arthro_sample_data$LC1_2018, method = "bh",list=TRUE)
dunn_arthro_richness_LC1_df <- as.data.frame(table_dunn_arthro_richness_LC1)
dunn_arthro_richness_LC1_df$Variable <- "LC1"

kruskal.test(Richness ~ LU2, data= a_arthro_sample_data)
kk_a_arthro_groups_LU2 <- kruskal(a_arthro_sample_data$Richness,a_arthro_sample_data$LU2,group=TRUE,p.adj="bonferroni")
kk_a_arthro_groups_LU2
dunn_Richness_a_arthro_LU2 <- kwAllPairsDunnTest(Richness ~ as.factor(LU2), data= a_arthro_sample_data, p.adjust.method = "BH")
dunn_Richness_a_arthro_LU2
table_dunn_a_arthro_richness_LU2 <- dunn.test(a_arthro_sample_data$Richness, g= a_arthro_sample_data$LU2, method = "bh",list=TRUE)
dunn_richness_a_arthro_LU2_df <- as.data.frame(table_dunn_a_arthro_richness_LU2)
dunn_richness_a_arthro_LU2_df$Variable <- "LU2"

kruskal.test(Richness ~ Sample_season, data= a_arthro_sample_data) 
#not significant

kruskal.test(Richness ~ pH_grouped, data= a_arthro_sample_data)
kk_a_arthro_groups_pH <- kruskal(a_arthro_sample_data$Richness,a_arthro_sample_data$pH_grouped,group=TRUE,p.adj="bonferroni")
kk_a_arthro_groups_pH
dunn_richness_arthro_pH <- kwAllPairsDunnTest(Richness ~ as.factor(pH_grouped), data= a_arthro_sample_data, f.adjust.method = "BH")
dunn_richness_arthro_pH
table_dunn_arthro_richness_pH <- dunn.test(a_arthro_sample_data$Richness, g= a_arthro_sample_data$pH_grouped, method = "bh",list=TRUE)
dunn_arthro_richness_pH_df <- as.data.frame(table_dunn_arthro_richness_pH)
dunn_arthro_richness_pH_df$Variable <- "pH"

kruskal.test(Richness ~ depth_grouped, data= a_arthro_sample_data)
kk_a_arthro_groups_depth <- kruskal(a_arthro_sample_data$Richness,a_arthro_sample_data$depth_grouped,group=TRUE,p.adj="bonferroni")
kk_a_arthro_groups_depth
dunn_richness_arthro_depth <- kwAllPairsDunnTest(Richness ~ as.factor(depth_grouped), data= a_arthro_sample_data, f.adjust.method = "BH")
dunn_richness_arthro_depth
table_dunn_arthro_richness_depth <- dunn.test(a_arthro_sample_data$Richness, g= a_arthro_sample_data$depth_grouped, method = "bh",list=TRUE)
dunn_arthro_richness_depth_df <- as.data.frame(table_dunn_arthro_richness_depth)
dunn_arthro_richness_depth_df$Variable <- "Depth"

kruskal.test(Richness ~ Erosion_grouped, data= a_arthro_sample_data)
kk_a_arthro_groups_Erosion <- kruskal(a_arthro_sample_data$Richness,a_arthro_sample_data$Erosion_grouped,group=TRUE,p.adj="bonferroni")
kk_a_arthro_groups_Erosion
dunn_Richness_arthro_erosion <- kwAllPairsDunnTest(Richness ~ as.factor(Erosion_grouped), data= p_sample_data, p.adjust.method = "BH")
dunn_Richness_arthro_erosion
table_dunn_arthro_richness_erosion <- dunn.test(a_arthro_sample_data$Richness, g= a_arthro_sample_data$Erosion_grouped, method = "bh",list=TRUE)
dunn_arthro_richness_erosion_df <- as.data.frame(table_dunn_arthro_richness_erosion)
dunn_arthro_richness_erosion_df$Variable <- "Erosion"

kruskal.test(Richness ~  cluster_temp_range, data= a_arthro_sample_data)
kk_arthro_groups_cluster_temp_range <- kruskal(a_arthro_sample_data$Richness,a_arthro_sample_data$cluster_temp_range,group=TRUE,p.adj="bonferroni")
kk_arthro_groups_cluster_temp_range
dunn_Richness_arthro_cluster_temp_range <- kwAllPairsDunnTest(Richness ~ as.factor(cluster_temp_range), data= a_arthro_sample_data, f.adjust.method = "BH")
dunn_Richness_arthro_cluster_temp_range
table_dunn_arthro_richness_cluster_temp_range <- dunn.test(a_arthro_sample_data$Richness, g= a_arthro_sample_data$cluster_temp_range, method = "bh",list=TRUE)
dunn_arthro_richness_cluster_temp_range_df <- as.data.frame(table_dunn_arthro_richness_cluster_temp_range)
dunn_arthro_richness_cluster_temp_range_df$Variable <- "Temp_range"

kruskal.test(Richness ~  cluster_temp, data= a_arthro_sample_data)
kk_arthro_groups_cluster_temp <- kruskal(a_arthro_sample_data$Richness,a_arthro_sample_data$cluster_temp,group=TRUE,p.adj="bonferroni")
kk_arthro_groups_cluster_temp
dunn_Richness_arthro_cluster_temp <- kwAllPairsDunnTest(Richness ~ as.factor(cluster_temp), data= a_arthro_sample_data, f.adjust.method = "BH")
dunn_Richness_arthro_cluster_temp
table_dunn_arthro_richness_cluster_temp <- dunn.test(a_arthro_sample_data$Richness, g= a_arthro_sample_data$cluster_temp, method = "bh",list=TRUE)
dunn_arthro_richness_cluster_temp_df <- as.data.frame(table_dunn_arthro_richness_cluster_temp)
dunn_arthro_richness_cluster_temp_df$Variable <- "Temp"

kruskal.test(Richness ~  cluster_prec, data= a_arthro_sample_data)
#not significant


#Shannon
kruskal.test(Shannon ~ LC1_2018, data= a_arthro_sample_data)
kk_H_a_arthro_groups_LC1 <- kruskal(a_arthro_sample_data$Shannon,a_arthro_sample_data$LC1_2018,group=TRUE,p.adj="bonferroni")
kk_H_a_arthro_groups_LC1 
dunn_Shannon_arthro_LC1 <- kwAllPairsDunnTest(Shannon ~ as.factor(LC1_2018), data= a_arthro_sample_data, f.adjust.method = "BH")
dunn_Shannon_arthro_LC1
table_dunn_arthro_shannon_LC1 <- dunn.test(a_arthro_sample_data$Shannon, g= a_arthro_sample_data$LC1_2018, method = "bh",list=TRUE)
dunn_arthro_shannon_LC1_df <- as.data.frame(table_dunn_arthro_shannon_LC1)
dunn_arthro_shannon_LC1_df$Variable <- "LC1"

kruskal.test(Shannon ~ LU2, data= a_arthro_sample_data)
kk_H_a_arthro_groups_LU2 <- kruskal(a_arthro_sample_data$Shannon,a_arthro_sample_data$LU2,group=TRUE,p.adj="bonferroni")
dunn_Shannon_a_arthro_LU2 <- kwAllPairsDunnTest(Shannon ~ as.factor(LU2), data= a_arthro_sample_data, p.adjust.method = "BH")
dunn_Shannon_a_arthro_LU2
table_dunn_a_arthro_Shannon_LU2 <- dunn.test(a_arthro_sample_data$Shannon, g= a_arthro_sample_data$LU2, method = "bh",list=TRUE)
dunn_Shannon_a_arthro_LU2_df <- as.data.frame(table_dunn_a_arthro_Shannon_LU2)
dunn_Shannon_a_arthro_LU2_df$Variable <- "LU2"

kruskal.test(Shannon ~ Sample_season, data= a_arthro_sample_data)
kk_H_a_arthro_groups_ss <- kruskal(a_arthro_sample_data$Shannon,a_arthro_sample_data$Sample_season,group=TRUE,p.adj="bonferroni")
kk_H_a_arthro_groups_ss
dunn_Shannon_arthro_season <- kwAllPairsDunnTest(Shannon ~ as.factor(Sample_season), data= a_arthro_sample_data, f.adjust.method = "BH")
dunn_Shannon_arthro_season
table_dunn_arthro_shannon_season <- dunn.test(a_arthro_sample_data$Shannon, g= a_arthro_sample_data$Sample_season, method = "bh",list=TRUE)
dunn_arthro_shannon_season_df <- as.data.frame(table_dunn_arthro_shannon_season)
dunn_arthro_shannon_season_df$Variable <- "Season"

kruskal.test(Shannon ~ pH_grouped, data= a_arthro_sample_data)
kk_H_a_arthro_groups_pH <- kruskal(a_arthro_sample_data$Shannon,a_arthro_sample_data$pH_grouped,group=TRUE,p.adj="bonferroni")
kk_H_a_arthro_groups_pH
dunn_Shannon_arthro_pH <- kwAllPairsDunnTest(Shannon ~ as.factor(pH_grouped), data= a_arthro_sample_data, f.adjust.method = "BH")
dunn_Shannon_arthro_pH
table_dunn_arthro_shannon_pH <- dunn.test(a_arthro_sample_data$Shannon, g= a_arthro_sample_data$pH_grouped, method = "bh",list=TRUE)
dunn_arthro_shannon_pH_df <- as.data.frame(table_dunn_arthro_shannon_pH)
dunn_arthro_shannon_pH_df$Variable <- "pH"

kruskal.test(Shannon ~ depth_grouped, data= a_arthro_sample_data) #not significant

kruskal.test(Shannon ~ Erosion_grouped, data= a_arthro_sample_data)
#not significant

kruskal.test(Shannon ~  cluster_temp, data= a_arthro_sample_data) 
#not sign

kruskal.test(Shannon ~  cluster_prec, data= a_arthro_sample_data) 
#not sign

kruskal.test(Shannon ~  cluster_temp_range, data= a_arthro_sample_data)
#not significant




#annelids
#Richness
kruskal.test(Richness ~ LC1_2018, data= a_anneli_sample_data)
kk_a_anneli_groups_LC1 <- kruskal(a_anneli_sample_data$Richness,a_anneli_sample_data$LC1_2018,group=TRUE,p.adj="bonferroni")
kk_a_anneli_groups_LC1
dunn_richness_anneli_LC1 <- kwAllPairsDunnTest(Richness ~ as.factor(LC1_2018), data= a_anneli_sample_data, f.adjust.method = "BH")
dunn_richness_anneli_LC1
table_dunn_anneli_richness_LC1 <- dunn.test(a_anneli_sample_data$Richness, g= a_anneli_sample_data$LC1_2018, method = "bh",list=TRUE)
dunn_anneli_richness_LC1_df <- as.data.frame(table_dunn_anneli_richness_LC1)
dunn_anneli_richness_LC1_df$Variable <- "LC1"

kruskal.test(Richness ~ LU2, data= a_anneli_sample_data)
kk_a_anneli_groups_LU2 <- kruskal(a_anneli_sample_data$Richness,a_anneli_sample_data$LU2,group=TRUE,p.adj="bonferroni")
kk_a_anneli_groups_LU2
dunn_Richness_a_anneli_LU2 <- kwAllPairsDunnTest(Richness ~ as.factor(LU2), data= a_anneli_sample_data, p.adjust.method = "BH")
dunn_Richness_a_anneli_LU2
table_dunn_a_anneli_richness_LU2 <- dunn.test(a_anneli_sample_data$Richness, g= a_anneli_sample_data$LU2, method = "bh",list=TRUE)
dunn_richness_a_anneli_LU2_df <- as.data.frame(table_dunn_a_anneli_richness_LU2)
dunn_richness_a_anneli_LU2_df$Variable <- "LU2"

kruskal.test(Richness ~ Sample_season, data= a_anneli_sample_data) 
#not significant

kruskal.test(Richness ~ pH_grouped, data= a_anneli_sample_data)
kk_a_anneli_groups_pH <- kruskal(a_anneli_sample_data$Richness,a_anneli_sample_data$pH_grouped,group=TRUE,p.adj="bonferroni")
kk_a_anneli_groups_pH
dunn_richness_anneli_pH <- kwAllPairsDunnTest(Richness ~ as.factor(pH_grouped), data= a_anneli_sample_data, f.adjust.method = "BH")
dunn_richness_anneli_pH
table_dunn_anneli_richness_pH <- dunn.test(a_anneli_sample_data$Richness, g= a_anneli_sample_data$pH_grouped, method = "bh",list=TRUE)
dunn_anneli_richness_pH_df <- as.data.frame(table_dunn_anneli_richness_pH)
dunn_anneli_richness_pH_df$Variable <- "pH"

kruskal.test(Richness ~ depth_grouped, data= a_anneli_sample_data)
#not significant

kruskal.test(Richness ~ Erosion_grouped, data= a_anneli_sample_data)
kk_a_anneli_groups_Erosion <- kruskal(a_anneli_sample_data$Richness,a_anneli_sample_data$Erosion_grouped,group=TRUE,p.adj="bonferroni")
kk_a_anneli_groups_Erosion
dunn_Richness_anneli_erosion <- kwAllPairsDunnTest(Richness ~ as.factor(Erosion_grouped), data= p_sample_data, p.adjust.method = "BH")
dunn_Richness_anneli_erosion
table_dunn_anneli_richness_erosion <- dunn.test(a_anneli_sample_data$Richness, g= a_anneli_sample_data$Erosion_grouped, method = "bh",list=TRUE)
dunn_anneli_richness_erosion_df <- as.data.frame(table_dunn_anneli_richness_erosion)
dunn_anneli_richness_erosion_df$Variable <- "Erosion"

kruskal.test(Richness ~  cluster_temp_range, data= a_anneli_sample_data)
#not significant

kruskal.test(Richness ~  cluster_temp, data= a_anneli_sample_data)
#not significant

kruskal.test(Richness ~  cluster_prec, data= a_anneli_sample_data)
kk_a_anneli_groups_cluster_prec <- kruskal(a_anneli_sample_data$Richness,a_anneli_sample_data$cluster_prec,group=TRUE,p.adj="bonferroni")
kk_a_anneli_groups_cluster_prec
dunn_Richness_anneli_cluster_prec <- kwAllPairsDunnTest(Richness ~ as.factor(cluster_prec), data= p_sample_data, p.adjust.method = "BH")
dunn_Richness_anneli_cluster_prec
table_dunn_anneli_richness_cluster_prec <- dunn.test(a_anneli_sample_data$Richness, g= a_anneli_sample_data$cluster_prec, method = "bh",list=TRUE)
dunn_anneli_richness_cluster_prec_df <- as.data.frame(table_dunn_anneli_richness_cluster_prec)
dunn_anneli_richness_cluster_prec_df$Variable <- "Prec"




#Shannon
kruskal.test(Shannon ~ LC1_2018, data= a_anneli_sample_data)
kk_H_a_anneli_groups_LC1 <- kruskal(a_anneli_sample_data$Shannon,a_anneli_sample_data$LC1_2018,group=TRUE,p.adj="bonferroni")
kk_H_a_anneli_groups_LC1 
dunn_Shannon_anneli_LC1 <- kwAllPairsDunnTest(Shannon ~ as.factor(LC1_2018), data= a_anneli_sample_data, f.adjust.method = "BH")
dunn_Shannon_anneli_LC1
table_dunn_anneli_shannon_LC1 <- dunn.test(a_anneli_sample_data$Shannon, g= a_anneli_sample_data$LC1_2018, method = "bh",list=TRUE)
dunn_anneli_shannon_LC1_df <- as.data.frame(table_dunn_anneli_shannon_LC1)
dunn_anneli_shannon_LC1_df$Variable <- "LC1"

kruskal.test(Shannon ~ LU2, data= a_anneli_sample_data)
kk_H_a_anneli_groups_LU2 <- kruskal(a_anneli_sample_data$Shannon,a_anneli_sample_data$LU2,group=TRUE,p.adj="bonferroni")
dunn_Shannon_a_anneli_LU2 <- kwAllPairsDunnTest(Shannon ~ as.factor(LU2), data= a_anneli_sample_data, p.adjust.method = "BH")
dunn_Shannon_a_anneli_LU2
table_dunn_a_anneli_Shannon_LU2 <- dunn.test(a_anneli_sample_data$Shannon, g= a_anneli_sample_data$LU2, method = "bh",list=TRUE)
dunn_Shannon_a_anneli_LU2_df <- as.data.frame(table_dunn_a_anneli_Shannon_LU2)
dunn_Shannon_a_anneli_LU2_df$Variable <- "LU2"

kruskal.test(Shannon ~ Sample_season, data= a_anneli_sample_data)
#not significant

kruskal.test(Shannon ~ pH_grouped, data= a_anneli_sample_data)
kk_H_a_anneli_groups_pH <- kruskal(a_anneli_sample_data$Shannon,a_anneli_sample_data$pH_grouped,group=TRUE,p.adj="bonferroni")
kk_H_a_anneli_groups_pH
dunn_Shannon_anneli_pH <- kwAllPairsDunnTest(Shannon ~ as.factor(pH_grouped), data= a_anneli_sample_data, f.adjust.method = "BH")
dunn_Shannon_anneli_pH
table_dunn_anneli_shannon_pH <- dunn.test(a_anneli_sample_data$Shannon, g= a_anneli_sample_data$pH_grouped, method = "bh",list=TRUE)
dunn_anneli_shannon_pH_df <- as.data.frame(table_dunn_anneli_shannon_pH)
dunn_anneli_shannon_pH_df$Variable <- "pH"

kruskal.test(Shannon ~ depth_grouped, data= a_anneli_sample_data) 
#not significant

kruskal.test(Shannon ~ Erosion_grouped, data= a_anneli_sample_data)
#not significant

kruskal.test(Shannon ~  cluster_temp, data= a_anneli_sample_data) 
#not sign

kruskal.test(Shannon ~  cluster_prec, data= a_anneli_sample_data) 
#not sign

kruskal.test(Shannon ~  cluster_temp_range, data= a_anneli_sample_data)
#not significant



#Common tables for all
tardi_shannon_dunn <- rbind(dunn_tardi_shannon_depth_df,dunn_tardi_shannon_LC1_df,dunn_Shannon_a_tardi_LU2_df,dunn_tardi_shannon_pH_df,dunn_tardi_Shannon_cluster_temp_range_df)
tardi_shannon_dunn$group <- "Tardigrades" 
nema_shannon_dunn <- rbind(dunn_nema_shannon_LC1_df,dunn_Shannon_a_nema_LU2_df,dunn_nema_Shannon_pH_df,dunn_nema_Shannon_cluster_temp_range_df)
nema_shannon_dunn$group <- "Nematodes" 
roti_shannon_dunn <- rbind(dunn_roti_shannon_erosion_df)
roti_shannon_dunn$group <- "Rotifers" 
anneli_shannon_dunn <- rbind(dunn_anneli_shannon_LC1_df,dunn_Shannon_a_anneli_LU2_df,dunn_anneli_shannon_pH_df)
anneli_shannon_dunn$group <- "Annelids" 
arthro_shannon_dunn <- rbind(dunn_arthro_shannon_season_df,dunn_arthro_shannon_LC1_df,dunn_Shannon_a_arthro_LU2_df,dunn_arthro_shannon_pH_df)
arthro_shannon_dunn$group <- "Arthropods" 
f_shannon_dunn <- rbind(dunn_f_shannon_depth_df,dunn_f_shannon_season_df,dunn_f_shannon_LC1_df,dunn_Shannon_f_LU2_df ,dunn_f_shannon_erosion_df,dunn_f_shannon_pH_df,dunn_f_Shannon_cluster_temp_df,dunn_f_Shannon_cluster_temp_range_df,dunn_f_Shannon_cluster_prec_df)
f_shannon_dunn$group <- "Fungi" 
p_shannon_dunn <- rbind(dunn_shannon_p_depth_df,dunn_shannon_p_season_df,dunn_shannon_p_LC1_df,dunn_Shannon_p_LU2_df ,dunn_shannon_p_erosion_df,dunn_shannon_p_pH_df,dunn_p_Shannon_cluster_temp_df,dunn_p_Shannon_cluster_temp_range_df)
p_shannon_dunn$group <- "Protists" 

dunn_shannon_complete <- rbind(f_shannon_dunn,p_shannon_dunn,roti_shannon_dunn,tardi_shannon_dunn,nema_shannon_dunn,arthro_shannon_dunn,anneli_shannon_dunn)
dunn_shannon_complete$Diversity_index <- "Shannon index" 

tardi_richness_dunn <- rbind(dunn_tardi_richness_depth_df,dunn_tardi_richness_LC1_df,dunn_richness_a_tardi_LU2_df,dunn_tardi_richness_pH_df,dunn_tardi_richness_cluster_temp_range_df)
tardi_richness_dunn$group <- "Tardigrades" 
nema_richness_dunn <- rbind(dunn_nema_Richness_LC1_df,dunn_nema_Richness_LU2_df,dunn_nema_Richness_pH_df,dunn_nema_Richness_cluster_temp_df,dunn_nema_Richness_cluster_temp_range_df)
nema_richness_dunn$group <- "Nematodes"
roti_richness_dunn <- rbind(dunn_roti_richness_LC1_df,dunn_richness_a_roti_LU2_df,dunn_roti_richness_cluster_temp_df)
roti_richness_dunn$group <- "Rotifers" 
anneli_richness_dunn <- rbind(dunn_anneli_richness_LC1_df,dunn_richness_a_anneli_LU2_df ,dunn_anneli_richness_erosion_df,dunn_anneli_richness_pH_df,dunn_anneli_richness_cluster_prec_df)
anneli_richness_dunn$group <- "Annelids" 
arthro_richness_dunn <- rbind(dunn_arthro_richness_depth_df,dunn_arthro_richness_LC1_df,dunn_richness_a_arthro_LU2_df ,dunn_arthro_richness_erosion_df,dunn_arthro_richness_pH_df,dunn_arthro_richness_cluster_temp_range_df,dunn_arthro_richness_cluster_temp_df)
arthro_richness_dunn$group <- "Arthropods" 
f_richness_dunn <- rbind(dunn_f_richness_depth_df,dunn_f_richness_season_df,dunn_f_richness_LC1_df,dunn_richness_f_LU2_df,dunn_f_richness_erosion_df,dunn_f_richness_pH_df,dunn_f_richness_cluster_temp_df,dunn_f_richness_cluster_temp_range_df,dunn_f_richness_cluster_prec_df)
f_richness_dunn$group <- "Fungi" 
p_richness_dunn <- rbind(dunn_richness_p_depth_df,dunn_richness_p_season_df,dunn_richness_p_LC1_df,dunn_richness_p_LU2_df,dunn_richness_p_erosion_df,dunn_richness_p_pH_df,dunn_p_richness_cluster_temp_df,dunn_p_richness_cluster_temp_range_df)
p_richness_dunn$group <- "Protists" 

dunn_richness_complete <- rbind(f_richness_dunn,p_richness_dunn,roti_richness_dunn,tardi_richness_dunn,nema_richness_dunn,arthro_richness_dunn,anneli_richness_dunn)
dunn_richness_complete$Diversity_index <- "Richness" 

dunn_complete <- rbind(dunn_shannon_complete,dunn_richness_complete)

write.csv(dunn_complete, file = "df_dunn_complete1805.csv", row.names = F)








#subgroup for median values of groups
sub_p_sample_data_CL_annual <- p_sample_data %>% filter(LC1_2018 == "CL_annual")
sum_p_sample_data_richness_CL_annual <- capture.output(summary(sub_p_sample_data_CL_annual$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_CL_annual <- as.data.frame(sum_p_sample_data_richness_CL_annual)
sum_p_sample_data_richness_CL_annual$group <- "Protists"
sum_p_sample_data_richness_CL_annual$diversity <- "Richness"
sum_p_sample_data_richness_CL_annual$class <- "Land Cover CL_annual"
sum_p_sample_data_richness_CL_annual <-sum_p_sample_data_richness_CL_annual[-1,]

sum_p_sample_data_shannon_CL_annual <- capture.output(summary(sub_p_sample_data_CL_annual$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_CL_annual <- as.data.frame(sum_p_sample_data_shannon_CL_annual)
sum_p_sample_data_shannon_CL_annual$group <- "Protists"
sum_p_sample_data_shannon_CL_annual$diversity <- "Shannon"
sum_p_sample_data_shannon_CL_annual$class <- "Land Cover CL_annual"
sum_p_sample_data_shannon_CL_annual <-sum_p_sample_data_shannon_CL_annual[-1,]

sub_p_sample_data_CL_permanent <- p_sample_data %>% filter(LC1_2018 == "CL_permanent")
sum_p_sample_data_richness_CL_permanent <- capture.output(summary(sub_p_sample_data_CL_permanent$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_CL_permanent <- as.data.frame(sum_p_sample_data_richness_CL_permanent)
sum_p_sample_data_richness_CL_permanent$group <- "Protists"
sum_p_sample_data_richness_CL_permanent$diversity <- "Richness"
sum_p_sample_data_richness_CL_permanent$class <- "Land Cover CL_permanent"
sum_p_sample_data_richness_CL_permanent <-sum_p_sample_data_richness_CL_permanent[-1,]

sum_p_sample_data_shannon_CL_permanent <- capture.output(summary(sub_p_sample_data_CL_permanent$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_CL_permanent <- as.data.frame(sum_p_sample_data_shannon_CL_permanent)
sum_p_sample_data_shannon_CL_permanent$group <- "Protists"
sum_p_sample_data_shannon_CL_permanent$diversity <- "Shannon"
sum_p_sample_data_shannon_CL_permanent$class <- "Land Cover CL_permanent"
sum_p_sample_data_shannon_CL_permanent <-sum_p_sample_data_shannon_CL_permanent[-1,]

sub_p_sample_data_GL_unmanaged <- p_sample_data %>% filter(LC1_2018 == "GL_unmanaged")
sum_p_sample_data_richness_GL_unmanaged <- capture.output(summary(sub_p_sample_data_GL_unmanaged$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_GL_unmanaged <- as.data.frame(sum_p_sample_data_richness_GL_unmanaged)
sum_p_sample_data_richness_GL_unmanaged$group <- "Protists"
sum_p_sample_data_richness_GL_unmanaged$diversity <- "Richness"
sum_p_sample_data_richness_GL_unmanaged$class <- "Land Cover GL_unmanaged"
sum_p_sample_data_richness_GL_unmanaged <-sum_p_sample_data_richness_GL_unmanaged[-1,]

sum_p_sample_data_shannon_GL_unmanaged <- capture.output(summary(sub_p_sample_data_GL_unmanaged$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_GL_unmanaged <- as.data.frame(sum_p_sample_data_shannon_GL_unmanaged)
sum_p_sample_data_shannon_GL_unmanaged$group <- "Protists"
sum_p_sample_data_shannon_GL_unmanaged$diversity <- "Shannon"
sum_p_sample_data_shannon_GL_unmanaged$class <- "Land Cover GL_unmanaged"
sum_p_sample_data_shannon_GL_unmanaged <-sum_p_sample_data_shannon_GL_unmanaged[-1,]

sub_p_sample_data_GL_managed <- p_sample_data %>% filter(LC1_2018 == "GL_managed")
sum_p_sample_data_richness_GL_managed <- capture.output(summary(sub_p_sample_data_GL_managed$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_GL_managed <- as.data.frame(sum_p_sample_data_richness_GL_managed)
sum_p_sample_data_richness_GL_managed$group <- "Protists"
sum_p_sample_data_richness_GL_managed$diversity <- "Richness"
sum_p_sample_data_richness_GL_managed$class <- "Land Cover GL_managed"
sum_p_sample_data_richness_GL_managed <-sum_p_sample_data_richness_GL_managed[-1,]

sum_p_sample_data_shannon_GL_managed <- capture.output(summary(sub_p_sample_data_GL_managed$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_GL_managed <- as.data.frame(sum_p_sample_data_shannon_GL_managed)
sum_p_sample_data_shannon_GL_managed$group <- "Protists"
sum_p_sample_data_shannon_GL_managed$diversity <- "Shannon"
sum_p_sample_data_shannon_GL_managed$class <- "Land Cover GL_managed"
sum_p_sample_data_shannon_GL_managed <-sum_p_sample_data_shannon_GL_managed[-1,]

sub_p_sample_data_WL_broadleaved <- p_sample_data %>% filter(LC1_2018 == "WL_broadleaved")
sum_p_sample_data_richness_WL_broadleaved <- capture.output(summary(sub_p_sample_data_WL_broadleaved$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_WL_broadleaved <- as.data.frame(sum_p_sample_data_richness_WL_broadleaved)
sum_p_sample_data_richness_WL_broadleaved$group <- "Protists"
sum_p_sample_data_richness_WL_broadleaved$diversity <- "Richness"
sum_p_sample_data_richness_WL_broadleaved$class <- "Land Cover WL_broadleaved"
sum_p_sample_data_richness_WL_broadleaved <-sum_p_sample_data_richness_WL_broadleaved[-1,]

sum_p_sample_data_shannon_WL_broadleaved <- capture.output(summary(sub_p_sample_data_WL_broadleaved$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_WL_broadleaved <- as.data.frame(sum_p_sample_data_shannon_WL_broadleaved)
sum_p_sample_data_shannon_WL_broadleaved$group <- "Protists"
sum_p_sample_data_shannon_WL_broadleaved$diversity <- "Shannon"
sum_p_sample_data_shannon_WL_broadleaved$class <- "Land Cover WL_broadleaved"
sum_p_sample_data_shannon_WL_broadleaved <-sum_p_sample_data_shannon_WL_broadleaved[-1,]

sub_p_sample_data_WL_coniferous <- p_sample_data %>% filter(LC1_2018 == "WL_coniferous")
sum_p_sample_data_richness_WL_coniferous <- capture.output(summary(sub_p_sample_data_WL_coniferous$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_WL_coniferous <- as.data.frame(sum_p_sample_data_richness_WL_coniferous)
sum_p_sample_data_richness_WL_coniferous$group <- "Protists"
sum_p_sample_data_richness_WL_coniferous$diversity <- "Richness"
sum_p_sample_data_richness_WL_coniferous$class <- "Land Cover WL_coniferous"
sum_p_sample_data_richness_WL_coniferous <-sum_p_sample_data_richness_WL_coniferous[-1,]

sum_p_sample_data_shannon_WL_coniferous <- capture.output(summary(sub_p_sample_data_WL_coniferous$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_WL_coniferous <- as.data.frame(sum_p_sample_data_shannon_WL_coniferous)
sum_p_sample_data_shannon_WL_coniferous$group <- "Protists"
sum_p_sample_data_shannon_WL_coniferous$diversity <- "Shannon"
sum_p_sample_data_shannon_WL_coniferous$class <- "Land Cover WL_coniferous"
sum_p_sample_data_shannon_WL_coniferous <-sum_p_sample_data_shannon_WL_coniferous[-1,]

sub_p_sample_data_temp_1 <- p_sample_data %>% filter(cluster_temp == "Very low")
sum_p_sample_data_shannon_temp_1 <- capture.output(summary(sub_p_sample_data_temp_1$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_temp_1 <- as.data.frame(sum_p_sample_data_shannon_temp_1)
sum_p_sample_data_shannon_temp_1$group <- "Protists"
sum_p_sample_data_shannon_temp_1$diversity <- "Shannon"
sum_p_sample_data_shannon_temp_1$class <- "Very low temperature"
sum_p_sample_data_shannon_temp_1 <-sum_p_sample_data_shannon_temp_1[-1,]

sum_p_sample_data_richness_temp_1 <- capture.output(summary(sub_p_sample_data_temp_1$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_temp_1 <- as.data.frame(sum_p_sample_data_richness_temp_1)
sum_p_sample_data_richness_temp_1$group <- "Protists"
sum_p_sample_data_richness_temp_1$diversity <- "Richness"
sum_p_sample_data_richness_temp_1$class <- "Very low temperature"
sum_p_sample_data_richness_temp_1 <-sum_p_sample_data_richness_temp_1[-1,]

sub_p_sample_data_temp_2 <- p_sample_data %>% filter(cluster_temp == "Low")
sum_p_sample_data_shannon_temp_2 <- capture.output(summary(sub_p_sample_data_temp_2$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_temp_2 <- as.data.frame(sum_p_sample_data_shannon_temp_2)
sum_p_sample_data_shannon_temp_2$group <- "Protists"
sum_p_sample_data_shannon_temp_2$diversity <- "Shannon"
sum_p_sample_data_shannon_temp_2$class <- "Low temperature"
sum_p_sample_data_shannon_temp_2 <-sum_p_sample_data_shannon_temp_2[-1,]

sum_p_sample_data_richness_temp_2 <- capture.output(summary(sub_p_sample_data_temp_2$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_temp_2 <- as.data.frame(sum_p_sample_data_richness_temp_2)
sum_p_sample_data_richness_temp_2$group <- "Protists"
sum_p_sample_data_richness_temp_2$diversity <- "Richness"
sum_p_sample_data_richness_temp_2$class <- "Low temperature"
sum_p_sample_data_richness_temp_2 <-sum_p_sample_data_richness_temp_2[-1,]

sub_p_sample_data_temp_3 <- p_sample_data %>% filter(cluster_temp == "Medium low")
sum_p_sample_data_shannon_temp_3 <- capture.output(summary(sub_p_sample_data_temp_3$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_temp_3 <- as.data.frame(sum_p_sample_data_shannon_temp_3)
sum_p_sample_data_shannon_temp_3$group <- "Protists"
sum_p_sample_data_shannon_temp_3$diversity <- "Shannon"
sum_p_sample_data_shannon_temp_3$class <- "Medium low temperature"
sum_p_sample_data_shannon_temp_3 <-sum_p_sample_data_shannon_temp_3[-1,]

sum_p_sample_data_richness_temp_3 <- capture.output(summary(sub_p_sample_data_temp_3$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_temp_3 <- as.data.frame(sum_p_sample_data_richness_temp_3)
sum_p_sample_data_richness_temp_3$group <- "Protists"
sum_p_sample_data_richness_temp_3$diversity <- "Richness"
sum_p_sample_data_richness_temp_3$class <- "Medium low temperature"
sum_p_sample_data_richness_temp_3 <-sum_p_sample_data_richness_temp_3[-1,]

sub_p_sample_data_temp_4 <- p_sample_data %>% filter(cluster_temp == "Medium high")
sum_p_sample_data_shannon_temp_4 <- capture.output(summary(sub_p_sample_data_temp_4$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_temp_4 <- as.data.frame(sum_p_sample_data_shannon_temp_4)
sum_p_sample_data_shannon_temp_4$group <- "Protists"
sum_p_sample_data_shannon_temp_4$diversity <- "Shannon"
sum_p_sample_data_shannon_temp_4$class <- "Medium high temperature"
sum_p_sample_data_shannon_temp_4 <-sum_p_sample_data_shannon_temp_4[-1,]

sum_p_sample_data_richness_temp_4 <- capture.output(summary(sub_p_sample_data_temp_4$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_temp_4 <- as.data.frame(sum_p_sample_data_richness_temp_4)
sum_p_sample_data_richness_temp_4$group <- "Protists"
sum_p_sample_data_richness_temp_4$diversity <- "Richness"
sum_p_sample_data_richness_temp_4$class <- "Medium high temperature"
sum_p_sample_data_richness_temp_4 <-sum_p_sample_data_richness_temp_4[-1,]

sub_p_sample_data_temp_5 <- p_sample_data %>% filter(cluster_temp == "High")
sum_p_sample_data_shannon_temp_5 <- capture.output(summary(sub_p_sample_data_temp_5$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_temp_5 <- as.data.frame(sum_p_sample_data_shannon_temp_5)
sum_p_sample_data_shannon_temp_5$group <- "Protists"
sum_p_sample_data_shannon_temp_5$diversity <- "Shannon"
sum_p_sample_data_shannon_temp_5$class <- "High temperature"
sum_p_sample_data_shannon_temp_5 <-sum_p_sample_data_shannon_temp_5[-1,]

sum_p_sample_data_richness_temp_5 <- capture.output(summary(sub_p_sample_data_temp_5$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_temp_5 <- as.data.frame(sum_p_sample_data_richness_temp_5)
sum_p_sample_data_richness_temp_5$group <- "Protists"
sum_p_sample_data_richness_temp_5$diversity <- "Richness"
sum_p_sample_data_richness_temp_5$class <- "High temperature"
sum_p_sample_data_richness_temp_5 <-sum_p_sample_data_richness_temp_5[-1,]

sub_p_sample_data_temp_6 <- p_sample_data %>% filter(cluster_temp == "Very high")
sum_p_sample_data_shannon_temp_6 <- capture.output(summary(sub_p_sample_data_temp_6$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_temp_6 <- as.data.frame(sum_p_sample_data_shannon_temp_6)
sum_p_sample_data_shannon_temp_6$group <- "Protists"
sum_p_sample_data_shannon_temp_6$diversity <- "Shannon"
sum_p_sample_data_shannon_temp_6$class <- "Very high temperature"
sum_p_sample_data_shannon_temp_6 <-sum_p_sample_data_shannon_temp_6[-1,]

sum_p_sample_data_richness_temp_6 <- capture.output(summary(sub_p_sample_data_temp_6$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_temp_6 <- as.data.frame(sum_p_sample_data_richness_temp_6)
sum_p_sample_data_richness_temp_6$group <- "Protists"
sum_p_sample_data_richness_temp_6$diversity <- "Richness"
sum_p_sample_data_richness_temp_6$class <- "Very high temperature"
sum_p_sample_data_richness_temp_6 <-sum_p_sample_data_richness_temp_6[-1,] 

sub_p_sample_data_temp_range_1 <- p_sample_data %>% filter(cluster_temp_range == "Very low")
sum_p_sample_data_shannon_temp_range_1 <- capture.output(summary(sub_p_sample_data_temp_range_1$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_temp_range_1 <- as.data.frame(sum_p_sample_data_shannon_temp_range_1)
sum_p_sample_data_shannon_temp_range_1$group <- "Protists"
sum_p_sample_data_shannon_temp_range_1$diversity <- "Shannon"
sum_p_sample_data_shannon_temp_range_1$class <- "Very low temperature range"
sum_p_sample_data_shannon_temp_range_1 <-sum_p_sample_data_shannon_temp_range_1[-1,]

sum_p_sample_data_richness_temp_range_1 <- capture.output(summary(sub_p_sample_data_temp_range_1$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_temp_range_1 <- as.data.frame(sum_p_sample_data_richness_temp_range_1)
sum_p_sample_data_richness_temp_range_1$group <- "Protists"
sum_p_sample_data_richness_temp_range_1$diversity <- "Richness"
sum_p_sample_data_richness_temp_range_1$class <- "Very low temperature range"
sum_p_sample_data_richness_temp_range_1 <-sum_p_sample_data_richness_temp_range_1[-1,]

sub_p_sample_data_temp_range_2 <- p_sample_data %>% filter(cluster_temp_range == "Low")
sum_p_sample_data_shannon_temp_range_2 <- capture.output(summary(sub_p_sample_data_temp_range_2$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_temp_range_2 <- as.data.frame(sum_p_sample_data_shannon_temp_range_2)
sum_p_sample_data_shannon_temp_range_2$group <- "Protists"
sum_p_sample_data_shannon_temp_range_2$diversity <- "Shannon"
sum_p_sample_data_shannon_temp_range_2$class <- "Low temperature range"
sum_p_sample_data_shannon_temp_range_2 <-sum_p_sample_data_shannon_temp_range_2[-1,]

sum_p_sample_data_richness_temp_range_2 <- capture.output(summary(sub_p_sample_data_temp_range_2$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_temp_range_2 <- as.data.frame(sum_p_sample_data_richness_temp_range_2)
sum_p_sample_data_richness_temp_range_2$group <- "Protists"
sum_p_sample_data_richness_temp_range_2$diversity <- "Richness"
sum_p_sample_data_richness_temp_range_2$class <- "Low temperature range"
sum_p_sample_data_richness_temp_range_2 <-sum_p_sample_data_richness_temp_range_2[-1,]

sub_p_sample_data_temp_range_3 <- p_sample_data %>% filter(cluster_temp_range == "Medium low")
sum_p_sample_data_shannon_temp_range_3 <- capture.output(summary(sub_p_sample_data_temp_range_3$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_temp_range_3 <- as.data.frame(sum_p_sample_data_shannon_temp_range_3)
sum_p_sample_data_shannon_temp_range_3$group <- "Protists"
sum_p_sample_data_shannon_temp_range_3$diversity <- "Shannon"
sum_p_sample_data_shannon_temp_range_3$class <- "Medium low temperature range"
sum_p_sample_data_shannon_temp_range_3 <-sum_p_sample_data_shannon_temp_range_3[-1,]

sum_p_sample_data_richness_temp_range_3 <- capture.output(summary(sub_p_sample_data_temp_range_3$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_temp_range_3 <- as.data.frame(sum_p_sample_data_richness_temp_range_3)
sum_p_sample_data_richness_temp_range_3$group <- "Protists"
sum_p_sample_data_richness_temp_range_3$diversity <- "Richness"
sum_p_sample_data_richness_temp_range_3$class <- "Medium low temperature range"
sum_p_sample_data_richness_temp_range_3 <-sum_p_sample_data_richness_temp_range_3[-1,]

sub_p_sample_data_temp_range_4 <- p_sample_data %>% filter(cluster_temp_range == "Medium high")
sum_p_sample_data_shannon_temp_range_4 <- capture.output(summary(sub_p_sample_data_temp_range_4$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_temp_range_4 <- as.data.frame(sum_p_sample_data_shannon_temp_range_4)
sum_p_sample_data_shannon_temp_range_4$group <- "Protists"
sum_p_sample_data_shannon_temp_range_4$diversity <- "Shannon"
sum_p_sample_data_shannon_temp_range_4$class <- "Medium high temperature range"
sum_p_sample_data_shannon_temp_range_4 <-sum_p_sample_data_shannon_temp_range_4[-1,]

sum_p_sample_data_richness_temp_range_4 <- capture.output(summary(sub_p_sample_data_temp_range_4$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_temp_range_4 <- as.data.frame(sum_p_sample_data_richness_temp_range_4)
sum_p_sample_data_richness_temp_range_4$group <- "Protists"
sum_p_sample_data_richness_temp_range_4$diversity <- "Richness"
sum_p_sample_data_richness_temp_range_4$class <- "Medium high temperature range"
sum_p_sample_data_richness_temp_range_4 <-sum_p_sample_data_richness_temp_range_4[-1,]

sub_p_sample_data_temp_range_5 <- p_sample_data %>% filter(cluster_temp_range == "High")
sum_p_sample_data_shannon_temp_range_5 <- capture.output(summary(sub_p_sample_data_temp_range_5$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_temp_range_5 <- as.data.frame(sum_p_sample_data_shannon_temp_range_5)
sum_p_sample_data_shannon_temp_range_5$group <- "Protists"
sum_p_sample_data_shannon_temp_range_5$diversity <- "Shannon"
sum_p_sample_data_shannon_temp_range_5$class <- "High temperature range"
sum_p_sample_data_shannon_temp_range_5 <-sum_p_sample_data_shannon_temp_range_5[-1,]

sum_p_sample_data_richness_temp_range_5 <- capture.output(summary(sub_p_sample_data_temp_range_5$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_temp_range_5 <- as.data.frame(sum_p_sample_data_richness_temp_range_5)
sum_p_sample_data_richness_temp_range_5$group <- "Protists"
sum_p_sample_data_richness_temp_range_5$diversity <- "Richness"
sum_p_sample_data_richness_temp_range_5$class <- "High temperature range"
sum_p_sample_data_richness_temp_range_5 <-sum_p_sample_data_richness_temp_range_5[-1,]

sub_p_sample_data_temp_range_6 <- p_sample_data %>% filter(cluster_temp_range == "Very high")
sum_p_sample_data_shannon_temp_range_6 <- capture.output(summary(sub_p_sample_data_temp_range_6$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_temp_range_6 <- as.data.frame(sum_p_sample_data_shannon_temp_range_6)
sum_p_sample_data_shannon_temp_range_6$group <- "Protists"
sum_p_sample_data_shannon_temp_range_6$diversity <- "Shannon"
sum_p_sample_data_shannon_temp_range_6$class <- "Very high temperature range"
sum_p_sample_data_shannon_temp_range_6 <-sum_p_sample_data_shannon_temp_range_6[-1,]

sum_p_sample_data_richness_temp_range_6 <- capture.output(summary(sub_p_sample_data_temp_range_6$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_temp_range_6 <- as.data.frame(sum_p_sample_data_richness_temp_range_6)
sum_p_sample_data_richness_temp_range_6$group <- "Protists"
sum_p_sample_data_richness_temp_range_6$diversity <- "Richness"
sum_p_sample_data_richness_temp_range_6$class <- "Very high temperature range"
sum_p_sample_data_richness_temp_range_6 <-sum_p_sample_data_richness_temp_range_6[-1,]

sub_p_sample_data_prec_1 <- p_sample_data %>% filter(cluster_prec == "Very low")
sum_p_sample_data_shannon_prec_1 <- capture.output(summary(sub_p_sample_data_prec_1$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_prec_1 <- as.data.frame(sum_p_sample_data_shannon_prec_1)
sum_p_sample_data_shannon_prec_1$group <- "Protists"
sum_p_sample_data_shannon_prec_1$diversity <- "Shannon"
sum_p_sample_data_shannon_prec_1$class <- "Very low precipitation "
sum_p_sample_data_shannon_prec_1 <-sum_p_sample_data_shannon_prec_1[-1,]

sum_p_sample_data_richness_prec_1 <- capture.output(summary(sub_p_sample_data_prec_1$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_prec_1 <- as.data.frame(sum_p_sample_data_richness_prec_1)
sum_p_sample_data_richness_prec_1$group <- "Protists"
sum_p_sample_data_richness_prec_1$diversity <- "Richness"
sum_p_sample_data_richness_prec_1$class <- "Very low precipitation "
sum_p_sample_data_richness_prec_1 <-sum_p_sample_data_richness_prec_1[-1,]

sub_p_sample_data_prec_2 <- p_sample_data %>% filter(cluster_prec == "Low")
sum_p_sample_data_shannon_prec_2 <- capture.output(summary(sub_p_sample_data_prec_2$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_prec_2 <- as.data.frame(sum_p_sample_data_shannon_prec_2)
sum_p_sample_data_shannon_prec_2$group <- "Protists"
sum_p_sample_data_shannon_prec_2$diversity <- "Shannon"
sum_p_sample_data_shannon_prec_2$class <- "Low precipitation"
sum_p_sample_data_shannon_prec_2 <-sum_p_sample_data_shannon_prec_2[-1,]

sum_p_sample_data_richness_prec_2 <- capture.output(summary(sub_p_sample_data_prec_2$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_prec_2 <- as.data.frame(sum_p_sample_data_richness_prec_2)
sum_p_sample_data_richness_prec_2$group <- "Protists"
sum_p_sample_data_richness_prec_2$diversity <- "Richness"
sum_p_sample_data_richness_prec_2$class <- "Low precipitation"
sum_p_sample_data_richness_prec_2 <-sum_p_sample_data_richness_prec_2[-1,]

sub_p_sample_data_prec_3 <- p_sample_data %>% filter(cluster_prec == "Medium low")
sum_p_sample_data_shannon_prec_3 <- capture.output(summary(sub_p_sample_data_prec_3$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_prec_3 <- as.data.frame(sum_p_sample_data_shannon_prec_3)
sum_p_sample_data_shannon_prec_3$group <- "Protists"
sum_p_sample_data_shannon_prec_3$diversity <- "Shannon"
sum_p_sample_data_shannon_prec_3$class <- "Medium low precipitation"
sum_p_sample_data_shannon_prec_3 <-sum_p_sample_data_shannon_prec_3[-1,]

sum_p_sample_data_richness_prec_3 <- capture.output(summary(sub_p_sample_data_prec_3$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_prec_3 <- as.data.frame(sum_p_sample_data_richness_prec_3)
sum_p_sample_data_richness_prec_3$group <- "Protists"
sum_p_sample_data_richness_prec_3$diversity <- "Richness"
sum_p_sample_data_richness_prec_3$class <- "Medium low precipitation"
sum_p_sample_data_richness_prec_3 <-sum_p_sample_data_richness_prec_3[-1,]

sub_p_sample_data_prec_4 <- p_sample_data %>% filter(cluster_prec == "Medium high")
sum_p_sample_data_shannon_prec_4 <- capture.output(summary(sub_p_sample_data_prec_4$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_prec_4 <- as.data.frame(sum_p_sample_data_shannon_prec_4)
sum_p_sample_data_shannon_prec_4$group <- "Protists"
sum_p_sample_data_shannon_prec_4$diversity <- "Shannon"
sum_p_sample_data_shannon_prec_4$class <- "Medium high precipitation"
sum_p_sample_data_shannon_prec_4 <-sum_p_sample_data_shannon_prec_4[-1,]

sum_p_sample_data_richness_prec_4 <- capture.output(summary(sub_p_sample_data_prec_4$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_prec_4 <- as.data.frame(sum_p_sample_data_richness_prec_4)
sum_p_sample_data_richness_prec_4$group <- "Protists"
sum_p_sample_data_richness_prec_4$diversity <- "Richness"
sum_p_sample_data_richness_prec_4$class <- "Medium high precipitation"
sum_p_sample_data_richness_prec_4 <-sum_p_sample_data_richness_prec_4[-1,]

sub_p_sample_data_prec_5 <- p_sample_data %>% filter(cluster_prec == "High")
sum_p_sample_data_shannon_prec_5 <- capture.output(summary(sub_p_sample_data_prec_5$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_prec_5 <- as.data.frame(sum_p_sample_data_shannon_prec_5)
sum_p_sample_data_shannon_prec_5$group <- "Protists"
sum_p_sample_data_shannon_prec_5$diversity <- "Shannon"
sum_p_sample_data_shannon_prec_5$class <- "High precipitation"
sum_p_sample_data_shannon_prec_5 <-sum_p_sample_data_shannon_prec_5[-1,]

sum_p_sample_data_richness_prec_5 <- capture.output(summary(sub_p_sample_data_prec_5$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_prec_5 <- as.data.frame(sum_p_sample_data_richness_prec_5)
sum_p_sample_data_richness_prec_5$group <- "Protists"
sum_p_sample_data_richness_prec_5$diversity <- "Richness"
sum_p_sample_data_richness_prec_5$class <- "High precipitation"
sum_p_sample_data_richness_prec_5 <-sum_p_sample_data_richness_prec_5[-1,]

sub_p_sample_data_prec_6 <- p_sample_data %>% filter(cluster_prec == "Very high")
sum_p_sample_data_shannon_prec_6 <- capture.output(summary(sub_p_sample_data_prec_6$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_prec_6 <- as.data.frame(sum_p_sample_data_shannon_prec_6)
sum_p_sample_data_shannon_prec_6$group <- "Protists"
sum_p_sample_data_shannon_prec_6$diversity <- "Shannon"
sum_p_sample_data_shannon_prec_6$class <- "Very high precipitation"
sum_p_sample_data_shannon_prec_6 <-sum_p_sample_data_shannon_prec_6[-1,]

sum_p_sample_data_richness_prec_6 <- capture.output(summary(sub_p_sample_data_prec_6$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_prec_6 <- as.data.frame(sum_p_sample_data_richness_prec_6)
sum_p_sample_data_richness_prec_6$group <- "Protists"
sum_p_sample_data_richness_prec_6$diversity <- "Richness"
sum_p_sample_data_richness_prec_6$class <- "Very high precipitation"
sum_p_sample_data_richness_prec_6 <-sum_p_sample_data_richness_prec_6[-1,] 

sub_f_sample_data_CL_annual <- f_sample_data %>% filter(LC1_2018 == "CL_annual")
sum_f_sample_data_richness_CL_annual <- capture.output(summary(sub_f_sample_data_CL_annual$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_CL_annual <- as.data.frame(sum_f_sample_data_richness_CL_annual)
sum_f_sample_data_richness_CL_annual$group <- "Fungi"
sum_f_sample_data_richness_CL_annual$diversity <- "Richness"
sum_f_sample_data_richness_CL_annual$class <- "Land Cover CL_annual"
sum_f_sample_data_richness_CL_annual <-sum_f_sample_data_richness_CL_annual[-1,]

sum_f_sample_data_shannon_CL_annual <- capture.output(summary(sub_f_sample_data_CL_annual$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_CL_annual <- as.data.frame(sum_f_sample_data_shannon_CL_annual)
sum_f_sample_data_shannon_CL_annual$group <- "Fungi"
sum_f_sample_data_shannon_CL_annual$diversity <- "Shannon"
sum_f_sample_data_shannon_CL_annual$class <- "Land Cover CL_annual"
sum_f_sample_data_shannon_CL_annual <-sum_f_sample_data_shannon_CL_annual[-1,]

sub_f_sample_data_CL_permanent <- f_sample_data %>% filter(LC1_2018 == "CL_permanent")
sum_f_sample_data_richness_CL_permanent <-capture.output(summary(sub_f_sample_data_CL_permanent$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_CL_permanent <- as.data.frame(sum_f_sample_data_richness_CL_permanent)
sum_f_sample_data_richness_CL_permanent$group <- "Fungi"
sum_f_sample_data_richness_CL_permanent$diversity <- "Richness"
sum_f_sample_data_richness_CL_permanent$class <- "Land Cover CL_permanent"
sum_f_sample_data_richness_CL_permanent <-sum_f_sample_data_richness_CL_permanent[-1,]

sum_f_sample_data_shannon_CL_permanent <- capture.output(summary(sub_f_sample_data_CL_permanent$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_CL_permanent <- as.data.frame(sum_f_sample_data_shannon_CL_permanent)
sum_f_sample_data_shannon_CL_permanent$group <- "Fungi"
sum_f_sample_data_shannon_CL_permanent$diversity <- "Shannon"
sum_f_sample_data_shannon_CL_permanent$class <- "Land Cover CL_permanent"
sum_f_sample_data_shannon_CL_permanent <-sum_f_sample_data_shannon_CL_permanent[-1,]

sub_f_sample_data_GL_unmanaged <- f_sample_data %>% filter(LC1_2018 == "GL_unmanaged")
sum_f_sample_data_richness_GL_unmanaged <-capture.output(summary(sub_f_sample_data_GL_unmanaged$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_GL_unmanaged <- as.data.frame(sum_f_sample_data_richness_GL_unmanaged)
sum_f_sample_data_richness_GL_unmanaged$group <- "Fungi"
sum_f_sample_data_richness_GL_unmanaged$diversity <- "Richness"
sum_f_sample_data_richness_GL_unmanaged$class <- "Land Cover GL_unmanaged"
sum_f_sample_data_richness_GL_unmanaged <-sum_f_sample_data_richness_GL_unmanaged[-1,]

sum_f_sample_data_shannon_GL_unmanaged <- capture.output(summary(sub_f_sample_data_GL_unmanaged$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_GL_unmanaged <- as.data.frame(sum_f_sample_data_shannon_GL_unmanaged)
sum_f_sample_data_shannon_GL_unmanaged$group <- "Fungi"
sum_f_sample_data_shannon_GL_unmanaged$diversity <- "Shannon"
sum_f_sample_data_shannon_GL_unmanaged$class <- "Land Cover GL_unmanaged"
sum_f_sample_data_shannon_GL_unmanaged <-sum_f_sample_data_shannon_GL_unmanaged[-1,]

sub_f_sample_data_GL_managed <- f_sample_data %>% filter(LC1_2018 == "GL_managed")
sum_f_sample_data_richness_GL_managed <- capture.output(summary(sub_f_sample_data_GL_managed$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_GL_managed <- as.data.frame(sum_f_sample_data_richness_GL_managed)
sum_f_sample_data_richness_GL_managed$group <- "Fungi"
sum_f_sample_data_richness_GL_managed$diversity <- "Richness"
sum_f_sample_data_richness_GL_managed$class <- "Land Cover GL_managed"
sum_f_sample_data_richness_GL_managed <-sum_f_sample_data_richness_GL_managed[-1,]

sum_f_sample_data_shannon_GL_managed <- capture.output(summary(sub_f_sample_data_GL_managed$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_GL_managed <- as.data.frame(sum_f_sample_data_shannon_GL_managed)
sum_f_sample_data_shannon_GL_managed$group <- "Fungi"
sum_f_sample_data_shannon_GL_managed$diversity <- "Shannon"
sum_f_sample_data_shannon_GL_managed$class <- "Land Cover GL_managed"
sum_f_sample_data_shannon_GL_managed <-sum_f_sample_data_shannon_GL_managed[-1,]

sub_f_sample_data_WL_broadleaved <- f_sample_data %>% filter(LC1_2018 == "WL_broadleaved")
sum_f_sample_data_richness_WL_broadleaved <- capture.output(summary(sub_f_sample_data_WL_broadleaved$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_WL_broadleaved <- as.data.frame(sum_f_sample_data_richness_WL_broadleaved)
sum_f_sample_data_richness_WL_broadleaved$group <- "Fungi"
sum_f_sample_data_richness_WL_broadleaved$diversity <- "Richness"
sum_f_sample_data_richness_WL_broadleaved$class <- "Land Cover WL_broadleaved"
sum_f_sample_data_richness_WL_broadleaved <-sum_f_sample_data_richness_WL_broadleaved[-1,]

sum_f_sample_data_shannon_WL_broadleaved <- capture.output(summary(sub_f_sample_data_WL_broadleaved$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_WL_broadleaved <- as.data.frame(sum_f_sample_data_shannon_WL_broadleaved)
sum_f_sample_data_shannon_WL_broadleaved$group <- "Fungi"
sum_f_sample_data_shannon_WL_broadleaved$diversity <- "Shannon"
sum_f_sample_data_shannon_WL_broadleaved$class <- "Land Cover WL_broadleaved"
sum_f_sample_data_shannon_WL_broadleaved <-sum_f_sample_data_shannon_WL_broadleaved[-1,]

sub_f_sample_data_WL_coniferous <- f_sample_data %>% filter(LC1_2018 == "WL_coniferous")
sum_f_sample_data_richness_WL_coniferous <- capture.output(summary(sub_f_sample_data_WL_coniferous$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_WL_coniferous <- as.data.frame(sum_f_sample_data_richness_WL_coniferous)
sum_f_sample_data_richness_WL_coniferous$group <- "Fungi"
sum_f_sample_data_richness_WL_coniferous$diversity <- "Richness"
sum_f_sample_data_richness_WL_coniferous$class <- "Land Cover WL_coniferous"
sum_f_sample_data_richness_WL_coniferous <-sum_f_sample_data_richness_WL_coniferous[-1,]

sum_f_sample_data_shannon_WL_coniferous <- capture.output(summary(sub_f_sample_data_WL_coniferous$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_WL_coniferous <- as.data.frame(sum_f_sample_data_shannon_WL_coniferous)
sum_f_sample_data_shannon_WL_coniferous$group <- "Fungi"
sum_f_sample_data_shannon_WL_coniferous$diversity <- "Shannon"
sum_f_sample_data_shannon_WL_coniferous$class <- "Land Cover WL_coniferous"
sum_f_sample_data_shannon_WL_coniferous <-sum_f_sample_data_shannon_WL_coniferous[-1,]

sub_f_sample_data_temp_1 <- f_sample_data %>% filter(cluster_temp == "Very low")
sum_f_sample_data_shannon_temp_1 <- capture.output(summary(sub_f_sample_data_temp_1$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_temp_1 <- as.data.frame(sum_f_sample_data_shannon_temp_1)
sum_f_sample_data_shannon_temp_1$group <- "Fungi"
sum_f_sample_data_shannon_temp_1$diversity <- "Shannon"
sum_f_sample_data_shannon_temp_1$class <- "Very low temperature"
sum_f_sample_data_shannon_temp_1 <-sum_f_sample_data_shannon_temp_1[-1,]

sum_f_sample_data_richness_temp_1 <- capture.output(summary(sub_f_sample_data_temp_1$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_temp_1 <- as.data.frame(sum_f_sample_data_richness_temp_1)
sum_f_sample_data_richness_temp_1$group <- "Fungi"
sum_f_sample_data_richness_temp_1$diversity <- "Richness"
sum_f_sample_data_richness_temp_1$class <- "Very low temperature"
sum_f_sample_data_richness_temp_1 <-sum_f_sample_data_richness_temp_1[-1,]

sub_f_sample_data_temp_2 <- p_sample_data %>% filter(cluster_temp == "Low")
sum_f_sample_data_shannon_temp_2 <- capture.output(summary(sub_f_sample_data_temp_2$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_temp_2 <- as.data.frame(sum_f_sample_data_shannon_temp_2)
sum_f_sample_data_shannon_temp_2$group <- "Fungi"
sum_f_sample_data_shannon_temp_2$diversity <- "Shannon"
sum_f_sample_data_shannon_temp_2$class <- "Low temperature"
sum_f_sample_data_shannon_temp_2 <-sum_f_sample_data_shannon_temp_2[-1,]

sum_f_sample_data_richness_temp_2 <- capture.output(summary(sub_f_sample_data_temp_2$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_temp_2 <- as.data.frame(sum_f_sample_data_richness_temp_2)
sum_f_sample_data_richness_temp_2$group <- "Fungi"
sum_f_sample_data_richness_temp_2$diversity <- "Richness"
sum_f_sample_data_richness_temp_2$class <- "Low temperature"
sum_f_sample_data_richness_temp_2 <-sum_f_sample_data_richness_temp_2[-1,]

sub_f_sample_data_temp_3 <- f_sample_data %>% filter(cluster_temp == "Medium low")
sum_f_sample_data_shannon_temp_3 <- capture.output(summary(sub_f_sample_data_temp_3$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_temp_3 <- as.data.frame(sum_f_sample_data_shannon_temp_3)
sum_f_sample_data_shannon_temp_3$group <- "Fungi"
sum_f_sample_data_shannon_temp_3$diversity <- "Shannon"
sum_f_sample_data_shannon_temp_3$class <- "Medium low temperature"
sum_f_sample_data_shannon_temp_3 <-sum_f_sample_data_shannon_temp_3[-1,]

sum_f_sample_data_richness_temp_3 <- capture.output(summary(sub_f_sample_data_temp_3$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_temp_3 <- as.data.frame(sum_f_sample_data_richness_temp_3)
sum_f_sample_data_richness_temp_3$group <- "Fungi"
sum_f_sample_data_richness_temp_3$diversity <- "Richness"
sum_f_sample_data_richness_temp_3$class <- "Medium low temperature"
sum_f_sample_data_richness_temp_3 <-sum_f_sample_data_richness_temp_3[-1,]

sub_f_sample_data_temp_4 <- f_sample_data %>% filter(cluster_temp == "Medium high")
sum_f_sample_data_shannon_temp_4 <- capture.output(summary(sub_f_sample_data_temp_4$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_temp_4 <- as.data.frame(sum_f_sample_data_shannon_temp_4)
sum_f_sample_data_shannon_temp_4$group <- "Fungi"
sum_f_sample_data_shannon_temp_4$diversity <- "Shannon"
sum_f_sample_data_shannon_temp_4$class <- "Medium high temperature"
sum_f_sample_data_shannon_temp_4 <-sum_f_sample_data_shannon_temp_4[-1,]

sum_f_sample_data_richness_temp_4 <- capture.output(summary(sub_f_sample_data_temp_4$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_temp_4 <- as.data.frame(sum_f_sample_data_richness_temp_4)
sum_f_sample_data_richness_temp_4$group <- "Fungi"
sum_f_sample_data_richness_temp_4$diversity <- "Richness"
sum_f_sample_data_richness_temp_4$class <- "Medium high temperature"
sum_f_sample_data_richness_temp_4 <-sum_f_sample_data_richness_temp_4[-1,]

sub_f_sample_data_temp_5 <- f_sample_data %>% filter(cluster_temp == "High")
sum_f_sample_data_shannon_temp_5 <- capture.output(summary(sub_f_sample_data_temp_5$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_temp_5 <- as.data.frame(sum_f_sample_data_shannon_temp_5)
sum_f_sample_data_shannon_temp_5$group <- "Fungi"
sum_f_sample_data_shannon_temp_5$diversity <- "Shannon"
sum_f_sample_data_shannon_temp_5$class <- "High temperature"
sum_f_sample_data_shannon_temp_5 <-sum_f_sample_data_shannon_temp_5[-1,]

sum_f_sample_data_richness_temp_5 <- capture.output(summary(sub_f_sample_data_temp_5$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_temp_5 <- as.data.frame(sum_f_sample_data_richness_temp_5)
sum_f_sample_data_richness_temp_5$group <- "Fungi"
sum_f_sample_data_richness_temp_5$diversity <- "Richness"
sum_f_sample_data_richness_temp_5$class <- "High temperature"
sum_f_sample_data_richness_temp_5 <-sum_f_sample_data_richness_temp_5[-1,]

sub_f_sample_data_temp_6 <- f_sample_data %>% filter(cluster_temp == "Very high")
sum_f_sample_data_shannon_temp_6 <- capture.output(summary(sub_f_sample_data_temp_6$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_temp_6 <- as.data.frame(sum_f_sample_data_shannon_temp_6)
sum_f_sample_data_shannon_temp_6$group <- "Fungi"
sum_f_sample_data_shannon_temp_6$diversity <- "Shannon"
sum_f_sample_data_shannon_temp_6$class <- "Very high temperature"
sum_f_sample_data_shannon_temp_6 <-sum_f_sample_data_shannon_temp_6[-1,]

sum_f_sample_data_richness_temp_6 <- capture.output(summary(sub_f_sample_data_temp_6$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_temp_6 <- as.data.frame(sum_f_sample_data_richness_temp_6)
sum_f_sample_data_richness_temp_6$group <- "Fungi"
sum_f_sample_data_richness_temp_6$diversity <- "Richness"
sum_f_sample_data_richness_temp_6$class <- "Very high temperature"
sum_f_sample_data_richness_temp_6 <-sum_f_sample_data_richness_temp_6[-1,] 

sub_f_sample_data_temp_range_1 <- f_sample_data %>% filter(cluster_temp_range == "Very low")
sum_f_sample_data_shannon_temp_range_1 <- capture.output(summary(sub_f_sample_data_temp_range_1$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_temp_range_1 <- as.data.frame(sum_f_sample_data_shannon_temp_range_1)
sum_f_sample_data_shannon_temp_range_1$group <- "Fungi"
sum_f_sample_data_shannon_temp_range_1$diversity <- "Shannon"
sum_f_sample_data_shannon_temp_range_1$class <- "Very low temperature range"
sum_f_sample_data_shannon_temp_range_1 <-sum_f_sample_data_shannon_temp_range_1[-1,]

sum_f_sample_data_richness_temp_range_1 <- capture.output(summary(sub_f_sample_data_temp_range_1$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_temp_range_1 <- as.data.frame(sum_f_sample_data_richness_temp_range_1)
sum_f_sample_data_richness_temp_range_1$group <- "Fungi"
sum_f_sample_data_richness_temp_range_1$diversity <- "Richness"
sum_f_sample_data_richness_temp_range_1$class <- "Very low temperature range"
sum_f_sample_data_richness_temp_range_1 <-sum_f_sample_data_richness_temp_range_1[-1,]

sub_f_sample_data_temp_range_2 <- p_sample_data %>% filter(cluster_temp_range == "Low")
sum_f_sample_data_shannon_temp_range_2 <- capture.output(summary(sub_f_sample_data_temp_range_2$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_temp_range_2 <- as.data.frame(sum_f_sample_data_shannon_temp_range_2)
sum_f_sample_data_shannon_temp_range_2$group <- "Fungi"
sum_f_sample_data_shannon_temp_range_2$diversity <- "Shannon"
sum_f_sample_data_shannon_temp_range_2$class <- "Low temperature range"
sum_f_sample_data_shannon_temp_range_2 <-sum_f_sample_data_shannon_temp_range_2[-1,]

sum_f_sample_data_richness_temp_range_2 <- capture.output(summary(sub_f_sample_data_temp_range_2$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_temp_range_2 <- as.data.frame(sum_f_sample_data_richness_temp_range_2)
sum_f_sample_data_richness_temp_range_2$group <- "Fungi"
sum_f_sample_data_richness_temp_range_2$diversity <- "Richness"
sum_f_sample_data_richness_temp_range_2$class <- "Low temperature range"
sum_f_sample_data_richness_temp_range_2 <-sum_f_sample_data_richness_temp_range_2[-1,]

sub_f_sample_data_temp_range_3 <- f_sample_data %>% filter(cluster_temp_range == "Medium low")
sum_f_sample_data_shannon_temp_range_3 <- capture.output(summary(sub_f_sample_data_temp_range_3$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_temp_range_3 <- as.data.frame(sum_f_sample_data_shannon_temp_range_3)
sum_f_sample_data_shannon_temp_range_3$group <- "Fungi"
sum_f_sample_data_shannon_temp_range_3$diversity <- "Shannon"
sum_f_sample_data_shannon_temp_range_3$class <- "Medium low temperature range"
sum_f_sample_data_shannon_temp_range_3 <-sum_f_sample_data_shannon_temp_range_3[-1,]

sum_f_sample_data_richness_temp_range_3 <- capture.output(summary(sub_f_sample_data_temp_range_3$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_temp_range_3 <- as.data.frame(sum_f_sample_data_richness_temp_range_3)
sum_f_sample_data_richness_temp_range_3$group <- "Fungi"
sum_f_sample_data_richness_temp_range_3$diversity <- "Richness"
sum_f_sample_data_richness_temp_range_3$class <- "Medium low temperature range"
sum_f_sample_data_richness_temp_range_3 <-sum_f_sample_data_richness_temp_range_3[-1,]

sub_f_sample_data_temp_range_4 <- f_sample_data %>% filter(cluster_temp_range == "Medium high")
sum_f_sample_data_shannon_temp_range_4 <- capture.output(summary(sub_f_sample_data_temp_range_4$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_temp_range_4 <- as.data.frame(sum_f_sample_data_shannon_temp_range_4)
sum_f_sample_data_shannon_temp_range_4$group <- "Fungi"
sum_f_sample_data_shannon_temp_range_4$diversity <- "Shannon"
sum_f_sample_data_shannon_temp_range_4$class <- "Medium high temperature range"
sum_f_sample_data_shannon_temp_range_4 <-sum_f_sample_data_shannon_temp_range_4[-1,]

sum_f_sample_data_richness_temp_range_4 <- capture.output(summary(sub_f_sample_data_temp_range_4$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_temp_range_4 <- as.data.frame(sum_f_sample_data_richness_temp_range_4)
sum_f_sample_data_richness_temp_range_4$group <- "Fungi"
sum_f_sample_data_richness_temp_range_4$diversity <- "Richness"
sum_f_sample_data_richness_temp_range_4$class <- "Medium high temperature range"
sum_f_sample_data_richness_temp_range_4 <-sum_f_sample_data_richness_temp_range_4[-1,]

sub_f_sample_data_temp_range_5 <- f_sample_data %>% filter(cluster_temp_range == "High")
sum_f_sample_data_shannon_temp_range_5 <- capture.output(summary(sub_f_sample_data_temp_range_5$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_temp_range_5 <- as.data.frame(sum_f_sample_data_shannon_temp_range_5)
sum_f_sample_data_shannon_temp_range_5$group <- "Fungi"
sum_f_sample_data_shannon_temp_range_5$diversity <- "Shannon"
sum_f_sample_data_shannon_temp_range_5$class <- "High temperature range"
sum_f_sample_data_shannon_temp_range_5 <-sum_f_sample_data_shannon_temp_range_5[-1,]

sum_f_sample_data_richness_temp_range_5 <- capture.output(summary(sub_f_sample_data_temp_range_5$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_temp_range_5 <- as.data.frame(sum_f_sample_data_richness_temp_range_5)
sum_f_sample_data_richness_temp_range_5$group <- "Fungi"
sum_f_sample_data_richness_temp_range_5$diversity <- "Richness"
sum_f_sample_data_richness_temp_range_5$class <- "High temperature range"
sum_f_sample_data_richness_temp_range_5 <-sum_f_sample_data_richness_temp_range_5[-1,]

sub_f_sample_data_temp_range_6 <- f_sample_data %>% filter(cluster_temp_range == "Very high")
sum_f_sample_data_shannon_temp_range_6 <- capture.output(summary(sub_f_sample_data_temp_range_6$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_temp_range_6 <- as.data.frame(sum_f_sample_data_shannon_temp_range_6)
sum_f_sample_data_shannon_temp_range_6$group <- "Fungi"
sum_f_sample_data_shannon_temp_range_6$diversity <- "Shannon"
sum_f_sample_data_shannon_temp_range_6$class <- "Very high temperature range"
sum_f_sample_data_shannon_temp_range_6 <-sum_f_sample_data_shannon_temp_range_6[-1,]

sum_f_sample_data_richness_temp_range_6 <- capture.output(summary(sub_f_sample_data_temp_range_6$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_temp_range_6 <- as.data.frame(sum_f_sample_data_richness_temp_range_6)
sum_f_sample_data_richness_temp_range_6$group <- "Fungi"
sum_f_sample_data_richness_temp_range_6$diversity <- "Richness"
sum_f_sample_data_richness_temp_range_6$class <- "Very high temperature range"
sum_f_sample_data_richness_temp_range_6 <-sum_f_sample_data_richness_temp_range_6[-1,]

sub_f_sample_data_prec_1 <- f_sample_data %>% filter(cluster_prec == "Very low")
sum_f_sample_data_shannon_prec_1 <- capture.output(summary(sub_f_sample_data_prec_1$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_prec_1 <- as.data.frame(sum_f_sample_data_shannon_prec_1)
sum_f_sample_data_shannon_prec_1$group <- "Fungi"
sum_f_sample_data_shannon_prec_1$diversity <- "Shannon"
sum_f_sample_data_shannon_prec_1$class <- "Very low precipitation"
sum_f_sample_data_shannon_prec_1 <-sum_f_sample_data_shannon_prec_1[-1,]

sum_f_sample_data_richness_prec_1 <- capture.output(summary(sub_f_sample_data_prec_1$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_prec_1 <- as.data.frame(sum_f_sample_data_richness_prec_1)
sum_f_sample_data_richness_prec_1$group <- "Fungi"
sum_f_sample_data_richness_prec_1$diversity <- "Richness"
sum_f_sample_data_richness_prec_1$class <- "Very low precipitation"
sum_f_sample_data_richness_prec_1 <-sum_f_sample_data_richness_prec_1[-1,]

sub_f_sample_data_prec_2 <- p_sample_data %>% filter(cluster_prec == "Low")
sum_f_sample_data_shannon_prec_2 <- capture.output(summary(sub_f_sample_data_prec_2$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_prec_2 <- as.data.frame(sum_f_sample_data_shannon_prec_2)
sum_f_sample_data_shannon_prec_2$group <- "Fungi"
sum_f_sample_data_shannon_prec_2$diversity <- "Shannon"
sum_f_sample_data_shannon_prec_2$class <- "Low precipitation"
sum_f_sample_data_shannon_prec_2 <-sum_f_sample_data_shannon_prec_2[-1,]

sum_f_sample_data_richness_prec_2 <- capture.output(summary(sub_f_sample_data_prec_2$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_prec_2 <- as.data.frame(sum_f_sample_data_richness_prec_2)
sum_f_sample_data_richness_prec_2$group <- "Fungi"
sum_f_sample_data_richness_prec_2$diversity <- "Richness"
sum_f_sample_data_richness_prec_2$class <- "Low precipitation"
sum_f_sample_data_richness_prec_2 <-sum_f_sample_data_richness_prec_2[-1,]

sub_f_sample_data_prec_3 <- f_sample_data %>% filter(cluster_prec == "Medium low")
sum_f_sample_data_shannon_prec_3 <- capture.output(summary(sub_f_sample_data_prec_3$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_prec_3 <- as.data.frame(sum_f_sample_data_shannon_prec_3)
sum_f_sample_data_shannon_prec_3$group <- "Fungi"
sum_f_sample_data_shannon_prec_3$diversity <- "Shannon"
sum_f_sample_data_shannon_prec_3$class <- "Medium low precipitation"
sum_f_sample_data_shannon_prec_3 <-sum_f_sample_data_shannon_prec_3[-1,]

sum_f_sample_data_richness_prec_3 <- capture.output(summary(sub_f_sample_data_prec_3$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_prec_3 <- as.data.frame(sum_f_sample_data_richness_prec_3)
sum_f_sample_data_richness_prec_3$group <- "Fungi"
sum_f_sample_data_richness_prec_3$diversity <- "Richness"
sum_f_sample_data_richness_prec_3$class <- "Medium low precipitation"
sum_f_sample_data_richness_prec_3 <-sum_f_sample_data_richness_prec_3[-1,]

sub_f_sample_data_prec_4 <- f_sample_data %>% filter(cluster_prec == "Medium high")
sum_f_sample_data_shannon_prec_4 <- capture.output(summary(sub_f_sample_data_prec_4$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_prec_4 <- as.data.frame(sum_f_sample_data_shannon_prec_4)
sum_f_sample_data_shannon_prec_4$group <- "Fungi"
sum_f_sample_data_shannon_prec_4$diversity <- "Shannon"
sum_f_sample_data_shannon_prec_4$class <- "Medium high precipitation"
sum_f_sample_data_shannon_prec_4 <-sum_f_sample_data_shannon_prec_4[-1,]

sum_f_sample_data_richness_prec_4 <- capture.output(summary(sub_f_sample_data_prec_4$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_prec_4 <- as.data.frame(sum_f_sample_data_richness_prec_4)
sum_f_sample_data_richness_prec_4$group <- "Fungi"
sum_f_sample_data_richness_prec_4$diversity <- "Richness"
sum_f_sample_data_richness_prec_4$class <- "Medium high precipitation"
sum_f_sample_data_richness_prec_4 <-sum_f_sample_data_richness_prec_4[-1,]

sub_f_sample_data_prec_5 <- f_sample_data %>% filter(cluster_prec == "High")
sum_f_sample_data_shannon_prec_5 <- capture.output(summary(sub_f_sample_data_prec_5$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_prec_5 <- as.data.frame(sum_f_sample_data_shannon_prec_5)
sum_f_sample_data_shannon_prec_5$group <- "Fungi"
sum_f_sample_data_shannon_prec_5$diversity <- "Shannon"
sum_f_sample_data_shannon_prec_5$class <- "High precipitation"
sum_f_sample_data_shannon_prec_5 <-sum_f_sample_data_shannon_prec_5[-1,]

sum_f_sample_data_richness_prec_5 <- capture.output(summary(sub_f_sample_data_prec_5$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_prec_5 <- as.data.frame(sum_f_sample_data_richness_prec_5)
sum_f_sample_data_richness_prec_5$group <- "Fungi"
sum_f_sample_data_richness_prec_5$diversity <- "Richness"
sum_f_sample_data_richness_prec_5$class <- "High precipitation"
sum_f_sample_data_richness_prec_5 <-sum_f_sample_data_richness_prec_5[-1,]

sub_f_sample_data_prec_6 <- f_sample_data %>% filter(cluster_prec == "Very high")
sum_f_sample_data_shannon_prec_6 <- capture.output(summary(sub_f_sample_data_prec_6$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_prec_6 <- as.data.frame(sum_f_sample_data_shannon_prec_6)
sum_f_sample_data_shannon_prec_6$group <- "Fungi"
sum_f_sample_data_shannon_prec_6$diversity <- "Shannon"
sum_f_sample_data_shannon_prec_6$class <- "Very high precipitation"
sum_f_sample_data_shannon_prec_6 <-sum_f_sample_data_shannon_prec_6[-1,]

sum_f_sample_data_richness_prec_6 <- capture.output(summary(sub_f_sample_data_prec_6$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_prec_6 <- as.data.frame(sum_f_sample_data_richness_prec_6)
sum_f_sample_data_richness_prec_6$group <- "Fungi"
sum_f_sample_data_richness_prec_6$diversity <- "Richness"
sum_f_sample_data_richness_prec_6$class <- "Very high precipitation"
sum_f_sample_data_richness_prec_6 <-sum_f_sample_data_richness_prec_6[-1,] 

sub_a_roti_sample_data_CL_annual <- a_roti_sample_data %>% filter(LC1_2018 == "CL_annual")
sum_a_roti_richness_CL_annual <- capture.output(summary(sub_a_roti_sample_data_CL_annual$Richness),file=NULL,append=FALSE)
sum_a_roti_richness_CL_annual <- as.data.frame(sum_a_roti_richness_CL_annual)
sum_a_roti_richness_CL_annual$group <- "Rotifers"
sum_a_roti_richness_CL_annual$diversity <- "Richness"
sum_a_roti_richness_CL_annual$class <- "Land Cover CL_annual"
sum_a_roti_richness_CL_annual <-sum_a_roti_richness_CL_annual[-1,]

sum_a_roti_shannon_CL_annual <- capture.output(summary(sub_a_roti_sample_data_CL_annual$Shannon),file=NULL,append=FALSE)
sum_a_roti_shannon_CL_annual <- as.data.frame(sum_a_roti_shannon_CL_annual)
sum_a_roti_shannon_CL_annual$group <- "Rotifers"
sum_a_roti_shannon_CL_annual$diversity <- "Shannon"
sum_a_roti_shannon_CL_annual$class <- "Land Cove CL_annualr"
sum_a_roti_shannon_CL_annual <-sum_a_roti_shannon_CL_annual[-1,]

sum_a_roti_sample_data_CL_permanent <- a_roti_sample_data %>% filter(LC1_2018 == "CL_permanent")
sum_a_roti_richness_CL_permanent <- capture.output(summary(sum_a_roti_sample_data_CL_permanent$Richness),file=NULL,append=FALSE)
sum_a_roti_richness_CL_permanent <- as.data.frame(sum_a_roti_richness_CL_permanent)
sum_a_roti_richness_CL_permanent$group <- "Rotifers"
sum_a_roti_richness_CL_permanent$diversity <- "Richness"
sum_a_roti_richness_CL_permanent$class <- "Land Cover CL_permanent"
sum_a_roti_richness_CL_permanent <-sum_a_roti_richness_CL_permanent[-1,]

sum_a_roti_shannon_CL_permanent <- capture.output(summary(sum_a_roti_sample_data_CL_permanent$Shannon),file=NULL,append=FALSE)
sum_a_roti_shannon_CL_permanent <- as.data.frame(sum_a_roti_shannon_CL_permanent)
sum_a_roti_shannon_CL_permanent$group <- "Rotifers"
sum_a_roti_shannon_CL_permanent$diversity <- "Shannon"
sum_a_roti_shannon_CL_permanent$class <- "Land Cover CL_permanent"
sum_a_roti_shannon_CL_permanent <-sum_a_roti_shannon_CL_permanent[-1,]

sub_a_roti_sample_data_GL_unmanaged <- a_roti_sample_data %>% filter(LC1_2018 == "GL_unmanaged")
sum_a_roti_richness_GL_unmanaged <- capture.output(summary(sub_a_roti_sample_data_GL_unmanaged$Richness),file=NULL,append=FALSE)
sum_a_roti_richness_GL_unmanaged <- as.data.frame(sum_a_roti_richness_GL_unmanaged)
sum_a_roti_richness_GL_unmanaged$group <- "Rotifers"
sum_a_roti_richness_GL_unmanaged$diversity <- "Richness"
sum_a_roti_richness_GL_unmanaged$class <- "Land Cover GL_unmanaged"
sum_a_roti_richness_GL_unmanaged <-sum_a_roti_richness_GL_unmanaged[-1,]

sum_a_roti_shannon_GL_unmanaged <- capture.output(summary(sub_a_roti_sample_data_GL_unmanaged$Shannon),file=NULL,append=FALSE)
sum_a_roti_shannon_GL_unmanaged <- as.data.frame(sum_a_roti_shannon_GL_unmanaged)
sum_a_roti_shannon_GL_unmanaged$group <- "Rotifers"
sum_a_roti_shannon_GL_unmanaged$diversity <- "Shannon"
sum_a_roti_shannon_GL_unmanaged$class <- "Land Cover GL_unmanaged"
sum_a_roti_shannon_GL_unmanaged <-sum_a_roti_shannon_GL_unmanaged[-1,]

sub_a_roti_sample_data_GL_managed <- a_roti_sample_data %>% filter(LC1_2018 == "GL_managed")
sum_a_roti_richness_GL_managed <- capture.output(summary(sub_a_roti_sample_data_GL_managed$Richness),file=NULL,append=FALSE)
sum_a_roti_richness_GL_managed <- as.data.frame(sum_a_roti_richness_GL_managed)
sum_a_roti_richness_GL_managed$group <- "Rotifers"
sum_a_roti_richness_GL_managed$diversity <- "Richness"
sum_a_roti_richness_GL_managed$class <- "Land Cover GL_managed"
sum_a_roti_richness_GL_managed <-sum_a_roti_richness_GL_managed[-1,]

sum_a_roti_shannon_GL_managed <- capture.output(summary(sub_a_roti_sample_data_GL_managed$Shannon),file=NULL,append=FALSE)
sum_a_roti_shannon_GL_managed <- as.data.frame(sum_a_roti_shannon_GL_managed)
sum_a_roti_shannon_GL_managed$group <- "Rotifers"
sum_a_roti_shannon_GL_managed$diversity <- "Shannon"
sum_a_roti_shannon_GL_managed$class <- "Land Cover GL_managed"
sum_a_roti_shannon_GL_managed <-sum_a_roti_shannon_GL_managed[-1,]

sub_a_roti_sample_data_WL_broadleaved <- a_roti_sample_data %>% filter(LC1_2018 == "WL_broadleaved")
sum_a_roti_richness_WL_broadleaved <- capture.output(summary(sub_a_roti_sample_data_WL_broadleaved$Richness),file=NULL,append=FALSE)
sum_a_roti_richness_WL_broadleaved <- as.data.frame(sum_a_roti_richness_WL_broadleaved)
sum_a_roti_richness_WL_broadleaved$group <- "Rotifers"
sum_a_roti_richness_WL_broadleaved$diversity <- "Richness"
sum_a_roti_richness_WL_broadleaved$class <- "Land Cover WL_broadleaved"
sum_a_roti_richness_WL_broadleaved <-sum_a_roti_richness_WL_broadleaved[-1,]

sum_a_roti_shannon_WL_broadleaved <- capture.output(summary(sub_a_roti_sample_data_WL_broadleaved$Shannon),file=NULL,append=FALSE)
sum_a_roti_shannon_WL_broadleaved <- as.data.frame(sum_a_roti_shannon_WL_broadleaved)
sum_a_roti_shannon_WL_broadleaved$group <- "Rotifers"
sum_a_roti_shannon_WL_broadleaved$diversity <- "Shannon"
sum_a_roti_shannon_WL_broadleaved$class <- "Land Cover WL_broadleaved"
sum_a_roti_shannon_WL_broadleaved <-sum_a_roti_shannon_WL_broadleaved[-1,]

sub_a_roti_sample_data_WL_coniferous <- a_roti_sample_data %>% filter(LC1_2018 == "WL_coniferous")
sum_a_roti_richness_WL_coniferous <- capture.output(summary(sub_a_roti_sample_data_WL_coniferous$Richness),file=NULL,append=FALSE)
sum_a_roti_richness_WL_coniferous <- as.data.frame(sum_a_roti_richness_WL_coniferous)
sum_a_roti_richness_WL_coniferous$group <- "Rotifers"
sum_a_roti_richness_WL_coniferous$diversity <- "Richness"
sum_a_roti_richness_WL_coniferous$class <- "Land Cover WL_coniferous"
sum_a_roti_richness_WL_coniferous <-sum_a_roti_richness_WL_coniferous[-1,]

sum_a_roti_shannon_WL_coniferous <- capture.output(summary(sub_a_roti_sample_data_WL_coniferous$Shannon),file=NULL,append=FALSE)
sum_a_roti_shannon_WL_coniferous <- as.data.frame(sum_a_roti_shannon_WL_coniferous)
sum_a_roti_shannon_WL_coniferous$group <- "Rotifers"
sum_a_roti_shannon_WL_coniferous$diversity <- "Shannon"
sum_a_roti_shannon_WL_coniferous$class <- "Land Cover WL_coniferous"
sum_a_roti_shannon_WL_coniferous <-sum_a_roti_shannon_WL_coniferous[-1,]

sub_a_roti_sample_data_temp_1 <- a_roti_sample_data %>% filter(cluster_temp == "Very low")
sum_a_roti_sample_data_shannon_temp_1 <- capture.output(summary(sub_a_roti_sample_data_temp_1$Shannon),file=NULL,append=FALSE)
sum_a_roti_sample_data_shannon_temp_1 <- as.data.frame(sum_a_roti_sample_data_shannon_temp_1)
sum_a_roti_sample_data_shannon_temp_1$group <- "Rotifers"
sum_a_roti_sample_data_shannon_temp_1$diversity <- "Shannon"
sum_a_roti_sample_data_shannon_temp_1$class <- "Very low temperature"
sum_a_roti_sample_data_shannon_temp_1 <-sum_a_roti_sample_data_shannon_temp_1[-1,]

sum_a_roti_sample_data_richness_temp_1 <- capture.output(summary(sub_a_roti_sample_data_temp_1$Richness),file=NULL,append=FALSE)
sum_a_roti_sample_data_richness_temp_1 <- as.data.frame(sum_a_roti_sample_data_richness_temp_1)
sum_a_roti_sample_data_richness_temp_1$group <- "Rotifers"
sum_a_roti_sample_data_richness_temp_1$diversity <- "Richness"
sum_a_roti_sample_data_richness_temp_1$class <- "Very low temperature"
sum_a_roti_sample_data_richness_temp_1 <-sum_a_roti_sample_data_richness_temp_1[-1,]

sub_a_roti_sample_data_temp_2 <- a_roti_sample_data %>% filter(cluster_temp == "Low")
sum_a_roti_sample_data_shannon_temp_2 <- capture.output(summary(sub_a_roti_sample_data_temp_2$Shannon),file=NULL,append=FALSE)
sum_a_roti_sample_data_shannon_temp_2 <- as.data.frame(sum_a_roti_sample_data_shannon_temp_2)
sum_a_roti_sample_data_shannon_temp_2$group <- "Rotifers"
sum_a_roti_sample_data_shannon_temp_2$diversity <- "Shannon"
sum_a_roti_sample_data_shannon_temp_2$class <- "Low temperature"
sum_a_roti_sample_data_shannon_temp_2 <-sum_a_roti_sample_data_shannon_temp_2[-1,]

sum_a_roti_sample_data_richness_temp_2 <- capture.output(summary(sub_a_roti_sample_data_temp_2$Richness),file=NULL,append=FALSE)
sum_a_roti_sample_data_richness_temp_2 <- as.data.frame(sum_a_roti_sample_data_richness_temp_2)
sum_a_roti_sample_data_richness_temp_2$group <- "Rotifers"
sum_a_roti_sample_data_richness_temp_2$diversity <- "Richness"
sum_a_roti_sample_data_richness_temp_2$class <- "Low temperature"
sum_a_roti_sample_data_richness_temp_2 <-sum_a_roti_sample_data_richness_temp_2[-1,]

sub_a_roti_sample_data_temp_3 <- a_roti_sample_data %>% filter(cluster_temp == "Medium low")
sum_a_roti_sample_data_shannon_temp_3 <- capture.output(summary(sub_a_roti_sample_data_temp_3$Shannon),file=NULL,append=FALSE)
sum_a_roti_sample_data_shannon_temp_3 <- as.data.frame(sum_a_roti_sample_data_shannon_temp_3)
sum_a_roti_sample_data_shannon_temp_3$group <- "Rotifers"
sum_a_roti_sample_data_shannon_temp_3$diversity <- "Shannon"
sum_a_roti_sample_data_shannon_temp_3$class <- "Medium low temperature"
sum_a_roti_sample_data_shannon_temp_3 <-sum_a_roti_sample_data_shannon_temp_3[-1,]

sum_a_roti_sample_data_richness_temp_3 <- capture.output(summary(sub_a_roti_sample_data_temp_3$Richness),file=NULL,append=FALSE)
sum_a_roti_sample_data_richness_temp_3 <- as.data.frame(sum_a_roti_sample_data_richness_temp_3)
sum_a_roti_sample_data_richness_temp_3$group <- "Rotifers"
sum_a_roti_sample_data_richness_temp_3$diversity <- "Richness"
sum_a_roti_sample_data_richness_temp_3$class <- "Medium low temperature"
sum_a_roti_sample_data_richness_temp_3 <-sum_a_roti_sample_data_richness_temp_3[-1,]

sub_a_roti_sample_data_temp_4 <- a_roti_sample_data %>% filter(cluster_temp == "Medium high")
sum_a_roti_sample_data_shannon_temp_4 <- capture.output(summary(sub_a_roti_sample_data_temp_4$Shannon),file=NULL,append=FALSE)
sum_a_roti_sample_data_shannon_temp_4 <- as.data.frame(sum_a_roti_sample_data_shannon_temp_4)
sum_a_roti_sample_data_shannon_temp_4$group <- "Rotifers"
sum_a_roti_sample_data_shannon_temp_4$diversity <- "Shannon"
sum_a_roti_sample_data_shannon_temp_4$class <- "Medium high temperature"
sum_a_roti_sample_data_shannon_temp_4 <-sum_a_roti_sample_data_shannon_temp_4[-1,]

sum_a_roti_sample_data_richness_temp_4 <- capture.output(summary(sub_a_roti_sample_data_temp_4$Richness),file=NULL,append=FALSE)
sum_a_roti_sample_data_richness_temp_4 <- as.data.frame(sum_a_roti_sample_data_richness_temp_4)
sum_a_roti_sample_data_richness_temp_4$group <- "Rotifers"
sum_a_roti_sample_data_richness_temp_4$diversity <- "Richness"
sum_a_roti_sample_data_richness_temp_4$class <- "Medium high temperature"
sum_a_roti_sample_data_richness_temp_4 <-sum_a_roti_sample_data_richness_temp_4[-1,]

sub_a_roti_sample_data_temp_5 <- a_roti_sample_data %>% filter(cluster_temp == "High")
sum_a_roti_sample_data_shannon_temp_5 <- capture.output(summary(sub_a_roti_sample_data_temp_5$Shannon),file=NULL,append=FALSE)
sum_a_roti_sample_data_shannon_temp_5 <- as.data.frame(sum_a_roti_sample_data_shannon_temp_5)
sum_a_roti_sample_data_shannon_temp_5$group <- "Rotifers"
sum_a_roti_sample_data_shannon_temp_5$diversity <- "Shannon"
sum_a_roti_sample_data_shannon_temp_5$class <- "High temperature"
sum_a_roti_sample_data_shannon_temp_5 <-sum_a_roti_sample_data_shannon_temp_5[-1,]

sum_a_roti_sample_data_richness_temp_5 <- capture.output(summary(sub_a_roti_sample_data_temp_5$Richness),file=NULL,append=FALSE)
sum_a_roti_sample_data_richness_temp_5 <- as.data.frame(sum_a_roti_sample_data_richness_temp_5)
sum_a_roti_sample_data_richness_temp_5$group <- "Rotifers"
sum_a_roti_sample_data_richness_temp_5$diversity <- "Richness"
sum_a_roti_sample_data_richness_temp_5$class <- "High temperature"
sum_a_roti_sample_data_richness_temp_5 <-sum_a_roti_sample_data_richness_temp_5[-1,]

sub_a_roti_sample_data_temp_6 <- a_roti_sample_data %>% filter(cluster_temp == "Very high")
sum_a_roti_sample_data_shannon_temp_6 <- capture.output(summary(sub_a_roti_sample_data_temp_6$Shannon),file=NULL,append=FALSE)
sum_a_roti_sample_data_shannon_temp_6 <- as.data.frame(sum_a_roti_sample_data_shannon_temp_6)
sum_a_roti_sample_data_shannon_temp_6$group <- "Rotifers"
sum_a_roti_sample_data_shannon_temp_6$diversity <- "Shannon"
sum_a_roti_sample_data_shannon_temp_6$class <- "Very high temperature"
sum_a_roti_sample_data_shannon_temp_6 <-sum_a_roti_sample_data_shannon_temp_6[-1,]

sum_a_roti_sample_data_richness_temp_6 <- capture.output(summary(sub_a_roti_sample_data_temp_6$Richness),file=NULL,append=FALSE)
sum_a_roti_sample_data_richness_temp_6 <- as.data.frame(sum_a_roti_sample_data_richness_temp_6)
sum_a_roti_sample_data_richness_temp_6$group <- "Rotifers"
sum_a_roti_sample_data_richness_temp_6$diversity <- "Richness"
sum_a_roti_sample_data_richness_temp_6$class <- "Very high temperature"
sum_a_roti_sample_data_richness_temp_6 <-sum_a_roti_sample_data_richness_temp_6[-1,] 

sub_a_roti_sample_data_temp_range_1 <- a_roti_sample_data %>% filter(cluster_temp_range == "Very low")
sum_a_roti_sample_data_shannon_temp_range_1 <- capture.output(summary(sub_a_roti_sample_data_temp_range_1$Shannon),file=NULL,append=FALSE)
sum_a_roti_sample_data_shannon_temp_range_1 <- as.data.frame(sum_a_roti_sample_data_shannon_temp_range_1)
sum_a_roti_sample_data_shannon_temp_range_1$group <- "Rotifers"
sum_a_roti_sample_data_shannon_temp_range_1$diversity <- "Shannon"
sum_a_roti_sample_data_shannon_temp_range_1$class <- "Very low temperature range"
sum_a_roti_sample_data_shannon_temp_range_1 <-sum_a_roti_sample_data_shannon_temp_range_1[-1,]

sum_a_roti_sample_data_richness_temp_range_1 <- capture.output(summary(sub_a_roti_sample_data_temp_range_1$Richness),file=NULL,append=FALSE)
sum_a_roti_sample_data_richness_temp_range_1 <- as.data.frame(sum_a_roti_sample_data_richness_temp_range_1)
sum_a_roti_sample_data_richness_temp_range_1$group <- "Rotifers"
sum_a_roti_sample_data_richness_temp_range_1$diversity <- "Richness"
sum_a_roti_sample_data_richness_temp_range_1$class <- "Very low temperature range"
sum_a_roti_sample_data_richness_temp_range_1 <-sum_a_roti_sample_data_richness_temp_range_1[-1,]

sub_a_roti_sample_data_temp_range_2 <- a_roti_sample_data %>% filter(cluster_temp_range == "Low")
sum_a_roti_sample_data_shannon_temp_range_2 <- capture.output(summary(sub_a_roti_sample_data_temp_range_2$Shannon),file=NULL,append=FALSE)
sum_a_roti_sample_data_shannon_temp_range_2 <- as.data.frame(sum_a_roti_sample_data_shannon_temp_range_2)
sum_a_roti_sample_data_shannon_temp_range_2$group <- "Rotifers"
sum_a_roti_sample_data_shannon_temp_range_2$diversity <- "Shannon"
sum_a_roti_sample_data_shannon_temp_range_2$class <- "Low temperature range"
sum_a_roti_sample_data_shannon_temp_range_2 <-sum_a_roti_sample_data_shannon_temp_range_2[-1,]

sum_a_roti_sample_data_richness_temp_range_2 <- capture.output(summary(sub_a_roti_sample_data_temp_range_2$Richness),file=NULL,append=FALSE)
sum_a_roti_sample_data_richness_temp_range_2 <- as.data.frame(sum_a_roti_sample_data_richness_temp_range_2)
sum_a_roti_sample_data_richness_temp_range_2$group <- "Rotifers"
sum_a_roti_sample_data_richness_temp_range_2$diversity <- "Richness"
sum_a_roti_sample_data_richness_temp_range_2$class <- "Low temperature range"
sum_a_roti_sample_data_richness_temp_range_2 <-sum_a_roti_sample_data_richness_temp_range_2[-1,]

sub_a_roti_sample_data_temp_range_3 <- a_roti_sample_data %>% filter(cluster_temp_range == "Medium low")
sum_a_roti_sample_data_shannon_temp_range_3 <- capture.output(summary(sub_a_roti_sample_data_temp_range_3$Shannon),file=NULL,append=FALSE)
sum_a_roti_sample_data_shannon_temp_range_3 <- as.data.frame(sum_a_roti_sample_data_shannon_temp_range_3)
sum_a_roti_sample_data_shannon_temp_range_3$group <- "Rotifers"
sum_a_roti_sample_data_shannon_temp_range_3$diversity <- "Shannon"
sum_a_roti_sample_data_shannon_temp_range_3$class <- "Medium low temperature range"
sum_a_roti_sample_data_shannon_temp_range_3 <-sum_a_roti_sample_data_shannon_temp_range_3[-1,]

sum_a_roti_sample_data_richness_temp_range_3 <- capture.output(summary(sub_a_roti_sample_data_temp_range_3$Richness),file=NULL,append=FALSE)
sum_a_roti_sample_data_richness_temp_range_3 <- as.data.frame(sum_a_roti_sample_data_richness_temp_range_3)
sum_a_roti_sample_data_richness_temp_range_3$group <- "Rotifers"
sum_a_roti_sample_data_richness_temp_range_3$diversity <- "Richness"
sum_a_roti_sample_data_richness_temp_range_3$class <- "Medium low temperature range"
sum_a_roti_sample_data_richness_temp_range_3 <-sum_a_roti_sample_data_richness_temp_range_3[-1,]

sub_a_roti_sample_data_temp_range_4 <- a_roti_sample_data %>% filter(cluster_temp_range == "Medium high")
sum_a_roti_sample_data_shannon_temp_range_4 <- capture.output(summary(sub_a_roti_sample_data_temp_range_4$Shannon),file=NULL,append=FALSE)
sum_a_roti_sample_data_shannon_temp_range_4 <- as.data.frame(sum_a_roti_sample_data_shannon_temp_range_4)
sum_a_roti_sample_data_shannon_temp_range_4$group <- "Rotifers"
sum_a_roti_sample_data_shannon_temp_range_4$diversity <- "Shannon"
sum_a_roti_sample_data_shannon_temp_range_4$class <- "Medium high temperature range"
sum_a_roti_sample_data_shannon_temp_range_4 <-sum_a_roti_sample_data_shannon_temp_range_4[-1,]

sum_a_roti_sample_data_richness_temp_range_4 <- capture.output(summary(sub_a_roti_sample_data_temp_range_4$Richness),file=NULL,append=FALSE)
sum_a_roti_sample_data_richness_temp_range_4 <- as.data.frame(sum_a_roti_sample_data_richness_temp_range_4)
sum_a_roti_sample_data_richness_temp_range_4$group <- "Rotifers"
sum_a_roti_sample_data_richness_temp_range_4$diversity <- "Richness"
sum_a_roti_sample_data_richness_temp_range_4$class <- "Medium high temperature range"
sum_a_roti_sample_data_richness_temp_range_4 <-sum_a_roti_sample_data_richness_temp_range_4[-1,]

sub_a_roti_sample_data_temp_range_5 <- a_roti_sample_data %>% filter(cluster_temp_range == "High")
sum_a_roti_sample_data_shannon_temp_range_5 <- capture.output(summary(sub_a_roti_sample_data_temp_range_5$Shannon),file=NULL,append=FALSE)
sum_a_roti_sample_data_shannon_temp_range_5 <- as.data.frame(sum_a_roti_sample_data_shannon_temp_range_5)
sum_a_roti_sample_data_shannon_temp_range_5$group <- "Rotifers"
sum_a_roti_sample_data_shannon_temp_range_5$diversity <- "Shannon"
sum_a_roti_sample_data_shannon_temp_range_5$class <- "High temperature range"
sum_a_roti_sample_data_shannon_temp_range_5 <-sum_a_roti_sample_data_shannon_temp_range_5[-1,]

sum_a_roti_sample_data_richness_temp_range_5 <- capture.output(summary(sub_a_roti_sample_data_temp_range_5$Richness),file=NULL,append=FALSE)
sum_a_roti_sample_data_richness_temp_range_5 <- as.data.frame(sum_a_roti_sample_data_richness_temp_range_5)
sum_a_roti_sample_data_richness_temp_range_5$group <- "Rotifers"
sum_a_roti_sample_data_richness_temp_range_5$diversity <- "Richness"
sum_a_roti_sample_data_richness_temp_range_5$class <- "High temperature range"
sum_a_roti_sample_data_richness_temp_range_5 <-sum_a_roti_sample_data_richness_temp_range_5[-1,]

sub_a_roti_sample_data_temp_range_6 <- a_roti_sample_data %>% filter(cluster_temp_range == "Very high")
sum_a_roti_sample_data_shannon_temp_range_6 <- capture.output(summary(sub_a_roti_sample_data_temp_range_6$Shannon),file=NULL,append=FALSE)
sum_a_roti_sample_data_shannon_temp_range_6 <- as.data.frame(sum_a_roti_sample_data_shannon_temp_range_6)
sum_a_roti_sample_data_shannon_temp_range_6$group <- "Rotifers"
sum_a_roti_sample_data_shannon_temp_range_6$diversity <- "Shannon"
sum_a_roti_sample_data_shannon_temp_range_6$class <- "Very high temperature range"
sum_a_roti_sample_data_shannon_temp_range_6 <-sum_a_roti_sample_data_shannon_temp_range_6[-1,]

sum_a_roti_sample_data_richness_temp_range_6 <- capture.output(summary(sub_a_roti_sample_data_temp_range_6$Richness),file=NULL,append=FALSE)
sum_a_roti_sample_data_richness_temp_range_6 <- as.data.frame(sum_a_roti_sample_data_richness_temp_range_6)
sum_a_roti_sample_data_richness_temp_range_6$group <- "Rotifers"
sum_a_roti_sample_data_richness_temp_range_6$diversity <- "Richness"
sum_a_roti_sample_data_richness_temp_range_6$class <- "Very high temperature range"
sum_a_roti_sample_data_richness_temp_range_6 <-sum_a_roti_sample_data_richness_temp_range_6[-1,]

sub_a_roti_sample_data_prec_1 <- a_roti_sample_data %>% filter(cluster_prec == "Very low")
sum_a_roti_sample_data_shannon_prec_1 <- capture.output(summary(sub_a_roti_sample_data_prec_1$Shannon),file=NULL,append=FALSE)
sum_a_roti_sample_data_shannon_prec_1 <- as.data.frame(sum_a_roti_sample_data_shannon_prec_1)
sum_a_roti_sample_data_shannon_prec_1$group <- "Rotifers"
sum_a_roti_sample_data_shannon_prec_1$diversity <- "Shannon"
sum_a_roti_sample_data_shannon_prec_1$class <- "Very low precipitation "
sum_a_roti_sample_data_shannon_prec_1 <-sum_a_roti_sample_data_shannon_prec_1[-1,]

sum_a_roti_sample_data_richness_prec_1 <- capture.output(summary(sub_a_roti_sample_data_prec_1$Richness),file=NULL,append=FALSE)
sum_a_roti_sample_data_richness_prec_1 <- as.data.frame(sum_a_roti_sample_data_richness_prec_1)
sum_a_roti_sample_data_richness_prec_1$group <- "Rotifers"
sum_a_roti_sample_data_richness_prec_1$diversity <- "Richness"
sum_a_roti_sample_data_richness_prec_1$class <- "Very low precipitation "
sum_a_roti_sample_data_richness_prec_1 <-sum_a_roti_sample_data_richness_prec_1[-1,]

sub_a_roti_sample_data_prec_2 <- a_roti_sample_data %>% filter(cluster_prec == "Low")
sum_a_roti_sample_data_shannon_prec_2 <- capture.output(summary(sub_a_roti_sample_data_prec_2$Shannon),file=NULL,append=FALSE)
sum_a_roti_sample_data_shannon_prec_2 <- as.data.frame(sum_a_roti_sample_data_shannon_prec_2)
sum_a_roti_sample_data_shannon_prec_2$group <- "Rotifers"
sum_a_roti_sample_data_shannon_prec_2$diversity <- "Shannon"
sum_a_roti_sample_data_shannon_prec_2$class <- "Low precipitation"
sum_a_roti_sample_data_shannon_prec_2 <-sum_a_roti_sample_data_shannon_prec_2[-1,]

sum_a_roti_sample_data_richness_prec_2 <- capture.output(summary(sub_a_roti_sample_data_prec_2$Richness),file=NULL,append=FALSE)
sum_a_roti_sample_data_richness_prec_2 <- as.data.frame(sum_a_roti_sample_data_richness_prec_2)
sum_a_roti_sample_data_richness_prec_2$group <- "Rotifers"
sum_a_roti_sample_data_richness_prec_2$diversity <- "Richness"
sum_a_roti_sample_data_richness_prec_2$class <- "Low precipitation"
sum_a_roti_sample_data_richness_prec_2 <-sum_a_roti_sample_data_richness_prec_2[-1,]

sub_a_roti_sample_data_prec_3 <- a_roti_sample_data %>% filter(cluster_prec == "Medium low")
sum_a_roti_sample_data_shannon_prec_3 <- capture.output(summary(sub_a_roti_sample_data_prec_3$Shannon),file=NULL,append=FALSE)
sum_a_roti_sample_data_shannon_prec_3 <- as.data.frame(sum_a_roti_sample_data_shannon_prec_3)
sum_a_roti_sample_data_shannon_prec_3$group <- "Rotifers"
sum_a_roti_sample_data_shannon_prec_3$diversity <- "Shannon"
sum_a_roti_sample_data_shannon_prec_3$class <- "Medium low precipitation"
sum_a_roti_sample_data_shannon_prec_3 <-sum_a_roti_sample_data_shannon_prec_3[-1,]

sum_a_roti_sample_data_richness_prec_3 <- capture.output(summary(sub_a_roti_sample_data_prec_3$Richness),file=NULL,append=FALSE)
sum_a_roti_sample_data_richness_prec_3 <- as.data.frame(sum_a_roti_sample_data_richness_prec_3)
sum_a_roti_sample_data_richness_prec_3$group <- "Rotifers"
sum_a_roti_sample_data_richness_prec_3$diversity <- "Richness"
sum_a_roti_sample_data_richness_prec_3$class <- "Medium low precipitation"
sum_a_roti_sample_data_richness_prec_3 <-sum_a_roti_sample_data_richness_prec_3[-1,]

sub_a_roti_sample_data_prec_4 <- a_roti_sample_data %>% filter(cluster_prec == "Medium high")
sum_a_roti_sample_data_shannon_prec_4 <- capture.output(summary(sub_a_roti_sample_data_prec_4$Shannon),file=NULL,append=FALSE)
sum_a_roti_sample_data_shannon_prec_4 <- as.data.frame(sum_a_roti_sample_data_shannon_prec_4)
sum_a_roti_sample_data_shannon_prec_4$group <- "Rotifers"
sum_a_roti_sample_data_shannon_prec_4$diversity <- "Shannon"
sum_a_roti_sample_data_shannon_prec_4$class <- "Medium high precipitation"
sum_a_roti_sample_data_shannon_prec_4 <-sum_a_roti_sample_data_shannon_prec_4[-1,]

sum_a_roti_sample_data_richness_prec_4 <- capture.output(summary(sub_a_roti_sample_data_prec_4$Richness),file=NULL,append=FALSE)
sum_a_roti_sample_data_richness_prec_4 <- as.data.frame(sum_a_roti_sample_data_richness_prec_4)
sum_a_roti_sample_data_richness_prec_4$group <- "Rotifers"
sum_a_roti_sample_data_richness_prec_4$diversity <- "Richness"
sum_a_roti_sample_data_richness_prec_4$class <- "Medium high precipitation"
sum_a_roti_sample_data_richness_prec_4 <-sum_a_roti_sample_data_richness_prec_4[-1,]

sub_a_roti_sample_data_prec_5 <- a_roti_sample_data %>% filter(cluster_prec == "High")
sum_a_roti_sample_data_shannon_prec_5 <- capture.output(summary(sub_a_roti_sample_data_prec_5$Shannon),file=NULL,append=FALSE)
sum_a_roti_sample_data_shannon_prec_5 <- as.data.frame(sum_a_roti_sample_data_shannon_prec_5)
sum_a_roti_sample_data_shannon_prec_5$group <- "Rotifers"
sum_a_roti_sample_data_shannon_prec_5$diversity <- "Shannon"
sum_a_roti_sample_data_shannon_prec_5$class <- "High precipitation"
sum_a_roti_sample_data_shannon_prec_5 <-sum_a_roti_sample_data_shannon_prec_5[-1,]

sum_a_roti_sample_data_richness_prec_5 <- capture.output(summary(sub_a_roti_sample_data_prec_5$Richness),file=NULL,append=FALSE)
sum_a_roti_sample_data_richness_prec_5 <- as.data.frame(sum_a_roti_sample_data_richness_prec_5)
sum_a_roti_sample_data_richness_prec_5$group <- "Rotifers"
sum_a_roti_sample_data_richness_prec_5$diversity <- "Richness"
sum_a_roti_sample_data_richness_prec_5$class <- "High precipitation"
sum_a_roti_sample_data_richness_prec_5 <-sum_a_roti_sample_data_richness_prec_5[-1,]

sub_a_roti_sample_data_prec_6 <- a_roti_sample_data %>% filter(cluster_prec == "Very high")
sum_a_roti_sample_data_shannon_prec_6 <- capture.output(summary(sub_a_roti_sample_data_prec_6$Shannon),file=NULL,append=FALSE)
sum_a_roti_sample_data_shannon_prec_6 <- as.data.frame(sum_a_roti_sample_data_shannon_prec_6)
sum_a_roti_sample_data_shannon_prec_6$group <- "Rotifers"
sum_a_roti_sample_data_shannon_prec_6$diversity <- "Shannon"
sum_a_roti_sample_data_shannon_prec_6$class <- "Very high precipitation"
sum_a_roti_sample_data_shannon_prec_6 <-sum_a_roti_sample_data_shannon_prec_6[-1,]

sum_a_roti_sample_data_richness_prec_6 <- capture.output(summary(sub_a_roti_sample_data_prec_6$Richness),file=NULL,append=FALSE)
sum_a_roti_sample_data_richness_prec_6 <- as.data.frame(sum_a_roti_sample_data_richness_prec_6)
sum_a_roti_sample_data_richness_prec_6$group <- "Rotifers"
sum_a_roti_sample_data_richness_prec_6$diversity <- "Richness"
sum_a_roti_sample_data_richness_prec_6$class <- "Very high precipitation"
sum_a_roti_sample_data_richness_prec_6 <-sum_a_roti_sample_data_richness_prec_6[-1,] 

sub_a_tardi_sample_data_CL_annual <- a_tardi_sample_data %>% filter(LC1_2018 == "CL_annual")
sum_a_tardi_richness_CL_annual <- capture.output(summary(sub_a_tardi_sample_data_CL_annual$Richness),file=NULL,append=FALSE)
sum_a_tardi_richness_CL_annual <- as.data.frame(sum_a_tardi_richness_CL_annual)
sum_a_tardi_richness_CL_annual$group <- "Tardigrades"
sum_a_tardi_richness_CL_annual$diversity <- "Richness"
sum_a_tardi_richness_CL_annual$class <- "Land Cover CL_annual"
sum_a_tardi_richness_CL_annual <-sum_a_tardi_richness_CL_annual[-1,]

sum_a_tardi_shannon_CL_annual <- capture.output(summary(sub_a_tardi_sample_data_CL_annual$Shannon),file=NULL,append=FALSE)
sum_a_tardi_shannon_CL_annual <- as.data.frame(sum_a_tardi_shannon_CL_annual)
sum_a_tardi_shannon_CL_annual$group <- "Tardigrades"
sum_a_tardi_shannon_CL_annual$diversity <- "Shannon"
sum_a_tardi_shannon_CL_annual$class <- "Land Cover CL_annual"
sum_a_tardi_shannon_CL_annual <-sum_a_tardi_shannon_CL_annual[-1,]

sub_a_tardi_sample_data_CL_permanent <- a_tardi_sample_data %>% filter(LC1_2018 == "CL_permanent")
sum_a_tardi_richness_CL_permanent <- capture.output(summary(sub_a_tardi_sample_data_CL_permanent$Richness),file=NULL,append=FALSE)
sum_a_tardi_richness_CL_permanent <- as.data.frame(sum_a_tardi_richness_CL_permanent)
sum_a_tardi_richness_CL_permanent$group <- "Tardigrades"
sum_a_tardi_richness_CL_permanent$diversity <- "Richness"
sum_a_tardi_richness_CL_permanent$class <- "Land Cover CL_permanent"
sum_a_tardi_richness_CL_permanent <-sum_a_tardi_richness_CL_permanent[-1,]

sum_a_tardi_shannon_CL_permanent <- capture.output(summary(sub_a_tardi_sample_data_CL_permanent$Shannon),file=NULL,append=FALSE)
sum_a_tardi_shannon_CL_permanent <- as.data.frame(sum_a_tardi_shannon_CL_permanent)
sum_a_tardi_shannon_CL_permanent$group <- "Tardigrades"
sum_a_tardi_shannon_CL_permanent$diversity <- "Shannon"
sum_a_tardi_shannon_CL_permanent$class <- "Land Cover CL_permanent"
sum_a_tardi_shannon_CL_permanent <-sum_a_tardi_shannon_CL_permanent[-1,]

sub_a_tardi_sample_data_GL_unmanaged <- a_tardi_sample_data %>% filter(LC1_2018 == "GL_unmanaged")
sum_a_tardi_richness_GL_unmanaged <- capture.output(summary(sub_a_tardi_sample_data_GL_unmanaged$Richness),file=NULL,append=FALSE)
sum_a_tardi_richness_GL_unmanaged <- as.data.frame(sum_a_tardi_richness_GL_unmanaged)
sum_a_tardi_richness_GL_unmanaged$group <- "Tardigrades"
sum_a_tardi_richness_GL_unmanaged$diversity <- "Richness"
sum_a_tardi_richness_GL_unmanaged$class <- "Land Cover GL_unmanaged"
sum_a_tardi_richness_GL_unmanaged <-sum_a_tardi_richness_GL_unmanaged[-1,]

sum_a_tardi_shannon_GL_unmanaged <- capture.output(summary(sub_a_tardi_sample_data_GL_unmanaged$Shannon),file=NULL,append=FALSE)
sum_a_tardi_shannon_GL_unmanaged <- as.data.frame(sum_a_tardi_shannon_GL_unmanaged)
sum_a_tardi_shannon_GL_unmanaged$group <- "Tardigrades"
sum_a_tardi_shannon_GL_unmanaged$diversity <- "Shannon"
sum_a_tardi_shannon_GL_unmanaged$class <- "Land Cover GL_unmanaged"
sum_a_tardi_shannon_GL_unmanaged <-sum_a_tardi_shannon_GL_unmanaged[-1,]

sub_a_tardi_sample_data_GL_managed <- a_tardi_sample_data %>% filter(LC1_2018 == "GL_managed")
sum_a_tardi_richness_GL_managed <- capture.output(summary(sub_a_tardi_sample_data_GL_managed$Richness),file=NULL,append=FALSE)
sum_a_tardi_richness_GL_managed <- as.data.frame(sum_a_tardi_richness_GL_managed)
sum_a_tardi_richness_GL_managed$group <- "Tardigrades"
sum_a_tardi_richness_GL_managed$diversity <- "Richness"
sum_a_tardi_richness_GL_managed$class <- "Land Cover GL_managed"
sum_a_tardi_richness_GL_managed <-sum_a_tardi_richness_GL_managed[-1,]

sum_a_tardi_shannon_GL_managed <- capture.output(summary(sub_a_tardi_sample_data_GL_managed$Shannon),file=NULL,append=FALSE)
sum_a_tardi_shannon_GL_managed <- as.data.frame(sum_a_tardi_shannon_GL_managed)
sum_a_tardi_shannon_GL_managed$group <- "Tardigrades"
sum_a_tardi_shannon_GL_managed$diversity <- "Shannon"
sum_a_tardi_shannon_GL_managed$class <- "Land Cover GL_managed"
sum_a_tardi_shannon_GL_managed <-sum_a_tardi_shannon_GL_managed[-1,]

sub_a_tardi_sample_data_WL_broadleaved <- a_tardi_sample_data %>% filter(LC1_2018 == "WL_broadleaved")
sum_a_tardi_richness_WL_broadleaved <- capture.output(summary(sub_a_tardi_sample_data_WL_broadleaved$Richness),file=NULL,append=FALSE)
sum_a_tardi_richness_WL_broadleaved <- as.data.frame(sum_a_tardi_richness_WL_broadleaved)
sum_a_tardi_richness_WL_broadleaved$group <- "Tardigrades"
sum_a_tardi_richness_WL_broadleaved$diversity <- "Richness"
sum_a_tardi_richness_WL_broadleaved$class <- "Land Cover WL_broadleaved"
sum_a_tardi_richness_WL_broadleaved <-sum_a_tardi_richness_WL_broadleaved[-1,]

sum_a_tardi_shannon_WL_broadleaved <- capture.output(summary(sub_a_tardi_sample_data_WL_broadleaved$Shannon),file=NULL,append=FALSE)
sum_a_tardi_shannon_WL_broadleaved <- as.data.frame(sum_a_tardi_shannon_WL_broadleaved)
sum_a_tardi_shannon_WL_broadleaved$group <- "Tardigrades"
sum_a_tardi_shannon_WL_broadleaved$diversity <- "Shannon"
sum_a_tardi_shannon_WL_broadleaved$class <- "Land Cover WL_broadleaved"
sum_a_tardi_shannon_WL_broadleaved <-sum_a_tardi_shannon_WL_broadleaved[-1,]

sub_a_tardi_sample_data_WL_coniferous <- a_tardi_sample_data %>% filter(LC1_2018 == "WL_coniferous")
sum_a_tardi_richness_WL_coniferous <- capture.output(summary(sub_a_tardi_sample_data_WL_coniferous$Richness),file=NULL,append=FALSE)
sum_a_tardi_richness_WL_coniferous <- as.data.frame(sum_a_tardi_richness_WL_coniferous)
sum_a_tardi_richness_WL_coniferous$group <- "Tardigrades"
sum_a_tardi_richness_WL_coniferous$diversity <- "Richness"
sum_a_tardi_richness_WL_coniferous$class <- "Land Cover WL_coniferous"
sum_a_tardi_richness_WL_coniferous <-sum_a_tardi_richness_WL_coniferous[-1,]

sum_a_tardi_shannon_WL_coniferous <- capture.output(summary(sub_a_tardi_sample_data_WL_coniferous$Shannon),file=NULL,append=FALSE)
sum_a_tardi_shannon_WL_coniferous <- as.data.frame(sum_a_tardi_shannon_WL_coniferous)
sum_a_tardi_shannon_WL_coniferous$group <- "Tardigrades"
sum_a_tardi_shannon_WL_coniferous$diversity <- "Shannon"
sum_a_tardi_shannon_WL_coniferous$class <- "Land Cover WL_coniferous"
sum_a_tardi_shannon_WL_coniferous <-sum_a_tardi_shannon_WL_coniferous[-1,]

sub_a_tardi_sample_data_temp_1 <- a_tardi_sample_data %>% filter(cluster_temp == "Very low")
sum_a_tardi_sample_data_shannon_temp_1 <- capture.output(summary(sub_a_tardi_sample_data_temp_1$Shannon),file=NULL,append=FALSE)
sum_a_tardi_sample_data_shannon_temp_1 <- as.data.frame(sum_a_tardi_sample_data_shannon_temp_1)
sum_a_tardi_sample_data_shannon_temp_1$group <- "Tardigrades"
sum_a_tardi_sample_data_shannon_temp_1$diversity <- "Shannon"
sum_a_tardi_sample_data_shannon_temp_1$class <- "Very low temperature"
sum_a_tardi_sample_data_shannon_temp_1 <-sum_a_tardi_sample_data_shannon_temp_1[-1,]

sum_a_tardi_sample_data_richness_temp_1 <- capture.output(summary(sub_a_tardi_sample_data_temp_1$Richness),file=NULL,append=FALSE)
sum_a_tardi_sample_data_richness_temp_1 <- as.data.frame(sum_a_tardi_sample_data_richness_temp_1)
sum_a_tardi_sample_data_richness_temp_1$group <- "Tardigrades"
sum_a_tardi_sample_data_richness_temp_1$diversity <- "Richness"
sum_a_tardi_sample_data_richness_temp_1$class <- "Very low temperature"
sum_a_tardi_sample_data_richness_temp_1 <-sum_a_tardi_sample_data_richness_temp_1[-1,]

sub_a_tardi_sample_data_temp_2 <- a_tardi_sample_data %>% filter(cluster_temp == "Low")
sum_a_tardi_sample_data_shannon_temp_2 <- capture.output(summary(sub_a_tardi_sample_data_temp_2$Shannon),file=NULL,append=FALSE)
sum_a_tardi_sample_data_shannon_temp_2 <- as.data.frame(sum_a_tardi_sample_data_shannon_temp_2)
sum_a_tardi_sample_data_shannon_temp_2$group <- "Tardigrades"
sum_a_tardi_sample_data_shannon_temp_2$diversity <- "Shannon"
sum_a_tardi_sample_data_shannon_temp_2$class <- "Low temperature"
sum_a_tardi_sample_data_shannon_temp_2 <-sum_a_tardi_sample_data_shannon_temp_2[-1,]

sum_a_tardi_sample_data_richness_temp_2 <- capture.output(summary(sub_a_tardi_sample_data_temp_2$Richness),file=NULL,append=FALSE)
sum_a_tardi_sample_data_richness_temp_2 <- as.data.frame(sum_a_tardi_sample_data_richness_temp_2)
sum_a_tardi_sample_data_richness_temp_2$group <- "Tardigrades"
sum_a_tardi_sample_data_richness_temp_2$diversity <- "Richness"
sum_a_tardi_sample_data_richness_temp_2$class <- "Low temperature"
sum_a_tardi_sample_data_richness_temp_2 <-sum_a_tardi_sample_data_richness_temp_2[-1,]

sub_a_tardi_sample_data_temp_3 <- a_tardi_sample_data %>% filter(cluster_temp == "Medium low")
sum_a_tardi_sample_data_shannon_temp_3 <- capture.output(summary(sub_a_tardi_sample_data_temp_3$Shannon),file=NULL,append=FALSE)
sum_a_tardi_sample_data_shannon_temp_3 <- as.data.frame(sum_a_tardi_sample_data_shannon_temp_3)
sum_a_tardi_sample_data_shannon_temp_3$group <- "Tardigrades"
sum_a_tardi_sample_data_shannon_temp_3$diversity <- "Shannon"
sum_a_tardi_sample_data_shannon_temp_3$class <- "Medium low temperature"
sum_a_tardi_sample_data_shannon_temp_3 <-sum_a_tardi_sample_data_shannon_temp_3[-1,]

sum_a_tardi_sample_data_richness_temp_3 <- capture.output(summary(sub_a_tardi_sample_data_temp_3$Richness),file=NULL,append=FALSE)
sum_a_tardi_sample_data_richness_temp_3 <- as.data.frame(sum_a_tardi_sample_data_richness_temp_3)
sum_a_tardi_sample_data_richness_temp_3$group <- "Tardigrades"
sum_a_tardi_sample_data_richness_temp_3$diversity <- "Richness"
sum_a_tardi_sample_data_richness_temp_3$class <- "Medium low temperature"
sum_a_tardi_sample_data_richness_temp_3 <-sum_a_tardi_sample_data_richness_temp_3[-1,]

sub_a_tardi_sample_data_temp_4 <- a_tardi_sample_data %>% filter(cluster_temp == "Medium high")
sum_a_tardi_sample_data_shannon_temp_4 <- capture.output(summary(sub_a_tardi_sample_data_temp_4$Shannon),file=NULL,append=FALSE)
sum_a_tardi_sample_data_shannon_temp_4 <- as.data.frame(sum_a_tardi_sample_data_shannon_temp_4)
sum_a_tardi_sample_data_shannon_temp_4$group <- "Tardigrades"
sum_a_tardi_sample_data_shannon_temp_4$diversity <- "Shannon"
sum_a_tardi_sample_data_shannon_temp_4$class <- "Medium high temperature"
sum_a_tardi_sample_data_shannon_temp_4 <-sum_a_tardi_sample_data_shannon_temp_4[-1,]

sum_a_tardi_sample_data_richness_temp_4 <- capture.output(summary(sub_a_tardi_sample_data_temp_4$Richness),file=NULL,append=FALSE)
sum_a_tardi_sample_data_richness_temp_4 <- as.data.frame(sum_a_tardi_sample_data_richness_temp_4)
sum_a_tardi_sample_data_richness_temp_4$group <- "Tardigrades"
sum_a_tardi_sample_data_richness_temp_4$diversity <- "Richness"
sum_a_tardi_sample_data_richness_temp_4$class <- "Medium high temperature"
sum_a_tardi_sample_data_richness_temp_4 <-sum_a_tardi_sample_data_richness_temp_4[-1,]

sub_a_tardi_sample_data_temp_5 <- a_tardi_sample_data %>% filter(cluster_temp == "High temperature")
sum_a_tardi_sample_data_shannon_temp_5 <- capture.output(summary(sub_a_tardi_sample_data_temp_5$Shannon),file=NULL,append=FALSE)
sum_a_tardi_sample_data_shannon_temp_5 <- as.data.frame(sum_a_tardi_sample_data_shannon_temp_5)
sum_a_tardi_sample_data_shannon_temp_5$group <- "Tardigrades"
sum_a_tardi_sample_data_shannon_temp_5$diversity <- "Shannon"
sum_a_tardi_sample_data_shannon_temp_5$class <- "High temperature"
sum_a_tardi_sample_data_shannon_temp_5 <-sum_a_tardi_sample_data_shannon_temp_5[-1,]

sum_a_tardi_sample_data_richness_temp_5 <- capture.output(summary(sub_a_tardi_sample_data_temp_5$Richness),file=NULL,append=FALSE)
sum_a_tardi_sample_data_richness_temp_5 <- as.data.frame(sum_a_tardi_sample_data_richness_temp_5)
sum_a_tardi_sample_data_richness_temp_5$group <- "Tardigrades"
sum_a_tardi_sample_data_richness_temp_5$diversity <- "Richness"
sum_a_tardi_sample_data_richness_temp_5$class <- "High temperature"
sum_a_tardi_sample_data_richness_temp_5 <-sum_a_tardi_sample_data_richness_temp_5[-1,]

sub_a_tardi_sample_data_temp_6 <- a_tardi_sample_data %>% filter(cluster_temp == "Very high")
sum_a_tardi_sample_data_shannon_temp_6 <- capture.output(summary(sub_a_tardi_sample_data_temp_6$Shannon),file=NULL,append=FALSE)
sum_a_tardi_sample_data_shannon_temp_6 <- as.data.frame(sum_a_tardi_sample_data_shannon_temp_6)
sum_a_tardi_sample_data_shannon_temp_6$group <- "Tardigrades"
sum_a_tardi_sample_data_shannon_temp_6$diversity <- "Shannon"
sum_a_tardi_sample_data_shannon_temp_6$class <- "Very high temperature"
sum_a_tardi_sample_data_shannon_temp_6 <-sum_a_tardi_sample_data_shannon_temp_6[-1,]

sum_a_tardi_sample_data_richness_temp_6 <- capture.output(summary(sub_a_tardi_sample_data_temp_6$Richness),file=NULL,append=FALSE)
sum_a_tardi_sample_data_richness_temp_6 <- as.data.frame(sum_a_tardi_sample_data_richness_temp_6)
sum_a_tardi_sample_data_richness_temp_6$group <- "Tardigrades"
sum_a_tardi_sample_data_richness_temp_6$diversity <- "Richness"
sum_a_tardi_sample_data_richness_temp_6$class <- "Very high temperature"
sum_a_tardi_sample_data_richness_temp_6 <-sum_a_tardi_sample_data_richness_temp_6[-1,] 

sub_a_tardi_sample_data_temp_range_1 <- a_tardi_sample_data %>% filter(cluster_temp_range == "Very low")
sum_a_tardi_sample_data_shannon_temp_range_1 <- capture.output(summary(sub_a_tardi_sample_data_temp_range_1$Shannon),file=NULL,append=FALSE)
sum_a_tardi_sample_data_shannon_temp_range_1 <- as.data.frame(sum_a_tardi_sample_data_shannon_temp_range_1)
sum_a_tardi_sample_data_shannon_temp_range_1$group <- "Tardigrades"
sum_a_tardi_sample_data_shannon_temp_range_1$diversity <- "Shannon"
sum_a_tardi_sample_data_shannon_temp_range_1$class <- "Very low temperature range"
sum_a_tardi_sample_data_shannon_temp_range_1 <-sum_a_tardi_sample_data_shannon_temp_range_1[-1,]

sum_a_tardi_sample_data_richness_temp_range_1 <- capture.output(summary(sub_a_tardi_sample_data_temp_range_1$Richness),file=NULL,append=FALSE)
sum_a_tardi_sample_data_richness_temp_range_1 <- as.data.frame(sum_a_tardi_sample_data_richness_temp_range_1)
sum_a_tardi_sample_data_richness_temp_range_1$group <- "Tardigrades"
sum_a_tardi_sample_data_richness_temp_range_1$diversity <- "Richness"
sum_a_tardi_sample_data_richness_temp_range_1$class <- "Very low temperature range"
sum_a_tardi_sample_data_richness_temp_range_1 <-sum_a_tardi_sample_data_richness_temp_range_1[-1,]

sub_a_tardi_sample_data_temp_range_2 <- a_tardi_sample_data %>% filter(cluster_temp_range == "Low")
sum_a_tardi_sample_data_shannon_temp_range_2 <- capture.output(summary(sub_a_tardi_sample_data_temp_range_2$Shannon),file=NULL,append=FALSE)
sum_a_tardi_sample_data_shannon_temp_range_2 <- as.data.frame(sum_a_tardi_sample_data_shannon_temp_range_2)
sum_a_tardi_sample_data_shannon_temp_range_2$group <- "Tardigrades"
sum_a_tardi_sample_data_shannon_temp_range_2$diversity <- "Shannon"
sum_a_tardi_sample_data_shannon_temp_range_2$class <- "Low temperature range"
sum_a_tardi_sample_data_shannon_temp_range_2 <-sum_a_tardi_sample_data_shannon_temp_range_2[-1,]

sum_a_tardi_sample_data_richness_temp_range_2 <- capture.output(summary(sub_a_tardi_sample_data_temp_range_2$Richness),file=NULL,append=FALSE)
sum_a_tardi_sample_data_richness_temp_range_2 <- as.data.frame(sum_a_tardi_sample_data_richness_temp_range_2)
sum_a_tardi_sample_data_richness_temp_range_2$group <- "Tardigrades"
sum_a_tardi_sample_data_richness_temp_range_2$diversity <- "Richness"
sum_a_tardi_sample_data_richness_temp_range_2$class <- "Low temperature range"
sum_a_tardi_sample_data_richness_temp_range_2 <-sum_a_tardi_sample_data_richness_temp_range_2[-1,]

sub_a_tardi_sample_data_temp_range_3 <- a_tardi_sample_data %>% filter(cluster_temp_range == "Medium low")
sum_a_tardi_sample_data_shannon_temp_range_3 <- capture.output(summary(sub_a_tardi_sample_data_temp_range_3$Shannon),file=NULL,append=FALSE)
sum_a_tardi_sample_data_shannon_temp_range_3 <- as.data.frame(sum_a_tardi_sample_data_shannon_temp_range_3)
sum_a_tardi_sample_data_shannon_temp_range_3$group <- "Tardigrades"
sum_a_tardi_sample_data_shannon_temp_range_3$diversity <- "Shannon"
sum_a_tardi_sample_data_shannon_temp_range_3$class <- "Medium low temperature range"
sum_a_tardi_sample_data_shannon_temp_range_3 <-sum_a_tardi_sample_data_shannon_temp_range_3[-1,]

sum_a_tardi_sample_data_richness_temp_range_3 <- capture.output(summary(sub_a_tardi_sample_data_temp_range_3$Richness),file=NULL,append=FALSE)
sum_a_tardi_sample_data_richness_temp_range_3 <- as.data.frame(sum_a_tardi_sample_data_richness_temp_range_3)
sum_a_tardi_sample_data_richness_temp_range_3$group <- "Tardigrades"
sum_a_tardi_sample_data_richness_temp_range_3$diversity <- "Richness"
sum_a_tardi_sample_data_richness_temp_range_3$class <- "Medium low temperature range"
sum_a_tardi_sample_data_richness_temp_range_3 <-sum_a_tardi_sample_data_richness_temp_range_3[-1,]

sub_a_tardi_sample_data_temp_range_4 <- a_tardi_sample_data %>% filter(cluster_temp_range == "Medium high")
sum_a_tardi_sample_data_shannon_temp_range_4 <- capture.output(summary(sub_a_tardi_sample_data_temp_range_4$Shannon),file=NULL,append=FALSE)
sum_a_tardi_sample_data_shannon_temp_range_4 <- as.data.frame(sum_a_tardi_sample_data_shannon_temp_range_4)
sum_a_tardi_sample_data_shannon_temp_range_4$group <- "Tardigrades"
sum_a_tardi_sample_data_shannon_temp_range_4$diversity <- "Shannon"
sum_a_tardi_sample_data_shannon_temp_range_4$class <- "Medium high temperature range"
sum_a_tardi_sample_data_shannon_temp_range_4 <-sum_a_tardi_sample_data_shannon_temp_range_4[-1,]

sum_a_tardi_sample_data_richness_temp_range_4 <- capture.output(summary(sub_a_tardi_sample_data_temp_range_4$Richness),file=NULL,append=FALSE)
sum_a_tardi_sample_data_richness_temp_range_4 <- as.data.frame(sum_a_tardi_sample_data_richness_temp_range_4)
sum_a_tardi_sample_data_richness_temp_range_4$group <- "Tardigrades"
sum_a_tardi_sample_data_richness_temp_range_4$diversity <- "Richness"
sum_a_tardi_sample_data_richness_temp_range_4$class <- "Medium high temperature range"
sum_a_tardi_sample_data_richness_temp_range_4 <-sum_a_tardi_sample_data_richness_temp_range_4[-1,]

sub_a_tardi_sample_data_temp_range_5 <- a_tardi_sample_data %>% filter(cluster_temp_range == "High")
sum_a_tardi_sample_data_shannon_temp_range_5 <- capture.output(summary(sub_a_tardi_sample_data_temp_range_5$Shannon),file=NULL,append=FALSE)
sum_a_tardi_sample_data_shannon_temp_range_5 <- as.data.frame(sum_a_tardi_sample_data_shannon_temp_range_5)
sum_a_tardi_sample_data_shannon_temp_range_5$group <- "Tardigrades"
sum_a_tardi_sample_data_shannon_temp_range_5$diversity <- "Shannon"
sum_a_tardi_sample_data_shannon_temp_range_5$class <- "High temperature range"
sum_a_tardi_sample_data_shannon_temp_range_5 <-sum_a_tardi_sample_data_shannon_temp_range_5[-1,]

sum_a_tardi_sample_data_richness_temp_range_5 <- capture.output(summary(sub_a_tardi_sample_data_temp_range_5$Richness),file=NULL,append=FALSE)
sum_a_tardi_sample_data_richness_temp_range_5 <- as.data.frame(sum_a_tardi_sample_data_richness_temp_range_5)
sum_a_tardi_sample_data_richness_temp_range_5$group <- "Tardigrades"
sum_a_tardi_sample_data_richness_temp_range_5$diversity <- "Richness"
sum_a_tardi_sample_data_richness_temp_range_5$class <- "High temperature range"
sum_a_tardi_sample_data_richness_temp_range_5 <-sum_a_tardi_sample_data_richness_temp_range_5[-1,]

sub_a_tardi_sample_data_temp_range_6 <- a_tardi_sample_data %>% filter(cluster_temp_range == "Very high")
sum_a_tardi_sample_data_shannon_temp_range_6 <- capture.output(summary(sub_a_tardi_sample_data_temp_range_6$Shannon),file=NULL,append=FALSE)
sum_a_tardi_sample_data_shannon_temp_range_6 <- as.data.frame(sum_a_tardi_sample_data_shannon_temp_range_6)
sum_a_tardi_sample_data_shannon_temp_range_6$group <- "Tardigrades"
sum_a_tardi_sample_data_shannon_temp_range_6$diversity <- "Shannon"
sum_a_tardi_sample_data_shannon_temp_range_6$class <- "Very high temperature range"
sum_a_tardi_sample_data_shannon_temp_range_6 <-sum_a_tardi_sample_data_shannon_temp_range_6[-1,]

sum_a_tardi_sample_data_richness_temp_range_6 <- capture.output(summary(sub_a_tardi_sample_data_temp_range_6$Richness),file=NULL,append=FALSE)
sum_a_tardi_sample_data_richness_temp_range_6 <- as.data.frame(sum_a_tardi_sample_data_richness_temp_range_6)
sum_a_tardi_sample_data_richness_temp_range_6$group <- "Tardigrades"
sum_a_tardi_sample_data_richness_temp_range_6$diversity <- "Richness"
sum_a_tardi_sample_data_richness_temp_range_6$class <- "Very high temperature range"
sum_a_tardi_sample_data_richness_temp_range_6 <-sum_a_tardi_sample_data_richness_temp_range_6[-1,]

sub_a_tardi_sample_data_prec_1 <- a_tardi_sample_data %>% filter(cluster_prec == "Very low")
sum_a_tardi_sample_data_shannon_prec_1 <- capture.output(summary(sub_a_tardi_sample_data_prec_1$Shannon),file=NULL,append=FALSE)
sum_a_tardi_sample_data_shannon_prec_1 <- as.data.frame(sum_a_tardi_sample_data_shannon_prec_1)
sum_a_tardi_sample_data_shannon_prec_1$group <- "Tardigrades"
sum_a_tardi_sample_data_shannon_prec_1$diversity <- "Shannon"
sum_a_tardi_sample_data_shannon_prec_1$class <- "Very low precipitation "
sum_a_tardi_sample_data_shannon_prec_1 <-sum_a_tardi_sample_data_shannon_prec_1[-1,]

sum_a_tardi_sample_data_richness_prec_1 <- capture.output(summary(sub_a_tardi_sample_data_prec_1$Richness),file=NULL,append=FALSE)
sum_a_tardi_sample_data_richness_prec_1 <- as.data.frame(sum_a_tardi_sample_data_richness_prec_1)
sum_a_tardi_sample_data_richness_prec_1$group <- "Tardigrades"
sum_a_tardi_sample_data_richness_prec_1$diversity <- "Richness"
sum_a_tardi_sample_data_richness_prec_1$class <- "Very low precipitation "
sum_a_tardi_sample_data_richness_prec_1 <-sum_a_tardi_sample_data_richness_prec_1[-1,]

sub_a_tardi_sample_data_prec_2 <- a_tardi_sample_data %>% filter(cluster_prec == "Low")
sum_a_tardi_sample_data_shannon_prec_2 <- capture.output(summary(sub_a_tardi_sample_data_prec_2$Shannon),file=NULL,append=FALSE)
sum_a_tardi_sample_data_shannon_prec_2 <- as.data.frame(sum_a_tardi_sample_data_shannon_prec_2)
sum_a_tardi_sample_data_shannon_prec_2$group <- "Tardigrades"
sum_a_tardi_sample_data_shannon_prec_2$diversity <- "Shannon"
sum_a_tardi_sample_data_shannon_prec_2$class <- "Low precipitation"
sum_a_tardi_sample_data_shannon_prec_2 <-sum_a_tardi_sample_data_shannon_prec_2[-1,]

sum_a_tardi_sample_data_richness_prec_2 <- capture.output(summary(sub_a_tardi_sample_data_prec_2$Richness),file=NULL,append=FALSE)
sum_a_tardi_sample_data_richness_prec_2 <- as.data.frame(sum_a_tardi_sample_data_richness_prec_2)
sum_a_tardi_sample_data_richness_prec_2$group <- "Tardigrades"
sum_a_tardi_sample_data_richness_prec_2$diversity <- "Richness"
sum_a_tardi_sample_data_richness_prec_2$class <- "Low precipitation"
sum_a_tardi_sample_data_richness_prec_2 <-sum_a_tardi_sample_data_richness_prec_2[-1,]

sub_a_tardi_sample_data_prec_3 <- a_tardi_sample_data %>% filter(cluster_prec == "Medium low")
sum_a_tardi_sample_data_shannon_prec_3 <- capture.output(summary(sub_a_tardi_sample_data_prec_3$Shannon),file=NULL,append=FALSE)
sum_a_tardi_sample_data_shannon_prec_3 <- as.data.frame(sum_a_tardi_sample_data_shannon_prec_3)
sum_a_tardi_sample_data_shannon_prec_3$group <- "Tardigrades"
sum_a_tardi_sample_data_shannon_prec_3$diversity <- "Shannon"
sum_a_tardi_sample_data_shannon_prec_3$class <- "Medium low precipitation"
sum_a_tardi_sample_data_shannon_prec_3 <-sum_a_tardi_sample_data_shannon_prec_3[-1,]

sum_a_tardi_sample_data_richness_prec_3 <- capture.output(summary(sub_a_tardi_sample_data_prec_3$Richness),file=NULL,append=FALSE)
sum_a_tardi_sample_data_richness_prec_3 <- as.data.frame(sum_a_tardi_sample_data_richness_prec_3)
sum_a_tardi_sample_data_richness_prec_3$group <- "Tardigrades"
sum_a_tardi_sample_data_richness_prec_3$diversity <- "Richness"
sum_a_tardi_sample_data_richness_prec_3$class <- "Medium low precipitation"
sum_a_tardi_sample_data_richness_prec_3 <-sum_a_tardi_sample_data_richness_prec_3[-1,]

sub_a_tardi_sample_data_prec_4 <- a_tardi_sample_data %>% filter(cluster_prec == "Medium high")
sum_a_tardi_sample_data_shannon_prec_4 <- capture.output(summary(sub_a_tardi_sample_data_prec_4$Shannon),file=NULL,append=FALSE)
sum_a_tardi_sample_data_shannon_prec_4 <- as.data.frame(sum_a_tardi_sample_data_shannon_prec_4)
sum_a_tardi_sample_data_shannon_prec_4$group <- "Tardigrades"
sum_a_tardi_sample_data_shannon_prec_4$diversity <- "Shannon"
sum_a_tardi_sample_data_shannon_prec_4$class <- "Medium high precipitation"
sum_a_tardi_sample_data_shannon_prec_4 <-sum_a_tardi_sample_data_shannon_prec_4[-1,]

sum_a_tardi_sample_data_richness_prec_4 <- capture.output(summary(sub_a_tardi_sample_data_prec_4$Richness),file=NULL,append=FALSE)
sum_a_tardi_sample_data_richness_prec_4 <- as.data.frame(sum_a_tardi_sample_data_richness_prec_4)
sum_a_tardi_sample_data_richness_prec_4$group <- "Tardigrades"
sum_a_tardi_sample_data_richness_prec_4$diversity <- "Richness"
sum_a_tardi_sample_data_richness_prec_4$class <- "Medium high precipitation"
sum_a_tardi_sample_data_richness_prec_4 <-sum_a_tardi_sample_data_richness_prec_4[-1,]

sub_a_tardi_sample_data_prec_5 <- a_tardi_sample_data %>% filter(cluster_prec == "High")
sum_a_tardi_sample_data_shannon_prec_5 <- capture.output(summary(sub_a_tardi_sample_data_prec_5$Shannon),file=NULL,append=FALSE)
sum_a_tardi_sample_data_shannon_prec_5 <- as.data.frame(sum_a_tardi_sample_data_shannon_prec_5)
sum_a_tardi_sample_data_shannon_prec_5$group <- "Tardigrades"
sum_a_tardi_sample_data_shannon_prec_5$diversity <- "Shannon"
sum_a_tardi_sample_data_shannon_prec_5$class <- "High precipitation"
sum_a_tardi_sample_data_shannon_prec_5 <-sum_a_tardi_sample_data_shannon_prec_5[-1,]

sum_a_tardi_sample_data_richness_prec_5 <- capture.output(summary(sub_a_tardi_sample_data_prec_5$Richness),file=NULL,append=FALSE)
sum_a_tardi_sample_data_richness_prec_5 <- as.data.frame(sum_a_tardi_sample_data_richness_prec_5)
sum_a_tardi_sample_data_richness_prec_5$group <- "Tardigrades"
sum_a_tardi_sample_data_richness_prec_5$diversity <- "Richness"
sum_a_tardi_sample_data_richness_prec_5$class <- "High precipitation"
sum_a_tardi_sample_data_richness_prec_5 <-sum_a_tardi_sample_data_richness_prec_5[-1,]

sub_a_tardi_sample_data_prec_6 <- a_tardi_sample_data %>% filter(cluster_prec == "Very high")
sum_a_tardi_sample_data_shannon_prec_6 <- capture.output(summary(sub_a_tardi_sample_data_prec_6$Shannon),file=NULL,append=FALSE)
sum_a_tardi_sample_data_shannon_prec_6 <- as.data.frame(sum_a_tardi_sample_data_shannon_prec_6)
sum_a_tardi_sample_data_shannon_prec_6$group <- "Tardigrades"
sum_a_tardi_sample_data_shannon_prec_6$diversity <- "Shannon"
sum_a_tardi_sample_data_shannon_prec_6$class <- "Very high precipitation"
sum_a_tardi_sample_data_shannon_prec_6 <-sum_a_tardi_sample_data_shannon_prec_6[-1,]

sum_a_tardi_sample_data_richness_prec_6 <- capture.output(summary(sub_a_tardi_sample_data_prec_6$Richness),file=NULL,append=FALSE)
sum_a_tardi_sample_data_richness_prec_6 <- as.data.frame(sum_a_tardi_sample_data_richness_prec_6)
sum_a_tardi_sample_data_richness_prec_6$group <- "Tardigrades"
sum_a_tardi_sample_data_richness_prec_6$diversity <- "Richness"
sum_a_tardi_sample_data_richness_prec_6$class <- "Very high precipitation"
sum_a_tardi_sample_data_richness_prec_6 <-sum_a_tardi_sample_data_richness_prec_6[-1,] 

sub_a_nema_sample_data_CL_annual <- a_nema_sample_data %>% filter(LC1_2018 == "CL_annual")
sum_a_nema_richness_CL_annual <- capture.output(summary(sub_a_nema_sample_data_CL_annual$Richness),file=NULL,append=FALSE)
sum_a_nema_richness_CL_annual <- as.data.frame(sum_a_nema_richness_CL_annual)
sum_a_nema_richness_CL_annual$group <- "Nematodes"
sum_a_nema_richness_CL_annual$diversity <- "Richness"
sum_a_nema_richness_CL_annual$class <- "Land Cover CL_annual"
sum_a_nema_richness_CL_annual <-sum_a_nema_richness_CL_annual[-1,]

sum_a_nema_sample_data_shannon_CL_annual <- capture.output(summary(sub_a_nema_sample_data_CL_annual$Shannon),file=NULL,append=FALSE)
sum_a_nema_shannon_CL_annual <- as.data.frame(sum_a_nema_sample_data_shannon_CL_annual)
sum_a_nema_shannon_CL_annual$group <- "Nematodes"
sum_a_nema_shannon_CL_annual$diversity <- "Shannon"
sum_a_nema_shannon_CL_annual$class <- "Land Cover CL_annual"
sum_a_nema_shannon_CL_annual <-sum_a_nema_shannon_CL_annual[-1,]

sub_a_nema_sample_data_CL_permanent <- a_nema_sample_data %>% filter(LC1_2018 == "CL_permanent")
sum_a_nema_sample_data_richness_CL_permanent <- capture.output(summary(sub_a_nema_sample_data_CL_permanent$Richness),file=NULL,append=FALSE)
sum_a_nema_richness_CL_permanent <- as.data.frame(sum_a_nema_sample_data_richness_CL_permanent)
sum_a_nema_richness_CL_permanent$group <- "Nematodes"
sum_a_nema_richness_CL_permanent$diversity <- "Richness"
sum_a_nema_richness_CL_permanent$class <- "Land Cover CL_permanent"
sum_a_nema_richness_CL_permanent <-sum_a_nema_richness_CL_permanent[-1,]

sum_a_nema_sample_data_shannon_CL_permanent <- capture.output(summary(sub_a_nema_sample_data_CL_permanent$Shannon),file=NULL,append=FALSE)
sum_a_nema_shannon_CL_permanent <- as.data.frame(sum_a_nema_sample_data_shannon_CL_permanent)
sum_a_nema_shannon_CL_permanent$group <- "Nematodes"
sum_a_nema_shannon_CL_permanent$diversity <- "Shannon"
sum_a_nema_shannon_CL_permanent$class <- "Land Cover CL_permanent"
sum_a_nema_shannon_CL_permanent <-sum_a_nema_shannon_CL_permanent[-1,]

sub_a_nema_sample_data_GL_unmanaged <- a_nema_sample_data %>% filter(LC1_2018 == "GL_unmanaged")
sum_a_nema_sample_data_richness_GL_unmanaged <- capture.output(summary(sub_a_nema_sample_data_GL_unmanaged$Richness),file=NULL,append=FALSE)
sum_a_nema_richness_GL_unmanaged <- as.data.frame(sum_a_nema_sample_data_richness_GL_unmanaged)
sum_a_nema_richness_GL_unmanaged$group <- "Nematodes"
sum_a_nema_richness_GL_unmanaged$diversity <- "Richness"
sum_a_nema_richness_GL_unmanaged$class <- "Land Cover GL_unmanaged"
sum_a_nema_richness_GL_unmanaged <-sum_a_nema_richness_GL_unmanaged[-1,]

sum_a_nema_sample_data_shannon_GL_unmanaged <- capture.output(summary(sub_a_nema_sample_data_GL_unmanaged$Shannon),file=NULL,append=FALSE)
sum_a_nema_shannon_GL_unmanaged <- as.data.frame(sum_a_nema_sample_data_shannon_GL_unmanaged)
sum_a_nema_shannon_GL_unmanaged$group <- "Nematodes"
sum_a_nema_shannon_GL_unmanaged$diversity <- "Shannon"
sum_a_nema_shannon_GL_unmanaged$class <- "Land Cover GL_unmanaged"
sum_a_nema_shannon_GL_unmanaged <-sum_a_nema_shannon_GL_unmanaged[-1,]

sub_a_nema_sample_data_GL_managed <- a_nema_sample_data %>% filter(LC1_2018 == "GL_managed")
sum_a_nema_sample_data_richness_GL_managed <- capture.output(summary(sub_a_nema_sample_data_GL_managed$Richness),file=NULL,append=FALSE)
sum_a_nema_richness_GL_managed <- as.data.frame(sum_a_nema_sample_data_richness_GL_managed)
sum_a_nema_richness_GL_managed$group <- "Nematodes"
sum_a_nema_richness_GL_managed$diversity <- "Richness"
sum_a_nema_richness_GL_managed$class <- "Land Cover GL_managed"
sum_a_nema_richness_GL_managed <-sum_a_nema_richness_GL_managed[-1,]

sum_a_nema_sample_data_shannon_GL_managed <- capture.output(summary(sub_a_nema_sample_data_GL_managed$Shannon),file=NULL,append=FALSE)
sum_a_nema_shannon_GL_managed <- as.data.frame(sum_a_nema_sample_data_shannon_GL_managed)
sum_a_nema_shannon_GL_managed$group <- "Nematodes"
sum_a_nema_shannon_GL_managed$diversity <- "Shannon"
sum_a_nema_shannon_GL_managed$class <- "Land Cover GL_managed"
sum_a_nema_shannon_GL_managed <-sum_a_nema_shannon_GL_managed[-1,]

sub_a_nema_sample_data_WL_broadleaved <- a_nema_sample_data %>% filter(LC1_2018 == "WL_broadleaved")
sum_a_nema_richness_WL_broadleaved <- capture.output(summary(sub_a_nema_sample_data_WL_broadleaved$Richness),file=NULL,append=FALSE)
sum_a_nema_richness_WL_broadleaved <- as.data.frame(sum_a_nema_richness_WL_broadleaved)
sum_a_nema_richness_WL_broadleaved$group <- "Nematodes"
sum_a_nema_richness_WL_broadleaved$diversity <- "Richness"
sum_a_nema_richness_WL_broadleaved$class <- "Land Cover WL_broadleaved"
sum_a_nema_richness_WL_broadleaved <-sum_a_nema_richness_WL_broadleaved[-1,]

sum_a_nema_sample_data_shannon_WL_broadleaved <- capture.output(summary(sub_a_nema_sample_data_WL_broadleaved$Shannon),file=NULL,append=FALSE)
sum_a_nema_shannon_WL_broadleaved <- as.data.frame(sum_a_nema_sample_data_shannon_WL_broadleaved)
sum_a_nema_shannon_WL_broadleaved$group <- "Nematodes"
sum_a_nema_shannon_WL_broadleaved$diversity <- "Shannon"
sum_a_nema_shannon_WL_broadleaved$class <- "Land Cover WL_broadleaved"
sum_a_nema_shannon_WL_broadleaved <-sum_a_nema_shannon_WL_broadleaved[-1,]

sub_a_nema_sample_data_WL_coniferous <- a_nema_sample_data %>% filter(LC1_2018 == "WL_coniferous")
sum_a_nema_sample_data_richness_WL_coniferous <- capture.output(summary(sub_a_nema_sample_data_WL_coniferous$Richness),file=NULL,append=FALSE)
sum_a_nema_richness_WL_coniferous <- as.data.frame(sum_a_nema_sample_data_richness_WL_coniferous)
sum_a_nema_richness_WL_coniferous$group <- "Nematodes"
sum_a_nema_richness_WL_coniferous$diversity <- "Richness"
sum_a_nema_richness_WL_coniferous$class <- "Land Cover WL_coniferous"
sum_a_nema_richness_WL_coniferous <-sum_a_nema_richness_WL_coniferous[-1,]

sum_a_nema_sample_data_shannon_WL_coniferous <- capture.output(summary(sub_a_nema_sample_data_WL_coniferous$Shannon),file=NULL,append=FALSE)
sum_a_nema_shannon_WL_coniferous <- as.data.frame(sum_a_nema_sample_data_shannon_WL_coniferous)
sum_a_nema_shannon_WL_coniferous$group <- "Nematodes"
sum_a_nema_shannon_WL_coniferous$diversity <- "Shannon"
sum_a_nema_shannon_WL_coniferous$class <- "Land Cover WL_coniferous"
sum_a_nema_shannon_WL_coniferous <-sum_a_nema_shannon_WL_coniferous[-1,]

sub_a_nema_sample_data_temp_1 <- a_nema_sample_data %>% filter(cluster_temp == "Very low")
sum_a_nema_sample_data_shannon_temp_1 <- capture.output(summary(sub_a_nema_sample_data_temp_1$Shannon),file=NULL,append=FALSE)
sum_a_nema_sample_data_shannon_temp_1 <- as.data.frame(sum_a_nema_sample_data_shannon_temp_1)
sum_a_nema_sample_data_shannon_temp_1$group <- "Nematodes"
sum_a_nema_sample_data_shannon_temp_1$diversity <- "Shannon"
sum_a_nema_sample_data_shannon_temp_1$class <- "Very low temperature"
sum_a_nema_sample_data_shannon_temp_1 <-sum_a_nema_sample_data_shannon_temp_1[-1,]

sum_a_nema_sample_data_richness_temp_1 <- capture.output(summary(sub_a_nema_sample_data_temp_1$Richness),file=NULL,append=FALSE)
sum_a_nema_sample_data_richness_temp_1 <- as.data.frame(sum_a_nema_sample_data_richness_temp_1)
sum_a_nema_sample_data_richness_temp_1$group <- "Nematodes"
sum_a_nema_sample_data_richness_temp_1$diversity <- "Richness"
sum_a_nema_sample_data_richness_temp_1$class <- "Very low temperature"
sum_a_nema_sample_data_richness_temp_1 <-sum_a_nema_sample_data_richness_temp_1[-1,]

sub_a_nema_sample_data_temp_2 <- a_nema_sample_data %>% filter(cluster_temp == "Low")
sum_a_nema_sample_data_shannon_temp_2 <- capture.output(summary(sub_a_nema_sample_data_temp_2$Shannon),file=NULL,append=FALSE)
sum_a_nema_sample_data_shannon_temp_2 <- as.data.frame(sum_a_nema_sample_data_shannon_temp_2)
sum_a_nema_sample_data_shannon_temp_2$group <- "Nematodes"
sum_a_nema_sample_data_shannon_temp_2$diversity <- "Shannon"
sum_a_nema_sample_data_shannon_temp_2$class <- "Low temperature"
sum_a_nema_sample_data_shannon_temp_2 <-sum_a_nema_sample_data_shannon_temp_2[-1,]

sum_a_nema_sample_data_richness_temp_2 <- capture.output(summary(sub_a_nema_sample_data_temp_2$Richness),file=NULL,append=FALSE)
sum_a_nema_sample_data_richness_temp_2 <- as.data.frame(sum_a_nema_sample_data_richness_temp_2)
sum_a_nema_sample_data_richness_temp_2$group <- "Nematodes"
sum_a_nema_sample_data_richness_temp_2$diversity <- "Richness"
sum_a_nema_sample_data_richness_temp_2$class <- "Low temperature"
sum_a_nema_sample_data_richness_temp_2 <-sum_a_nema_sample_data_richness_temp_2[-1,]

sub_a_nema_sample_data_temp_3 <- a_nema_sample_data %>% filter(cluster_temp == "Medium low")
sum_a_nema_sample_data_shannon_temp_3 <- capture.output(summary(sub_a_nema_sample_data_temp_3$Shannon),file=NULL,append=FALSE)
sum_a_nema_sample_data_shannon_temp_3 <- as.data.frame(sum_a_nema_sample_data_shannon_temp_3)
sum_a_nema_sample_data_shannon_temp_3$group <- "Nematodes"
sum_a_nema_sample_data_shannon_temp_3$diversity <- "Shannon"
sum_a_nema_sample_data_shannon_temp_3$class <- "Medium low temperature"
sum_a_nema_sample_data_shannon_temp_3 <-sum_a_nema_sample_data_shannon_temp_3[-1,]

sum_a_nema_sample_data_richness_temp_3 <- capture.output(summary(sub_a_nema_sample_data_temp_3$Richness),file=NULL,append=FALSE)
sum_a_nema_sample_data_richness_temp_3 <- as.data.frame(sum_a_nema_sample_data_richness_temp_3)
sum_a_nema_sample_data_richness_temp_3$group <- "Nematodes"
sum_a_nema_sample_data_richness_temp_3$diversity <- "Richness"
sum_a_nema_sample_data_richness_temp_3$class <- "Medium low temperature"
sum_a_nema_sample_data_richness_temp_3 <-sum_a_nema_sample_data_richness_temp_3[-1,]

sub_a_nema_sample_data_temp_4 <- a_nema_sample_data %>% filter(cluster_temp == "Medium high")
sum_a_nema_sample_data_shannon_temp_4 <- capture.output(summary(sub_a_nema_sample_data_temp_4$Shannon),file=NULL,append=FALSE)
sum_a_nema_sample_data_shannon_temp_4 <- as.data.frame(sum_a_nema_sample_data_shannon_temp_4)
sum_a_nema_sample_data_shannon_temp_4$group <- "Nematodes"
sum_a_nema_sample_data_shannon_temp_4$diversity <- "Shannon"
sum_a_nema_sample_data_shannon_temp_4$class <- "Medium high temperature"
sum_a_nema_sample_data_shannon_temp_4 <-sum_a_nema_sample_data_shannon_temp_4[-1,]

sum_a_nema_sample_data_richness_temp_4 <- capture.output(summary(sub_a_nema_sample_data_temp_4$Richness),file=NULL,append=FALSE)
sum_a_nema_sample_data_richness_temp_4 <- as.data.frame(sum_a_nema_sample_data_richness_temp_4)
sum_a_nema_sample_data_richness_temp_4$group <- "Nematodes"
sum_a_nema_sample_data_richness_temp_4$diversity <- "Richness"
sum_a_nema_sample_data_richness_temp_4$class <- "Medium high temperature"
sum_a_nema_sample_data_richness_temp_4 <-sum_a_nema_sample_data_richness_temp_4[-1,]

sub_a_nema_sample_data_temp_5 <- a_nema_sample_data %>% filter(cluster_temp == "High")
sum_a_nema_sample_data_shannon_temp_5 <- capture.output(summary(sub_a_nema_sample_data_temp_5$Shannon),file=NULL,append=FALSE)
sum_a_nema_sample_data_shannon_temp_5 <- as.data.frame(sum_a_nema_sample_data_shannon_temp_5)
sum_a_nema_sample_data_shannon_temp_5$group <- "Nematodes"
sum_a_nema_sample_data_shannon_temp_5$diversity <- "Shannon"
sum_a_nema_sample_data_shannon_temp_5$class <- "High temperature"
sum_a_nema_sample_data_shannon_temp_5 <-sum_a_nema_sample_data_shannon_temp_5[-1,]

sum_a_nema_sample_data_richness_temp_5 <- capture.output(summary(sub_a_nema_sample_data_temp_5$Richness),file=NULL,append=FALSE)
sum_a_nema_sample_data_richness_temp_5 <- as.data.frame(sum_a_nema_sample_data_richness_temp_5)
sum_a_nema_sample_data_richness_temp_5$group <- "Nematodes"
sum_a_nema_sample_data_richness_temp_5$diversity <- "Richness"
sum_a_nema_sample_data_richness_temp_5$class <- "High temperature"
sum_a_nema_sample_data_richness_temp_5 <-sum_a_nema_sample_data_richness_temp_5[-1,]

sub_a_nema_sample_data_temp_6 <- a_nema_sample_data %>% filter(cluster_temp == "Very high")
sum_a_nema_sample_data_shannon_temp_6 <- capture.output(summary(sub_a_nema_sample_data_temp_6$Shannon),file=NULL,append=FALSE)
sum_a_nema_sample_data_shannon_temp_6 <- as.data.frame(sum_a_nema_sample_data_shannon_temp_6)
sum_a_nema_sample_data_shannon_temp_6$group <- "Nematodes"
sum_a_nema_sample_data_shannon_temp_6$diversity <- "Shannon"
sum_a_nema_sample_data_shannon_temp_6$class <- "Very high temperature"
sum_a_nema_sample_data_shannon_temp_6 <-sum_a_nema_sample_data_shannon_temp_6[-1,]

sum_a_nema_sample_data_richness_temp_6 <- capture.output(summary(sub_a_nema_sample_data_temp_6$Richness),file=NULL,append=FALSE)
sum_a_nema_sample_data_richness_temp_6 <- as.data.frame(sum_a_nema_sample_data_richness_temp_6)
sum_a_nema_sample_data_richness_temp_6$group <- "Nematodes"
sum_a_nema_sample_data_richness_temp_6$diversity <- "Richness"
sum_a_nema_sample_data_richness_temp_6$class <- "Very high temperature"
sum_a_nema_sample_data_richness_temp_6 <-sum_a_nema_sample_data_richness_temp_6[-1,] 

sub_a_nema_sample_data_temp_range_1 <- a_nema_sample_data %>% filter(cluster_temp_range == "Very low")
sum_a_nema_sample_data_shannon_temp_range_1 <- capture.output(summary(sub_a_nema_sample_data_temp_range_1$Shannon),file=NULL,append=FALSE)
sum_a_nema_sample_data_shannon_temp_range_1 <- as.data.frame(sum_a_nema_sample_data_shannon_temp_range_1)
sum_a_nema_sample_data_shannon_temp_range_1$group <- "Nematodes"
sum_a_nema_sample_data_shannon_temp_range_1$diversity <- "Shannon"
sum_a_nema_sample_data_shannon_temp_range_1$class <- "Very low temperature range"
sum_a_nema_sample_data_shannon_temp_range_1 <-sum_a_nema_sample_data_shannon_temp_range_1[-1,]

sum_a_nema_sample_data_richness_temp_range_1 <- capture.output(summary(sub_a_nema_sample_data_temp_range_1$Richness),file=NULL,append=FALSE)
sum_a_nema_sample_data_richness_temp_range_1 <- as.data.frame(sum_a_nema_sample_data_richness_temp_range_1)
sum_a_nema_sample_data_richness_temp_range_1$group <- "Nematodes"
sum_a_nema_sample_data_richness_temp_range_1$diversity <- "Richness"
sum_a_nema_sample_data_richness_temp_range_1$class <- "Very low temperature range"
sum_a_nema_sample_data_richness_temp_range_1 <-sum_a_nema_sample_data_richness_temp_range_1[-1,]

sub_a_nema_sample_data_temp_range_2 <- a_nema_sample_data %>% filter(cluster_temp_range == "Low")
sum_a_nema_sample_data_shannon_temp_range_2 <- capture.output(summary(sub_a_nema_sample_data_temp_range_2$Shannon),file=NULL,append=FALSE)
sum_a_nema_sample_data_shannon_temp_range_2 <- as.data.frame(sum_a_nema_sample_data_shannon_temp_range_2)
sum_a_nema_sample_data_shannon_temp_range_2$group <- "Nematodes"
sum_a_nema_sample_data_shannon_temp_range_2$diversity <- "Shannon"
sum_a_nema_sample_data_shannon_temp_range_2$class <- "Low temperature range"
sum_a_nema_sample_data_shannon_temp_range_2 <-sum_a_nema_sample_data_shannon_temp_range_2[-1,]

sum_a_nema_sample_data_richness_temp_range_2 <- capture.output(summary(sub_a_nema_sample_data_temp_range_2$Richness),file=NULL,append=FALSE)
sum_a_nema_sample_data_richness_temp_range_2 <- as.data.frame(sum_a_nema_sample_data_richness_temp_range_2)
sum_a_nema_sample_data_richness_temp_range_2$group <- "Nematodes"
sum_a_nema_sample_data_richness_temp_range_2$diversity <- "Richness"
sum_a_nema_sample_data_richness_temp_range_2$class <- "Low temperature range"
sum_a_nema_sample_data_richness_temp_range_2 <-sum_a_nema_sample_data_richness_temp_range_2[-1,]

sub_a_nema_sample_data_temp_range_3 <- a_nema_sample_data %>% filter(cluster_temp_range == "Medium low")
sum_a_nema_sample_data_shannon_temp_range_3 <- capture.output(summary(sub_a_nema_sample_data_temp_range_3$Shannon),file=NULL,append=FALSE)
sum_a_nema_sample_data_shannon_temp_range_3 <- as.data.frame(sum_a_nema_sample_data_shannon_temp_range_3)
sum_a_nema_sample_data_shannon_temp_range_3$group <- "Nematodes"
sum_a_nema_sample_data_shannon_temp_range_3$diversity <- "Shannon"
sum_a_nema_sample_data_shannon_temp_range_3$class <- "Medium low temperature range"
sum_a_nema_sample_data_shannon_temp_range_3 <-sum_a_nema_sample_data_shannon_temp_range_3[-1,]

sum_a_nema_sample_data_richness_temp_range_3 <- capture.output(summary(sub_a_nema_sample_data_temp_range_3$Richness),file=NULL,append=FALSE)
sum_a_nema_sample_data_richness_temp_range_3 <- as.data.frame(sum_a_nema_sample_data_richness_temp_range_3)
sum_a_nema_sample_data_richness_temp_range_3$group <- "Nematodes"
sum_a_nema_sample_data_richness_temp_range_3$diversity <- "Richness"
sum_a_nema_sample_data_richness_temp_range_3$class <- "Medium low temperature range"
sum_a_nema_sample_data_richness_temp_range_3 <-sum_a_nema_sample_data_richness_temp_range_3[-1,]

sub_a_nema_sample_data_temp_range_4 <- a_nema_sample_data %>% filter(cluster_temp_range == "Medium high")
sum_a_nema_sample_data_shannon_temp_range_4 <- capture.output(summary(sub_a_nema_sample_data_temp_range_4$Shannon),file=NULL,append=FALSE)
sum_a_nema_sample_data_shannon_temp_range_4 <- as.data.frame(sum_a_nema_sample_data_shannon_temp_range_4)
sum_a_nema_sample_data_shannon_temp_range_4$group <- "Nematodes"
sum_a_nema_sample_data_shannon_temp_range_4$diversity <- "Shannon"
sum_a_nema_sample_data_shannon_temp_range_4$class <- "Medium high temperature range"
sum_a_nema_sample_data_shannon_temp_range_4 <-sum_a_nema_sample_data_shannon_temp_range_4[-1,]

sum_a_nema_sample_data_richness_temp_range_4 <- capture.output(summary(sub_a_nema_sample_data_temp_range_4$Richness),file=NULL,append=FALSE)
sum_a_nema_sample_data_richness_temp_range_4 <- as.data.frame(sum_a_nema_sample_data_richness_temp_range_4)
sum_a_nema_sample_data_richness_temp_range_4$group <- "Nematodes"
sum_a_nema_sample_data_richness_temp_range_4$diversity <- "Richness"
sum_a_nema_sample_data_richness_temp_range_4$class <- "Medium high temperature range"
sum_a_nema_sample_data_richness_temp_range_4 <-sum_a_nema_sample_data_richness_temp_range_4[-1,]

sub_a_nema_sample_data_temp_range_5 <- a_nema_sample_data %>% filter(cluster_temp_range == "High")
sum_a_nema_sample_data_shannon_temp_range_5 <- capture.output(summary(sub_a_nema_sample_data_temp_range_5$Shannon),file=NULL,append=FALSE)
sum_a_nema_sample_data_shannon_temp_range_5 <- as.data.frame(sum_a_nema_sample_data_shannon_temp_range_5)
sum_a_nema_sample_data_shannon_temp_range_5$group <- "Nematodes"
sum_a_nema_sample_data_shannon_temp_range_5$diversity <- "Shannon"
sum_a_nema_sample_data_shannon_temp_range_5$class <- "High temperature range"
sum_a_nema_sample_data_shannon_temp_range_5 <-sum_a_nema_sample_data_shannon_temp_range_5[-1,]

sum_a_nema_sample_data_richness_temp_range_5 <- capture.output(summary(sub_a_nema_sample_data_temp_range_5$Richness),file=NULL,append=FALSE)
sum_a_nema_sample_data_richness_temp_range_5 <- as.data.frame(sum_a_nema_sample_data_richness_temp_range_5)
sum_a_nema_sample_data_richness_temp_range_5$group <- "Nematodes"
sum_a_nema_sample_data_richness_temp_range_5$diversity <- "Richness"
sum_a_nema_sample_data_richness_temp_range_5$class <- "High temperature range"
sum_a_nema_sample_data_richness_temp_range_5 <-sum_a_nema_sample_data_richness_temp_range_5[-1,]

sub_a_nema_sample_data_temp_range_6 <- a_nema_sample_data %>% filter(cluster_temp_range == "Very high")
sum_a_nema_sample_data_shannon_temp_range_6 <- capture.output(summary(sub_a_nema_sample_data_temp_range_6$Shannon),file=NULL,append=FALSE)
sum_a_nema_sample_data_shannon_temp_range_6 <- as.data.frame(sum_a_nema_sample_data_shannon_temp_range_6)
sum_a_nema_sample_data_shannon_temp_range_6$group <- "Nematodes"
sum_a_nema_sample_data_shannon_temp_range_6$diversity <- "Shannon"
sum_a_nema_sample_data_shannon_temp_range_6$class <- "Very high temperature range"
sum_a_nema_sample_data_shannon_temp_range_6 <-sum_a_nema_sample_data_shannon_temp_range_6[-1,]

sum_a_nema_sample_data_richness_temp_range_6 <- capture.output(summary(sub_a_nema_sample_data_temp_range_6$Richness),file=NULL,append=FALSE)
sum_a_nema_sample_data_richness_temp_range_6 <- as.data.frame(sum_a_nema_sample_data_richness_temp_range_6)
sum_a_nema_sample_data_richness_temp_range_6$group <- "Nematodes"
sum_a_nema_sample_data_richness_temp_range_6$diversity <- "Richness"
sum_a_nema_sample_data_richness_temp_range_6$class <- "Very high temperature range"
sum_a_nema_sample_data_richness_temp_range_6 <-sum_a_nema_sample_data_richness_temp_range_6[-1,]

sub_a_nema_sample_data_prec_1 <- a_nema_sample_data %>% filter(cluster_prec == "Very low")
sum_a_nema_sample_data_shannon_prec_1 <- capture.output(summary(sub_a_nema_sample_data_prec_1$Shannon),file=NULL,append=FALSE)
sum_a_nema_sample_data_shannon_prec_1 <- as.data.frame(sum_a_nema_sample_data_shannon_prec_1)
sum_a_nema_sample_data_shannon_prec_1$group <- "Nematodes"
sum_a_nema_sample_data_shannon_prec_1$diversity <- "Shannon"
sum_a_nema_sample_data_shannon_prec_1$class <- "Very low precipitation "
sum_a_nema_sample_data_shannon_prec_1 <-sum_a_nema_sample_data_shannon_prec_1[-1,]

sum_a_nema_sample_data_richness_prec_1 <- capture.output(summary(sub_a_nema_sample_data_prec_1$Richness),file=NULL,append=FALSE)
sum_a_nema_sample_data_richness_prec_1 <- as.data.frame(sum_a_nema_sample_data_richness_prec_1)
sum_a_nema_sample_data_richness_prec_1$group <- "Nematodes"
sum_a_nema_sample_data_richness_prec_1$diversity <- "Richness"
sum_a_nema_sample_data_richness_prec_1$class <- "Very low precipitation "
sum_a_nema_sample_data_richness_prec_1 <-sum_a_nema_sample_data_richness_prec_1[-1,]

sub_a_nema_sample_data_prec_2 <- a_nema_sample_data %>% filter(cluster_prec == "Low")
sum_a_nema_sample_data_shannon_prec_2 <- capture.output(summary(sub_a_nema_sample_data_prec_2$Shannon),file=NULL,append=FALSE)
sum_a_nema_sample_data_shannon_prec_2 <- as.data.frame(sum_a_nema_sample_data_shannon_prec_2)
sum_a_nema_sample_data_shannon_prec_2$group <- "Nematodes"
sum_a_nema_sample_data_shannon_prec_2$diversity <- "Shannon"
sum_a_nema_sample_data_shannon_prec_2$class <- "Low precipitation"
sum_a_nema_sample_data_shannon_prec_2 <-sum_a_nema_sample_data_shannon_prec_2[-1,]

sum_a_nema_sample_data_richness_prec_2 <- capture.output(summary(sub_a_nema_sample_data_prec_2$Richness),file=NULL,append=FALSE)
sum_a_nema_sample_data_richness_prec_2 <- as.data.frame(sum_a_nema_sample_data_richness_prec_2)
sum_a_nema_sample_data_richness_prec_2$group <- "Nematodes"
sum_a_nema_sample_data_richness_prec_2$diversity <- "Richness"
sum_a_nema_sample_data_richness_prec_2$class <- "Low precipitation"
sum_a_nema_sample_data_richness_prec_2 <-sum_a_nema_sample_data_richness_prec_2[-1,]

sub_a_nema_sample_data_prec_3 <- a_nema_sample_data %>% filter(cluster_prec == "Medium low")
sum_a_nema_sample_data_shannon_prec_3 <- capture.output(summary(sub_a_nema_sample_data_prec_3$Shannon),file=NULL,append=FALSE)
sum_a_nema_sample_data_shannon_prec_3 <- as.data.frame(sum_a_nema_sample_data_shannon_prec_3)
sum_a_nema_sample_data_shannon_prec_3$group <- "Nematodes"
sum_a_nema_sample_data_shannon_prec_3$diversity <- "Shannon"
sum_a_nema_sample_data_shannon_prec_3$class <- "Medium low precipitation"
sum_a_nema_sample_data_shannon_prec_3 <-sum_a_nema_sample_data_shannon_prec_3[-1,]

sum_a_nema_sample_data_richness_prec_3 <- capture.output(summary(sub_a_nema_sample_data_prec_3$Richness),file=NULL,append=FALSE)
sum_a_nema_sample_data_richness_prec_3 <- as.data.frame(sum_a_nema_sample_data_richness_prec_3)
sum_a_nema_sample_data_richness_prec_3$group <- "Nematodes"
sum_a_nema_sample_data_richness_prec_3$diversity <- "Richness"
sum_a_nema_sample_data_richness_prec_3$class <- "Medium low precipitation"
sum_a_nema_sample_data_richness_prec_3 <-sum_a_nema_sample_data_richness_prec_3[-1,]

sub_a_nema_sample_data_prec_4 <- a_nema_sample_data %>% filter(cluster_prec == "Medium high")
sum_a_nema_sample_data_shannon_prec_4 <- capture.output(summary(sub_a_nema_sample_data_prec_4$Shannon),file=NULL,append=FALSE)
sum_a_nema_sample_data_shannon_prec_4 <- as.data.frame(sum_a_nema_sample_data_shannon_prec_4)
sum_a_nema_sample_data_shannon_prec_4$group <- "Nematodes"
sum_a_nema_sample_data_shannon_prec_4$diversity <- "Shannon"
sum_a_nema_sample_data_shannon_prec_4$class <- "Medium high precipitation"
sum_a_nema_sample_data_shannon_prec_4 <-sum_a_nema_sample_data_shannon_prec_4[-1,]

sum_a_nema_sample_data_richness_prec_4 <- capture.output(summary(sub_a_nema_sample_data_prec_4$Richness),file=NULL,append=FALSE)
sum_a_nema_sample_data_richness_prec_4 <- as.data.frame(sum_a_nema_sample_data_richness_prec_4)
sum_a_nema_sample_data_richness_prec_4$group <- "Nematodes"
sum_a_nema_sample_data_richness_prec_4$diversity <- "Richness"
sum_a_nema_sample_data_richness_prec_4$class <- "Medium high precipitation"
sum_a_nema_sample_data_richness_prec_4 <-sum_a_nema_sample_data_richness_prec_4[-1,]

sub_a_nema_sample_data_prec_5 <- a_nema_sample_data %>% filter(cluster_prec == "High")
sum_a_nema_sample_data_shannon_prec_5 <- capture.output(summary(sub_a_nema_sample_data_prec_5$Shannon),file=NULL,append=FALSE)
sum_a_nema_sample_data_shannon_prec_5 <- as.data.frame(sum_a_nema_sample_data_shannon_prec_5)
sum_a_nema_sample_data_shannon_prec_5$group <- "Nematodes"
sum_a_nema_sample_data_shannon_prec_5$diversity <- "Shannon"
sum_a_nema_sample_data_shannon_prec_5$class <- "High precipitation"
sum_a_nema_sample_data_shannon_prec_5 <-sum_a_nema_sample_data_shannon_prec_5[-1,]

sum_a_nema_sample_data_richness_prec_5 <- capture.output(summary(sub_a_nema_sample_data_prec_5$Richness),file=NULL,append=FALSE)
sum_a_nema_sample_data_richness_prec_5 <- as.data.frame(sum_a_nema_sample_data_richness_prec_5)
sum_a_nema_sample_data_richness_prec_5$group <- "Nematodes"
sum_a_nema_sample_data_richness_prec_5$diversity <- "Richness"
sum_a_nema_sample_data_richness_prec_5$class <- "High precipitation"
sum_a_nema_sample_data_richness_prec_5 <-sum_a_nema_sample_data_richness_prec_5[-1,]

sub_a_nema_sample_data_prec_6 <- a_nema_sample_data %>% filter(cluster_prec == "Very high")
sum_a_nema_sample_data_shannon_prec_6 <- capture.output(summary(sub_a_nema_sample_data_prec_6$Shannon),file=NULL,append=FALSE)
sum_a_nema_sample_data_shannon_prec_6 <- as.data.frame(sum_a_nema_sample_data_shannon_prec_6)
sum_a_nema_sample_data_shannon_prec_6$group <- "Nematodes"
sum_a_nema_sample_data_shannon_prec_6$diversity <- "Shannon"
sum_a_nema_sample_data_shannon_prec_6$class <- "Very high precipitation"
sum_a_nema_sample_data_shannon_prec_6 <-sum_a_nema_sample_data_shannon_prec_6[-1,]

sum_a_nema_sample_data_richness_prec_6 <- capture.output(summary(sub_a_nema_sample_data_prec_6$Richness),file=NULL,append=FALSE)
sum_a_nema_sample_data_richness_prec_6 <- as.data.frame(sum_a_nema_sample_data_richness_prec_6)
sum_a_nema_sample_data_richness_prec_6$group <- "Nematodes"
sum_a_nema_sample_data_richness_prec_6$diversity <- "Richness"
sum_a_nema_sample_data_richness_prec_6$class <- "Very high precipitation"
sum_a_nema_sample_data_richness_prec_6 <-sum_a_nema_sample_data_richness_prec_6[-1,] 

sub_p_sample_data_CL1 <- p_sample_data %>% filter(LU2 == "CL1")
sum_p_sample_data_richness_CL1 <- capture.output(summary(sub_p_sample_data_CL1$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_CL1 <- as.data.frame(sum_p_sample_data_richness_CL1)
sum_p_sample_data_richness_CL1$group <- "Protists"
sum_p_sample_data_richness_CL1$diversity <- "Richness"
sum_p_sample_data_richness_CL1$class <- "Land Cover CL1"
sum_p_sample_data_richness_CL1 <- sum_p_sample_data_richness_CL1[-1,]

sum_p_sample_data_shannon_CL1 <- capture.output(summary(sub_p_sample_data_CL1$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_CL1 <- as.data.frame(sum_p_sample_data_shannon_CL1)
sum_p_sample_data_shannon_CL1$group <- "Protists"
sum_p_sample_data_shannon_CL1$diversity <- "Shannon"
sum_p_sample_data_shannon_CL1$class <- "Land Cover CL1"
sum_p_sample_data_shannon_CL1 <-sum_p_sample_data_shannon_CL1[-1,]

sub_p_sample_data_CL2 <- p_sample_data %>% filter(LU2 == "CL2")
sum_p_sample_data_richness_CL2 <- capture.output(summary(sub_p_sample_data_CL2$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_CL2 <- as.data.frame(sum_p_sample_data_richness_CL2)
sum_p_sample_data_richness_CL2$group <- "Protists"
sum_p_sample_data_richness_CL2$diversity <- "Richness"
sum_p_sample_data_richness_CL2$class <- "Land Cover CL2"
sum_p_sample_data_richness_CL2 <-sum_p_sample_data_richness_CL2[-1,]

sum_p_sample_data_shannon_CL2 <- capture.output(summary(sub_p_sample_data_CL2$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_CL2 <- as.data.frame(sum_p_sample_data_shannon_CL2)
sum_p_sample_data_shannon_CL2$group <- "Protists"
sum_p_sample_data_shannon_CL2$diversity <- "Shannon"
sum_p_sample_data_shannon_CL2$class <- "Land Cover CL2"
sum_p_sample_data_shannon_CL2 <-sum_p_sample_data_shannon_CL2[-1,]

sub_p_sample_data_CL3 <- p_sample_data %>% filter(LU2 == "CL3")
sum_p_sample_data_richness_CL3 <- capture.output(summary(sub_p_sample_data_CL3$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_CL3 <- as.data.frame(sum_p_sample_data_richness_CL3)
sum_p_sample_data_richness_CL3$group <- "Protists"
sum_p_sample_data_richness_CL3$diversity <- "Richness"
sum_p_sample_data_richness_CL3$class <- "Land Cover CL3"
sum_p_sample_data_richness_CL3 <- sum_p_sample_data_richness_CL3[-1,]

sum_p_sample_data_shannon_CL3 <- capture.output(summary(sub_p_sample_data_CL3$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_CL3 <- as.data.frame(sum_p_sample_data_shannon_CL3)
sum_p_sample_data_shannon_CL3$group <- "Protists"
sum_p_sample_data_shannon_CL3$diversity <- "Shannon"
sum_p_sample_data_shannon_CL3$class <- "Land Cover CL3"
sum_p_sample_data_shannon_CL3 <-sum_p_sample_data_shannon_CL3[-1,]

sub_p_sample_data_CL4 <- p_sample_data %>% filter(LU2 == "CL4")
sum_p_sample_data_richness_CL4 <- capture.output(summary(sub_p_sample_data_CL4$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_CL4 <- as.data.frame(sum_p_sample_data_richness_CL4)
sum_p_sample_data_richness_CL4$group <- "Protists"
sum_p_sample_data_richness_CL4$diversity <- "Richness"
sum_p_sample_data_richness_CL4$class <- "Land Cover CL4"
sum_p_sample_data_richness_CL4 <-sum_p_sample_data_richness_CL4[-1,]

sum_p_sample_data_shannon_CL4 <- capture.output(summary(sub_p_sample_data_CL4$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_CL4 <- as.data.frame(sum_p_sample_data_shannon_CL4)
sum_p_sample_data_shannon_CL4$group <- "Protists"
sum_p_sample_data_shannon_CL4$diversity <- "Shannon"
sum_p_sample_data_shannon_CL4$class <- "Land Cover CL4"
sum_p_sample_data_shannon_CL4 <-sum_p_sample_data_shannon_CL4[-1,]

sub_p_sample_data_GL1 <- p_sample_data %>% filter(LU2 == "GL1")
sum_p_sample_data_richness_GL1 <- capture.output(summary(sub_p_sample_data_GL1$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_GL1 <- as.data.frame(sum_p_sample_data_richness_GL1)
sum_p_sample_data_richness_GL1$group <- "Protists"
sum_p_sample_data_richness_GL1$diversity <- "Richness"
sum_p_sample_data_richness_GL1$class <- "Land Cover GL1"
sum_p_sample_data_richness_GL1 <- sum_p_sample_data_richness_GL1[-1,]

sum_p_sample_data_shannon_GL1 <- capture.output(summary(sub_p_sample_data_GL1$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_GL1 <- as.data.frame(sum_p_sample_data_shannon_GL1)
sum_p_sample_data_shannon_GL1$group <- "Protists"
sum_p_sample_data_shannon_GL1$diversity <- "Shannon"
sum_p_sample_data_shannon_GL1$class <- "Land Cover GL1"
sum_p_sample_data_shannon_GL1 <-sum_p_sample_data_shannon_GL1[-1,]

sub_p_sample_data_GL2 <- p_sample_data %>% filter(LU2 == "GL2")
sum_p_sample_data_richness_GL2 <- capture.output(summary(sub_p_sample_data_GL2$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_GL2 <- as.data.frame(sum_p_sample_data_richness_GL2)
sum_p_sample_data_richness_GL2$group <- "Protists"
sum_p_sample_data_richness_GL2$diversity <- "Richness"
sum_p_sample_data_richness_GL2$class <- "Land Cover GL2"
sum_p_sample_data_richness_GL2 <- sum_p_sample_data_richness_GL2[-1,]

sum_p_sample_data_shannon_GL2 <- capture.output(summary(sub_p_sample_data_GL2$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_GL2 <- as.data.frame(sum_p_sample_data_shannon_GL2)
sum_p_sample_data_shannon_GL2$group <- "Protists"
sum_p_sample_data_shannon_GL2$diversity <- "Shannon"
sum_p_sample_data_shannon_GL2$class <- "Land Cover GL2"
sum_p_sample_data_shannon_GL2 <-sum_p_sample_data_shannon_GL2[-1,]

sub_p_sample_data_GL3 <- p_sample_data %>% filter(LU2 == "GL3")
sum_p_sample_data_richness_GL3 <- capture.output(summary(sub_p_sample_data_GL3$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_GL3 <- as.data.frame(sum_p_sample_data_richness_GL3)
sum_p_sample_data_richness_GL3$group <- "Protists"
sum_p_sample_data_richness_GL3$diversity <- "Richness"
sum_p_sample_data_richness_GL3$class <- "Land Cover GL3"
sum_p_sample_data_richness_GL3 <- sum_p_sample_data_richness_GL3[-1,]

sum_p_sample_data_shannon_GL3 <- capture.output(summary(sub_p_sample_data_GL3$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_GL3 <- as.data.frame(sum_p_sample_data_shannon_GL3)
sum_p_sample_data_shannon_GL3$group <- "Protists"
sum_p_sample_data_shannon_GL3$diversity <- "Shannon"
sum_p_sample_data_shannon_GL3$class <- "Land Cover GL3"
sum_p_sample_data_shannon_GL3 <-sum_p_sample_data_shannon_GL3[-1,]

sub_p_sample_data_WL1 <- p_sample_data %>% filter(LU2 == "WL1")
sum_p_sample_data_richness_WL1 <- capture.output(summary(sub_p_sample_data_WL1$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_WL1 <- as.data.frame(sum_p_sample_data_richness_WL1)
sum_p_sample_data_richness_WL1$group <- "Protists"
sum_p_sample_data_richness_WL1$diversity <- "Richness"
sum_p_sample_data_richness_WL1$class <- "Land Cover WL1"
sum_p_sample_data_richness_WL1 <- sum_p_sample_data_richness_WL1[-1,]

sum_p_sample_data_shannon_WL1 <- capture.output(summary(sub_p_sample_data_WL1$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_WL1 <- as.data.frame(sum_p_sample_data_shannon_WL1)
sum_p_sample_data_shannon_WL1$group <- "Protists"
sum_p_sample_data_shannon_WL1$diversity <- "Shannon"
sum_p_sample_data_shannon_WL1$class <- "Land Cover WL1"
sum_p_sample_data_shannon_WL1 <-sum_p_sample_data_shannon_WL1[-1,]

sub_p_sample_data_WL2 <- p_sample_data %>% filter(LU2 == "WL2")
sum_p_sample_data_richness_WL2 <- capture.output(summary(sub_p_sample_data_WL2$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_WL2 <- as.data.frame(sum_p_sample_data_richness_WL2)
sum_p_sample_data_richness_WL2$group <- "Protists"
sum_p_sample_data_richness_WL2$diversity <- "Richness"
sum_p_sample_data_richness_WL2$class <- "Land Cover WL2"
sum_p_sample_data_richness_WL2 <- sum_p_sample_data_richness_WL2[-1,]

sum_p_sample_data_shannon_WL2 <- capture.output(summary(sub_p_sample_data_WL2$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_WL2 <- as.data.frame(sum_p_sample_data_shannon_WL2)
sum_p_sample_data_shannon_WL2$group <- "Protists"
sum_p_sample_data_shannon_WL2$diversity <- "Shannon"
sum_p_sample_data_shannon_WL2$class <- "Land Cover WL2"
sum_p_sample_data_shannon_WL2 <-sum_p_sample_data_shannon_WL2[-1,]

sub_p_sample_data_WL3 <- p_sample_data %>% filter(LU2 == "WL3")
sum_p_sample_data_richness_WL3 <- capture.output(summary(sub_p_sample_data_WL3$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_WL3 <- as.data.frame(sum_p_sample_data_richness_WL3)
sum_p_sample_data_richness_WL3$group <- "Protists"
sum_p_sample_data_richness_WL3$diversity <- "Richness"
sum_p_sample_data_richness_WL3$class <- "Land Cover WL3"
sum_p_sample_data_richness_WL3 <- sum_p_sample_data_richness_WL3[-1,]

sum_p_sample_data_shannon_WL3 <- capture.output(summary(sub_p_sample_data_WL3$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_WL3 <- as.data.frame(sum_p_sample_data_shannon_WL3)
sum_p_sample_data_shannon_WL3$group <- "Protists"
sum_p_sample_data_shannon_WL3$diversity <- "Shannon"
sum_p_sample_data_shannon_WL3$class <- "Land Cover WL3"
sum_p_sample_data_shannon_WL3 <-sum_p_sample_data_shannon_WL3[-1,]

sub_p_sample_data_WL4 <- p_sample_data %>% filter(LU2 == "WL4")
sum_p_sample_data_richness_WL4 <- capture.output(summary(sub_p_sample_data_WL4$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_WL4 <- as.data.frame(sum_p_sample_data_richness_WL4)
sum_p_sample_data_richness_WL4$group <- "Protists"
sum_p_sample_data_richness_WL4$diversity <- "Richness"
sum_p_sample_data_richness_WL4$class <- "Land Cover WL4"
sum_p_sample_data_richness_WL4 <- sum_p_sample_data_richness_WL4[-1,]

sum_p_sample_data_shannon_WL4 <- capture.output(summary(sub_p_sample_data_WL4$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_WL4 <- as.data.frame(sum_p_sample_data_shannon_WL4)
sum_p_sample_data_shannon_WL4$group <- "Protists"
sum_p_sample_data_shannon_WL4$diversity <- "Shannon"
sum_p_sample_data_shannon_WL4$class <- "Land Cover WL4"
sum_p_sample_data_shannon_WL4 <-sum_p_sample_data_shannon_WL4[-1,]

sub_p_sample_data_No_C <- p_sample_data %>% filter(LU2 == "No_C")
sum_p_sample_data_richness_No_C <- capture.output(summary(sub_p_sample_data_No_C$Richness),file=NULL,append=FALSE)
sum_p_sample_data_richness_No_C <- as.data.frame(sum_p_sample_data_richness_No_C)
sum_p_sample_data_richness_No_C$group <- "Protists"
sum_p_sample_data_richness_No_C$diversity <- "Richness"
sum_p_sample_data_richness_No_C$class <- "Land Cover No_C"
sum_p_sample_data_richness_No_C <- sum_p_sample_data_richness_No_C[-1,]

sum_p_sample_data_shannon_No_C <- capture.output(summary(sub_p_sample_data_No_C$Shannon),file=NULL,append=FALSE)
sum_p_sample_data_shannon_No_C <- as.data.frame(sum_p_sample_data_shannon_No_C)
sum_p_sample_data_shannon_No_C$group <- "Protists"
sum_p_sample_data_shannon_No_C$diversity <- "Shannon"
sum_p_sample_data_shannon_No_C$class <- "Land Cover No_C"
sum_p_sample_data_shannon_No_C <-sum_p_sample_data_shannon_No_C[-1,]

sub_f_sample_data_CL1 <- f_sample_data %>% filter(LU2 == "CL1")
sum_f_sample_data_richness_CL1 <- capture.output(summary(sub_f_sample_data_CL1$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_CL1 <- as.data.frame(sum_f_sample_data_richness_CL1)
sum_f_sample_data_richness_CL1$group <- "Fungi"
sum_f_sample_data_richness_CL1$diversity <- "Richness"
sum_f_sample_data_richness_CL1$class <- "Land Cover CL1"
sum_f_sample_data_richness_CL1 <- sum_f_sample_data_richness_CL1[-1,]

sum_f_sample_data_shannon_CL1 <- capture.output(summary(sub_f_sample_data_CL1$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_CL1 <- as.data.frame(sum_f_sample_data_shannon_CL1)
sum_f_sample_data_shannon_CL1$group <- "Fungi"
sum_f_sample_data_shannon_CL1$diversity <- "Shannon"
sum_f_sample_data_shannon_CL1$class <- "Land Cover CL1"
sum_f_sample_data_shannon_CL1 <-sum_f_sample_data_shannon_CL1[-1,]

sub_f_sample_data_CL2 <- p_sample_data %>% filter(LU2 == "CL2")
sum_f_sample_data_richness_CL2 <- capture.output(summary(sub_f_sample_data_CL2$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_CL2 <- as.data.frame(sum_f_sample_data_richness_CL2)
sum_f_sample_data_richness_CL2$group <- "Fungi"
sum_f_sample_data_richness_CL2$diversity <- "Richness"
sum_f_sample_data_richness_CL2$class <- "Land Cover CL2"
sum_f_sample_data_richness_CL2 <- sum_f_sample_data_richness_CL2[-1,]

sum_f_sample_data_shannon_CL2 <- capture.output(summary(sub_f_sample_data_CL2$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_CL2 <- as.data.frame(sum_f_sample_data_shannon_CL2)
sum_f_sample_data_shannon_CL2$group <- "Fungi"
sum_f_sample_data_shannon_CL2$diversity <- "Shannon"
sum_f_sample_data_shannon_CL2$class <- "Land Cover CL2"
sum_f_sample_data_shannon_CL2 <-sum_f_sample_data_shannon_CL2[-1,]

sub_f_sample_data_CL3 <- p_sample_data %>% filter(LU2 == "CL3")
sum_f_sample_data_richness_CL3 <- capture.output(summary(sub_f_sample_data_CL3$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_CL3 <- as.data.frame(sum_f_sample_data_richness_CL3)
sum_f_sample_data_richness_CL3$group <- "Fungi"
sum_f_sample_data_richness_CL3$diversity <- "Richness"
sum_f_sample_data_richness_CL3$class <- "Land Cover CL3"
sum_f_sample_data_richness_CL3 <- sum_f_sample_data_richness_CL3[-1,]

sum_f_sample_data_shannon_CL3 <- capture.output(summary(sub_f_sample_data_CL3$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_CL3 <- as.data.frame(sum_f_sample_data_shannon_CL3)
sum_f_sample_data_shannon_CL3$group <- "Fungi"
sum_f_sample_data_shannon_CL3$diversity <- "Shannon"
sum_f_sample_data_shannon_CL3$class <- "Land Cover CL3"
sum_f_sample_data_shannon_CL3 <-sum_f_sample_data_shannon_CL3[-1,]

sub_f_sample_data_CL4 <- p_sample_data %>% filter(LU2 == "CL4")
sum_f_sample_data_richness_CL4 <- capture.output(summary(sub_f_sample_data_CL4$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_CL4 <- as.data.frame(sum_f_sample_data_richness_CL4)
sum_f_sample_data_richness_CL4$group <- "Fungi"
sum_f_sample_data_richness_CL4$diversity <- "Richness"
sum_f_sample_data_richness_CL4$class <- "Land Cover CL4"
sum_f_sample_data_richness_CL4 <- sum_f_sample_data_richness_CL4[-1,]

sum_f_sample_data_shannon_CL4 <- capture.output(summary(sub_f_sample_data_CL4$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_CL4 <- as.data.frame(sum_f_sample_data_shannon_CL4)
sum_f_sample_data_shannon_CL4$group <- "Fungi"
sum_f_sample_data_shannon_CL4$diversity <- "Shannon"
sum_f_sample_data_shannon_CL4$class <- "Land Cover CL4"
sum_f_sample_data_shannon_CL4 <-sum_f_sample_data_shannon_CL4[-1,]

sub_f_sample_data_GL1 <- p_sample_data %>% filter(LU2 == "GL1")
sum_f_sample_data_richness_GL1 <- capture.output(summary(sub_f_sample_data_GL1$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_GL1 <- as.data.frame(sum_f_sample_data_richness_GL1)
sum_f_sample_data_richness_GL1$group <- "Fungi"
sum_f_sample_data_richness_GL1$diversity <- "Richness"
sum_f_sample_data_richness_GL1$class <- "Land Cover GL1"
sum_f_sample_data_richness_GL1 <- sum_f_sample_data_richness_GL1[-1,]

sum_f_sample_data_shannon_GL1 <- capture.output(summary(sub_f_sample_data_GL1$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_GL1 <- as.data.frame(sum_f_sample_data_shannon_GL1)
sum_f_sample_data_shannon_GL1$group <- "Fungi"
sum_f_sample_data_shannon_GL1$diversity <- "Shannon"
sum_f_sample_data_shannon_GL1$class <- "Land Cover GL1"
sum_f_sample_data_shannon_GL1 <-sum_f_sample_data_shannon_GL1[-1,]

sub_f_sample_data_GL2 <- p_sample_data %>% filter(LU2 == "GL2")
sum_f_sample_data_richness_GL2 <- capture.output(summary(sub_f_sample_data_GL2$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_GL2 <- as.data.frame(sum_f_sample_data_richness_GL2)
sum_f_sample_data_richness_GL2$group <- "Fungi"
sum_f_sample_data_richness_GL2$diversity <- "Richness"
sum_f_sample_data_richness_GL2$class <- "Land Cover GL2"
sum_f_sample_data_richness_GL2 <- sum_f_sample_data_richness_GL2[-1,]

sum_f_sample_data_shannon_GL2 <- capture.output(summary(sub_f_sample_data_GL2$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_GL2 <- as.data.frame(sum_f_sample_data_shannon_GL2)
sum_f_sample_data_shannon_GL2$group <- "Fungi"
sum_f_sample_data_shannon_GL2$diversity <- "Shannon"
sum_f_sample_data_shannon_GL2$class <- "Land Cover GL2"
sum_f_sample_data_shannon_GL2 <-sum_f_sample_data_shannon_GL2[-1,]

sub_f_sample_data_GL3 <- p_sample_data %>% filter(LU2 == "GL3")
sum_f_sample_data_richness_GL3 <- capture.output(summary(sub_f_sample_data_GL3$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_GL3 <- as.data.frame(sum_f_sample_data_richness_GL3)
sum_f_sample_data_richness_GL3$group <- "Fungi"
sum_f_sample_data_richness_GL3$diversity <- "Richness"
sum_f_sample_data_richness_GL3$class <- "Land Cover GL3"
sum_f_sample_data_richness_GL3 <- sum_f_sample_data_richness_GL3[-1,]

sum_f_sample_data_shannon_GL3 <- capture.output(summary(sub_f_sample_data_GL3$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_GL3 <- as.data.frame(sum_f_sample_data_shannon_GL3)
sum_f_sample_data_shannon_GL3$group <- "Fungi"
sum_f_sample_data_shannon_GL3$diversity <- "Shannon"
sum_f_sample_data_shannon_GL3$class <- "Land Cover GL3"
sum_f_sample_data_shannon_GL3 <-sum_f_sample_data_shannon_GL3[-1,]

sub_f_sample_data_WL1 <- p_sample_data %>% filter(LU2 == "WL1")
sum_f_sample_data_richness_WL1 <- capture.output(summary(sub_f_sample_data_WL1$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_WL1 <- as.data.frame(sum_f_sample_data_richness_WL1)
sum_f_sample_data_richness_WL1$group <- "Fungi"
sum_f_sample_data_richness_WL1$diversity <- "Richness"
sum_f_sample_data_richness_WL1$class <- "Land Cover WL1"
sum_f_sample_data_richness_WL1 <- sum_f_sample_data_richness_WL1[-1,]

sum_f_sample_data_shannon_WL1 <- capture.output(summary(sub_f_sample_data_WL1$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_WL1 <- as.data.frame(sum_f_sample_data_shannon_WL1)
sum_f_sample_data_shannon_WL1$group <- "Fungi"
sum_f_sample_data_shannon_WL1$diversity <- "Shannon"
sum_f_sample_data_shannon_WL1$class <- "Land Cover WL1"
sum_f_sample_data_shannon_WL1 <-sum_f_sample_data_shannon_WL1[-1,]

sub_f_sample_data_WL2 <- p_sample_data %>% filter(LU2 == "WL2")
sum_f_sample_data_richness_WL2 <- capture.output(summary(sub_f_sample_data_WL2$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_WL2 <- as.data.frame(sum_f_sample_data_richness_WL2)
sum_f_sample_data_richness_WL2$group <- "Fungi"
sum_f_sample_data_richness_WL2$diversity <- "Richness"
sum_f_sample_data_richness_WL2$class <- "Land Cover WL2"
sum_f_sample_data_richness_WL2 <- sum_f_sample_data_richness_WL2[-1,]

sum_f_sample_data_shannon_WL2 <- capture.output(summary(sub_f_sample_data_WL2$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_WL2 <- as.data.frame(sum_f_sample_data_shannon_WL2)
sum_f_sample_data_shannon_WL2$group <- "Fungi"
sum_f_sample_data_shannon_WL2$diversity <- "Shannon"
sum_f_sample_data_shannon_WL2$class <- "Land Cover WL2"
sum_f_sample_data_shannon_WL2 <-sum_f_sample_data_shannon_WL2[-1,]

sub_f_sample_data_WL3 <- p_sample_data %>% filter(LU2 == "WL3")
sum_f_sample_data_richness_WL3 <- capture.output(summary(sub_f_sample_data_WL3$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_WL3 <- as.data.frame(sum_f_sample_data_richness_WL3)
sum_f_sample_data_richness_WL3$group <- "Fungi"
sum_f_sample_data_richness_WL3$diversity <- "Richness"
sum_f_sample_data_richness_WL3$class <- "Land Cover WL3"
sum_f_sample_data_richness_WL3 <- sum_f_sample_data_richness_WL3[-1,]

sum_f_sample_data_shannon_WL3 <- capture.output(summary(sub_f_sample_data_WL3$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_WL3 <- as.data.frame(sum_f_sample_data_shannon_WL3)
sum_f_sample_data_shannon_WL3$group <- "Fungi"
sum_f_sample_data_shannon_WL3$diversity <- "Shannon"
sum_f_sample_data_shannon_WL3$class <- "Land Cover WL3"
sum_f_sample_data_shannon_WL3 <-sum_f_sample_data_shannon_WL3[-1,]

sub_f_sample_data_WL4 <- p_sample_data %>% filter(LU2 == "WL4")
sum_f_sample_data_richness_WL4 <- capture.output(summary(sub_f_sample_data_WL4$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_WL4 <- as.data.frame(sum_f_sample_data_richness_WL4)
sum_f_sample_data_richness_WL4$group <- "Fungi"
sum_f_sample_data_richness_WL4$diversity <- "Richness"
sum_f_sample_data_richness_WL4$class <- "Land Cover WL4"
sum_f_sample_data_richness_WL4 <- sum_f_sample_data_richness_WL4[-1,]

sum_f_sample_data_shannon_WL4 <- capture.output(summary(sub_f_sample_data_WL4$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_WL4 <- as.data.frame(sum_f_sample_data_shannon_WL4)
sum_f_sample_data_shannon_WL4$group <- "Fungi"
sum_f_sample_data_shannon_WL4$diversity <- "Shannon"
sum_f_sample_data_shannon_WL4$class <- "Land Cover WL4"
sum_f_sample_data_shannon_WL4 <-sum_f_sample_data_shannon_WL4[-1,]

sub_f_sample_data_No_C <- p_sample_data %>% filter(LU2 == "No_C")
sum_f_sample_data_richness_No_C <- capture.output(summary(sub_f_sample_data_No_C$Richness),file=NULL,append=FALSE)
sum_f_sample_data_richness_No_C <- as.data.frame(sum_f_sample_data_richness_No_C)
sum_f_sample_data_richness_No_C$group <- "Fungi"
sum_f_sample_data_richness_No_C$diversity <- "Richness"
sum_f_sample_data_richness_No_C$class <- "Land Cover No_C"
sum_f_sample_data_richness_No_C <- sum_f_sample_data_richness_No_C[-1,]

sum_f_sample_data_shannon_No_C <- capture.output(summary(sub_f_sample_data_No_C$Shannon),file=NULL,append=FALSE)
sum_f_sample_data_shannon_No_C <- as.data.frame(sum_f_sample_data_shannon_No_C)
sum_f_sample_data_shannon_No_C$group <- "Fungi"
sum_f_sample_data_shannon_No_C$diversity <- "Shannon"
sum_f_sample_data_shannon_No_C$class <- "Land Cover No_C"
sum_f_sample_data_shannon_No_C <-sum_f_sample_data_shannon_No_C[-1,]

sub_roti_sample_data_CL1 <- a_roti_sample_data %>% filter(LU2 == "CL1")
sum_roti_sample_data_richness_CL1 <- capture.output(summary(sub_roti_sample_data_CL1$Richness),file=NULL,append=FALSE)
sum_roti_sample_data_richness_CL1 <- as.data.frame(sum_roti_sample_data_richness_CL1)
sum_roti_sample_data_richness_CL1$group <- "Rotifers"
sum_roti_sample_data_richness_CL1$diversity <- "Richness"
sum_roti_sample_data_richness_CL1$class <- "Land Cover CL1"
sum_roti_sample_data_richness_CL1 <- sum_roti_sample_data_richness_CL1[-1,]

sum_roti_sample_data_shannon_CL1 <- capture.output(summary(sub_roti_sample_data_CL1$Shannon),file=NULL,append=FALSE)
sum_roti_sample_data_shannon_CL1 <- as.data.frame(sum_roti_sample_data_shannon_CL1)
sum_roti_sample_data_shannon_CL1$group <- "Rotifers"
sum_roti_sample_data_shannon_CL1$diversity <- "Shannon"
sum_roti_sample_data_shannon_CL1$class <- "Land Cover CL1"
sum_roti_sample_data_shannon_CL1 <-sum_roti_sample_data_shannon_CL1[-1,]

sub_roti_sample_data_CL2 <- a_roti_sample_data %>% filter(LU2 == "CL2")
sum_roti_sample_data_richness_CL2 <- capture.output(summary(sub_roti_sample_data_CL2$Richness),file=NULL,append=FALSE)
sum_roti_sample_data_richness_CL2 <- as.data.frame(sum_roti_sample_data_richness_CL2)
sum_roti_sample_data_richness_CL2$group <- "Rotifers"
sum_roti_sample_data_richness_CL2$diversity <- "Richness"
sum_roti_sample_data_richness_CL2$class <- "Land Cover CL2"
sum_roti_sample_data_richness_CL2 <- sum_roti_sample_data_richness_CL2[-1,]

sum_roti_sample_data_shannon_CL2 <- capture.output(summary(sub_roti_sample_data_CL2$Shannon),file=NULL,append=FALSE)
sum_roti_sample_data_shannon_CL2 <- as.data.frame(sum_roti_sample_data_shannon_CL2)
sum_roti_sample_data_shannon_CL2$group <- "Rotifers"
sum_roti_sample_data_shannon_CL2$diversity <- "Shannon"
sum_roti_sample_data_shannon_CL2$class <- "Land Cover CL2"
sum_roti_sample_data_shannon_CL2 <-sum_roti_sample_data_shannon_CL2[-1,]

sub_roti_sample_data_CL3 <- a_roti_sample_data %>% filter(LU2 == "CL3")
sum_roti_sample_data_richness_CL3 <- capture.output(summary(sub_roti_sample_data_CL3$Richness),file=NULL,append=FALSE)
sum_roti_sample_data_richness_CL3 <- as.data.frame(sum_roti_sample_data_richness_CL3)
sum_roti_sample_data_richness_CL3$group <- "Rotifers"
sum_roti_sample_data_richness_CL3$diversity <- "Richness"
sum_roti_sample_data_richness_CL3$class <- "Land Cover CL3"
sum_roti_sample_data_richness_CL3 <- sum_roti_sample_data_richness_CL3[-1,]

sum_roti_sample_data_shannon_CL3 <- capture.output(summary(sub_roti_sample_data_CL3$Shannon),file=NULL,append=FALSE)
sum_roti_sample_data_shannon_CL3 <- as.data.frame(sum_roti_sample_data_shannon_CL3)
sum_roti_sample_data_shannon_CL3$group <- "Rotifers"
sum_roti_sample_data_shannon_CL3$diversity <- "Shannon"
sum_roti_sample_data_shannon_CL3$class <- "Land Cover CL3"
sum_roti_sample_data_shannon_CL3 <-sum_roti_sample_data_shannon_CL3[-1,]

sub_roti_sample_data_CL4 <- a_roti_sample_data %>% filter(LU2 == "CL4")
sum_roti_sample_data_richness_CL4 <- capture.output(summary(sub_roti_sample_data_CL4$Richness),file=NULL,append=FALSE)
sum_roti_sample_data_richness_CL4 <- as.data.frame(sum_roti_sample_data_richness_CL4)
sum_roti_sample_data_richness_CL4$group <- "Rotifers"
sum_roti_sample_data_richness_CL4$diversity <- "Richness"
sum_roti_sample_data_richness_CL4$class <- "Land Cover CL4"
sum_roti_sample_data_richness_CL4 <- sum_roti_sample_data_richness_CL4[-1,]

sum_roti_sample_data_shannon_CL4 <- capture.output(summary(sub_roti_sample_data_CL4$Shannon),file=NULL,append=FALSE)
sum_roti_sample_data_shannon_CL4 <- as.data.frame(sum_roti_sample_data_shannon_CL4)
sum_roti_sample_data_shannon_CL4$group <- "Rotifers"
sum_roti_sample_data_shannon_CL4$diversity <- "Shannon"
sum_roti_sample_data_shannon_CL4$class <- "Land Cover CL4"
sum_roti_sample_data_shannon_CL4 <-sum_roti_sample_data_shannon_CL4[-1,]

sub_roti_sample_data_GL1 <- a_roti_sample_data %>% filter(LU2 == "GL1")
sum_roti_sample_data_richness_GL1 <- capture.output(summary(sub_roti_sample_data_GL1$Richness),file=NULL,append=FALSE)
sum_roti_sample_data_richness_GL1 <- as.data.frame(sum_roti_sample_data_richness_GL1)
sum_roti_sample_data_richness_GL1$group <- "Rotifers"
sum_roti_sample_data_richness_GL1$diversity <- "Richness"
sum_roti_sample_data_richness_GL1$class <- "Land Cover GL1"
sum_roti_sample_data_richness_GL1 <- sum_roti_sample_data_richness_GL1[-1,]

sum_roti_sample_data_shannon_GL1 <- capture.output(summary(sub_roti_sample_data_GL1$Shannon),file=NULL,append=FALSE)
sum_roti_sample_data_shannon_GL1 <- as.data.frame(sum_roti_sample_data_shannon_GL1)
sum_roti_sample_data_shannon_GL1$group <- "Rotifers"
sum_roti_sample_data_shannon_GL1$diversity <- "Shannon"
sum_roti_sample_data_shannon_GL1$class <- "Land Cover GL1"
sum_roti_sample_data_shannon_GL1 <-sum_roti_sample_data_shannon_GL1[-1,]

sub_roti_sample_data_GL2 <- a_roti_sample_data %>% filter(LU2 == "GL2")
sum_roti_sample_data_richness_GL2 <- capture.output(summary(sub_roti_sample_data_GL2$Richness),file=NULL,append=FALSE)
sum_roti_sample_data_richness_GL2 <- as.data.frame(sum_roti_sample_data_richness_GL2)
sum_roti_sample_data_richness_GL2$group <- "Rotifers"
sum_roti_sample_data_richness_GL2$diversity <- "Richness"
sum_roti_sample_data_richness_GL2$class <- "Land Cover GL2"
sum_roti_sample_data_richness_GL2 <- sum_roti_sample_data_richness_GL2[-1,]

sum_roti_sample_data_shannon_GL2 <- capture.output(summary(sub_roti_sample_data_GL2$Shannon),file=NULL,append=FALSE)
sum_roti_sample_data_shannon_GL2 <- as.data.frame(sum_roti_sample_data_shannon_GL2)
sum_roti_sample_data_shannon_GL2$group <- "Rotifers"
sum_roti_sample_data_shannon_GL2$diversity <- "Shannon"
sum_roti_sample_data_shannon_GL2$class <- "Land Cover GL2"
sum_roti_sample_data_shannon_GL2 <-sum_roti_sample_data_shannon_GL2[-1,]

sub_roti_sample_data_GL3 <- a_roti_sample_data %>% filter(LU2 == "GL3")
sum_roti_sample_data_richness_GL3 <- capture.output(summary(sub_roti_sample_data_GL3$Richness),file=NULL,append=FALSE)
sum_roti_sample_data_richness_GL3 <- as.data.frame(sum_roti_sample_data_richness_GL3)
sum_roti_sample_data_richness_GL3$group <- "Rotifers"
sum_roti_sample_data_richness_GL3$diversity <- "Richness"
sum_roti_sample_data_richness_GL3$class <- "Land Cover GL3"
sum_roti_sample_data_richness_GL3 <- sum_roti_sample_data_richness_GL3[-1,]

sum_roti_sample_data_shannon_GL3 <- capture.output(summary(sub_roti_sample_data_GL3$Shannon),file=NULL,append=FALSE)
sum_roti_sample_data_shannon_GL3 <- as.data.frame(sum_roti_sample_data_shannon_GL3)
sum_roti_sample_data_shannon_GL3$group <- "Rotifers"
sum_roti_sample_data_shannon_GL3$diversity <- "Shannon"
sum_roti_sample_data_shannon_GL3$class <- "Land Cover GL3"
sum_roti_sample_data_shannon_GL3 <-sum_roti_sample_data_shannon_GL3[-1,]

sub_roti_sample_data_WL1 <- a_roti_sample_data %>% filter(LU2 == "WL1")
sum_roti_sample_data_richness_WL1 <- capture.output(summary(sub_roti_sample_data_WL1$Richness),file=NULL,append=FALSE)
sum_roti_sample_data_richness_WL1 <- as.data.frame(sum_roti_sample_data_richness_WL1)
sum_roti_sample_data_richness_WL1$group <- "Rotifers"
sum_roti_sample_data_richness_WL1$diversity <- "Richness"
sum_roti_sample_data_richness_WL1$class <- "Land Cover WL1"
sum_roti_sample_data_richness_WL1 <- sum_roti_sample_data_richness_WL1[-1,]

sum_roti_sample_data_shannon_WL1 <- capture.output(summary(sub_roti_sample_data_WL1$Shannon),file=NULL,append=FALSE)
sum_roti_sample_data_shannon_WL1 <- as.data.frame(sum_roti_sample_data_shannon_WL1)
sum_roti_sample_data_shannon_WL1$group <- "Rotifers"
sum_roti_sample_data_shannon_WL1$diversity <- "Shannon"
sum_roti_sample_data_shannon_WL1$class <- "Land Cover WL1"
sum_roti_sample_data_shannon_WL1 <-sum_roti_sample_data_shannon_WL1[-1,]

sub_roti_sample_data_WL2 <- a_roti_sample_data %>% filter(LU2 == "WL2")
sum_roti_sample_data_richness_WL2 <- capture.output(summary(sub_roti_sample_data_WL2$Richness),file=NULL,append=FALSE)
sum_roti_sample_data_richness_WL2 <- as.data.frame(sum_roti_sample_data_richness_WL2)
sum_roti_sample_data_richness_WL2$group <- "Rotifers"
sum_roti_sample_data_richness_WL2$diversity <- "Richness"
sum_roti_sample_data_richness_WL2$class <- "Land Cover WL2"
sum_roti_sample_data_richness_WL2 <- sum_roti_sample_data_richness_WL2[-1,]

sum_roti_sample_data_shannon_WL2 <- capture.output(summary(sub_roti_sample_data_WL2$Shannon),file=NULL,append=FALSE)
sum_roti_sample_data_shannon_WL2 <- as.data.frame(sum_roti_sample_data_shannon_WL2)
sum_roti_sample_data_shannon_WL2$group <- "Rotifers"
sum_roti_sample_data_shannon_WL2$diversity <- "Shannon"
sum_roti_sample_data_shannon_WL2$class <- "Land Cover WL2"
sum_roti_sample_data_shannon_WL2 <-sum_roti_sample_data_shannon_WL2[-1,]

sub_roti_sample_data_WL3 <- a_roti_sample_data %>% filter(LU2 == "WL3")
sum_roti_sample_data_richness_WL3 <- capture.output(summary(sub_roti_sample_data_WL3$Richness),file=NULL,append=FALSE)
sum_roti_sample_data_richness_WL3 <- as.data.frame(sum_roti_sample_data_richness_WL3)
sum_roti_sample_data_richness_WL3$group <- "Rotifers"
sum_roti_sample_data_richness_WL3$diversity <- "Richness"
sum_roti_sample_data_richness_WL3$class <- "Land Cover WL3"
sum_roti_sample_data_richness_WL3 <- sum_roti_sample_data_richness_WL3[-1,]

sum_roti_sample_data_shannon_WL3 <- capture.output(summary(sub_roti_sample_data_WL3$Shannon),file=NULL,append=FALSE)
sum_roti_sample_data_shannon_WL3 <- as.data.frame(sum_roti_sample_data_shannon_WL3)
sum_roti_sample_data_shannon_WL3$group <- "Rotifers"
sum_roti_sample_data_shannon_WL3$diversity <- "Shannon"
sum_roti_sample_data_shannon_WL3$class <- "Land Cover WL3"
sum_roti_sample_data_shannon_WL3 <-sum_roti_sample_data_shannon_WL3[-1,]

sub_roti_sample_data_WL4 <- a_roti_sample_data %>% filter(LU2 == "WL4")
sum_roti_sample_data_richness_WL4 <- capture.output(summary(sub_roti_sample_data_WL4$Richness),file=NULL,append=FALSE)
sum_roti_sample_data_richness_WL4 <- as.data.frame(sum_roti_sample_data_richness_WL4)
sum_roti_sample_data_richness_WL4$group <- "Rotifers"
sum_roti_sample_data_richness_WL4$diversity <- "Richness"
sum_roti_sample_data_richness_WL4$class <- "Land Cover WL4"
sum_roti_sample_data_richness_WL4 <- sum_roti_sample_data_richness_WL4[-1,]

sum_roti_sample_data_shannon_WL4 <- capture.output(summary(sub_roti_sample_data_WL4$Shannon),file=NULL,append=FALSE)
sum_roti_sample_data_shannon_WL4 <- as.data.frame(sum_roti_sample_data_shannon_WL4)
sum_roti_sample_data_shannon_WL4$group <- "Rotifers"
sum_roti_sample_data_shannon_WL4$diversity <- "Shannon"
sum_roti_sample_data_shannon_WL4$class <- "Land Cover WL4"
sum_roti_sample_data_shannon_WL4 <-sum_roti_sample_data_shannon_WL4[-1,]

sub_roti_sample_data_No_C <- a_roti_sample_data %>% filter(LU2 == "No_C")
sum_roti_sample_data_richness_No_C <- capture.output(summary(sub_roti_sample_data_No_C$Richness),file=NULL,append=FALSE)
sum_roti_sample_data_richness_No_C <- as.data.frame(sum_roti_sample_data_richness_No_C)
sum_roti_sample_data_richness_No_C$group <- "Rotifers"
sum_roti_sample_data_richness_No_C$diversity <- "Richness"
sum_roti_sample_data_richness_No_C$class <- "Land Cover No_C"
sum_roti_sample_data_richness_No_C <- sum_roti_sample_data_richness_No_C[-1,]

sum_roti_sample_data_shannon_No_C <- capture.output(summary(sub_roti_sample_data_No_C$Shannon),file=NULL,append=FALSE)
sum_roti_sample_data_shannon_No_C <- as.data.frame(sum_roti_sample_data_shannon_No_C)
sum_roti_sample_data_shannon_No_C$group <- "Rotifers"
sum_roti_sample_data_shannon_No_C$diversity <- "Shannon"
sum_roti_sample_data_shannon_No_C$class <- "Land Cover No_C"
sum_roti_sample_data_shannon_No_C <-sum_roti_sample_data_shannon_No_C[-1,]


sub_nema_sample_data_CL1 <- a_nema_sample_data %>% filter(LU2 == "CL1")
sum_nema_sample_data_richness_CL1 <- capture.output(summary(sub_nema_sample_data_CL1$Richness),file=NULL,append=FALSE)
sum_nema_sample_data_richness_CL1 <- as.data.frame(sum_nema_sample_data_richness_CL1)
sum_nema_sample_data_richness_CL1$group <- "Nematodes"
sum_nema_sample_data_richness_CL1$diversity <- "Richness"
sum_nema_sample_data_richness_CL1$class <- "Land Cover CL1"
sum_nema_sample_data_richness_CL1 <- sum_nema_sample_data_richness_CL1[-1,]

sum_nema_sample_data_shannon_CL1 <- capture.output(summary(sub_nema_sample_data_CL1$Shannon),file=NULL,append=FALSE)
sum_nema_sample_data_shannon_CL1 <- as.data.frame(sum_nema_sample_data_shannon_CL1)
sum_nema_sample_data_shannon_CL1$group <- "Nematodes"
sum_nema_sample_data_shannon_CL1$diversity <- "Shannon"
sum_nema_sample_data_shannon_CL1$class <- "Land Cover CL1"
sum_nema_sample_data_shannon_CL1 <-sum_nema_sample_data_shannon_CL1[-1,]

sub_nema_sample_data_CL2 <- a_nema_sample_data %>% filter(LU2 == "CL2")
sum_nema_sample_data_richness_CL2 <- capture.output(summary(sub_nema_sample_data_CL2$Richness),file=NULL,append=FALSE)
sum_nema_sample_data_richness_CL2 <- as.data.frame(sum_nema_sample_data_richness_CL2)
sum_nema_sample_data_richness_CL2$group <- "Nematodes"
sum_nema_sample_data_richness_CL2$diversity <- "Richness"
sum_nema_sample_data_richness_CL2$class <- "Land Cover CL2"
sum_nema_sample_data_richness_CL2 <- sum_nema_sample_data_richness_CL2[-1,]

sum_nema_sample_data_shannon_CL2 <- capture.output(summary(sub_nema_sample_data_CL2$Shannon),file=NULL,append=FALSE)
sum_nema_sample_data_shannon_CL2 <- as.data.frame(sum_nema_sample_data_shannon_CL2)
sum_nema_sample_data_shannon_CL2$group <- "Nematodes"
sum_nema_sample_data_shannon_CL2$diversity <- "Shannon"
sum_nema_sample_data_shannon_CL2$class <- "Land Cover CL2"
sum_nema_sample_data_shannon_CL2 <-sum_nema_sample_data_shannon_CL2[-1,]

sub_nema_sample_data_CL3 <- a_nema_sample_data %>% filter(LU2 == "CL3")
sum_nema_sample_data_richness_CL3 <- capture.output(summary(sub_nema_sample_data_CL3$Richness),file=NULL,append=FALSE)
sum_nema_sample_data_richness_CL3 <- as.data.frame(sum_nema_sample_data_richness_CL3)
sum_nema_sample_data_richness_CL3$group <- "Nematodes"
sum_nema_sample_data_richness_CL3$diversity <- "Richness"
sum_nema_sample_data_richness_CL3$class <- "Land Cover CL3"
sum_nema_sample_data_richness_CL3 <- sum_nema_sample_data_richness_CL3[-1,]

sum_nema_sample_data_shannon_CL3 <- capture.output(summary(sub_nema_sample_data_CL3$Shannon),file=NULL,append=FALSE)
sum_nema_sample_data_shannon_CL3 <- as.data.frame(sum_nema_sample_data_shannon_CL3)
sum_nema_sample_data_shannon_CL3$group <- "Nematodes"
sum_nema_sample_data_shannon_CL3$diversity <- "Shannon"
sum_nema_sample_data_shannon_CL3$class <- "Land Cover CL3"
sum_nema_sample_data_shannon_CL3 <-sum_nema_sample_data_shannon_CL3[-1,]

sub_nema_sample_data_CL4 <- a_nema_sample_data %>% filter(LU2 == "CL4")
sum_nema_sample_data_richness_CL4 <- capture.output(summary(sub_nema_sample_data_CL4$Richness),file=NULL,append=FALSE)
sum_nema_sample_data_richness_CL4 <- as.data.frame(sum_nema_sample_data_richness_CL4)
sum_nema_sample_data_richness_CL4$group <- "Nematodes"
sum_nema_sample_data_richness_CL4$diversity <- "Richness"
sum_nema_sample_data_richness_CL4$class <- "Land Cover CL4"

sum_nema_sample_data_shannon_CL4 <- capture.output(summary(sub_nema_sample_data_CL4$Shannon),file=NULL,append=FALSE)
sum_nema_sample_data_shannon_CL4 <- as.data.frame(sum_nema_sample_data_shannon_CL4)
sum_nema_sample_data_shannon_CL4$group <- "Nematodes"
sum_nema_sample_data_shannon_CL4$diversity <- "Shannon"
sum_nema_sample_data_shannon_CL4$class <- "Land Cover CL4"
sum_nema_sample_data_shannon_CL4 <-sum_nema_sample_data_shannon_CL4[-1,]
sum_nema_sample_data_richness_CL4 <- sum_nema_sample_data_richness_CL4[-1,]

sub_nema_sample_data_GL1 <- a_nema_sample_data %>% filter(LU2 == "GL1")
sum_nema_sample_data_richness_GL1 <- capture.output(summary(sub_nema_sample_data_GL1$Richness),file=NULL,append=FALSE)
sum_nema_sample_data_richness_GL1 <- as.data.frame(sum_nema_sample_data_richness_GL1)
sum_nema_sample_data_richness_GL1$group <- "Nematodes"
sum_nema_sample_data_richness_GL1$diversity <- "Richness"
sum_nema_sample_data_richness_GL1$class <- "Land Cover GL1"
sum_nema_sample_data_richness_GL1 <- sum_nema_sample_data_richness_GL1[-1,]

sum_nema_sample_data_shannon_GL1 <- capture.output(summary(sub_nema_sample_data_GL1$Shannon),file=NULL,append=FALSE)
sum_nema_sample_data_shannon_GL1 <- as.data.frame(sum_nema_sample_data_shannon_GL1)
sum_nema_sample_data_shannon_GL1$group <- "Nematodes"
sum_nema_sample_data_shannon_GL1$diversity <- "Shannon"
sum_nema_sample_data_shannon_GL1$class <- "Land Cover GL1"
sum_nema_sample_data_shannon_GL1 <-sum_nema_sample_data_shannon_GL1[-1,]

sub_nema_sample_data_GL2 <- a_nema_sample_data %>% filter(LU2 == "GL2")
sum_nema_sample_data_richness_GL2 <- capture.output(summary(sub_nema_sample_data_GL2$Richness),file=NULL,append=FALSE)
sum_nema_sample_data_richness_GL2 <- as.data.frame(sum_nema_sample_data_richness_GL2)
sum_nema_sample_data_richness_GL2$group <- "Nematodes"
sum_nema_sample_data_richness_GL2$diversity <- "Richness"
sum_nema_sample_data_richness_GL2$class <- "Land Cover GL2"
sum_nema_sample_data_richness_GL2 <- sum_nema_sample_data_richness_GL2[-1,]

sum_nema_sample_data_shannon_GL2 <- capture.output(summary(sub_nema_sample_data_GL2$Shannon),file=NULL,append=FALSE)
sum_nema_sample_data_shannon_GL2 <- as.data.frame(sum_nema_sample_data_shannon_GL2)
sum_nema_sample_data_shannon_GL2$group <- "Nematodes"
sum_nema_sample_data_shannon_GL2$diversity <- "Shannon"
sum_nema_sample_data_shannon_GL2$class <- "Land Cover GL2"
sum_nema_sample_data_shannon_GL2 <-sum_nema_sample_data_shannon_GL2[-1,]

sub_nema_sample_data_GL3 <- a_nema_sample_data %>% filter(LU2 == "GL3")
sum_nema_sample_data_richness_GL3 <- capture.output(summary(sub_nema_sample_data_GL3$Richness),file=NULL,append=FALSE)
sum_nema_sample_data_richness_GL3 <- as.data.frame(sum_nema_sample_data_richness_GL3)
sum_nema_sample_data_richness_GL3$group <- "Nematodes"
sum_nema_sample_data_richness_GL3$diversity <- "Richness"
sum_nema_sample_data_richness_GL3$class <- "Land Cover GL3"
sum_nema_sample_data_richness_GL3 <- sum_nema_sample_data_richness_GL3[-1,]

sum_nema_sample_data_shannon_GL3 <- capture.output(summary(sub_nema_sample_data_GL3$Shannon),file=NULL,append=FALSE)
sum_nema_sample_data_shannon_GL3 <- as.data.frame(sum_nema_sample_data_shannon_GL3)
sum_nema_sample_data_shannon_GL3$group <- "Nematodes"
sum_nema_sample_data_shannon_GL3$diversity <- "Shannon"
sum_nema_sample_data_shannon_GL3$class <- "Land Cover GL3"
sum_nema_sample_data_shannon_GL3 <-sum_nema_sample_data_shannon_GL3[-1,]

sub_nema_sample_data_WL1 <- a_nema_sample_data %>% filter(LU2 == "WL1")
sum_nema_sample_data_richness_WL1 <- capture.output(summary(sub_nema_sample_data_WL1$Richness),file=NULL,append=FALSE)
sum_nema_sample_data_richness_WL1 <- as.data.frame(sum_nema_sample_data_richness_WL1)
sum_nema_sample_data_richness_WL1$group <- "Nematodes"
sum_nema_sample_data_richness_WL1$diversity <- "Richness"
sum_nema_sample_data_richness_WL1$class <- "Land Cover WL1"
sum_nema_sample_data_richness_WL1 <- sum_nema_sample_data_richness_WL1[-1,]

sum_nema_sample_data_shannon_WL1 <- capture.output(summary(sub_nema_sample_data_WL1$Shannon),file=NULL,append=FALSE)
sum_nema_sample_data_shannon_WL1 <- as.data.frame(sum_nema_sample_data_shannon_WL1)
sum_nema_sample_data_shannon_WL1$group <- "Nematodes"
sum_nema_sample_data_shannon_WL1$diversity <- "Shannon"
sum_nema_sample_data_shannon_WL1$class <- "Land Cover WL1"
sum_nema_sample_data_shannon_WL1 <-sum_nema_sample_data_shannon_WL1[-1,]

sub_nema_sample_data_WL2 <- a_nema_sample_data %>% filter(LU2 == "WL2")
sum_nema_sample_data_richness_WL2 <- capture.output(summary(sub_nema_sample_data_WL2$Richness),file=NULL,append=FALSE)
sum_nema_sample_data_richness_WL2 <- as.data.frame(sum_nema_sample_data_richness_WL2)
sum_nema_sample_data_richness_WL2$group <- "Nematodes"
sum_nema_sample_data_richness_WL2$diversity <- "Richness"
sum_nema_sample_data_richness_WL2$class <- "Land Cover WL2"
sum_nema_sample_data_richness_WL2 <- sum_nema_sample_data_richness_WL2[-1,]

sum_nema_sample_data_shannon_WL2 <- capture.output(summary(sub_nema_sample_data_WL2$Shannon),file=NULL,append=FALSE)
sum_nema_sample_data_shannon_WL2 <- as.data.frame(sum_nema_sample_data_shannon_WL2)
sum_nema_sample_data_shannon_WL2$group <- "Nematodes"
sum_nema_sample_data_shannon_WL2$diversity <- "Shannon"
sum_nema_sample_data_shannon_WL2$class <- "Land Cover WL2"
sum_nema_sample_data_shannon_WL2 <-sum_nema_sample_data_shannon_WL2[-1,]

sub_nema_sample_data_WL3 <- a_nema_sample_data %>% filter(LU2 == "WL3")
sum_nema_sample_data_richness_WL3 <- capture.output(summary(sub_nema_sample_data_WL3$Richness),file=NULL,append=FALSE)
sum_nema_sample_data_richness_WL3 <- as.data.frame(sum_nema_sample_data_richness_WL3)
sum_nema_sample_data_richness_WL3$group <- "Nematodes"
sum_nema_sample_data_richness_WL3$diversity <- "Richness"
sum_nema_sample_data_richness_WL3$class <- "Land Cover WL3"
sum_nema_sample_data_richness_WL3 <- sum_nema_sample_data_richness_WL3[-1,]

sum_nema_sample_data_shannon_WL3 <- capture.output(summary(sub_nema_sample_data_WL3$Shannon),file=NULL,append=FALSE)
sum_nema_sample_data_shannon_WL3 <- as.data.frame(sum_nema_sample_data_shannon_WL3)
sum_nema_sample_data_shannon_WL3$group <- "Nematodes"
sum_nema_sample_data_shannon_WL3$diversity <- "Shannon"
sum_nema_sample_data_shannon_WL3$class <- "Land Cover WL3"
sum_nema_sample_data_shannon_WL3 <-sum_nema_sample_data_shannon_WL3[-1,]

sub_nema_sample_data_WL4 <- a_nema_sample_data %>% filter(LU2 == "WL4")
sum_nema_sample_data_richness_WL4 <- capture.output(summary(sub_nema_sample_data_WL4$Richness),file=NULL,append=FALSE)
sum_nema_sample_data_richness_WL4 <- as.data.frame(sum_nema_sample_data_richness_WL4)
sum_nema_sample_data_richness_WL4$group <- "Nematodes"
sum_nema_sample_data_richness_WL4$diversity <- "Richness"
sum_nema_sample_data_richness_WL4$class <- "Land Cover WL4"
sum_nema_sample_data_richness_WL4 <- sum_nema_sample_data_richness_WL4[-1,]

sum_nema_sample_data_shannon_WL4 <- capture.output(summary(sub_nema_sample_data_WL4$Shannon),file=NULL,append=FALSE)
sum_nema_sample_data_shannon_WL4 <- as.data.frame(sum_nema_sample_data_shannon_WL4)
sum_nema_sample_data_shannon_WL4$group <- "Nematodes"
sum_nema_sample_data_shannon_WL4$diversity <- "Shannon"
sum_nema_sample_data_shannon_WL4$class <- "Land Cover WL4"
sum_nema_sample_data_shannon_WL4 <-sum_nema_sample_data_shannon_WL4[-1,]

sub_nema_sample_data_No_C <- a_nema_sample_data %>% filter(LU2 == "No_C")
sum_nema_sample_data_richness_No_C <- capture.output(summary(sub_nema_sample_data_No_C$Richness),file=NULL,append=FALSE)
sum_nema_sample_data_richness_No_C <- as.data.frame(sum_nema_sample_data_richness_No_C)
sum_nema_sample_data_richness_No_C$group <- "Nematodes"
sum_nema_sample_data_richness_No_C$diversity <- "Richness"
sum_nema_sample_data_richness_No_C$class <- "Land Cover No_C"
sum_nema_sample_data_richness_No_C <- sum_nema_sample_data_richness_No_C[-1,]

sum_nema_sample_data_shannon_No_C <- capture.output(summary(sub_nema_sample_data_No_C$Shannon),file=NULL,append=FALSE)
sum_nema_sample_data_shannon_No_C <- as.data.frame(sum_nema_sample_data_shannon_No_C)
sum_nema_sample_data_shannon_No_C$group <- "Nematodes"
sum_nema_sample_data_shannon_No_C$diversity <- "Shannon"
sum_nema_sample_data_shannon_No_C$class <- "Land Cover No_C"
sum_nema_sample_data_shannon_No_C <-sum_nema_sample_data_shannon_No_C[-1,]

sub_tardi_sample_data_CL1 <- a_tardi_sample_data %>% filter(LU2 == "CL1")
sum_tardi_sample_data_richness_CL1 <- capture.output(summary(sub_tardi_sample_data_CL1$Richness),file=NULL,append=FALSE)
sum_tardi_sample_data_richness_CL1 <- as.data.frame(sum_tardi_sample_data_richness_CL1)
sum_tardi_sample_data_richness_CL1$group <- "Tardigrades"
sum_tardi_sample_data_richness_CL1$diversity <- "Richness"
sum_tardi_sample_data_richness_CL1$class <- "Land Cover CL1"
sum_tardi_sample_data_richness_CL1 <- sum_tardi_sample_data_richness_CL1[-1,]

sum_tardi_sample_data_shannon_CL1 <- capture.output(summary(sub_tardi_sample_data_CL1$Shannon),file=NULL,append=FALSE)
sum_tardi_sample_data_shannon_CL1 <- as.data.frame(sum_tardi_sample_data_shannon_CL1)
sum_tardi_sample_data_shannon_CL1$group <- "Tardigrades"
sum_tardi_sample_data_shannon_CL1$diversity <- "Shannon"
sum_tardi_sample_data_shannon_CL1$class <- "Land Cover CL1"
sum_tardi_sample_data_shannon_CL1 <-sum_tardi_sample_data_shannon_CL1[-1,]

sub_tardi_sample_data_CL2 <- a_tardi_sample_data %>% filter(LU2 == "CL2")
sum_tardi_sample_data_richness_CL2 <- capture.output(summary(sub_tardi_sample_data_CL2$Richness),file=NULL,append=FALSE)
sum_tardi_sample_data_richness_CL2 <- as.data.frame(sum_tardi_sample_data_richness_CL2)
sum_tardi_sample_data_richness_CL2$group <- "Tardigrades"
sum_tardi_sample_data_richness_CL2$diversity <- "Richness"
sum_tardi_sample_data_richness_CL2$class <- "Land Cover CL2"
sum_tardi_sample_data_richness_CL2 <- sum_tardi_sample_data_richness_CL2[-1,]

sum_tardi_sample_data_shannon_CL2 <- capture.output(summary(sub_tardi_sample_data_CL2$Shannon),file=NULL,append=FALSE)
sum_tardi_sample_data_shannon_CL2 <- as.data.frame(sum_tardi_sample_data_shannon_CL2)
sum_tardi_sample_data_shannon_CL2$group <- "Tardigrades"
sum_tardi_sample_data_shannon_CL2$diversity <- "Shannon"
sum_tardi_sample_data_shannon_CL2$class <- "Land Cover CL2"
sum_tardi_sample_data_shannon_CL2 <-sum_tardi_sample_data_shannon_CL2[-1,]

sub_tardi_sample_data_CL3 <- a_tardi_sample_data %>% filter(LU2 == "CL3")
sum_tardi_sample_data_richness_CL3 <- capture.output(summary(sub_tardi_sample_data_CL3$Richness),file=NULL,append=FALSE)
sum_tardi_sample_data_richness_CL3 <- as.data.frame(sum_tardi_sample_data_richness_CL3)
sum_tardi_sample_data_richness_CL3$group <- "Tardigrades"
sum_tardi_sample_data_richness_CL3$diversity <- "Richness"
sum_tardi_sample_data_richness_CL3$class <- "Land Cover CL3"
sum_tardi_sample_data_richness_CL3 <- sum_tardi_sample_data_richness_CL3[-1,]

sum_tardi_sample_data_shannon_CL3 <- capture.output(summary(sub_tardi_sample_data_CL3$Shannon),file=NULL,append=FALSE)
sum_tardi_sample_data_shannon_CL3 <- as.data.frame(sum_tardi_sample_data_shannon_CL3)
sum_tardi_sample_data_shannon_CL3$group <- "Tardigrades"
sum_tardi_sample_data_shannon_CL3$diversity <- "Shannon"
sum_tardi_sample_data_shannon_CL3$class <- "Land Cover CL3"
sum_tardi_sample_data_shannon_CL3 <-sum_tardi_sample_data_shannon_CL3[-1,]

sub_tardi_sample_data_CL4 <- a_tardi_sample_data %>% filter(LU2 == "CL4")
sum_tardi_sample_data_richness_CL4 <- capture.output(summary(sub_tardi_sample_data_CL4$Richness),file=NULL,append=FALSE)
sum_tardi_sample_data_richness_CL4 <- as.data.frame(sum_tardi_sample_data_richness_CL4)
sum_tardi_sample_data_richness_CL4$group <- "Tardigrades"
sum_tardi_sample_data_richness_CL4$diversity <- "Richness"
sum_tardi_sample_data_richness_CL4$class <- "Land Cover CL4"
sum_tardi_sample_data_richness_CL4 <- sum_tardi_sample_data_richness_CL4[-1,]

sum_tardi_sample_data_shannon_CL4 <- capture.output(summary(sub_tardi_sample_data_CL4$Shannon),file=NULL,append=FALSE)
sum_tardi_sample_data_shannon_CL4 <- as.data.frame(sum_tardi_sample_data_shannon_CL4)
sum_tardi_sample_data_shannon_CL4$group <- "Tardigrades"
sum_tardi_sample_data_shannon_CL4$diversity <- "Shannon"
sum_tardi_sample_data_shannon_CL4$class <- "Land Cover CL4"
sum_tardi_sample_data_shannon_CL4 <-sum_tardi_sample_data_shannon_CL4[-1,]

sub_tardi_sample_data_GL1 <- a_tardi_sample_data %>% filter(LU2 == "GL1")
sum_tardi_sample_data_richness_GL1 <- capture.output(summary(sub_tardi_sample_data_GL1$Richness),file=NULL,append=FALSE)
sum_tardi_sample_data_richness_GL1 <- as.data.frame(sum_tardi_sample_data_richness_GL1)
sum_tardi_sample_data_richness_GL1$group <- "Tardigrades"
sum_tardi_sample_data_richness_GL1$diversity <- "Richness"
sum_tardi_sample_data_richness_GL1$class <- "Land Cover GL1"
sum_tardi_sample_data_richness_GL1 <- sum_tardi_sample_data_richness_GL1[-1,]

sum_tardi_sample_data_shannon_GL1 <- capture.output(summary(sub_tardi_sample_data_GL1$Shannon),file=NULL,append=FALSE)
sum_tardi_sample_data_shannon_GL1 <- as.data.frame(sum_tardi_sample_data_shannon_GL1)
sum_tardi_sample_data_shannon_GL1$group <- "Tardigrades"
sum_tardi_sample_data_shannon_GL1$diversity <- "Shannon"
sum_tardi_sample_data_shannon_GL1$class <- "Land Cover GL1"
sum_tardi_sample_data_shannon_GL1 <-sum_tardi_sample_data_shannon_GL1[-1,]

sub_tardi_sample_data_GL2 <- a_tardi_sample_data %>% filter(LU2 == "GL2")
sum_tardi_sample_data_richness_GL2 <- capture.output(summary(sub_tardi_sample_data_GL2$Richness),file=NULL,append=FALSE)
sum_tardi_sample_data_richness_GL2 <- as.data.frame(sum_tardi_sample_data_richness_GL2)
sum_tardi_sample_data_richness_GL2$group <- "Tardigrades"
sum_tardi_sample_data_richness_GL2$diversity <- "Richness"
sum_tardi_sample_data_richness_GL2$class <- "Land Cover GL2"
sum_tardi_sample_data_richness_GL2 <- sum_tardi_sample_data_richness_GL2[-1,]

sum_tardi_sample_data_shannon_GL2 <- capture.output(summary(sub_tardi_sample_data_GL2$Shannon),file=NULL,append=FALSE)
sum_tardi_sample_data_shannon_GL2 <- as.data.frame(sum_tardi_sample_data_shannon_GL2)
sum_tardi_sample_data_shannon_GL2$group <- "Tardigrades"
sum_tardi_sample_data_shannon_GL2$diversity <- "Shannon"
sum_tardi_sample_data_shannon_GL2$class <- "Land Cover GL2"
sum_tardi_sample_data_shannon_GL2 <-sum_tardi_sample_data_shannon_GL2[-1,]

sub_tardi_sample_data_GL3 <- a_tardi_sample_data %>% filter(LU2 == "GL3")
sum_tardi_sample_data_richness_GL3 <- capture.output(summary(sub_tardi_sample_data_GL3$Richness),file=NULL,append=FALSE)
sum_tardi_sample_data_richness_GL3 <- as.data.frame(sum_tardi_sample_data_richness_GL3)
sum_tardi_sample_data_richness_GL3$group <- "Tardigrades"
sum_tardi_sample_data_richness_GL3$diversity <- "Richness"
sum_tardi_sample_data_richness_GL3$class <- "Land Cover GL3"
sum_tardi_sample_data_richness_GL3 <- sum_tardi_sample_data_richness_GL3[-1,]

sum_tardi_sample_data_shannon_GL3 <- capture.output(summary(sub_tardi_sample_data_GL3$Shannon),file=NULL,append=FALSE)
sum_tardi_sample_data_shannon_GL3 <- as.data.frame(sum_tardi_sample_data_shannon_GL3)
sum_tardi_sample_data_shannon_GL3$group <- "Tardigrades"
sum_tardi_sample_data_shannon_GL3$diversity <- "Shannon"
sum_tardi_sample_data_shannon_GL3$class <- "Land Cover GL3"
sum_tardi_sample_data_shannon_GL3 <-sum_tardi_sample_data_shannon_GL3[-1,]

sub_tardi_sample_data_WL1 <- a_tardi_sample_data %>% filter(LU2 == "WL1")
sum_tardi_sample_data_richness_WL1 <- capture.output(summary(sub_tardi_sample_data_WL1$Richness),file=NULL,append=FALSE)
sum_tardi_sample_data_richness_WL1 <- as.data.frame(sum_tardi_sample_data_richness_WL1)
sum_tardi_sample_data_richness_WL1$group <- "Tardigrades"
sum_tardi_sample_data_richness_WL1$diversity <- "Richness"
sum_tardi_sample_data_richness_WL1$class <- "Land Cover WL1"
sum_tardi_sample_data_richness_WL1 <- sum_tardi_sample_data_richness_WL1[-1,]

sum_tardi_sample_data_shannon_WL1 <- capture.output(summary(sub_tardi_sample_data_WL1$Shannon),file=NULL,append=FALSE)
sum_tardi_sample_data_shannon_WL1 <- as.data.frame(sum_tardi_sample_data_shannon_WL1)
sum_tardi_sample_data_shannon_WL1$group <- "Tardigrades"
sum_tardi_sample_data_shannon_WL1$diversity <- "Shannon"
sum_tardi_sample_data_shannon_WL1$class <- "Land Cover WL1"
sum_tardi_sample_data_shannon_WL1 <-sum_tardi_sample_data_shannon_WL1[-1,]

sub_tardi_sample_data_WL2 <- a_tardi_sample_data %>% filter(LU2 == "WL2")
sum_tardi_sample_data_richness_WL2 <- capture.output(summary(sub_tardi_sample_data_WL2$Richness),file=NULL,append=FALSE)
sum_tardi_sample_data_richness_WL2 <- as.data.frame(sum_tardi_sample_data_richness_WL2)
sum_tardi_sample_data_richness_WL2$group <- "Tardigrades"
sum_tardi_sample_data_richness_WL2$diversity <- "Richness"
sum_tardi_sample_data_richness_WL2$class <- "Land Cover WL2"
sum_tardi_sample_data_richness_WL2 <- sum_tardi_sample_data_richness_WL2[-1,]

sum_tardi_sample_data_shannon_WL2 <- capture.output(summary(sub_tardi_sample_data_WL2$Shannon),file=NULL,append=FALSE)
sum_tardi_sample_data_shannon_WL2 <- as.data.frame(sum_tardi_sample_data_shannon_WL2)
sum_tardi_sample_data_shannon_WL2$group <- "Tardigrades"
sum_tardi_sample_data_shannon_WL2$diversity <- "Shannon"
sum_tardi_sample_data_shannon_WL2$class <- "Land Cover WL2"
sum_tardi_sample_data_shannon_WL2 <-sum_tardi_sample_data_shannon_WL2[-1,]

sub_tardi_sample_data_WL3 <- a_tardi_sample_data %>% filter(LU2 == "WL3")
sum_tardi_sample_data_richness_WL3 <- capture.output(summary(sub_tardi_sample_data_WL3$Richness),file=NULL,append=FALSE)
sum_tardi_sample_data_richness_WL3 <- as.data.frame(sum_tardi_sample_data_richness_WL3)
sum_tardi_sample_data_richness_WL3$group <- "Tardigrades"
sum_tardi_sample_data_richness_WL3$diversity <- "Richness"
sum_tardi_sample_data_richness_WL3$class <- "Land Cover WL3"
sum_tardi_sample_data_richness_WL3 <- sum_tardi_sample_data_richness_WL3[-1,]

sum_tardi_sample_data_shannon_WL3 <- capture.output(summary(sub_tardi_sample_data_WL3$Shannon),file=NULL,append=FALSE)
sum_tardi_sample_data_shannon_WL3 <- as.data.frame(sum_tardi_sample_data_shannon_WL3)
sum_tardi_sample_data_shannon_WL3$group <- "Tardigrades"
sum_tardi_sample_data_shannon_WL3$diversity <- "Shannon"
sum_tardi_sample_data_shannon_WL3$class <- "Land Cover WL3"
sum_tardi_sample_data_shannon_WL3 <-sum_tardi_sample_data_shannon_WL3[-1,]

sub_tardi_sample_data_WL4 <- a_tardi_sample_data %>% filter(LU2 == "WL4")
sum_tardi_sample_data_richness_WL4 <- capture.output(summary(sub_tardi_sample_data_WL4$Richness),file=NULL,append=FALSE)
sum_tardi_sample_data_richness_WL4 <- as.data.frame(sum_tardi_sample_data_richness_WL4)
sum_tardi_sample_data_richness_WL4$group <- "Tardigrades"
sum_tardi_sample_data_richness_WL4$diversity <- "Richness"
sum_tardi_sample_data_richness_WL4$class <- "Land Cover WL4"
sum_tardi_sample_data_richness_WL4 <- sum_tardi_sample_data_richness_WL4[-1,]

sum_tardi_sample_data_shannon_WL4 <- capture.output(summary(sub_tardi_sample_data_WL4$Shannon),file=NULL,append=FALSE)
sum_tardi_sample_data_shannon_WL4 <- as.data.frame(sum_tardi_sample_data_shannon_WL4)
sum_tardi_sample_data_shannon_WL4$group <- "Tardigrades"
sum_tardi_sample_data_shannon_WL4$diversity <- "Shannon"
sum_tardi_sample_data_shannon_WL4$class <- "Land Cover WL4"
sum_tardi_sample_data_shannon_WL4 <-sum_tardi_sample_data_shannon_WL4[-1,]

sub_tardi_sample_data_No_C <- a_tardi_sample_data %>% filter(LU2 == "No_C")
sum_tardi_sample_data_richness_No_C <- capture.output(summary(sub_tardi_sample_data_No_C$Richness),file=NULL,append=FALSE)
sum_tardi_sample_data_richness_No_C <- as.data.frame(sum_tardi_sample_data_richness_No_C)
sum_tardi_sample_data_richness_No_C$group <- "Tardigrades"
sum_tardi_sample_data_richness_No_C$diversity <- "Richness"
sum_tardi_sample_data_richness_No_C$class <- "Land Cover No_C"
sum_tardi_sample_data_richness_No_C <- sum_tardi_sample_data_richness_No_C[-1,]

sum_tardi_sample_data_shannon_No_C <- capture.output(summary(sub_tardi_sample_data_No_C$Shannon),file=NULL,append=FALSE)
sum_tardi_sample_data_shannon_No_C <- as.data.frame(sum_tardi_sample_data_shannon_No_C)
sum_tardi_sample_data_shannon_No_C$group <- "Tardigrades"
sum_tardi_sample_data_shannon_No_C$diversity <- "Shannon"
sum_tardi_sample_data_shannon_No_C$class <- "Land Cover No_C"
sum_tardi_sample_data_shannon_No_C <-sum_tardi_sample_data_shannon_No_C[-1,]

sub_a_nema_sample_data_erosion1 <- a_nema_sample_data %>% filter(Erosion_grouped == "Low")
sum_a_nema_richness_erosion1 <- capture.output(summary(sub_a_nema_sample_data_erosion1$Richness),file=NULL,append=FALSE)
sum_a_nema_richness_erosion1 <- as.data.frame(sum_a_nema_richness_erosion1)
sum_a_nema_richness_erosion1$group <- "Nematodes"
sum_a_nema_richness_erosion1$diversity <- "Richness"
sum_a_nema_richness_erosion1$class <- "Erosion 0-1"
sum_a_nema_richness_erosion1 <-sum_a_nema_richness_erosion1[-1,]

sum_a_nema_shannon_erosion1 <- capture.output(summary(sub_a_nema_sample_data_erosion1$Shannon),file=NULL,append=FALSE)
sum_a_nema_shannon_erosion1 <- as.data.frame(sum_a_nema_shannon_erosion1)
sum_a_nema_shannon_erosion1$group <- "Nematodes"
sum_a_nema_shannon_erosion1$diversity <- "Shannon"
sum_a_nema_shannon_erosion1$class <- "Erosion 0-1"
sum_a_nema_shannon_erosion1 <-sum_a_nema_shannon_erosion1[-1,]

sub_a_nema_sample_data_erosion2 <- a_nema_sample_data %>% filter(Erosion_grouped == "Medium high")
sum_a_nema_richness_erosion2 <- capture.output(summary(sub_a_nema_sample_data_erosion2$Richness),file=NULL,append=FALSE)
sum_a_nema_richness_erosion2 <- as.data.frame(sum_a_nema_richness_erosion2)
sum_a_nema_richness_erosion2$group <- "Nematodes"
sum_a_nema_richness_erosion2$diversity <- "Richness"
sum_a_nema_richness_erosion2$class <- "Erosion 1-5"
sum_a_nema_richness_erosion2 <-sum_a_nema_richness_erosion2[-1,]

sum_a_nema_shannon_erosion2 <- capture.output(summary(sub_a_nema_sample_data_erosion2$Shannon),file=NULL,append=FALSE)
sum_a_nema_shannon_erosion2 <- as.data.frame(sum_a_nema_shannon_erosion2)
sum_a_nema_shannon_erosion2$group <- "Nematodes"
sum_a_nema_shannon_erosion2$diversity <- "Shannon"
sum_a_nema_shannon_erosion2$class <- "Erosion 1-5"
sum_a_nema_shannon_erosion2 <-sum_a_nema_shannon_erosion2[-1,]

sub_a_nema_sample_data_erosion3 <- a_nema_sample_data %>% filter(Erosion_grouped == "High")
sum_a_nema_richness_erosion3 <- capture.output(summary(sub_a_nema_sample_data_erosion3$Richness),file=NULL,append=FALSE)
sum_a_nema_richness_erosion3 <- as.data.frame(sum_a_nema_richness_erosion3)
sum_a_nema_richness_erosion3$group <- "Nematodes"
sum_a_nema_richness_erosion3$diversity <- "Richness"
sum_a_nema_richness_erosion3$class <- "Erosion 5-10"
sum_a_nema_richness_erosion3 <-sum_a_nema_richness_erosion3[-1,]

sum_a_nema_sample_data_shannon_erosion3 <- capture.output(summary(sub_a_nema_sample_data_erosion3$Shannon),file=NULL,append=FALSE)
sum_a_nema_shannon_erosion3 <- as.data.frame(sum_a_nema_sample_data_shannon_erosion3)
sum_a_nema_shannon_erosion3$group <- "Nematodes"
sum_a_nema_shannon_erosion3$diversity <- "Shannon"
sum_a_nema_shannon_erosion3$class <- "Erosion 5-10"
sum_a_nema_shannon_erosion3 <-sum_a_nema_shannon_erosion3[-1,]

sub_a_nema_sample_data_erosion4 <- a_nema_sample_data %>% filter(Erosion_grouped == "Very high")
sum_a_nema_richness_erosion4 <- capture.output(summary(sub_a_nema_sample_data_erosion4$Richness),file=NULL,append=FALSE)
sum_a_nema_richness_erosion4 <- as.data.frame(sum_a_nema_richness_erosion4)
sum_a_nema_richness_erosion4$group <- "Nematodes"
sum_a_nema_richness_erosion4$diversity <- "Richness"
sum_a_nema_richness_erosion4$class <- "Erosion 10+"
sum_a_nema_richness_erosion4 <-sum_a_nema_richness_erosion4[-1,]

sum_a_nema_shannon_erosion4 <- capture.output(summary(sub_a_nema_sample_data_erosion4$Shannon),file=NULL,append=FALSE)
sum_a_nema_shannon_erosion4 <- as.data.frame(sum_a_nema_shannon_erosion4)
sum_a_nema_shannon_erosion4$group <- "Nematodes"
sum_a_nema_shannon_erosion4$diversity <- "Shannon"
sum_a_nema_shannon_erosion4$class <- "Erosion 10+"
sum_a_nema_shannon_erosion4 <-sum_a_nema_shannon_erosion4[-1,]

sub_a_tardi_sample_data_erosion1 <- a_tardi_sample_data %>% filter(Erosion_grouped == "Low")
sum_a_tardi_richness_erosion1 <- capture.output(summary(sub_a_tardi_sample_data_erosion1$Richness),file=NULL,append=FALSE)
sum_a_tardi_richness_erosion1 <- as.data.frame(sum_a_tardi_richness_erosion1)
sum_a_tardi_richness_erosion1$group <- "Tardigrades"
sum_a_tardi_richness_erosion1$diversity <- "Richness"
sum_a_tardi_richness_erosion1$class <- "Erosion 0-1"
sum_a_tardi_richness_erosion1 <-sum_a_tardi_richness_erosion1[-1,]

sum_a_tardi_shannon_erosion1 <- capture.output(summary(sub_a_tardi_sample_data_erosion1$Shannon),file=NULL,append=FALSE)
sum_a_tardi_shannon_erosion1 <- as.data.frame(sum_a_tardi_shannon_erosion1)
sum_a_tardi_shannon_erosion1$group <- "Tardigrades"
sum_a_tardi_shannon_erosion1$diversity <- "Shannon"
sum_a_tardi_shannon_erosion1$class <- "Erosion 0-1"
sum_a_tardi_shannon_erosion1 <-sum_a_tardi_shannon_erosion1[-1,]

sub_a_tardi_sample_data_erosion2 <- a_tardi_sample_data %>% filter(Erosion_grouped == "Medium high")
sum_a_tardi_richness_erosion2 <- capture.output(summary(sub_a_tardi_sample_data_erosion2$Richness),file=NULL,append=FALSE)
sum_a_tardi_richness_erosion2 <- as.data.frame(sum_a_tardi_richness_erosion2)
sum_a_tardi_richness_erosion2$group <- "Tardigrades"
sum_a_tardi_richness_erosion2$diversity <- "Richness"
sum_a_tardi_richness_erosion2$class <- "Erosion 1-5"
sum_a_tardi_richness_erosion2 <-sum_a_tardi_richness_erosion2[-1,]

sum_a_tardi_shannon_erosion2 <- capture.output(summary(sub_a_tardi_sample_data_erosion2$Shannon),file=NULL,append=FALSE)
sum_a_tardi_shannon_erosion2 <- as.data.frame(sum_a_tardi_shannon_erosion2)
sum_a_tardi_shannon_erosion2$group <- "Tardigrades"
sum_a_tardi_shannon_erosion2$diversity <- "Shannon"
sum_a_tardi_shannon_erosion2$class <- "Erosion 1-5"
sum_a_tardi_shannon_erosion2 <-sum_a_tardi_shannon_erosion2[-1,]

sub_a_tardi_sample_data_erosion3 <- a_tardi_sample_data %>% filter(Erosion_grouped == "High")
sum_a_tardi_richness_erosion3 <-capture.output(summary(sub_a_tardi_sample_data_erosion3$Richness),file=NULL,append=FALSE)
sum_a_tardi_richness_erosion3 <- as.data.frame(sum_a_tardi_richness_erosion3)
sum_a_tardi_richness_erosion3$group <- "Tardigrades"
sum_a_tardi_richness_erosion3$diversity <- "Richness"
sum_a_tardi_richness_erosion3$class <- "Erosion 5-10"
sum_a_tardi_richness_erosion3 <-sum_a_tardi_richness_erosion3[-1,]

sum_a_tardi_shannon_erosion3 <- capture.output(summary(sub_a_tardi_sample_data_erosion3$Shannon),file=NULL,append=FALSE)
sum_a_tardi_shannon_erosion3 <- as.data.frame(sum_a_tardi_shannon_erosion3)
sum_a_tardi_shannon_erosion3$group <- "Tardigrades"
sum_a_tardi_shannon_erosion3$diversity <- "Shannon"
sum_a_tardi_shannon_erosion3$class <- "Erosion 5-10"
sum_a_tardi_shannon_erosion3 <-sum_a_tardi_shannon_erosion3[-1,]

sub_a_tardi_sample_data_erosion4 <- a_tardi_sample_data %>% filter(Erosion_grouped == "Very high")
sum_a_tardi_richness_erosion4 <-capture.output(summary(sub_a_tardi_sample_data_erosion4$Richness),file=NULL,append=FALSE)
sum_a_tardi_richness_erosion4 <- as.data.frame(sum_a_tardi_richness_erosion4)
sum_a_tardi_richness_erosion4$group <- "Tardigrades"
sum_a_tardi_richness_erosion4$diversity <- "Richness"
sum_a_tardi_richness_erosion4$class <- "Erosion 10+"
sum_a_tardi_richness_erosion4 <-sum_a_tardi_richness_erosion4[-1,]

sum_a_tardi_shannon_erosion4 <- capture.output(summary(sub_a_tardi_sample_data_erosion4$Shannon),file=NULL,append=FALSE)
sum_a_tardi_shannon_erosion4 <- as.data.frame(sum_a_tardi_shannon_erosion4)
sum_a_tardi_shannon_erosion4$group <- "Tardigrades"
sum_a_tardi_shannon_erosion4$diversity <- "Shannon"
sum_a_tardi_shannon_erosion4$class <- "Erosion 10+"
sum_a_tardi_shannon_erosion4 <-sum_a_tardi_shannon_erosion4[-1,]

sub_a_roti_sample_data_erosion1 <- a_roti_sample_data %>% filter(Erosion_grouped == "Low")
sum_a_roti_richness_erosion1 <- capture.output(summary(sub_a_roti_sample_data_erosion1$Richness),file=NULL,append=FALSE)
sum_a_roti_richness_erosion1 <- as.data.frame(sum_a_roti_richness_erosion1)
sum_a_roti_richness_erosion1$group <- "Rotifers"
sum_a_roti_richness_erosion1$diversity <- "Richness"
sum_a_roti_richness_erosion1$class <- "Erosion 0-1"
sum_a_roti_richness_erosion1 <-sum_a_roti_richness_erosion1[-1,]

sum_a_roti_shannon_erosion1 <- capture.output(summary(sub_a_roti_sample_data_erosion1$Shannon),file=NULL,append=FALSE)
sum_a_roti_shannon_erosion1 <- as.data.frame(sum_a_roti_shannon_erosion1)
sum_a_roti_shannon_erosion1$group <- "Rotifers"
sum_a_roti_shannon_erosion1$diversity <- "Shannon"
sum_a_roti_shannon_erosion1$class <- "Erosion 0-1"
sum_a_roti_shannon_erosion1 <-sum_a_roti_shannon_erosion1[-1,]

sub_a_roti_sample_data_erosion2 <- a_roti_sample_data %>% filter(Erosion_grouped == "Medium high")
sum_a_roti_richness_erosion2 <- capture.output(summary(sub_a_roti_sample_data_erosion2$Richness),file=NULL,append=FALSE)
sum_a_roti_richness_erosion2 <- as.data.frame(sum_a_roti_richness_erosion2)
sum_a_roti_richness_erosion2$group <- "Rotifers"
sum_a_roti_richness_erosion2$diversity <- "Richness"
sum_a_roti_richness_erosion2$class <- "Erosion 1-5"
sum_a_roti_richness_erosion2 <-sum_a_roti_richness_erosion2[-1,]

sum_a_roti_shannon_erosion2 <- capture.output(summary(sub_a_roti_sample_data_erosion2$Shannon),file=NULL,append=FALSE)
sum_a_roti_shannon_erosion2 <- as.data.frame(sum_a_roti_shannon_erosion2)
sum_a_roti_shannon_erosion2$group <- "Rotifers"
sum_a_roti_shannon_erosion2$diversity <- "Shannon"
sum_a_roti_shannon_erosion2$class <- "Erosion 1-5"
sum_a_roti_shannon_erosion2 <-sum_a_roti_shannon_erosion2[-1,]

sub_a_roti_sample_data_erosion3 <- a_roti_sample_data %>% filter(Erosion_grouped == "High")
sum_a_roti_richness_erosion3 <- capture.output(summary(sub_a_roti_sample_data_erosion3$Richness),file=NULL,append=FALSE)
sum_a_roti_richness_erosion3 <- as.data.frame(sum_a_roti_richness_erosion3)
sum_a_roti_richness_erosion3$group <- "Rotifers"
sum_a_roti_richness_erosion3$diversity <- "Richness"
sum_a_roti_richness_erosion3$class <- "Erosion 5-10"
sum_a_roti_richness_erosion3 <-sum_a_roti_richness_erosion3[-1,]

sum_a_roti_shannon_erosion3 <- capture.output(summary(sub_a_roti_sample_data_erosion3$Shannon),file=NULL,append=FALSE)
sum_a_roti_shannon_erosion3 <- as.data.frame(sum_a_roti_shannon_erosion3)
sum_a_roti_shannon_erosion3$group <- "Rotifers"
sum_a_roti_shannon_erosion3$diversity <- "Shannon"
sum_a_roti_shannon_erosion3$class <- "Erosion 5-10"
sum_a_roti_shannon_erosion3 <-sum_a_roti_shannon_erosion3[-1,]

sub_a_roti_sample_data_erosion4 <- a_roti_sample_data %>% filter(Erosion_grouped == "Very high")
sum_a_roti_richness_erosion4 <- capture.output(summary(sub_a_roti_sample_data_erosion4$Richness),file=NULL,append=FALSE)
sum_a_roti_richness_erosion4 <- as.data.frame(sum_a_roti_richness_erosion4)
sum_a_roti_richness_erosion4$group <- "Rotifers"
sum_a_roti_richness_erosion4$diversity <- "Richness"
sum_a_roti_richness_erosion4$class <- "Erosion 10+"
sum_a_roti_richness_erosion4 <-sum_a_roti_richness_erosion4[-1,]

sum_a_roti_shannon_erosion4 <- capture.output(summary(sub_a_roti_sample_data_erosion4$Shannon),file=NULL,append=FALSE)
sum_a_roti_shannon_erosion4 <- as.data.frame(sum_a_roti_shannon_erosion4)
sum_a_roti_shannon_erosion4$group <- "Rotifers"
sum_a_roti_shannon_erosion4$diversity <- "Shannon"
sum_a_roti_shannon_erosion4$class <- "Erosion 10+"
sum_a_roti_shannon_erosion4 <-sum_a_roti_shannon_erosion4[-1,]

sub_p_sample_data_erosion1 <- p_sample_data %>% filter(Erosion_grouped == "Low")
sum_p_richness_erosion1 <- capture.output(summary(sub_p_sample_data_erosion1$Richness),file=NULL,append=FALSE)
sum_p_richness_erosion1 <- as.data.frame(sum_p_richness_erosion1)
sum_p_richness_erosion1$group <- "Protists"
sum_p_richness_erosion1$diversity <- "Richness"
sum_p_richness_erosion1$class <- "Erosion 0-1"
sum_p_richness_erosion1 <-sum_p_richness_erosion1[-1,]

sum_p_shannon_erosion1 <- capture.output(summary(sub_p_sample_data_erosion1$Shannon),file=NULL,append=FALSE)
sum_p_shannon_erosion1 <- as.data.frame(sum_p_shannon_erosion1)
sum_p_shannon_erosion1$group <- "Protists"
sum_p_shannon_erosion1$diversity <- "Shannon"
sum_p_shannon_erosion1$class <- "Erosion 0-1"
sum_p_shannon_erosion1 <-sum_p_shannon_erosion1[-1,]

sub_p_sample_data_erosion2 <- p_sample_data %>% filter(Erosion_grouped == "Medium high")
sum_p_richness_erosion2 <- capture.output(summary(sub_p_sample_data_erosion2$Richness),file=NULL,append=FALSE)
sum_p_richness_erosion2 <- as.data.frame(sum_p_richness_erosion2)
sum_p_richness_erosion2$group <- "Protists"
sum_p_richness_erosion2$diversity <- "Richness"
sum_p_richness_erosion2$class <- "Erosion 1-5"
sum_p_richness_erosion2 <-sum_p_richness_erosion2[-1,]

sum_p_shannon_erosion2 <- capture.output(summary(sub_p_sample_data_erosion2$Shannon),file=NULL,append=FALSE)
sum_p_shannon_erosion2 <- as.data.frame(sum_p_shannon_erosion2)
sum_p_shannon_erosion2$group <- "Protists"
sum_p_shannon_erosion2$diversity <- "Shannon"
sum_p_shannon_erosion2$class <- "Erosion 1-5"
sum_p_shannon_erosion2 <-sum_p_shannon_erosion2[-1,]

sub_p_sample_data_erosion3 <- p_sample_data %>% filter(Erosion_grouped == "High")
sum_p_richness_erosion3 <- capture.output(summary(sub_p_sample_data_erosion3$Richness),file=NULL,append=FALSE)
sum_p_richness_erosion3 <- as.data.frame(sum_p_richness_erosion3)
sum_p_richness_erosion3$group <- "Protists"
sum_p_richness_erosion3$diversity <- "Richness"
sum_p_richness_erosion3$class <- "Erosion 5-10"
sum_p_richness_erosion3 <-sum_p_richness_erosion3[-1,]

sum_p_shannon_erosion3 <- capture.output(summary(sub_p_sample_data_erosion3$Shannon),file=NULL,append=FALSE)
sum_p_shannon_erosion3 <- as.data.frame(sum_p_shannon_erosion3)
sum_p_shannon_erosion3$group <- "Protists"
sum_p_shannon_erosion3$diversity <- "Shannon"
sum_p_shannon_erosion3$class <- "Erosion 5-10"
sum_p_shannon_erosion3 <-sum_p_shannon_erosion3[-1,]

sub_p_sample_data_erosion4 <- p_sample_data %>% filter(Erosion_grouped == "Very high")
sum_p_richness_erosion4 <- capture.output(summary(sub_p_sample_data_erosion4$Richness),file=NULL,append=FALSE)
sum_p_richness_erosion4 <- as.data.frame(sum_p_richness_erosion4)
sum_p_richness_erosion4$group <- "Protists"
sum_p_richness_erosion4$diversity <- "Richness"
sum_p_richness_erosion4$class <- "Erosion 10+"
sum_p_richness_erosion4 <-sum_p_richness_erosion4[-1,]

sum_p_shannon_erosion4 <- capture.output(summary(sub_p_sample_data_erosion4$Shannon),file=NULL,append=FALSE)
sum_p_shannon_erosion4 <- as.data.frame(sum_p_shannon_erosion4)
sum_p_shannon_erosion4$group <- "Protists"
sum_p_shannon_erosion4$diversity <- "Shannon"
sum_p_shannon_erosion4$class <- "Erosion 10+"
sum_p_shannon_erosion4 <-sum_p_shannon_erosion4[-1,]

sub_f_sample_data_erosion1 <- f_sample_data %>% filter(Erosion_grouped == "Low")
sum_f_richness_erosion1 <- capture.output(summary(sub_f_sample_data_erosion1$Richness),file=NULL,append=FALSE)
sum_f_richness_erosion1 <- as.data.frame(sum_f_richness_erosion1)
sum_f_richness_erosion1$group <- "Fungi"
sum_f_richness_erosion1$diversity <- "Richness"
sum_f_richness_erosion1$class <- "Erosion 0-1"
sum_f_richness_erosion1 <-sum_f_richness_erosion1[-1,]

sum_f_shannon_erosion1 <- capture.output(summary(sub_f_sample_data_erosion1$Shannon),file=NULL,append=FALSE)
sum_f_shannon_erosion1 <- as.data.frame(sum_f_shannon_erosion1)
sum_f_shannon_erosion1$group <- "Fungi"
sum_f_shannon_erosion1$diversity <- "Shannon"
sum_f_shannon_erosion1$class <- "Erosion 0-1"
sum_f_shannon_erosion1 <-sum_f_shannon_erosion1[-1,]

sub_f_sample_data_erosion2 <- f_sample_data %>% filter(Erosion_grouped == "Medium high")
sum_f_richness_erosion2 <- capture.output(summary(sub_f_sample_data_erosion2$Richness),file=NULL,append=FALSE)
sum_f_richness_erosion2 <- as.data.frame(sum_f_richness_erosion2)
sum_f_richness_erosion2$group <- "Fungi"
sum_f_richness_erosion2$diversity <- "Richness"
sum_f_richness_erosion2$class <- "Erosion 1-5"
sum_f_richness_erosion2 <-sum_f_richness_erosion2[-1,]

sum_f_shannon_erosion2 <- capture.output(summary(sub_f_sample_data_erosion2$Shannon),file=NULL,append=FALSE)
sum_f_shannon_erosion2 <- as.data.frame(sum_f_shannon_erosion2)
sum_f_shannon_erosion2$group <- "Fungi"
sum_f_shannon_erosion2$diversity <- "Shannon"
sum_f_shannon_erosion2$class <- "Erosion 1-5"
sum_f_shannon_erosion2 <-sum_f_shannon_erosion2[-1,]

sub_f_sample_data_erosion3 <- f_sample_data %>% filter(Erosion_grouped == "High")
sum_f_richness_erosion3 <- capture.output(summary(sub_f_sample_data_erosion3$Richness),file=NULL,append=FALSE)
sum_f_richness_erosion3 <- as.data.frame(sum_f_richness_erosion3)
sum_f_richness_erosion3$group <- "Fungi"
sum_f_richness_erosion3$diversity <- "Richness"
sum_f_richness_erosion3$class <- "Erosion 5-10"
sum_f_richness_erosion3 <-sum_f_richness_erosion3[-1,]

sum_f_shannon_erosion3 <- capture.output(summary(sub_f_sample_data_erosion3$Shannon),file=NULL,append=FALSE)
sum_f_shannon_erosion3 <- as.data.frame(sum_f_shannon_erosion3)
sum_f_shannon_erosion3$group <- "Fungi"
sum_f_shannon_erosion3$diversity <- "Shannon"
sum_f_shannon_erosion3$class <- "Erosion 5-10"
sum_f_shannon_erosion3 <-sum_f_shannon_erosion3[-1,]

sub_f_sample_data_erosion4 <- f_sample_data %>% filter(Erosion_grouped == "Very high")
sum_f_richness_erosion4 <-capture.output(summary(sub_f_sample_data_erosion4$Richness),file=NULL,append=FALSE)
sum_f_richness_erosion4 <- as.data.frame(sum_f_richness_erosion4)
sum_f_richness_erosion4$group <- "Fungi"
sum_f_richness_erosion4$diversity <- "Richness"
sum_f_richness_erosion4$class <- "Erosion 10+"
sum_f_richness_erosion4 <-sum_f_richness_erosion4[-1,]

sum_f_shannon_erosion4 <- capture.output(summary(sub_f_sample_data_erosion4$Shannon),file=NULL,append=FALSE)
sum_f_shannon_erosion4 <- as.data.frame(sum_f_shannon_erosion4)
sum_f_shannon_erosion4$group <- "Fungi"
sum_f_shannon_erosion4$diversity <- "Shannon"
sum_f_shannon_erosion4$class <- "Erosion 10+"
sum_f_shannon_erosion4 <-sum_f_shannon_erosion4[-1,]

#Season
sub_f_sample_data_season1 <- f_sample_data %>% filter(Sample_season == "Spring")
sum_f_richness_season1 <- capture.output(summary(sub_f_sample_data_season1$Richness),file=NULL,append=FALSE)
sum_f_richness_season1 <- as.data.frame(sum_f_richness_season1)
sum_f_richness_season1$group <- "Fungi"
sum_f_richness_season1$diversity <- "Richness"
sum_f_richness_season1$class <- "Season - Spring"
sum_f_richness_season1 <-sum_f_richness_season1[-1,]

sum_f_shannon_season1 <- capture.output(summary(sub_f_sample_data_season1$Shannon),file=NULL,append=FALSE)
sum_f_shannon_season1 <- as.data.frame(sum_f_shannon_season1)
sum_f_shannon_season1$group <- "Fungi"
sum_f_shannon_season1$diversity <- "Shannon"
sum_f_shannon_season1$class <- "Season - Spring"
sum_f_shannon_season1 <-sum_f_shannon_season1[-1,]

sub_f_sample_data_season2 <- f_sample_data %>% filter(Sample_season == "Summer")
sum_f_richness_season2 <- capture.output(summary(sub_f_sample_data_season2$Richness),file=NULL,append=FALSE)
sum_f_richness_season2 <- as.data.frame(sum_f_richness_season2)
sum_f_richness_season2$group <- "Fungi"
sum_f_richness_season2$diversity <- "Richness"
sum_f_richness_season2$class <- "Season - Summer"
sum_f_richness_season2 <-sum_f_richness_season2[-1,]

sum_f_shannon_season2 <- capture.output(summary(sub_f_sample_data_season2$Shannon),file=NULL,append=FALSE)
sum_f_shannon_season2 <- as.data.frame(sum_f_shannon_season2)
sum_f_shannon_season2$group <- "Fungi"
sum_f_shannon_season2$diversity <- "Shannon"
sum_f_shannon_season2$class <- "Season - Summer"
sum_f_shannon_season2 <-sum_f_shannon_season2[-1,]

sub_f_sample_data_season3 <- f_sample_data %>% filter(Sample_season == "Autumn")
sum_f_richness_season3 <- capture.output(summary(sub_f_sample_data_season3$Richness),file=NULL,append=FALSE)
sum_f_richness_season3 <- as.data.frame(sum_f_richness_season3)
sum_f_richness_season3$group <- "Fungi"
sum_f_richness_season3$diversity <- "Richness"
sum_f_richness_season3$class <- "Season - Autumn"
sum_f_richness_season3 <-sum_f_richness_season3[-1,]

sum_f_shannon_season3 <- capture.output(summary(sub_f_sample_data_season3$Shannon),file=NULL,append=FALSE)
sum_f_shannon_season3 <- as.data.frame(sum_f_shannon_season3)
sum_f_shannon_season3$group <- "Fungi"
sum_f_shannon_season3$diversity <- "Shannon"
sum_f_shannon_season3$class <- "Season - Autumn"
sum_f_shannon_season3 <-sum_f_shannon_season3[-1,]

sub_p_sample_data_season1 <- p_sample_data %>% filter(Sample_season == "Spring")
sum_p_richness_season1 <- capture.output(summary(sub_p_sample_data_season1$Richness),file=NULL,append=FALSE)
sum_p_richness_season1 <- as.data.frame(sum_p_richness_season1)
sum_p_richness_season1$group <- "Protists"
sum_p_richness_season1$diversity <- "Richness"
sum_p_richness_season1$class <- "Season - Spring"
sum_p_richness_season1 <-sum_p_richness_season1[-1,]

sum_p_shannon_season1 <- capture.output(summary(sub_p_sample_data_season1$Shannon),file=NULL,append=FALSE)
sum_p_shannon_season1 <- as.data.frame(sum_p_shannon_season1)
sum_p_shannon_season1$group <- "Protists"
sum_p_shannon_season1$diversity <- "Shannon"
sum_p_shannon_season1$class <- "Season - Spring"
sum_p_shannon_season1 <-sum_p_shannon_season1[-1,]

sub_p_sample_data_season2 <- p_sample_data %>% filter(Sample_season == "Summer")
sum_p_richness_season2 <- capture.output(summary(sub_p_sample_data_season2$Richness),file=NULL,append=FALSE)
sum_p_richness_season2 <- as.data.frame(sum_p_richness_season2)
sum_p_richness_season2$group <- "Protists"
sum_p_richness_season2$diversity <- "Richness"
sum_p_richness_season2$class <- "Season - Summer"
sum_p_richness_season2 <-sum_p_richness_season2[-1,]

sum_p_shannon_season2 <- capture.output(summary(sub_p_sample_data_season2$Shannon),file=NULL,append=FALSE)
sum_p_shannon_season2 <- as.data.frame(sum_p_shannon_season2)
sum_p_shannon_season2$group <- "Protists"
sum_p_shannon_season2$diversity <- "Shannon"
sum_p_shannon_season2$class <- "Season - Summer"
sum_p_shannon_season2 <-sum_p_shannon_season2[-1,]

sub_p_sample_data_season3 <- p_sample_data %>% filter(Sample_season == "Autumn")
sum_p_richness_season3 <- capture.output(summary(sub_p_sample_data_season3$Richness),file=NULL,append=FALSE)
sum_p_richness_season3 <- as.data.frame(sum_p_richness_season3)
sum_p_richness_season3$group <- "Protists"
sum_p_richness_season3$diversity <- "Richness"
sum_p_richness_season3$class <- "Season - Autumn"
sum_p_richness_season3 <-sum_p_richness_season3[-1,]

sum_p_shannon_season3 <- capture.output(summary(sub_p_sample_data_season3$Shannon),file=NULL,append=FALSE)
sum_p_shannon_season3 <- as.data.frame(sum_p_shannon_season3)
sum_p_shannon_season3$group <- "Protists"
sum_p_shannon_season3$diversity <- "Shannon"
sum_p_shannon_season3$class <- "Season - Autumn"
sum_p_shannon_season3 <-sum_p_shannon_season3[-1,]

sub_a_roti_sample_data_season1 <- a_roti_sample_data %>% filter(Sample_season == "Spring")
sum_a_roti_richness_season1 <- capture.output(summary(sub_a_roti_sample_data_season1$Richness),file=NULL,append=FALSE)
sum_a_roti_richness_season1 <- as.data.frame(sum_a_roti_richness_season1)
sum_a_roti_richness_season1$group <- "Rotifers"
sum_a_roti_richness_season1$diversity <- "Richness"
sum_a_roti_richness_season1$class <- "Season - Spring"
sum_a_roti_richness_season1 <-sum_a_roti_richness_season1[-1,]

sum_a_roti_shannon_season1 <- capture.output(summary(sub_a_roti_sample_data_season1$Shannon),file=NULL,append=FALSE)
sum_a_roti_shannon_season1 <- as.data.frame(sum_a_roti_shannon_season1)
sum_a_roti_shannon_season1$group <- "Rotifers"
sum_a_roti_shannon_season1$diversity <- "Shannon"
sum_a_roti_shannon_season1$class <- "Season - Spring"
sum_a_roti_shannon_season1 <-sum_a_roti_shannon_season1[-1,]

sub_a_roti_sample_data_season2 <- a_roti_sample_data %>% filter(Sample_season == "Summer")
sum_a_roti_richness_season2 <- capture.output(summary(sub_a_roti_sample_data_season2$Richness),file=NULL,append=FALSE)
sum_a_roti_richness_season2 <- as.data.frame(sum_a_roti_richness_season2)
sum_a_roti_richness_season2$group <- "Rotifers"
sum_a_roti_richness_season2$diversity <- "Richness"
sum_a_roti_richness_season2$class <- "Season - Summer"
sum_a_roti_richness_season2 <-sum_a_roti_richness_season2[-1,]

sum_a_roti_shannon_season2 <- capture.output(summary(sub_a_roti_sample_data_season2$Shannon),file=NULL,append=FALSE)
sum_a_roti_shannon_season2 <- as.data.frame(sum_a_roti_shannon_season2)
sum_a_roti_shannon_season2$group <- "Rotifers"
sum_a_roti_shannon_season2$diversity <- "Shannon"
sum_a_roti_shannon_season2$class <- "Season - Summer"
sum_a_roti_shannon_season2 <-sum_a_roti_shannon_season2[-1,]

sub_a_roti_sample_data_season3 <- a_roti_sample_data %>% filter(Sample_season == "Autumn")
sum_a_roti_richness_season3 <- capture.output(summary(sub_a_roti_sample_data_season3$Richness),file=NULL,append=FALSE)
sum_a_roti_richness_season3 <- as.data.frame(sum_a_roti_richness_season3)
sum_a_roti_richness_season3$group <- "Rotifers"
sum_a_roti_richness_season3$diversity <- "Richness"
sum_a_roti_richness_season3$class <- "Season - Autumn"
sum_a_roti_richness_season3 <-sum_a_roti_richness_season3[-1,]

sum_a_roti_shannon_season3 <- capture.output(summary(sub_a_roti_sample_data_season3$Shannon),file=NULL,append=FALSE)
sum_a_roti_shannon_season3 <- as.data.frame(sum_a_roti_shannon_season3)
sum_a_roti_shannon_season3$group <- "Rotifers"
sum_a_roti_shannon_season3$diversity <- "Shannon"
sum_a_roti_shannon_season3$class <- "Season - Autumn"
sum_a_roti_shannon_season3 <-sum_a_roti_shannon_season3[-1,]

sub_a_nema_sample_data_season1 <- a_nema_sample_data %>% filter(Sample_season == "Spring")
sum_a_nema_richness_season1 <- capture.output(summary(sub_a_nema_sample_data_season1$Richness),file=NULL,append=FALSE)
sum_a_nema_richness_season1 <- as.data.frame(sum_a_nema_richness_season1)
sum_a_nema_richness_season1$group <- "Nematodes"
sum_a_nema_richness_season1$diversity <- "Richness"
sum_a_nema_richness_season1$class <- "Season - Spring"
sum_a_nema_richness_season1 <-sum_a_nema_richness_season1[-1,]

sum_a_nema_shannon_season1 <- capture.output(summary(sub_a_nema_sample_data_season1$Shannon),file=NULL,append=FALSE)
sum_a_nema_shannon_season1 <- as.data.frame(sum_a_nema_shannon_season1)
sum_a_nema_shannon_season1$group <- "Nematodes"
sum_a_nema_shannon_season1$diversity <- "Shannon"
sum_a_nema_shannon_season1$class <- "Season - Spring"
sum_a_nema_shannon_season1 <-sum_a_nema_shannon_season1[-1,]

sub_a_nema_sample_data_season2 <- a_nema_sample_data %>% filter(Sample_season == "Summer")
sum_a_nema_richness_season2 <- capture.output(summary(sub_a_nema_sample_data_season2$Richness),file=NULL,append=FALSE)
sum_a_nema_richness_season2 <- as.data.frame(sum_a_nema_richness_season2)
sum_a_nema_richness_season2$group <- "Nematodes"
sum_a_nema_richness_season2$diversity <- "Richness"
sum_a_nema_richness_season2$class <- "Season - Summer"
sum_a_nema_richness_season2 <-sum_a_nema_richness_season2[-1,]

sum_a_nema_shannon_season2 <- capture.output(summary(sub_a_nema_sample_data_season2$Shannon),file=NULL,append=FALSE)
sum_a_nema_shannon_season2 <- as.data.frame(sum_a_nema_shannon_season2)
sum_a_nema_shannon_season2$group <- "Nematodes"
sum_a_nema_shannon_season2$diversity <- "Shannon"
sum_a_nema_shannon_season2$class <- "Season - Summer"
sum_a_nema_shannon_season2 <-sum_a_nema_shannon_season2[-1,]

sub_a_nema_sample_data_season3 <- a_nema_sample_data %>% filter(Sample_season == "Autumn")
sum_a_nema_richness_season3 <- capture.output(summary(sub_a_nema_sample_data_season3$Richness),file=NULL,append=FALSE)
sum_a_nema_richness_season3 <- as.data.frame(sum_a_nema_richness_season3)
sum_a_nema_richness_season3$group <- "Nematodes"
sum_a_nema_richness_season3$diversity <- "Richness"
sum_a_nema_richness_season3$class <- "Season - Autumn"
sum_a_nema_richness_season3 <-sum_a_nema_richness_season3[-1,]

sum_a_nema_shannon_season3 <- capture.output(summary(sub_a_nema_sample_data_season3$Shannon),file=NULL,append=FALSE)
sum_a_nema_shannon_season3 <- as.data.frame(sum_a_nema_shannon_season3)
sum_a_nema_shannon_season3$group <- "Nematodes"
sum_a_nema_shannon_season3$diversity <- "Shannon"
sum_a_nema_shannon_season3$class <- "Season - Autumn"
sum_a_nema_shannon_season3 <-sum_a_nema_shannon_season3[-1,]

sub_a_tardi_sample_data_season1 <- a_tardi_sample_data %>% filter(Sample_season == "Spring")
sum_a_tardi_richness_season1 <- capture.output(summary(sub_a_tardi_sample_data_season1$Richness),file=NULL,append=FALSE)
sum_a_tardi_richness_season1 <- as.data.frame(sum_a_tardi_richness_season1)
sum_a_tardi_richness_season1$group <- "Tardigrades"
sum_a_tardi_richness_season1$diversity <- "Richness"
sum_a_tardi_richness_season1$class <- "Season - Spring"
sum_a_tardi_richness_season1 <-sum_a_tardi_richness_season1[-1,]

sum_a_tardi_shannon_season1 <- capture.output(summary(sub_a_tardi_sample_data_season1$Shannon),file=NULL,append=FALSE)
sum_a_tardi_shannon_season1 <- as.data.frame(sum_a_tardi_shannon_season1)
sum_a_tardi_shannon_season1$group <- "Tardigrades"
sum_a_tardi_shannon_season1$diversity <- "Shannon"
sum_a_tardi_shannon_season1$class <- "Season - Spring"
sum_a_tardi_shannon_season1 <-sum_a_tardi_shannon_season1[-1,]

sub_a_tardi_sample_data_season2 <- a_tardi_sample_data %>% filter(Sample_season == "Summer")
sum_a_tardi_richness_season2 <- capture.output(summary(sub_a_tardi_sample_data_season2$Richness),file=NULL,append=FALSE)
sum_a_tardi_richness_season2 <- as.data.frame(sum_a_tardi_richness_season2)
sum_a_tardi_richness_season2$group <- "Tardigrades"
sum_a_tardi_richness_season2$diversity <- "Richness"
sum_a_tardi_richness_season2$class <- "Season - Summer"
sum_a_tardi_richness_season2 <-sum_a_tardi_richness_season2[-1,]

sum_a_tardi_shannon_season2 <- capture.output(summary(sub_a_tardi_sample_data_season2$Shannon),file=NULL,append=FALSE)
sum_a_tardi_shannon_season2 <- as.data.frame(sum_a_tardi_shannon_season2)
sum_a_tardi_shannon_season2$group <- "Tardigrades"
sum_a_tardi_shannon_season2$diversity <- "Shannon"
sum_a_tardi_shannon_season2$class <- "Season - Summer"
sum_a_tardi_shannon_season2 <-sum_a_tardi_shannon_season2[-1,]

sub_a_tardi_sample_data_season3 <- a_tardi_sample_data %>% filter(Sample_season == "Autumn")
sum_a_tardi_richness_season3 <- capture.output(summary(sub_a_tardi_sample_data_season3$Richness),file=NULL,append=FALSE)
sum_a_tardi_richness_season3 <- as.data.frame(sum_a_tardi_richness_season3)
sum_a_tardi_richness_season3$group <- "Tardigrades"
sum_a_tardi_richness_season3$diversity <- "Richness"
sum_a_tardi_richness_season3$class <- "Season - Autumn"
sum_a_tardi_richness_season3 <-sum_a_tardi_richness_season3[-1,]

sum_a_tardi_shannon_season3 <- capture.output(summary(sub_a_tardi_sample_data_season3$Shannon),file=NULL,append=FALSE)
sum_a_tardi_shannon_season3 <- as.data.frame(sum_a_tardi_shannon_season3)
sum_a_tardi_shannon_season3$group <- "Tardigrades"
sum_a_tardi_shannon_season3$diversity <- "Shannon"
sum_a_tardi_shannon_season3$class <- "Season - Autumn"
sum_a_tardi_shannon_season3 <-sum_a_tardi_shannon_season3[-1,]

sub_a_tardi_sample_data_depth1 <- a_tardi_sample_data %>% filter(depth_grouped == "Shallow")
sum_a_tardi_richness_depth1 <-capture.output(summary(sub_a_tardi_sample_data_depth1$Richness),file=NULL,append=FALSE)
sum_a_tardi_richness_depth1 <- as.data.frame(sum_a_tardi_richness_depth1)
sum_a_tardi_richness_depth1$group <- "Tardigrades"
sum_a_tardi_richness_depth1$diversity <- "Richness"
sum_a_tardi_richness_depth1$class <- "Depth <1m"
sum_a_tardi_richness_depth1 <-sum_a_tardi_richness_depth1[-1,]

sum_a_tardi_shannon_depth1 <-capture.output(summary(sub_a_tardi_sample_data_depth1$Shannon),file=NULL,append=FALSE)
sum_a_tardi_shannon_depth1 <- as.data.frame(sum_a_tardi_shannon_depth1)
sum_a_tardi_shannon_depth1$group <- "Tardigrades"
sum_a_tardi_shannon_depth1$diversity <- "Shannon"
sum_a_tardi_shannon_depth1$class <- "Depth <1m"
sum_a_tardi_shannon_depth1 <-sum_a_tardi_shannon_depth1[-1,]

sub_a_tardi_sample_data_depth2 <- a_tardi_sample_data %>% filter(depth_grouped == "Moderately deep")
sum_a_tardi_richness_depth2 <-capture.output(summary(sub_a_tardi_sample_data_depth2$Richness),file=NULL,append=FALSE)
sum_a_tardi_richness_depth2 <- as.data.frame(sum_a_tardi_richness_depth2)
sum_a_tardi_richness_depth2$group <- "Tardigrades"
sum_a_tardi_richness_depth2$diversity <- "Richness"
sum_a_tardi_richness_depth2$class <- "Depth 1-2m"
sum_a_tardi_richness_depth2 <-sum_a_tardi_richness_depth2[-1,]

sum_a_tardi_shannon_depth2 <- capture.output(summary(sub_a_tardi_sample_data_depth2$Shannon),file=NULL,append=FALSE)
sum_a_tardi_shannon_depth2 <- as.data.frame(sum_a_tardi_shannon_depth2)
sum_a_tardi_shannon_depth2$group <- "Tardigrades"
sum_a_tardi_shannon_depth2$diversity <- "Shannon"
sum_a_tardi_shannon_depth2$class <- "Depth 1-2m"
sum_a_tardi_shannon_depth2 <-sum_a_tardi_shannon_depth2[-1,]

sub_a_tardi_sample_data_depth3 <- a_tardi_sample_data %>% filter(depth_grouped == "Deep")
sum_a_tardi_richness_depth3 <- capture.output(summary(sub_a_tardi_sample_data_depth3$Richness),file=NULL,append=FALSE)
sum_a_tardi_richness_depth3 <- as.data.frame(sum_a_tardi_richness_depth3)
sum_a_tardi_richness_depth3$group <- "Tardigrades"
sum_a_tardi_richness_depth3$diversity <- "Richness"
sum_a_tardi_richness_depth3$class <- "Depth >2m"
sum_a_tardi_richness_depth3 <-sum_a_tardi_richness_depth3[-1,]

sum_a_tardi_shannon_depth3 <- capture.output(summary(sub_a_tardi_sample_data_depth3$Shannon),file=NULL,append=FALSE)
sum_a_tardi_shannon_depth3 <- as.data.frame(sum_a_tardi_shannon_depth3)
sum_a_tardi_shannon_depth3$group <- "Tardigrades"
sum_a_tardi_shannon_depth3$diversity <- "Shannon"
sum_a_tardi_shannon_depth3$class <- "Depth >2m"
sum_a_tardi_shannon_depth3 <-sum_a_tardi_shannon_depth3[-1,]

sub_a_nema_sample_data_depth1 <- a_nema_sample_data %>% filter(depth_grouped == "Shallow")
sum_a_nema_richness_depth1 <- capture.output(summary(sub_a_nema_sample_data_depth1$Richness),file=NULL,append=FALSE)
sum_a_nema_richness_depth1 <- as.data.frame(sum_a_nema_richness_depth1)
sum_a_nema_richness_depth1$group <- "Nematodes"
sum_a_nema_richness_depth1$diversity <- "Richness"
sum_a_nema_richness_depth1$class <- "Depth <1m"
sum_a_nema_richness_depth1 <-sum_a_nema_richness_depth1[-1,]

sum_a_nema_shannon_depth1 <- capture.output(summary(sub_a_nema_sample_data_depth1$Shannon),file=NULL,append=FALSE)
sum_a_nema_shannon_depth1 <- as.data.frame(sum_a_nema_shannon_depth1)
sum_a_nema_shannon_depth1$group <- "Nematodes"
sum_a_nema_shannon_depth1$diversity <- "Shannon"
sum_a_nema_shannon_depth1$class <- "Depth <1m"
sum_a_nema_shannon_depth1 <-sum_a_nema_shannon_depth1[-1,]

sub_a_nema_sample_data_depth2 <- a_nema_sample_data %>% filter(depth_grouped == "Moderately deep")
sum_a_nema_richness_depth2 <- capture.output(summary(sub_a_nema_sample_data_depth2$Richness),file=NULL,append=FALSE)
sum_a_nema_richness_depth2 <- as.data.frame(sum_a_nema_richness_depth2)
sum_a_nema_richness_depth2$group <- "Nematodes"
sum_a_nema_richness_depth2$diversity <- "Richness"
sum_a_nema_richness_depth2$class <- "Depth 1-2m"
sum_a_nema_richness_depth2 <-sum_a_nema_richness_depth2[-1,]

sum_a_nema_shannon_depth2 <- capture.output(summary(sub_a_nema_sample_data_depth2$Shannon),file=NULL,append=FALSE)
sum_a_nema_shannon_depth2 <- as.data.frame(sum_a_nema_shannon_depth2)
sum_a_nema_shannon_depth2$group <- "Nematodes"
sum_a_nema_shannon_depth2$diversity <- "Shannon"
sum_a_nema_shannon_depth2$class <- "Depth 1-2m"
sum_a_nema_shannon_depth2 <-sum_a_nema_shannon_depth2[-1,]

sub_a_nema_sample_data_depth3 <- a_nema_sample_data %>% filter(depth_grouped == "Deep")
sum_a_nema_richness_depth3 <- capture.output(summary(sub_a_nema_sample_data_depth3$Richness),file=NULL,append=FALSE)
sum_a_nema_richness_depth3 <- as.data.frame(sum_a_nema_richness_depth3)
sum_a_nema_richness_depth3$group <- "Nematodes"
sum_a_nema_richness_depth3$diversity <- "Richness"
sum_a_nema_richness_depth3$class <- "Depth >2m"
sum_a_nema_richness_depth3 <-sum_a_nema_richness_depth3[-1,]

sum_a_nema_shannon_depth3 <- capture.output(summary(sub_a_nema_sample_data_depth3$Shannon),file=NULL,append=FALSE)
sum_a_nema_shannon_depth3 <- as.data.frame(sum_a_nema_shannon_depth3)
sum_a_nema_shannon_depth3$group <- "Nematodes"
sum_a_nema_shannon_depth3$diversity <- "Shannon"
sum_a_nema_shannon_depth3$class <- "Depth >2m"
sum_a_nema_shannon_depth3 <-sum_a_nema_shannon_depth3[-1,]

sub_a_roti_sample_data_depth1 <- a_roti_sample_data %>% filter(depth_grouped == "Shallow")
sum_a_roti_richness_depth1 <- capture.output(summary(sub_a_roti_sample_data_depth1$Richness),file=NULL,append=FALSE)
sum_a_roti_richness_depth1 <- as.data.frame(sum_a_roti_richness_depth1)
sum_a_roti_richness_depth1$group <- "Rotifers"
sum_a_roti_richness_depth1$diversity <- "Richness"
sum_a_roti_richness_depth1$class <- "Depth <1m"
sum_a_roti_richness_depth1 <-sum_a_roti_richness_depth1[-1,]

sum_a_roti_shannon_depth1 <- capture.output(summary(sub_a_roti_sample_data_depth1$Shannon),file=NULL,append=FALSE)
sum_a_roti_shannon_depth1 <- as.data.frame(sum_a_roti_shannon_depth1)
sum_a_roti_shannon_depth1$group <- "Rotifers"
sum_a_roti_shannon_depth1$diversity <- "Shannon"
sum_a_roti_shannon_depth1$class <- "Depth >1m"
sum_a_roti_shannon_depth1 <-sum_a_roti_shannon_depth1[-1,]

sub_a_roti_sample_data_depth2 <- a_roti_sample_data %>% filter(depth_grouped == "Moderately deep")
sum_a_roti_richness_depth2 <- capture.output(summary(sub_a_roti_sample_data_depth2$Richness),file=NULL,append=FALSE)
sum_a_roti_richness_depth2 <- as.data.frame(sum_a_roti_richness_depth2)
sum_a_roti_richness_depth2$group <- "Rotifers"
sum_a_roti_richness_depth2$diversity <- "Richness"
sum_a_roti_richness_depth2$class <- "Depth 1-2m"
sum_a_roti_richness_depth2 <-sum_a_roti_richness_depth2[-1,]

sum_a_roti_shannon_depth2 <- capture.output(summary(sub_a_roti_sample_data_depth2$Shannon),file=NULL,append=FALSE)
sum_a_roti_shannon_depth2 <- as.data.frame(sum_a_roti_shannon_depth2)
sum_a_roti_shannon_depth2$group <- "Rotifers"
sum_a_roti_shannon_depth2$diversity <- "Shannon"
sum_a_roti_shannon_depth2$class <- "Depth 1-2m"
sum_a_roti_shannon_depth2 <-sum_a_roti_shannon_depth2[-1,]

sub_a_roti_sample_data_depth3 <- a_roti_sample_data %>% filter(depth_grouped == "Deep")
sum_a_roti_richness_depth3 <- capture.output(summary(sub_a_roti_sample_data_depth3$Richness),file=NULL,append=FALSE)
sum_a_roti_richness_depth3 <- as.data.frame(sum_a_roti_richness_depth3)
sum_a_roti_richness_depth3$group <- "Rotifers"
sum_a_roti_richness_depth3$diversity <- "Richness"
sum_a_roti_richness_depth3$class <- "Depth >2m"
sum_a_roti_richness_depth3 <-sum_a_roti_richness_depth3[-1,]

sum_a_roti_shannon_depth3 <- capture.output(summary(sub_a_roti_sample_data_depth3$Shannon),file=NULL,append=FALSE)
sum_a_roti_shannon_depth3 <- as.data.frame(sum_a_roti_shannon_depth3)
sum_a_roti_shannon_depth3$group <- "Rotifers"
sum_a_roti_shannon_depth3$diversity <- "Shannon"
sum_a_roti_shannon_depth3$class <- "Depth >2m"
sum_a_roti_shannon_depth3 <-sum_a_roti_shannon_depth3[-1,]

sub_p_sample_data_depth1 <- p_sample_data %>% filter(depth_grouped == "Shallow")
sum_p_richness_depth1 <- capture.output(summary(sub_p_sample_data_depth1$Richness),file=NULL,append=FALSE)
sum_p_richness_depth1 <- as.data.frame(sum_p_richness_depth1)
sum_p_richness_depth1$group <- "Protists"
sum_p_richness_depth1$diversity <- "Richness"
sum_p_richness_depth1$class <- "Depth <1m"
sum_p_richness_depth1 <-sum_p_richness_depth1[-1,]

sum_p_shannon_depth1 <- capture.output(summary(sub_p_sample_data_depth1$Shannon),file=NULL,append=FALSE)
sum_p_shannon_depth1 <- as.data.frame(sum_p_shannon_depth1)
sum_p_shannon_depth1$group <- "Protists"
sum_p_shannon_depth1$diversity <- "Shannon"
sum_p_shannon_depth1$class <- "Depth <1m"
sum_p_shannon_depth1 <-sum_p_shannon_depth1[-1,]

sub_p_sample_data_depth2 <- p_sample_data %>% filter(depth_grouped == "Moderately deep")
sum_p_richness_depth2 <- capture.output(summary(sub_p_sample_data_depth2$Richness),file=NULL,append=FALSE)
sum_p_richness_depth2 <- as.data.frame(sum_p_richness_depth2)
sum_p_richness_depth2$group <- "Protists"
sum_p_richness_depth2$diversity <- "Richness"
sum_p_richness_depth2$class <- "Depth 1-2m"
sum_p_richness_depth2 <-sum_p_richness_depth2[-1,]

sum_p_shannon_depth2 <-capture.output(summary(sub_p_sample_data_depth2$Shannon),file=NULL,append=FALSE)
sum_p_shannon_depth2 <- as.data.frame(sum_p_shannon_depth2)
sum_p_shannon_depth2$group <- "Protists"
sum_p_shannon_depth2$diversity <- "Shannon"
sum_p_shannon_depth2$class <- "Depth 1-2m"
sum_p_shannon_depth2 <-sum_p_shannon_depth2[-1,]

sub_p_sample_data_depth3 <- p_sample_data %>% filter(depth_grouped == "Deep")
sum_p_richness_depth3 <- capture.output(summary(sub_p_sample_data_depth3$Richness),file=NULL,append=FALSE)
sum_p_richness_depth3 <- as.data.frame(sum_p_richness_depth3)
sum_p_richness_depth3$group <- "Protists"
sum_p_richness_depth3$diversity <- "Richness"
sum_p_richness_depth3$class <- "Depth >2m"
sum_p_richness_depth3 <-sum_p_richness_depth3[-1,]

sum_p_shannon_depth3 <-capture.output(summary(sub_p_sample_data_depth3$Shannon),file=NULL,append=FALSE)
sum_p_shannon_depth3 <- as.data.frame(sum_p_shannon_depth3)
sum_p_shannon_depth3$group <- "Protists"
sum_p_shannon_depth3$diversity <- "Shannon"
sum_p_shannon_depth3$class <- "Depth >2m"
sum_p_shannon_depth3 <-sum_p_shannon_depth3[-1,]

sub_f_sample_data_depth1 <- f_sample_data %>% filter(depth_grouped == "Shallow")
sum_f_richness_depth1 <- capture.output(summary(sub_f_sample_data_depth1$Richness),file=NULL,append=FALSE)
sum_f_richness_depth1 <- as.data.frame(sum_f_richness_depth1)
sum_f_richness_depth1$group <- "Fungi"
sum_f_richness_depth1$diversity <- "Richness"
sum_f_richness_depth1$class <- "Depth <1m"
sum_f_richness_depth1 <-sum_f_richness_depth1[-1,]

sum_f_shannon_depth1 <- capture.output(summary(sub_f_sample_data_depth1$Shannon),file=NULL,append=FALSE)
sum_f_shannon_depth1 <- as.data.frame(sum_f_shannon_depth1)
sum_f_shannon_depth1$group <- "Fungi"
sum_f_shannon_depth1$diversity <- "Shannon"
sum_f_shannon_depth1$class <- "Depth <1m"
sum_f_shannon_depth1 <-sum_f_shannon_depth1[-1,]

sub_f_sample_data_depth2 <- f_sample_data %>% filter(depth_grouped == "Moderately deep")
sum_f_richness_depth2 <- capture.output(summary(sub_f_sample_data_depth2$Richness),file=NULL,append=FALSE)
sum_f_richness_depth2 <- as.data.frame(sum_f_richness_depth2)
sum_f_richness_depth2$group <- "Fungi"
sum_f_richness_depth2$diversity <- "Richness"
sum_f_richness_depth2$class <- "Depth 1-2m"
sum_f_richness_depth2 <-sum_f_richness_depth2[-1,]

sum_f_shannon_depth2 <-capture.output(summary(sub_f_sample_data_depth2$Shannon),file=NULL,append=FALSE)
sum_f_shannon_depth2 <- as.data.frame(sum_f_shannon_depth2)
sum_f_shannon_depth2$group <- "Fungi"
sum_f_shannon_depth2$diversity <- "Shannon"
sum_f_shannon_depth2$class <- "Depth 1-2m"
sum_f_shannon_depth2 <-sum_f_shannon_depth2[-1,]

sub_f_sample_data_depth3 <- f_sample_data %>% filter(depth_grouped == "Deep")
sum_f_richness_depth3 <- capture.output(summary(sub_f_sample_data_depth3$Richness),file=NULL,append=FALSE)
sum_f_richness_depth3 <- as.data.frame(sum_f_richness_depth3)
sum_f_richness_depth3$group <- "Fungi"
sum_f_richness_depth3$diversity <- "Richness"
sum_f_richness_depth3$class <- "Depth >2m"
sum_f_richness_depth3 <-sum_f_richness_depth3[-1,]

sum_f_shannon_depth3 <- capture.output(summary(sub_f_sample_data_depth3$Shannon),file=NULL,append=FALSE)
sum_f_shannon_depth3 <- as.data.frame(sum_f_shannon_depth3)
sum_f_shannon_depth3$group <- "Fungi"
sum_f_shannon_depth3$diversity <- "Shannon"
sum_f_shannon_depth3$class <- "Depth >2m"
sum_f_shannon_depth3 <-sum_f_shannon_depth3[-1,]


#pH
sub_a_tardi_sample_data_pH1 <- a_tardi_sample_data %>% filter(pH_grouped == "Acidic")
sum_a_tardi_richness_pH1 <-capture.output(summary(sub_a_tardi_sample_data_pH1$Richness),file=NULL,append=FALSE)
sum_a_tardi_richness_pH1 <- as.data.frame(sum_a_tardi_richness_pH1)
sum_a_tardi_richness_pH1$group <- "Tardigrades"
sum_a_tardi_richness_pH1$diversity <- "Richness"
sum_a_tardi_richness_pH1$class <- "pH 3.3-6.5"
sum_a_tardi_richness_pH1 <-sum_a_tardi_richness_pH1[-1,]

sum_a_tardi_shannon_pH1 <-capture.output(summary(sub_a_tardi_sample_data_pH1$Shannon),file=NULL,append=FALSE)
sum_a_tardi_shannon_pH1 <- as.data.frame(sum_a_tardi_shannon_pH1)
sum_a_tardi_shannon_pH1$group <- "Tardigrades"
sum_a_tardi_shannon_pH1$diversity <- "Shannon"
sum_a_tardi_shannon_pH1$class <- "pH 3.3-6.5"
sum_a_tardi_shannon_pH1 <-sum_a_tardi_shannon_pH1[-1,]

sub_a_tardi_sample_data_pH2 <- a_tardi_sample_data %>% filter(pH_grouped == "Neutral")
sum_a_tardi_richness_pH2 <-capture.output(summary(sub_a_tardi_sample_data_pH2$Richness),file=NULL,append=FALSE)
sum_a_tardi_richness_pH2 <- as.data.frame(sum_a_tardi_richness_pH2)
sum_a_tardi_richness_pH2$group <- "Tardigrades"
sum_a_tardi_richness_pH2$diversity <- "Richness"
sum_a_tardi_richness_pH2$class <- "pH 6.5-7.5"
sum_a_tardi_richness_pH2 <-sum_a_tardi_richness_pH2[-1,]

sum_a_tardi_shannon_pH2 <- capture.output(summary(sub_a_tardi_sample_data_pH2$Shannon),file=NULL,append=FALSE)
sum_a_tardi_shannon_pH2 <- as.data.frame(sum_a_tardi_shannon_pH2)
sum_a_tardi_shannon_pH2$group <- "Tardigrades"
sum_a_tardi_shannon_pH2$diversity <- "Shannon"
sum_a_tardi_shannon_pH2$class <- "pH 6.5-7.5"
sum_a_tardi_shannon_pH2 <-sum_a_tardi_shannon_pH2[-1,]

sub_a_tardi_sample_data_pH3 <- a_tardi_sample_data %>% filter(pH_grouped == "Alkaline")
sum_a_tardi_richness_pH3 <- capture.output(summary(sub_a_tardi_sample_data_pH3$Richness),file=NULL,append=FALSE)
sum_a_tardi_richness_pH3 <- as.data.frame(sum_a_tardi_richness_pH3)
sum_a_tardi_richness_pH3$group <- "Tardigrades"
sum_a_tardi_richness_pH3$diversity <- "Richness"
sum_a_tardi_richness_pH3$class <- "pH 7.5+"
sum_a_tardi_richness_pH3 <-sum_a_tardi_richness_pH3[-1,]

sum_a_tardi_shannon_pH3 <- capture.output(summary(sub_a_tardi_sample_data_pH3$Shannon),file=NULL,append=FALSE)
sum_a_tardi_shannon_pH3 <- as.data.frame(sum_a_tardi_shannon_pH3)
sum_a_tardi_shannon_pH3$group <- "Tardigrades"
sum_a_tardi_shannon_pH3$diversity <- "Shannon"
sum_a_tardi_shannon_pH3$class <- "pH 7.5+"
sum_a_tardi_shannon_pH3 <-sum_a_tardi_shannon_pH3[-1,]

sub_a_nema_sample_data_pH1 <- a_nema_sample_data %>% filter(pH_grouped == "Acidic")
sum_a_nema_richness_pH1 <- capture.output(summary(sub_a_nema_sample_data_pH1$Richness),file=NULL,append=FALSE)
sum_a_nema_richness_pH1 <- as.data.frame(sum_a_nema_richness_pH1)
sum_a_nema_richness_pH1$group <- "Nematodes"
sum_a_nema_richness_pH1$diversity <- "Richness"
sum_a_nema_richness_pH1$class <- "pH 3.3-6.5"
sum_a_nema_richness_pH1 <-sum_a_nema_richness_pH1[-1,]

sum_a_nema_shannon_pH1 <- capture.output(summary(sub_a_nema_sample_data_pH1$Shannon),file=NULL,append=FALSE)
sum_a_nema_shannon_pH1 <- as.data.frame(sum_a_nema_shannon_pH1)
sum_a_nema_shannon_pH1$group <- "Nematodes"
sum_a_nema_shannon_pH1$diversity <- "Shannon"
sum_a_nema_shannon_pH1$class <- "pH 3.3-6.5"
sum_a_nema_shannon_pH1 <-sum_a_nema_shannon_pH1[-1,]

sub_a_nema_sample_data_pH2 <- a_nema_sample_data %>% filter(pH_grouped == "Neutral")
sum_a_nema_richness_pH2 <- capture.output(summary(sub_a_nema_sample_data_pH2$Richness),file=NULL,append=FALSE)
sum_a_nema_richness_pH2 <- as.data.frame(sum_a_nema_richness_pH2)
sum_a_nema_richness_pH2$group <- "Nematodes"
sum_a_nema_richness_pH2$diversity <- "Richness"
sum_a_nema_richness_pH2$class <- "pH 6.5-7.5"
sum_a_nema_richness_pH2 <-sum_a_nema_richness_pH2[-1,]

sum_a_nema_shannon_pH2 <- capture.output(summary(sub_a_nema_sample_data_pH2$Shannon),file=NULL,append=FALSE)
sum_a_nema_shannon_pH2 <- as.data.frame(sum_a_nema_shannon_pH2)
sum_a_nema_shannon_pH2$group <- "Nematodes"
sum_a_nema_shannon_pH2$diversity <- "Shannon"
sum_a_nema_shannon_pH2$class <- "pH 6.5-7.5"
sum_a_nema_shannon_pH2 <-sum_a_nema_shannon_pH2[-1,]

sub_a_nema_sample_data_pH3 <- a_nema_sample_data %>% filter(pH_grouped == "Alkaline")
sum_a_nema_richness_pH3 <- capture.output(summary(sub_a_nema_sample_data_pH3$Richness),file=NULL,append=FALSE)
sum_a_nema_richness_pH3 <- as.data.frame(sum_a_nema_richness_pH3)
sum_a_nema_richness_pH3$group <- "Nematodes"
sum_a_nema_richness_pH3$diversity <- "Richness"
sum_a_nema_richness_pH3$class <- "pH 7.5+"
sum_a_nema_richness_pH3 <-sum_a_nema_richness_pH3[-1,]

sum_a_nema_shannon_pH3 <- capture.output(summary(sub_a_nema_sample_data_pH3$Shannon),file=NULL,append=FALSE)
sum_a_nema_shannon_pH3 <- as.data.frame(sum_a_nema_shannon_pH3)
sum_a_nema_shannon_pH3$group <- "Nematodes"
sum_a_nema_shannon_pH3$diversity <- "Shannon"
sum_a_nema_shannon_pH3$class <- "pH 7.5+"
sum_a_nema_shannon_pH3 <-sum_a_nema_shannon_pH3[-1,]

sub_a_roti_sample_data_pH1 <- a_roti_sample_data %>% filter(pH_grouped == "Acidic")
sum_a_roti_richness_pH1 <- capture.output(summary(sub_a_roti_sample_data_pH1$Richness),file=NULL,append=FALSE)
sum_a_roti_richness_pH1 <- as.data.frame(sum_a_roti_richness_pH1)
sum_a_roti_richness_pH1$group <- "Rotifers"
sum_a_roti_richness_pH1$diversity <- "Richness"
sum_a_roti_richness_pH1$class <- "pH 3.3-6.5"
sum_a_roti_richness_pH1 <-sum_a_roti_richness_pH1[-1,]

sum_a_roti_shannon_pH1 <- capture.output(summary(sub_a_roti_sample_data_pH1$Shannon),file=NULL,append=FALSE)
sum_a_roti_shannon_pH1 <- as.data.frame(sum_a_roti_shannon_pH1)
sum_a_roti_shannon_pH1$group <- "Rotifers"
sum_a_roti_shannon_pH1$diversity <- "Shannon"
sum_a_roti_shannon_pH1$class <- "pH 3.3-6.5"
sum_a_roti_shannon_pH1 <-sum_a_roti_shannon_pH1[-1,]

sub_a_roti_sample_data_pH2 <- a_roti_sample_data %>% filter(pH_grouped == "Neutral")
sum_a_roti_richness_pH2 <- capture.output(summary(sub_a_roti_sample_data_pH2$Richness),file=NULL,append=FALSE)
sum_a_roti_richness_pH2 <- as.data.frame(sum_a_roti_richness_pH2)
sum_a_roti_richness_pH2$group <- "Rotifers"
sum_a_roti_richness_pH2$diversity <- "Richness"
sum_a_roti_richness_pH2$class <- "pH 6.5-7.5"
sum_a_roti_richness_pH2 <-sum_a_roti_richness_pH2[-1,]

sum_a_roti_shannon_pH2 <- capture.output(summary(sub_a_roti_sample_data_pH2$Shannon),file=NULL,append=FALSE)
sum_a_roti_shannon_pH2 <- as.data.frame(sum_a_roti_shannon_pH2)
sum_a_roti_shannon_pH2$group <- "Rotifers"
sum_a_roti_shannon_pH2$diversity <- "Shannon"
sum_a_roti_shannon_pH2$class <- "pH 6.5-7.5"
sum_a_roti_shannon_pH2 <-sum_a_roti_shannon_pH2[-1,]

sub_a_roti_sample_data_pH3 <- a_roti_sample_data %>% filter(pH_grouped == "Alkaline")
sum_a_roti_richness_pH3 <- capture.output(summary(sub_a_roti_sample_data_pH3$Richness),file=NULL,append=FALSE)
sum_a_roti_richness_pH3 <- as.data.frame(sum_a_roti_richness_pH3)
sum_a_roti_richness_pH3$group <- "Rotifers"
sum_a_roti_richness_pH3$diversity <- "Richness"
sum_a_roti_richness_pH3$class <- "pH 7.5+"
sum_a_roti_richness_pH3 <-sum_a_roti_richness_pH3[-1,]

sum_a_roti_shannon_pH3 <- capture.output(summary(sub_a_roti_sample_data_pH3$Shannon),file=NULL,append=FALSE)
sum_a_roti_shannon_pH3 <- as.data.frame(sum_a_roti_shannon_pH3)
sum_a_roti_shannon_pH3$group <- "Rotifers"
sum_a_roti_shannon_pH3$diversity <- "Shannon"
sum_a_roti_shannon_pH3$class <- "pH 7.5+"
sum_a_roti_shannon_pH3 <-sum_a_roti_shannon_pH3[-1,]

sub_p_sample_data_pH1 <- p_sample_data %>% filter(pH_grouped == "Acidic")
sum_p_richness_pH1 <- capture.output(summary(sub_p_sample_data_pH1$Richness),file=NULL,append=FALSE)
sum_p_richness_pH1 <- as.data.frame(sum_p_richness_pH1)
sum_p_richness_pH1$group <- "Protists"
sum_p_richness_pH1$diversity <- "Richness"
sum_p_richness_pH1$class <- "pH 3.3-6.5"
sum_p_richness_pH1 <-sum_p_richness_pH1[-1,]

sum_p_shannon_pH1 <- capture.output(summary(sub_p_sample_data_pH1$Shannon),file=NULL,append=FALSE)
sum_p_shannon_pH1 <- as.data.frame(sum_p_shannon_pH1)
sum_p_shannon_pH1$group <- "Protists"
sum_p_shannon_pH1$diversity <- "Shannon"
sum_p_shannon_pH1$class <- "pH 3.3-6.5"
sum_p_shannon_pH1 <-sum_p_shannon_pH1[-1,]

sub_p_sample_data_pH2 <- p_sample_data %>% filter(pH_grouped == "Neutral")
sum_p_richness_pH2 <- capture.output(summary(sub_p_sample_data_pH2$Richness),file=NULL,append=FALSE)
sum_p_richness_pH2 <- as.data.frame(sum_p_richness_pH2)
sum_p_richness_pH2$group <- "Protists"
sum_p_richness_pH2$diversity <- "Richness"
sum_p_richness_pH2$class <- "pH 6.5-7.5"
sum_p_richness_pH2 <-sum_p_richness_pH2[-1,]

sum_p_shannon_pH2 <-capture.output(summary(sub_p_sample_data_pH2$Shannon),file=NULL,append=FALSE)
sum_p_shannon_pH2 <- as.data.frame(sum_p_shannon_pH2)
sum_p_shannon_pH2$group <- "Protists"
sum_p_shannon_pH2$diversity <- "Shannon"
sum_p_shannon_pH2$class <- "pH 6.5-7.5"
sum_p_shannon_pH2 <-sum_p_shannon_pH2[-1,]

sub_p_sample_data_pH3 <- p_sample_data %>% filter(pH_grouped == "Alkaline")
sum_p_richness_pH3 <- capture.output(summary(sub_p_sample_data_pH3$Richness),file=NULL,append=FALSE)
sum_p_richness_pH3 <- as.data.frame(sum_p_richness_pH3)
sum_p_richness_pH3$group <- "Protists"
sum_p_richness_pH3$diversity <- "Richness"
sum_p_richness_pH3$class <- "pH 7.5+"
sum_p_richness_pH3 <-sum_p_richness_pH3[-1,]

sum_p_shannon_pH3 <-capture.output(summary(sub_p_sample_data_pH3$Shannon),file=NULL,append=FALSE)
sum_p_shannon_pH3 <- as.data.frame(sum_p_shannon_pH3)
sum_p_shannon_pH3$group <- "Protists"
sum_p_shannon_pH3$diversity <- "Shannon"
sum_p_shannon_pH3$class <- "pH 7.5+"
sum_p_shannon_pH3 <-sum_p_shannon_pH3[-1,]

sub_f_sample_data_pH1 <- f_sample_data %>% filter(pH_grouped == "Acidic")
sum_f_richness_pH1 <- capture.output(summary(sub_f_sample_data_pH1$Richness),file=NULL,append=FALSE)
sum_f_richness_pH1 <- as.data.frame(sum_f_richness_pH1)
sum_f_richness_pH1$group <- "Fungi"
sum_f_richness_pH1$diversity <- "Richness"
sum_f_richness_pH1$class <- "pH 3.3-6.5"
sum_f_richness_pH1 <-sum_f_richness_pH1[-1,]

sum_f_shannon_pH1 <- capture.output(summary(sub_f_sample_data_pH1$Shannon),file=NULL,append=FALSE)
sum_f_shannon_pH1 <- as.data.frame(sum_f_shannon_pH1)
sum_f_shannon_pH1$group <- "Fungi"
sum_f_shannon_pH1$diversity <- "Shannon"
sum_f_shannon_pH1$class <- "pH 3.3-6.5"
sum_f_shannon_pH1 <-sum_f_shannon_pH1[-1,]

sub_f_sample_data_pH2 <- f_sample_data %>% filter(pH_grouped == "Neutral")
sum_f_richness_pH2 <- capture.output(summary(sub_f_sample_data_pH2$Richness),file=NULL,append=FALSE)
sum_f_richness_pH2 <- as.data.frame(sum_f_richness_pH2)
sum_f_richness_pH2$group <- "Fungi"
sum_f_richness_pH2$diversity <- "Richness"
sum_f_richness_pH2$class <- "pH 6.5-7.5"
sum_f_richness_pH2 <-sum_f_richness_pH2[-1,]

sum_f_shannon_pH2 <-capture.output(summary(sub_f_sample_data_pH2$Shannon),file=NULL,append=FALSE)
sum_f_shannon_pH2 <- as.data.frame(sum_f_shannon_pH2)
sum_f_shannon_pH2$group <- "Fungi"
sum_f_shannon_pH2$diversity <- "Shannon"
sum_f_shannon_pH2$class <- "pH 6.5-7.5"
sum_f_shannon_pH2 <-sum_f_shannon_pH2[-1,]

sub_f_sample_data_pH3 <- f_sample_data %>% filter(pH_grouped == "Alkaline")
sum_f_richness_pH3 <- capture.output(summary(sub_f_sample_data_pH3$Richness),file=NULL,append=FALSE)
sum_f_richness_pH3 <- as.data.frame(sum_f_richness_pH3)
sum_f_richness_pH3$group <- "Fungi"
sum_f_richness_pH3$diversity <- "Richness"
sum_f_richness_pH3$class <- "pH 7.5+"
sum_f_richness_pH3 <-sum_f_richness_pH3[-1,]

sum_f_shannon_pH3 <- capture.output(summary(sub_f_sample_data_pH3$Shannon),file=NULL,append=FALSE)
sum_f_shannon_pH3 <- as.data.frame(sum_f_shannon_pH3)
sum_f_shannon_pH3$group <- "Fungi"
sum_f_shannon_pH3$diversity <- "Shannon"
sum_f_shannon_pH3$class <- "pH 7.5+"
sum_f_shannon_pH3 <-sum_f_shannon_pH3[-1,]


#Arthropods

sub_a_arthro_sample_data_season1 <- a_arthro_sample_data %>% filter(Sample_season == "Spring")
sum_a_arthro_richness_season1 <- capture.output(summary(sub_a_arthro_sample_data_season1$Richness),file=NULL,append=FALSE)
sum_a_arthro_richness_season1 <- as.data.frame(sum_a_arthro_richness_season1)
sum_a_arthro_richness_season1$group <- "Arthropods"
sum_a_arthro_richness_season1$diversity <- "Richness"
sum_a_arthro_richness_season1$class <- "Season - Spring"
sum_a_arthro_richness_season1 <-sum_a_arthro_richness_season1[-1,]

sum_a_arthro_shannon_season1 <- capture.output(summary(sub_a_arthro_sample_data_season1$Shannon),file=NULL,append=FALSE)
sum_a_arthro_shannon_season1 <- as.data.frame(sum_a_arthro_shannon_season1)
sum_a_arthro_shannon_season1$group <- "Arthropods"
sum_a_arthro_shannon_season1$diversity <- "Shannon"
sum_a_arthro_shannon_season1$class <- "Season - Spring"
sum_a_arthro_shannon_season1 <-sum_a_arthro_shannon_season1[-1,]

sub_a_arthro_sample_data_season2 <- a_arthro_sample_data %>% filter(Sample_season == "Summer")
sum_a_arthro_richness_season2 <- capture.output(summary(sub_a_arthro_sample_data_season2$Richness),file=NULL,append=FALSE)
sum_a_arthro_richness_season2 <- as.data.frame(sum_a_arthro_richness_season2)
sum_a_arthro_richness_season2$group <- "Arthropods"
sum_a_arthro_richness_season2$diversity <- "Richness"
sum_a_arthro_richness_season2$class <- "Season - Summer"
sum_a_arthro_richness_season2 <-sum_a_arthro_richness_season2[-1,]

sum_a_arthro_shannon_season2 <- capture.output(summary(sub_a_arthro_sample_data_season2$Shannon),file=NULL,append=FALSE)
sum_a_arthro_shannon_season2 <- as.data.frame(sum_a_arthro_shannon_season2)
sum_a_arthro_shannon_season2$group <- "Arthropods"
sum_a_arthro_shannon_season2$diversity <- "Shannon"
sum_a_arthro_shannon_season2$class <- "Season - Summer"
sum_a_arthro_shannon_season2 <-sum_a_arthro_shannon_season2[-1,]

sub_a_arthro_sample_data_season3 <- a_arthro_sample_data %>% filter(Sample_season == "Autumn")
sum_a_arthro_richness_season3 <- capture.output(summary(sub_a_arthro_sample_data_season3$Richness),file=NULL,append=FALSE)
sum_a_arthro_richness_season3 <- as.data.frame(sum_a_arthro_richness_season3)
sum_a_arthro_richness_season3$group <- "Arthropods"
sum_a_arthro_richness_season3$diversity <- "Richness"
sum_a_arthro_richness_season3$class <- "Season - Autumn"
sum_a_arthro_richness_season3 <-sum_a_arthro_richness_season3[-1,]

sum_a_arthro_shannon_season3 <- capture.output(summary(sub_a_arthro_sample_data_season3$Shannon),file=NULL,append=FALSE)
sum_a_arthro_shannon_season3 <- as.data.frame(sum_a_arthro_shannon_season3)
sum_a_arthro_shannon_season3$group <- "Arthropods"
sum_a_arthro_shannon_season3$diversity <- "Shannon"
sum_a_arthro_shannon_season3$class <- "Season - Autumn"
sum_a_arthro_shannon_season3 <-sum_a_arthro_shannon_season3[-1,]

sub_a_arthro_sample_data_depth1 <- a_arthro_sample_data %>% filter(depth_grouped == "Shallow")
sum_a_arthro_richness_depth1 <- capture.output(summary(sub_a_arthro_sample_data_depth1$Richness),file=NULL,append=FALSE)
sum_a_arthro_richness_depth1 <- as.data.frame(sum_a_arthro_richness_depth1)
sum_a_arthro_richness_depth1$group <- "Arthropods"
sum_a_arthro_richness_depth1$diversity <- "Richness"
sum_a_arthro_richness_depth1$class <- "Depth <1m"
sum_a_arthro_richness_depth1 <-sum_a_arthro_richness_depth1[-1,]

sum_a_arthro_shannon_depth1 <- capture.output(summary(sub_a_arthro_sample_data_depth1$Shannon),file=NULL,append=FALSE)
sum_a_arthro_shannon_depth1 <- as.data.frame(sum_a_arthro_shannon_depth1)
sum_a_arthro_shannon_depth1$group <- "Arthropods"
sum_a_arthro_shannon_depth1$diversity <- "Shannon"
sum_a_arthro_shannon_depth1$class <- "Depth >1m"
sum_a_arthro_shannon_depth1 <-sum_a_arthro_shannon_depth1[-1,]

sub_a_arthro_sample_data_depth2 <- a_arthro_sample_data %>% filter(depth_grouped == "Moderately deep")
sum_a_arthro_richness_depth2 <- capture.output(summary(sub_a_arthro_sample_data_depth2$Richness),file=NULL,append=FALSE)
sum_a_arthro_richness_depth2 <- as.data.frame(sum_a_arthro_richness_depth2)
sum_a_arthro_richness_depth2$group <- "Arthropods"
sum_a_arthro_richness_depth2$diversity <- "Richness"
sum_a_arthro_richness_depth2$class <- "Depth 1-2m"
sum_a_arthro_richness_depth2 <-sum_a_arthro_richness_depth2[-1,]

sum_a_arthro_shannon_depth2 <- capture.output(summary(sub_a_arthro_sample_data_depth2$Shannon),file=NULL,append=FALSE)
sum_a_arthro_shannon_depth2 <- as.data.frame(sum_a_arthro_shannon_depth2)
sum_a_arthro_shannon_depth2$group <- "Arthropods"
sum_a_arthro_shannon_depth2$diversity <- "Shannon"
sum_a_arthro_shannon_depth2$class <- "Depth 1-2m"
sum_a_arthro_shannon_depth2 <-sum_a_arthro_shannon_depth2[-1,]

sub_a_arthro_sample_data_depth3 <- a_arthro_sample_data %>% filter(depth_grouped == "Deep")
sum_a_arthro_richness_depth3 <- capture.output(summary(sub_a_arthro_sample_data_depth3$Richness),file=NULL,append=FALSE)
sum_a_arthro_richness_depth3 <- as.data.frame(sum_a_arthro_richness_depth3)
sum_a_arthro_richness_depth3$group <- "Arthropods"
sum_a_arthro_richness_depth3$diversity <- "Richness"
sum_a_arthro_richness_depth3$class <- "Depth >2m"
sum_a_arthro_richness_depth3 <-sum_a_arthro_richness_depth3[-1,]

sum_a_arthro_shannon_depth3 <- capture.output(summary(sub_a_arthro_sample_data_depth3$Shannon),file=NULL,append=FALSE)
sum_a_arthro_shannon_depth3 <- as.data.frame(sum_a_arthro_shannon_depth3)
sum_a_arthro_shannon_depth3$group <- "Arthropods"
sum_a_arthro_shannon_depth3$diversity <- "Shannon"
sum_a_arthro_shannon_depth3$class <- "Depth >2m"
sum_a_arthro_shannon_depth3 <-sum_a_arthro_shannon_depth3[-1,]

sub_a_arthro_sample_data_pH1 <- a_arthro_sample_data %>% filter(pH_grouped == "Acidic")
sum_a_arthro_richness_pH1 <- capture.output(summary(sub_a_arthro_sample_data_pH1$Richness),file=NULL,append=FALSE)
sum_a_arthro_richness_pH1 <- as.data.frame(sum_a_arthro_richness_pH1)
sum_a_arthro_richness_pH1$group <- "Arthropods"
sum_a_arthro_richness_pH1$diversity <- "Richness"
sum_a_arthro_richness_pH1$class <- "pH 3.3-6.5"
sum_a_arthro_richness_pH1 <-sum_a_arthro_richness_pH1[-1,]

sum_a_arthro_shannon_pH1 <- capture.output(summary(sub_a_arthro_sample_data_pH1$Shannon),file=NULL,append=FALSE)
sum_a_arthro_shannon_pH1 <- as.data.frame(sum_a_arthro_shannon_pH1)
sum_a_arthro_shannon_pH1$group <- "Arthropods"
sum_a_arthro_shannon_pH1$diversity <- "Shannon"
sum_a_arthro_shannon_pH1$class <- "pH 3.3-6.5"
sum_a_arthro_shannon_pH1 <-sum_a_arthro_shannon_pH1[-1,]

sub_a_arthro_sample_data_pH2 <- a_arthro_sample_data %>% filter(pH_grouped == "Neutral")
sum_a_arthro_richness_pH2 <- capture.output(summary(sub_a_arthro_sample_data_pH2$Richness),file=NULL,append=FALSE)
sum_a_arthro_richness_pH2 <- as.data.frame(sum_a_arthro_richness_pH2)
sum_a_arthro_richness_pH2$group <- "Arthropods"
sum_a_arthro_richness_pH2$diversity <- "Richness"
sum_a_arthro_richness_pH2$class <- "pH 6.5-7.5"
sum_a_arthro_richness_pH2 <-sum_a_arthro_richness_pH2[-1,]

sum_a_arthro_shannon_pH2 <- capture.output(summary(sub_a_arthro_sample_data_pH2$Shannon),file=NULL,append=FALSE)
sum_a_arthro_shannon_pH2 <- as.data.frame(sum_a_arthro_shannon_pH2)
sum_a_arthro_shannon_pH2$group <- "Arthropods"
sum_a_arthro_shannon_pH2$diversity <- "Shannon"
sum_a_arthro_shannon_pH2$class <- "pH 6.5-7.5"
sum_a_arthro_shannon_pH2 <-sum_a_arthro_shannon_pH2[-1,]

sub_a_arthro_sample_data_pH3 <- a_arthro_sample_data %>% filter(pH_grouped == "Alkaline")
sum_a_arthro_richness_pH3 <- capture.output(summary(sub_a_arthro_sample_data_pH3$Richness),file=NULL,append=FALSE)
sum_a_arthro_richness_pH3 <- as.data.frame(sum_a_arthro_richness_pH3)
sum_a_arthro_richness_pH3$group <- "Arthropods"
sum_a_arthro_richness_pH3$diversity <- "Richness"
sum_a_arthro_richness_pH3$class <- "pH 7.5+"
sum_a_arthro_richness_pH3 <-sum_a_arthro_richness_pH3[-1,]

sum_a_arthro_shannon_pH3 <- capture.output(summary(sub_a_arthro_sample_data_pH3$Shannon),file=NULL,append=FALSE)
sum_a_arthro_shannon_pH3 <- as.data.frame(sum_a_arthro_shannon_pH3)
sum_a_arthro_shannon_pH3$group <- "Arthropods"
sum_a_arthro_shannon_pH3$diversity <- "Shannon"
sum_a_arthro_shannon_pH3$class <- "pH 7.5+"
sum_a_arthro_shannon_pH3 <-sum_a_arthro_shannon_pH3[-1,]

sub_a_arthro_sample_data_erosion1 <- a_arthro_sample_data %>% filter(Erosion_grouped == "Low")
sum_a_arthro_richness_erosion1 <- capture.output(summary(sub_a_arthro_sample_data_erosion1$Richness),file=NULL,append=FALSE)
sum_a_arthro_richness_erosion1 <- as.data.frame(sum_a_arthro_richness_erosion1)
sum_a_arthro_richness_erosion1$group <- "Arthropods"
sum_a_arthro_richness_erosion1$diversity <- "Richness"
sum_a_arthro_richness_erosion1$class <- "Erosion 0-1"
sum_a_arthro_richness_erosion1 <-sum_a_arthro_richness_erosion1[-1,]

sum_a_arthro_shannon_erosion1 <- capture.output(summary(sub_a_arthro_sample_data_erosion1$Shannon),file=NULL,append=FALSE)
sum_a_arthro_shannon_erosion1 <- as.data.frame(sum_a_arthro_shannon_erosion1)
sum_a_arthro_shannon_erosion1$group <- "Arthropods"
sum_a_arthro_shannon_erosion1$diversity <- "Shannon"
sum_a_arthro_shannon_erosion1$class <- "Erosion 0-1"
sum_a_arthro_shannon_erosion1 <-sum_a_arthro_shannon_erosion1[-1,]

sub_a_arthro_sample_data_erosion2 <- a_arthro_sample_data %>% filter(Erosion_grouped == "Medium high")
sum_a_arthro_richness_erosion2 <- capture.output(summary(sub_a_arthro_sample_data_erosion2$Richness),file=NULL,append=FALSE)
sum_a_arthro_richness_erosion2 <- as.data.frame(sum_a_arthro_richness_erosion2)
sum_a_arthro_richness_erosion2$group <- "Arthropods"
sum_a_arthro_richness_erosion2$diversity <- "Richness"
sum_a_arthro_richness_erosion2$class <- "Erosion 1-5"
sum_a_arthro_richness_erosion2 <-sum_a_arthro_richness_erosion2[-1,]

sum_a_arthro_shannon_erosion2 <- capture.output(summary(sub_a_arthro_sample_data_erosion2$Shannon),file=NULL,append=FALSE)
sum_a_arthro_shannon_erosion2 <- as.data.frame(sum_a_arthro_shannon_erosion2)
sum_a_arthro_shannon_erosion2$group <- "Arthropods"
sum_a_arthro_shannon_erosion2$diversity <- "Shannon"
sum_a_arthro_shannon_erosion2$class <- "Erosion 1-5"
sum_a_arthro_shannon_erosion2 <-sum_a_arthro_shannon_erosion2[-1,]

sub_a_arthro_sample_data_erosion3 <- a_arthro_sample_data %>% filter(Erosion_grouped == "High")
sum_a_arthro_richness_erosion3 <- capture.output(summary(sub_a_arthro_sample_data_erosion3$Richness),file=NULL,append=FALSE)
sum_a_arthro_richness_erosion3 <- as.data.frame(sum_a_arthro_richness_erosion3)
sum_a_arthro_richness_erosion3$group <- "Arthropods"
sum_a_arthro_richness_erosion3$diversity <- "Richness"
sum_a_arthro_richness_erosion3$class <- "Erosion 5-10"
sum_a_arthro_richness_erosion3 <-sum_a_arthro_richness_erosion3[-1,]

sum_a_arthro_shannon_erosion3 <- capture.output(summary(sub_a_arthro_sample_data_erosion3$Shannon),file=NULL,append=FALSE)
sum_a_arthro_shannon_erosion3 <- as.data.frame(sum_a_arthro_shannon_erosion3)
sum_a_arthro_shannon_erosion3$group <- "Arthropods"
sum_a_arthro_shannon_erosion3$diversity <- "Shannon"
sum_a_arthro_shannon_erosion3$class <- "Erosion 5-10"
sum_a_arthro_shannon_erosion3 <-sum_a_arthro_shannon_erosion3[-1,]

sub_a_arthro_sample_data_erosion4 <- a_arthro_sample_data %>% filter(Erosion_grouped == "Very high")
sum_a_arthro_richness_erosion4 <- capture.output(summary(sub_a_arthro_sample_data_erosion4$Richness),file=NULL,append=FALSE)
sum_a_arthro_richness_erosion4 <- as.data.frame(sum_a_arthro_richness_erosion4)
sum_a_arthro_richness_erosion4$group <- "Arthropods"
sum_a_arthro_richness_erosion4$diversity <- "Richness"
sum_a_arthro_richness_erosion4$class <- "Erosion 10+"
sum_a_arthro_richness_erosion4 <-sum_a_arthro_richness_erosion4[-1,]

sum_a_arthro_shannon_erosion4 <- capture.output(summary(sub_a_arthro_sample_data_erosion4$Shannon),file=NULL,append=FALSE)
sum_a_arthro_shannon_erosion4 <- as.data.frame(sum_a_arthro_shannon_erosion4)
sum_a_arthro_shannon_erosion4$group <- "Arthropods"
sum_a_arthro_shannon_erosion4$diversity <- "Shannon"
sum_a_arthro_shannon_erosion4$class <- "Erosion 10+"
sum_a_arthro_shannon_erosion4 <-sum_a_arthro_shannon_erosion4[-1,]


sub_arthro_sample_data_CL1 <- a_arthro_sample_data %>% filter(LU2 == "CL1")
sum_arthro_sample_data_richness_CL1 <- capture.output(summary(sub_arthro_sample_data_CL1$Richness),file=NULL,append=FALSE)
sum_arthro_sample_data_richness_CL1 <- as.data.frame(sum_arthro_sample_data_richness_CL1)
sum_arthro_sample_data_richness_CL1$group <- "Arthropods"
sum_arthro_sample_data_richness_CL1$diversity <- "Richness"
sum_arthro_sample_data_richness_CL1$class <- "Land Cover CL1"
sum_arthro_sample_data_richness_CL1 <- sum_arthro_sample_data_richness_CL1[-1,]

sum_arthro_sample_data_shannon_CL1 <- capture.output(summary(sub_arthro_sample_data_CL1$Shannon),file=NULL,append=FALSE)
sum_arthro_sample_data_shannon_CL1 <- as.data.frame(sum_arthro_sample_data_shannon_CL1)
sum_arthro_sample_data_shannon_CL1$group <- "Arthropods"
sum_arthro_sample_data_shannon_CL1$diversity <- "Shannon"
sum_arthro_sample_data_shannon_CL1$class <- "Land Cover CL1"
sum_arthro_sample_data_shannon_CL1 <-sum_arthro_sample_data_shannon_CL1[-1,]

sub_arthro_sample_data_CL2 <- a_arthro_sample_data %>% filter(LU2 == "CL2")
sum_arthro_sample_data_richness_CL2 <- capture.output(summary(sub_arthro_sample_data_CL2$Richness),file=NULL,append=FALSE)
sum_arthro_sample_data_richness_CL2 <- as.data.frame(sum_arthro_sample_data_richness_CL2)
sum_arthro_sample_data_richness_CL2$group <- "Arthropods"
sum_arthro_sample_data_richness_CL2$diversity <- "Richness"
sum_arthro_sample_data_richness_CL2$class <- "Land Cover CL2"
sum_arthro_sample_data_richness_CL2 <- sum_arthro_sample_data_richness_CL2[-1,]

sum_arthro_sample_data_shannon_CL2 <- capture.output(summary(sub_arthro_sample_data_CL2$Shannon),file=NULL,append=FALSE)
sum_arthro_sample_data_shannon_CL2 <- as.data.frame(sum_arthro_sample_data_shannon_CL2)
sum_arthro_sample_data_shannon_CL2$group <- "Arthropods"
sum_arthro_sample_data_shannon_CL2$diversity <- "Shannon"
sum_arthro_sample_data_shannon_CL2$class <- "Land Cover CL2"
sum_arthro_sample_data_shannon_CL2 <-sum_arthro_sample_data_shannon_CL2[-1,]

sub_arthro_sample_data_CL3 <- a_arthro_sample_data %>% filter(LU2 == "CL3")
sum_arthro_sample_data_richness_CL3 <- capture.output(summary(sub_arthro_sample_data_CL3$Richness),file=NULL,append=FALSE)
sum_arthro_sample_data_richness_CL3 <- as.data.frame(sum_arthro_sample_data_richness_CL3)
sum_arthro_sample_data_richness_CL3$group <- "Arthropods"
sum_arthro_sample_data_richness_CL3$diversity <- "Richness"
sum_arthro_sample_data_richness_CL3$class <- "Land Cover CL3"
sum_arthro_sample_data_richness_CL3 <- sum_arthro_sample_data_richness_CL3[-1,]

sum_arthro_sample_data_shannon_CL3 <- capture.output(summary(sub_arthro_sample_data_CL3$Shannon),file=NULL,append=FALSE)
sum_arthro_sample_data_shannon_CL3 <- as.data.frame(sum_arthro_sample_data_shannon_CL3)
sum_arthro_sample_data_shannon_CL3$group <- "Arthropods"
sum_arthro_sample_data_shannon_CL3$diversity <- "Shannon"
sum_arthro_sample_data_shannon_CL3$class <- "Land Cover CL3"
sum_arthro_sample_data_shannon_CL3 <-sum_arthro_sample_data_shannon_CL3[-1,]

sub_arthro_sample_data_CL4 <- a_arthro_sample_data %>% filter(LU2 == "CL4")
sum_arthro_sample_data_richness_CL4 <- capture.output(summary(sub_arthro_sample_data_CL4$Richness),file=NULL,append=FALSE)
sum_arthro_sample_data_richness_CL4 <- as.data.frame(sum_arthro_sample_data_richness_CL4)
sum_arthro_sample_data_richness_CL4$group <- "Arthropods"
sum_arthro_sample_data_richness_CL4$diversity <- "Richness"
sum_arthro_sample_data_richness_CL4$class <- "Land Cover CL4"
sum_arthro_sample_data_richness_CL4 <- sum_arthro_sample_data_richness_CL4[-1,]

sum_arthro_sample_data_shannon_CL4 <- capture.output(summary(sub_arthro_sample_data_CL4$Shannon),file=NULL,append=FALSE)
sum_arthro_sample_data_shannon_CL4 <- as.data.frame(sum_arthro_sample_data_shannon_CL4)
sum_arthro_sample_data_shannon_CL4$group <- "Arthropods"
sum_arthro_sample_data_shannon_CL4$diversity <- "Shannon"
sum_arthro_sample_data_shannon_CL4$class <- "Land Cover CL4"
sum_arthro_sample_data_shannon_CL4 <-sum_arthro_sample_data_shannon_CL4[-1,]

sub_arthro_sample_data_GL1 <- a_arthro_sample_data %>% filter(LU2 == "GL1")
sum_arthro_sample_data_richness_GL1 <- capture.output(summary(sub_arthro_sample_data_GL1$Richness),file=NULL,append=FALSE)
sum_arthro_sample_data_richness_GL1 <- as.data.frame(sum_arthro_sample_data_richness_GL1)
sum_arthro_sample_data_richness_GL1$group <- "Arthropods"
sum_arthro_sample_data_richness_GL1$diversity <- "Richness"
sum_arthro_sample_data_richness_GL1$class <- "Land Cover GL1"
sum_arthro_sample_data_richness_GL1 <- sum_arthro_sample_data_richness_GL1[-1,]

sum_arthro_sample_data_shannon_GL1 <- capture.output(summary(sub_arthro_sample_data_GL1$Shannon),file=NULL,append=FALSE)
sum_arthro_sample_data_shannon_GL1 <- as.data.frame(sum_arthro_sample_data_shannon_GL1)
sum_arthro_sample_data_shannon_GL1$group <- "Arthropods"
sum_arthro_sample_data_shannon_GL1$diversity <- "Shannon"
sum_arthro_sample_data_shannon_GL1$class <- "Land Cover GL1"
sum_arthro_sample_data_shannon_GL1 <-sum_arthro_sample_data_shannon_GL1[-1,]

sub_arthro_sample_data_GL2 <- a_arthro_sample_data %>% filter(LU2 == "GL2")
sum_arthro_sample_data_richness_GL2 <- capture.output(summary(sub_arthro_sample_data_GL2$Richness),file=NULL,append=FALSE)
sum_arthro_sample_data_richness_GL2 <- as.data.frame(sum_arthro_sample_data_richness_GL2)
sum_arthro_sample_data_richness_GL2$group <- "Arthropods"
sum_arthro_sample_data_richness_GL2$diversity <- "Richness"
sum_arthro_sample_data_richness_GL2$class <- "Land Cover GL2"
sum_arthro_sample_data_richness_GL2 <- sum_arthro_sample_data_richness_GL2[-1,]

sum_arthro_sample_data_shannon_GL2 <- capture.output(summary(sub_arthro_sample_data_GL2$Shannon),file=NULL,append=FALSE)
sum_arthro_sample_data_shannon_GL2 <- as.data.frame(sum_arthro_sample_data_shannon_GL2)
sum_arthro_sample_data_shannon_GL2$group <- "Arthropods"
sum_arthro_sample_data_shannon_GL2$diversity <- "Shannon"
sum_arthro_sample_data_shannon_GL2$class <- "Land Cover GL2"
sum_arthro_sample_data_shannon_GL2 <-sum_arthro_sample_data_shannon_GL2[-1,]

sub_arthro_sample_data_GL3 <- a_arthro_sample_data %>% filter(LU2 == "GL3")
sum_arthro_sample_data_richness_GL3 <- capture.output(summary(sub_arthro_sample_data_GL3$Richness),file=NULL,append=FALSE)
sum_arthro_sample_data_richness_GL3 <- as.data.frame(sum_arthro_sample_data_richness_GL3)
sum_arthro_sample_data_richness_GL3$group <- "Arthropods"
sum_arthro_sample_data_richness_GL3$diversity <- "Richness"
sum_arthro_sample_data_richness_GL3$class <- "Land Cover GL3"
sum_arthro_sample_data_richness_GL3 <- sum_arthro_sample_data_richness_GL3[-1,]

sum_arthro_sample_data_shannon_GL3 <- capture.output(summary(sub_arthro_sample_data_GL3$Shannon),file=NULL,append=FALSE)
sum_arthro_sample_data_shannon_GL3 <- as.data.frame(sum_arthro_sample_data_shannon_GL3)
sum_arthro_sample_data_shannon_GL3$group <- "Arthropods"
sum_arthro_sample_data_shannon_GL3$diversity <- "Shannon"
sum_arthro_sample_data_shannon_GL3$class <- "Land Cover GL3"
sum_arthro_sample_data_shannon_GL3 <-sum_arthro_sample_data_shannon_GL3[-1,]

sub_arthro_sample_data_WL1 <- a_arthro_sample_data %>% filter(LU2 == "WL1")
sum_arthro_sample_data_richness_WL1 <- capture.output(summary(sub_arthro_sample_data_WL1$Richness),file=NULL,append=FALSE)
sum_arthro_sample_data_richness_WL1 <- as.data.frame(sum_arthro_sample_data_richness_WL1)
sum_arthro_sample_data_richness_WL1$group <- "Arthropods"
sum_arthro_sample_data_richness_WL1$diversity <- "Richness"
sum_arthro_sample_data_richness_WL1$class <- "Land Cover WL1"
sum_arthro_sample_data_richness_WL1 <- sum_arthro_sample_data_richness_WL1[-1,]

sum_arthro_sample_data_shannon_WL1 <- capture.output(summary(sub_arthro_sample_data_WL1$Shannon),file=NULL,append=FALSE)
sum_arthro_sample_data_shannon_WL1 <- as.data.frame(sum_arthro_sample_data_shannon_WL1)
sum_arthro_sample_data_shannon_WL1$group <- "Arthropods"
sum_arthro_sample_data_shannon_WL1$diversity <- "Shannon"
sum_arthro_sample_data_shannon_WL1$class <- "Land Cover WL1"
sum_arthro_sample_data_shannon_WL1 <-sum_arthro_sample_data_shannon_WL1[-1,]

sub_arthro_sample_data_WL2 <- a_arthro_sample_data %>% filter(LU2 == "WL2")
sum_arthro_sample_data_richness_WL2 <- capture.output(summary(sub_arthro_sample_data_WL2$Richness),file=NULL,append=FALSE)
sum_arthro_sample_data_richness_WL2 <- as.data.frame(sum_arthro_sample_data_richness_WL2)
sum_arthro_sample_data_richness_WL2$group <- "Arthropods"
sum_arthro_sample_data_richness_WL2$diversity <- "Richness"
sum_arthro_sample_data_richness_WL2$class <- "Land Cover WL2"
sum_arthro_sample_data_richness_WL2 <- sum_arthro_sample_data_richness_WL2[-1,]

sum_arthro_sample_data_shannon_WL2 <- capture.output(summary(sub_arthro_sample_data_WL2$Shannon),file=NULL,append=FALSE)
sum_arthro_sample_data_shannon_WL2 <- as.data.frame(sum_arthro_sample_data_shannon_WL2)
sum_arthro_sample_data_shannon_WL2$group <- "Arthropods"
sum_arthro_sample_data_shannon_WL2$diversity <- "Shannon"
sum_arthro_sample_data_shannon_WL2$class <- "Land Cover WL2"
sum_arthro_sample_data_shannon_WL2 <-sum_arthro_sample_data_shannon_WL2[-1,]

sub_arthro_sample_data_WL3 <- a_arthro_sample_data %>% filter(LU2 == "WL3")
sum_arthro_sample_data_richness_WL3 <- capture.output(summary(sub_arthro_sample_data_WL3$Richness),file=NULL,append=FALSE)
sum_arthro_sample_data_richness_WL3 <- as.data.frame(sum_arthro_sample_data_richness_WL3)
sum_arthro_sample_data_richness_WL3$group <- "Arthropods"
sum_arthro_sample_data_richness_WL3$diversity <- "Richness"
sum_arthro_sample_data_richness_WL3$class <- "Land Cover WL3"
sum_arthro_sample_data_richness_WL3 <- sum_arthro_sample_data_richness_WL3[-1,]

sum_arthro_sample_data_shannon_WL3 <- capture.output(summary(sub_arthro_sample_data_WL3$Shannon),file=NULL,append=FALSE)
sum_arthro_sample_data_shannon_WL3 <- as.data.frame(sum_arthro_sample_data_shannon_WL3)
sum_arthro_sample_data_shannon_WL3$group <- "Arthropods"
sum_arthro_sample_data_shannon_WL3$diversity <- "Shannon"
sum_arthro_sample_data_shannon_WL3$class <- "Land Cover WL3"
sum_arthro_sample_data_shannon_WL3 <-sum_arthro_sample_data_shannon_WL3[-1,]

sub_arthro_sample_data_WL4 <- a_arthro_sample_data %>% filter(LU2 == "WL4")
sum_arthro_sample_data_richness_WL4 <- capture.output(summary(sub_arthro_sample_data_WL4$Richness),file=NULL,append=FALSE)
sum_arthro_sample_data_richness_WL4 <- as.data.frame(sum_arthro_sample_data_richness_WL4)
sum_arthro_sample_data_richness_WL4$group <- "Arthropods"
sum_arthro_sample_data_richness_WL4$diversity <- "Richness"
sum_arthro_sample_data_richness_WL4$class <- "Land Cover WL4"
sum_arthro_sample_data_richness_WL4 <- sum_arthro_sample_data_richness_WL4[-1,]

sum_arthro_sample_data_shannon_WL4 <- capture.output(summary(sub_arthro_sample_data_WL4$Shannon),file=NULL,append=FALSE)
sum_arthro_sample_data_shannon_WL4 <- as.data.frame(sum_arthro_sample_data_shannon_WL4)
sum_arthro_sample_data_shannon_WL4$group <- "Arthropods"
sum_arthro_sample_data_shannon_WL4$diversity <- "Shannon"
sum_arthro_sample_data_shannon_WL4$class <- "Land Cover WL4"
sum_arthro_sample_data_shannon_WL4 <-sum_arthro_sample_data_shannon_WL4[-1,]

sub_arthro_sample_data_No_C <- a_arthro_sample_data %>% filter(LU2 == "No_C")
sum_arthro_sample_data_richness_No_C <- capture.output(summary(sub_arthro_sample_data_No_C$Richness),file=NULL,append=FALSE)
sum_arthro_sample_data_richness_No_C <- as.data.frame(sum_arthro_sample_data_richness_No_C)
sum_arthro_sample_data_richness_No_C$group <- "Arthropods"
sum_arthro_sample_data_richness_No_C$diversity <- "Richness"
sum_arthro_sample_data_richness_No_C$class <- "Land Cover No_C"
sum_arthro_sample_data_richness_No_C <- sum_arthro_sample_data_richness_No_C[-1,]

sum_arthro_sample_data_shannon_No_C <- capture.output(summary(sub_arthro_sample_data_No_C$Shannon),file=NULL,append=FALSE)
sum_arthro_sample_data_shannon_No_C <- as.data.frame(sum_arthro_sample_data_shannon_No_C)
sum_arthro_sample_data_shannon_No_C$group <- "Arthropods"
sum_arthro_sample_data_shannon_No_C$diversity <- "Shannon"
sum_arthro_sample_data_shannon_No_C$class <- "Land Cover No_C"
sum_arthro_sample_data_shannon_No_C <-sum_arthro_sample_data_shannon_No_C[-1,]



sub_a_arthro_sample_data_CL_annual <- a_arthro_sample_data %>% filter(LC1_2018 == "CL_annual")
sum_a_arthro_richness_CL_annual <- capture.output(summary(sub_a_arthro_sample_data_CL_annual$Richness),file=NULL,append=FALSE)
sum_a_arthro_richness_CL_annual <- as.data.frame(sum_a_arthro_richness_CL_annual)
sum_a_arthro_richness_CL_annual$group <- "Arthropods"
sum_a_arthro_richness_CL_annual$diversity <- "Richness"
sum_a_arthro_richness_CL_annual$class <- "Land Cover CL_annual"
sum_a_arthro_richness_CL_annual <-sum_a_arthro_richness_CL_annual[-1,]

sum_a_arthro_shannon_CL_annual <- capture.output(summary(sub_a_arthro_sample_data_CL_annual$Shannon),file=NULL,append=FALSE)
sum_a_arthro_shannon_CL_annual <- as.data.frame(sum_a_arthro_shannon_CL_annual)
sum_a_arthro_shannon_CL_annual$group <- "Arthropods"
sum_a_arthro_shannon_CL_annual$diversity <- "Shannon"
sum_a_arthro_shannon_CL_annual$class <- "Land Cove CL_annualr"
sum_a_arthro_shannon_CL_annual <-sum_a_arthro_shannon_CL_annual[-1,]

sum_a_arthro_sample_data_CL_permanent <- a_arthro_sample_data %>% filter(LC1_2018 == "CL_permanent")
sum_a_arthro_richness_CL_permanent <- capture.output(summary(sum_a_arthro_sample_data_CL_permanent$Richness),file=NULL,append=FALSE)
sum_a_arthro_richness_CL_permanent <- as.data.frame(sum_a_arthro_richness_CL_permanent)
sum_a_arthro_richness_CL_permanent$group <- "Arthropods"
sum_a_arthro_richness_CL_permanent$diversity <- "Richness"
sum_a_arthro_richness_CL_permanent$class <- "Land Cover CL_permanent"
sum_a_arthro_richness_CL_permanent <-sum_a_arthro_richness_CL_permanent[-1,]

sum_a_arthro_shannon_CL_permanent <- capture.output(summary(sum_a_arthro_sample_data_CL_permanent$Shannon),file=NULL,append=FALSE)
sum_a_arthro_shannon_CL_permanent <- as.data.frame(sum_a_arthro_shannon_CL_permanent)
sum_a_arthro_shannon_CL_permanent$group <- "Arthropods"
sum_a_arthro_shannon_CL_permanent$diversity <- "Shannon"
sum_a_arthro_shannon_CL_permanent$class <- "Land Cover CL_permanent"
sum_a_arthro_shannon_CL_permanent <-sum_a_arthro_shannon_CL_permanent[-1,]

sub_a_arthro_sample_data_GL_unmanaged <- a_arthro_sample_data %>% filter(LC1_2018 == "GL_unmanaged")
sum_a_arthro_richness_GL_unmanaged <- capture.output(summary(sub_a_arthro_sample_data_GL_unmanaged$Richness),file=NULL,append=FALSE)
sum_a_arthro_richness_GL_unmanaged <- as.data.frame(sum_a_arthro_richness_GL_unmanaged)
sum_a_arthro_richness_GL_unmanaged$group <- "Arthropods"
sum_a_arthro_richness_GL_unmanaged$diversity <- "Richness"
sum_a_arthro_richness_GL_unmanaged$class <- "Land Cover GL_unmanaged"
sum_a_arthro_richness_GL_unmanaged <-sum_a_arthro_richness_GL_unmanaged[-1,]

sum_a_arthro_shannon_GL_unmanaged <- capture.output(summary(sub_a_arthro_sample_data_GL_unmanaged$Shannon),file=NULL,append=FALSE)
sum_a_arthro_shannon_GL_unmanaged <- as.data.frame(sum_a_arthro_shannon_GL_unmanaged)
sum_a_arthro_shannon_GL_unmanaged$group <- "Arthropods"
sum_a_arthro_shannon_GL_unmanaged$diversity <- "Shannon"
sum_a_arthro_shannon_GL_unmanaged$class <- "Land Cover GL_unmanaged"
sum_a_arthro_shannon_GL_unmanaged <-sum_a_arthro_shannon_GL_unmanaged[-1,]

sub_a_arthro_sample_data_GL_managed <- a_arthro_sample_data %>% filter(LC1_2018 == "GL_managed")
sum_a_arthro_richness_GL_managed <- capture.output(summary(sub_a_arthro_sample_data_GL_managed$Richness),file=NULL,append=FALSE)
sum_a_arthro_richness_GL_managed <- as.data.frame(sum_a_arthro_richness_GL_managed)
sum_a_arthro_richness_GL_managed$group <- "Arthropods"
sum_a_arthro_richness_GL_managed$diversity <- "Richness"
sum_a_arthro_richness_GL_managed$class <- "Land Cover GL_managed"
sum_a_arthro_richness_GL_managed <-sum_a_arthro_richness_GL_managed[-1,]

sum_a_arthro_shannon_GL_managed <- capture.output(summary(sub_a_arthro_sample_data_GL_managed$Shannon),file=NULL,append=FALSE)
sum_a_arthro_shannon_GL_managed <- as.data.frame(sum_a_arthro_shannon_GL_managed)
sum_a_arthro_shannon_GL_managed$group <- "Arthropods"
sum_a_arthro_shannon_GL_managed$diversity <- "Shannon"
sum_a_arthro_shannon_GL_managed$class <- "Land Cover GL_managed"
sum_a_arthro_shannon_GL_managed <-sum_a_arthro_shannon_GL_managed[-1,]

sub_a_arthro_sample_data_WL_broadleaved <- a_arthro_sample_data %>% filter(LC1_2018 == "WL_broadleaved")
sum_a_arthro_richness_WL_broadleaved <- capture.output(summary(sub_a_arthro_sample_data_WL_broadleaved$Richness),file=NULL,append=FALSE)
sum_a_arthro_richness_WL_broadleaved <- as.data.frame(sum_a_arthro_richness_WL_broadleaved)
sum_a_arthro_richness_WL_broadleaved$group <- "Arthropods"
sum_a_arthro_richness_WL_broadleaved$diversity <- "Richness"
sum_a_arthro_richness_WL_broadleaved$class <- "Land Cover WL_broadleaved"
sum_a_arthro_richness_WL_broadleaved <-sum_a_arthro_richness_WL_broadleaved[-1,]

sum_a_arthro_shannon_WL_broadleaved <- capture.output(summary(sub_a_arthro_sample_data_WL_broadleaved$Shannon),file=NULL,append=FALSE)
sum_a_arthro_shannon_WL_broadleaved <- as.data.frame(sum_a_arthro_shannon_WL_broadleaved)
sum_a_arthro_shannon_WL_broadleaved$group <- "Arthropods"
sum_a_arthro_shannon_WL_broadleaved$diversity <- "Shannon"
sum_a_arthro_shannon_WL_broadleaved$class <- "Land Cover WL_broadleaved"
sum_a_arthro_shannon_WL_broadleaved <-sum_a_arthro_shannon_WL_broadleaved[-1,]

sub_a_arthro_sample_data_WL_coniferous <- a_arthro_sample_data %>% filter(LC1_2018 == "WL_coniferous")
sum_a_arthro_richness_WL_coniferous <- capture.output(summary(sub_a_arthro_sample_data_WL_coniferous$Richness),file=NULL,append=FALSE)
sum_a_arthro_richness_WL_coniferous <- as.data.frame(sum_a_arthro_richness_WL_coniferous)
sum_a_arthro_richness_WL_coniferous$group <- "Arthropods"
sum_a_arthro_richness_WL_coniferous$diversity <- "Richness"
sum_a_arthro_richness_WL_coniferous$class <- "Land Cover WL_coniferous"
sum_a_arthro_richness_WL_coniferous <-sum_a_arthro_richness_WL_coniferous[-1,]

sum_a_arthro_shannon_WL_coniferous <- capture.output(summary(sub_a_arthro_sample_data_WL_coniferous$Shannon),file=NULL,append=FALSE)
sum_a_arthro_shannon_WL_coniferous <- as.data.frame(sum_a_arthro_shannon_WL_coniferous)
sum_a_arthro_shannon_WL_coniferous$group <- "Arthropods"
sum_a_arthro_shannon_WL_coniferous$diversity <- "Shannon"
sum_a_arthro_shannon_WL_coniferous$class <- "Land Cover WL_coniferous"
sum_a_arthro_shannon_WL_coniferous <-sum_a_arthro_shannon_WL_coniferous[-1,]

sub_a_arthro_sample_data_temp_1 <- a_arthro_sample_data %>% filter(cluster_temp == "Very low")
sum_a_arthro_sample_data_shannon_temp_1 <- capture.output(summary(sub_a_arthro_sample_data_temp_1$Shannon),file=NULL,append=FALSE)
sum_a_arthro_sample_data_shannon_temp_1 <- as.data.frame(sum_a_arthro_sample_data_shannon_temp_1)
sum_a_arthro_sample_data_shannon_temp_1$group <- "Arthropods"
sum_a_arthro_sample_data_shannon_temp_1$diversity <- "Shannon"
sum_a_arthro_sample_data_shannon_temp_1$class <- "Very low temperature"
sum_a_arthro_sample_data_shannon_temp_1 <-sum_a_arthro_sample_data_shannon_temp_1[-1,]

sum_a_arthro_sample_data_richness_temp_1 <- capture.output(summary(sub_a_arthro_sample_data_temp_1$Richness),file=NULL,append=FALSE)
sum_a_arthro_sample_data_richness_temp_1 <- as.data.frame(sum_a_arthro_sample_data_richness_temp_1)
sum_a_arthro_sample_data_richness_temp_1$group <- "Arthropods"
sum_a_arthro_sample_data_richness_temp_1$diversity <- "Richness"
sum_a_arthro_sample_data_richness_temp_1$class <- "Very low temperature"
sum_a_arthro_sample_data_richness_temp_1 <-sum_a_arthro_sample_data_richness_temp_1[-1,]

sub_a_arthro_sample_data_temp_2 <- a_arthro_sample_data %>% filter(cluster_temp == "Low")
sum_a_arthro_sample_data_shannon_temp_2 <- capture.output(summary(sub_a_arthro_sample_data_temp_2$Shannon),file=NULL,append=FALSE)
sum_a_arthro_sample_data_shannon_temp_2 <- as.data.frame(sum_a_arthro_sample_data_shannon_temp_2)
sum_a_arthro_sample_data_shannon_temp_2$group <- "Arthropods"
sum_a_arthro_sample_data_shannon_temp_2$diversity <- "Shannon"
sum_a_arthro_sample_data_shannon_temp_2$class <- "Low temperature"
sum_a_arthro_sample_data_shannon_temp_2 <-sum_a_arthro_sample_data_shannon_temp_2[-1,]

sum_a_arthro_sample_data_richness_temp_2 <- capture.output(summary(sub_a_arthro_sample_data_temp_2$Richness),file=NULL,append=FALSE)
sum_a_arthro_sample_data_richness_temp_2 <- as.data.frame(sum_a_arthro_sample_data_richness_temp_2)
sum_a_arthro_sample_data_richness_temp_2$group <- "Arthropods"
sum_a_arthro_sample_data_richness_temp_2$diversity <- "Richness"
sum_a_arthro_sample_data_richness_temp_2$class <- "Low temperature"
sum_a_arthro_sample_data_richness_temp_2 <-sum_a_arthro_sample_data_richness_temp_2[-1,]

sub_a_arthro_sample_data_temp_3 <- a_arthro_sample_data %>% filter(cluster_temp == "Medium low")
sum_a_arthro_sample_data_shannon_temp_3 <- capture.output(summary(sub_a_arthro_sample_data_temp_3$Shannon),file=NULL,append=FALSE)
sum_a_arthro_sample_data_shannon_temp_3 <- as.data.frame(sum_a_arthro_sample_data_shannon_temp_3)
sum_a_arthro_sample_data_shannon_temp_3$group <- "Arthropods"
sum_a_arthro_sample_data_shannon_temp_3$diversity <- "Shannon"
sum_a_arthro_sample_data_shannon_temp_3$class <- "Medium low temperature"
sum_a_arthro_sample_data_shannon_temp_3 <-sum_a_arthro_sample_data_shannon_temp_3[-1,]

sum_a_arthro_sample_data_richness_temp_3 <- capture.output(summary(sub_a_arthro_sample_data_temp_3$Richness),file=NULL,append=FALSE)
sum_a_arthro_sample_data_richness_temp_3 <- as.data.frame(sum_a_arthro_sample_data_richness_temp_3)
sum_a_arthro_sample_data_richness_temp_3$group <- "Arthropods"
sum_a_arthro_sample_data_richness_temp_3$diversity <- "Richness"
sum_a_arthro_sample_data_richness_temp_3$class <- "Medium low temperature"
sum_a_arthro_sample_data_richness_temp_3 <-sum_a_arthro_sample_data_richness_temp_3[-1,]

sub_a_arthro_sample_data_temp_4 <- a_arthro_sample_data %>% filter(cluster_temp == "Medium high")
sum_a_arthro_sample_data_shannon_temp_4 <- capture.output(summary(sub_a_arthro_sample_data_temp_4$Shannon),file=NULL,append=FALSE)
sum_a_arthro_sample_data_shannon_temp_4 <- as.data.frame(sum_a_arthro_sample_data_shannon_temp_4)
sum_a_arthro_sample_data_shannon_temp_4$group <- "Arthropods"
sum_a_arthro_sample_data_shannon_temp_4$diversity <- "Shannon"
sum_a_arthro_sample_data_shannon_temp_4$class <- "Medium high temperature"
sum_a_arthro_sample_data_shannon_temp_4 <-sum_a_arthro_sample_data_shannon_temp_4[-1,]

sum_a_arthro_sample_data_richness_temp_4 <- capture.output(summary(sub_a_arthro_sample_data_temp_4$Richness),file=NULL,append=FALSE)
sum_a_arthro_sample_data_richness_temp_4 <- as.data.frame(sum_a_arthro_sample_data_richness_temp_4)
sum_a_arthro_sample_data_richness_temp_4$group <- "Arthropods"
sum_a_arthro_sample_data_richness_temp_4$diversity <- "Richness"
sum_a_arthro_sample_data_richness_temp_4$class <- "Medium high temperature"
sum_a_arthro_sample_data_richness_temp_4 <-sum_a_arthro_sample_data_richness_temp_4[-1,]

sub_a_arthro_sample_data_temp_5 <- a_arthro_sample_data %>% filter(cluster_temp == "High")
sum_a_arthro_sample_data_shannon_temp_5 <- capture.output(summary(sub_a_arthro_sample_data_temp_5$Shannon),file=NULL,append=FALSE)
sum_a_arthro_sample_data_shannon_temp_5 <- as.data.frame(sum_a_arthro_sample_data_shannon_temp_5)
sum_a_arthro_sample_data_shannon_temp_5$group <- "Arthropods"
sum_a_arthro_sample_data_shannon_temp_5$diversity <- "Shannon"
sum_a_arthro_sample_data_shannon_temp_5$class <- "High temperature"
sum_a_arthro_sample_data_shannon_temp_5 <-sum_a_arthro_sample_data_shannon_temp_5[-1,]

sum_a_arthro_sample_data_richness_temp_5 <- capture.output(summary(sub_a_arthro_sample_data_temp_5$Richness),file=NULL,append=FALSE)
sum_a_arthro_sample_data_richness_temp_5 <- as.data.frame(sum_a_arthro_sample_data_richness_temp_5)
sum_a_arthro_sample_data_richness_temp_5$group <- "Arthropods"
sum_a_arthro_sample_data_richness_temp_5$diversity <- "Richness"
sum_a_arthro_sample_data_richness_temp_5$class <- "High temperature"
sum_a_arthro_sample_data_richness_temp_5 <-sum_a_arthro_sample_data_richness_temp_5[-1,]

sub_a_arthro_sample_data_temp_6 <- a_arthro_sample_data %>% filter(cluster_temp == "Very high")
sum_a_arthro_sample_data_shannon_temp_6 <- capture.output(summary(sub_a_arthro_sample_data_temp_6$Shannon),file=NULL,append=FALSE)
sum_a_arthro_sample_data_shannon_temp_6 <- as.data.frame(sum_a_arthro_sample_data_shannon_temp_6)
sum_a_arthro_sample_data_shannon_temp_6$group <- "Arthropods"
sum_a_arthro_sample_data_shannon_temp_6$diversity <- "Shannon"
sum_a_arthro_sample_data_shannon_temp_6$class <- "Very high temperature"
sum_a_arthro_sample_data_shannon_temp_6 <-sum_a_arthro_sample_data_shannon_temp_6[-1,]

sum_a_arthro_sample_data_richness_temp_6 <- capture.output(summary(sub_a_arthro_sample_data_temp_6$Richness),file=NULL,append=FALSE)
sum_a_arthro_sample_data_richness_temp_6 <- as.data.frame(sum_a_arthro_sample_data_richness_temp_6)
sum_a_arthro_sample_data_richness_temp_6$group <- "Arthropods"
sum_a_arthro_sample_data_richness_temp_6$diversity <- "Richness"
sum_a_arthro_sample_data_richness_temp_6$class <- "Very high temperature"
sum_a_arthro_sample_data_richness_temp_6 <-sum_a_arthro_sample_data_richness_temp_6[-1,] 

sub_a_arthro_sample_data_temp_range_1 <- a_arthro_sample_data %>% filter(cluster_temp_range == "Very low")
sum_a_arthro_sample_data_shannon_temp_range_1 <- capture.output(summary(sub_a_arthro_sample_data_temp_range_1$Shannon),file=NULL,append=FALSE)
sum_a_arthro_sample_data_shannon_temp_range_1 <- as.data.frame(sum_a_arthro_sample_data_shannon_temp_range_1)
sum_a_arthro_sample_data_shannon_temp_range_1$group <- "Arthropods"
sum_a_arthro_sample_data_shannon_temp_range_1$diversity <- "Shannon"
sum_a_arthro_sample_data_shannon_temp_range_1$class <- "Very low temperature range"
sum_a_arthro_sample_data_shannon_temp_range_1 <-sum_a_arthro_sample_data_shannon_temp_range_1[-1,]

sum_a_arthro_sample_data_richness_temp_range_1 <- capture.output(summary(sub_a_arthro_sample_data_temp_range_1$Richness),file=NULL,append=FALSE)
sum_a_arthro_sample_data_richness_temp_range_1 <- as.data.frame(sum_a_arthro_sample_data_richness_temp_range_1)
sum_a_arthro_sample_data_richness_temp_range_1$group <- "Arthropods"
sum_a_arthro_sample_data_richness_temp_range_1$diversity <- "Richness"
sum_a_arthro_sample_data_richness_temp_range_1$class <- "Very low temperature range"
sum_a_arthro_sample_data_richness_temp_range_1 <-sum_a_arthro_sample_data_richness_temp_range_1[-1,]

sub_a_arthro_sample_data_temp_range_2 <- a_arthro_sample_data %>% filter(cluster_temp_range == "Low")
sum_a_arthro_sample_data_shannon_temp_range_2 <- capture.output(summary(sub_a_arthro_sample_data_temp_range_2$Shannon),file=NULL,append=FALSE)
sum_a_arthro_sample_data_shannon_temp_range_2 <- as.data.frame(sum_a_arthro_sample_data_shannon_temp_range_2)
sum_a_arthro_sample_data_shannon_temp_range_2$group <- "Arthropods"
sum_a_arthro_sample_data_shannon_temp_range_2$diversity <- "Shannon"
sum_a_arthro_sample_data_shannon_temp_range_2$class <- "Low temperature range"
sum_a_arthro_sample_data_shannon_temp_range_2 <-sum_a_arthro_sample_data_shannon_temp_range_2[-1,]

sum_a_arthro_sample_data_richness_temp_range_2 <- capture.output(summary(sub_a_arthro_sample_data_temp_range_2$Richness),file=NULL,append=FALSE)
sum_a_arthro_sample_data_richness_temp_range_2 <- as.data.frame(sum_a_arthro_sample_data_richness_temp_range_2)
sum_a_arthro_sample_data_richness_temp_range_2$group <- "Arthropods"
sum_a_arthro_sample_data_richness_temp_range_2$diversity <- "Richness"
sum_a_arthro_sample_data_richness_temp_range_2$class <- "Low temperature range"
sum_a_arthro_sample_data_richness_temp_range_2 <-sum_a_arthro_sample_data_richness_temp_range_2[-1,]

sub_a_arthro_sample_data_temp_range_3 <- a_arthro_sample_data %>% filter(cluster_temp_range == "Medium low")
sum_a_arthro_sample_data_shannon_temp_range_3 <- capture.output(summary(sub_a_arthro_sample_data_temp_range_3$Shannon),file=NULL,append=FALSE)
sum_a_arthro_sample_data_shannon_temp_range_3 <- as.data.frame(sum_a_arthro_sample_data_shannon_temp_range_3)
sum_a_arthro_sample_data_shannon_temp_range_3$group <- "Arthropods"
sum_a_arthro_sample_data_shannon_temp_range_3$diversity <- "Shannon"
sum_a_arthro_sample_data_shannon_temp_range_3$class <- "Medium low temperature range"
sum_a_arthro_sample_data_shannon_temp_range_3 <-sum_a_arthro_sample_data_shannon_temp_range_3[-1,]

sum_a_arthro_sample_data_richness_temp_range_3 <- capture.output(summary(sub_a_arthro_sample_data_temp_range_3$Richness),file=NULL,append=FALSE)
sum_a_arthro_sample_data_richness_temp_range_3 <- as.data.frame(sum_a_arthro_sample_data_richness_temp_range_3)
sum_a_arthro_sample_data_richness_temp_range_3$group <- "Arthropods"
sum_a_arthro_sample_data_richness_temp_range_3$diversity <- "Richness"
sum_a_arthro_sample_data_richness_temp_range_3$class <- "Medium low temperature range"
sum_a_arthro_sample_data_richness_temp_range_3 <-sum_a_arthro_sample_data_richness_temp_range_3[-1,]

sub_a_arthro_sample_data_temp_range_4 <- a_arthro_sample_data %>% filter(cluster_temp_range == "Medium high")
sum_a_arthro_sample_data_shannon_temp_range_4 <- capture.output(summary(sub_a_arthro_sample_data_temp_range_4$Shannon),file=NULL,append=FALSE)
sum_a_arthro_sample_data_shannon_temp_range_4 <- as.data.frame(sum_a_arthro_sample_data_shannon_temp_range_4)
sum_a_arthro_sample_data_shannon_temp_range_4$group <- "Arthropods"
sum_a_arthro_sample_data_shannon_temp_range_4$diversity <- "Shannon"
sum_a_arthro_sample_data_shannon_temp_range_4$class <- "Medium high temperature range"
sum_a_arthro_sample_data_shannon_temp_range_4 <-sum_a_arthro_sample_data_shannon_temp_range_4[-1,]

sum_a_arthro_sample_data_richness_temp_range_4 <- capture.output(summary(sub_a_arthro_sample_data_temp_range_4$Richness),file=NULL,append=FALSE)
sum_a_arthro_sample_data_richness_temp_range_4 <- as.data.frame(sum_a_arthro_sample_data_richness_temp_range_4)
sum_a_arthro_sample_data_richness_temp_range_4$group <- "Arthropods"
sum_a_arthro_sample_data_richness_temp_range_4$diversity <- "Richness"
sum_a_arthro_sample_data_richness_temp_range_4$class <- "Medium high temperature range"
sum_a_arthro_sample_data_richness_temp_range_4 <-sum_a_arthro_sample_data_richness_temp_range_4[-1,]

sub_a_arthro_sample_data_temp_range_5 <- a_arthro_sample_data %>% filter(cluster_temp_range == "High")
sum_a_arthro_sample_data_shannon_temp_range_5 <- capture.output(summary(sub_a_arthro_sample_data_temp_range_5$Shannon),file=NULL,append=FALSE)
sum_a_arthro_sample_data_shannon_temp_range_5 <- as.data.frame(sum_a_arthro_sample_data_shannon_temp_range_5)
sum_a_arthro_sample_data_shannon_temp_range_5$group <- "Arthropods"
sum_a_arthro_sample_data_shannon_temp_range_5$diversity <- "Shannon"
sum_a_arthro_sample_data_shannon_temp_range_5$class <- "High temperature range"
sum_a_arthro_sample_data_shannon_temp_range_5 <-sum_a_arthro_sample_data_shannon_temp_range_5[-1,]

sum_a_arthro_sample_data_richness_temp_range_5 <- capture.output(summary(sub_a_arthro_sample_data_temp_range_5$Richness),file=NULL,append=FALSE)
sum_a_arthro_sample_data_richness_temp_range_5 <- as.data.frame(sum_a_arthro_sample_data_richness_temp_range_5)
sum_a_arthro_sample_data_richness_temp_range_5$group <- "Arthropods"
sum_a_arthro_sample_data_richness_temp_range_5$diversity <- "Richness"
sum_a_arthro_sample_data_richness_temp_range_5$class <- "High temperature range"
sum_a_arthro_sample_data_richness_temp_range_5 <-sum_a_arthro_sample_data_richness_temp_range_5[-1,]

sub_a_arthro_sample_data_temp_range_6 <- a_arthro_sample_data %>% filter(cluster_temp_range == "Very high")
sum_a_arthro_sample_data_shannon_temp_range_6 <- capture.output(summary(sub_a_arthro_sample_data_temp_range_6$Shannon),file=NULL,append=FALSE)
sum_a_arthro_sample_data_shannon_temp_range_6 <- as.data.frame(sum_a_arthro_sample_data_shannon_temp_range_6)
sum_a_arthro_sample_data_shannon_temp_range_6$group <- "Arthropods"
sum_a_arthro_sample_data_shannon_temp_range_6$diversity <- "Shannon"
sum_a_arthro_sample_data_shannon_temp_range_6$class <- "Very high temperature range"
sum_a_arthro_sample_data_shannon_temp_range_6 <-sum_a_arthro_sample_data_shannon_temp_range_6[-1,]

sum_a_arthro_sample_data_richness_temp_range_6 <- capture.output(summary(sub_a_arthro_sample_data_temp_range_6$Richness),file=NULL,append=FALSE)
sum_a_arthro_sample_data_richness_temp_range_6 <- as.data.frame(sum_a_arthro_sample_data_richness_temp_range_6)
sum_a_arthro_sample_data_richness_temp_range_6$group <- "Arthropods"
sum_a_arthro_sample_data_richness_temp_range_6$diversity <- "Richness"
sum_a_arthro_sample_data_richness_temp_range_6$class <- "Very high temperature range"
sum_a_arthro_sample_data_richness_temp_range_6 <-sum_a_arthro_sample_data_richness_temp_range_6[-1,]

sub_a_arthro_sample_data_prec_1 <- a_arthro_sample_data %>% filter(cluster_prec == "Very low")
sum_a_arthro_sample_data_shannon_prec_1 <- capture.output(summary(sub_a_arthro_sample_data_prec_1$Shannon),file=NULL,append=FALSE)
sum_a_arthro_sample_data_shannon_prec_1 <- as.data.frame(sum_a_arthro_sample_data_shannon_prec_1)
sum_a_arthro_sample_data_shannon_prec_1$group <- "Arthropods"
sum_a_arthro_sample_data_shannon_prec_1$diversity <- "Shannon"
sum_a_arthro_sample_data_shannon_prec_1$class <- "Very low precipitation "
sum_a_arthro_sample_data_shannon_prec_1 <-sum_a_arthro_sample_data_shannon_prec_1[-1,]

sum_a_arthro_sample_data_richness_prec_1 <- capture.output(summary(sub_a_arthro_sample_data_prec_1$Richness),file=NULL,append=FALSE)
sum_a_arthro_sample_data_richness_prec_1 <- as.data.frame(sum_a_arthro_sample_data_richness_prec_1)
sum_a_arthro_sample_data_richness_prec_1$group <- "Arthropods"
sum_a_arthro_sample_data_richness_prec_1$diversity <- "Richness"
sum_a_arthro_sample_data_richness_prec_1$class <- "Very low precipitation "
sum_a_arthro_sample_data_richness_prec_1 <-sum_a_arthro_sample_data_richness_prec_1[-1,]

sub_a_arthro_sample_data_prec_2 <- a_arthro_sample_data %>% filter(cluster_prec == "Low")
sum_a_arthro_sample_data_shannon_prec_2 <- capture.output(summary(sub_a_arthro_sample_data_prec_2$Shannon),file=NULL,append=FALSE)
sum_a_arthro_sample_data_shannon_prec_2 <- as.data.frame(sum_a_arthro_sample_data_shannon_prec_2)
sum_a_arthro_sample_data_shannon_prec_2$group <- "Arthropods"
sum_a_arthro_sample_data_shannon_prec_2$diversity <- "Shannon"
sum_a_arthro_sample_data_shannon_prec_2$class <- "Low precipitation"
sum_a_arthro_sample_data_shannon_prec_2 <-sum_a_arthro_sample_data_shannon_prec_2[-1,]

sum_a_arthro_sample_data_richness_prec_2 <- capture.output(summary(sub_a_arthro_sample_data_prec_2$Richness),file=NULL,append=FALSE)
sum_a_arthro_sample_data_richness_prec_2 <- as.data.frame(sum_a_arthro_sample_data_richness_prec_2)
sum_a_arthro_sample_data_richness_prec_2$group <- "Arthropods"
sum_a_arthro_sample_data_richness_prec_2$diversity <- "Richness"
sum_a_arthro_sample_data_richness_prec_2$class <- "Low precipitation"
sum_a_arthro_sample_data_richness_prec_2 <-sum_a_arthro_sample_data_richness_prec_2[-1,]

sub_a_arthro_sample_data_prec_3 <- a_arthro_sample_data %>% filter(cluster_prec == "Medium low")
sum_a_arthro_sample_data_shannon_prec_3 <- capture.output(summary(sub_a_arthro_sample_data_prec_3$Shannon),file=NULL,append=FALSE)
sum_a_arthro_sample_data_shannon_prec_3 <- as.data.frame(sum_a_arthro_sample_data_shannon_prec_3)
sum_a_arthro_sample_data_shannon_prec_3$group <- "Arthropods"
sum_a_arthro_sample_data_shannon_prec_3$diversity <- "Shannon"
sum_a_arthro_sample_data_shannon_prec_3$class <- "Medium low precipitation"
sum_a_arthro_sample_data_shannon_prec_3 <-sum_a_arthro_sample_data_shannon_prec_3[-1,]

sum_a_arthro_sample_data_richness_prec_3 <- capture.output(summary(sub_a_arthro_sample_data_prec_3$Richness),file=NULL,append=FALSE)
sum_a_arthro_sample_data_richness_prec_3 <- as.data.frame(sum_a_arthro_sample_data_richness_prec_3)
sum_a_arthro_sample_data_richness_prec_3$group <- "Arthropods"
sum_a_arthro_sample_data_richness_prec_3$diversity <- "Richness"
sum_a_arthro_sample_data_richness_prec_3$class <- "Medium low precipitation"
sum_a_arthro_sample_data_richness_prec_3 <-sum_a_arthro_sample_data_richness_prec_3[-1,]

sub_a_arthro_sample_data_prec_4 <- a_arthro_sample_data %>% filter(cluster_prec == "Medium high")
sum_a_arthro_sample_data_shannon_prec_4 <- capture.output(summary(sub_a_arthro_sample_data_prec_4$Shannon),file=NULL,append=FALSE)
sum_a_arthro_sample_data_shannon_prec_4 <- as.data.frame(sum_a_arthro_sample_data_shannon_prec_4)
sum_a_arthro_sample_data_shannon_prec_4$group <- "Arthropods"
sum_a_arthro_sample_data_shannon_prec_4$diversity <- "Shannon"
sum_a_arthro_sample_data_shannon_prec_4$class <- "Medium high precipitation"
sum_a_arthro_sample_data_shannon_prec_4 <-sum_a_arthro_sample_data_shannon_prec_4[-1,]

sum_a_arthro_sample_data_richness_prec_4 <- capture.output(summary(sub_a_arthro_sample_data_prec_4$Richness),file=NULL,append=FALSE)
sum_a_arthro_sample_data_richness_prec_4 <- as.data.frame(sum_a_arthro_sample_data_richness_prec_4)
sum_a_arthro_sample_data_richness_prec_4$group <- "Arthropods"
sum_a_arthro_sample_data_richness_prec_4$diversity <- "Richness"
sum_a_arthro_sample_data_richness_prec_4$class <- "Medium high precipitation"
sum_a_arthro_sample_data_richness_prec_4 <-sum_a_arthro_sample_data_richness_prec_4[-1,]

sub_a_arthro_sample_data_prec_5 <- a_arthro_sample_data %>% filter(cluster_prec == "High")
sum_a_arthro_sample_data_shannon_prec_5 <- capture.output(summary(sub_a_arthro_sample_data_prec_5$Shannon),file=NULL,append=FALSE)
sum_a_arthro_sample_data_shannon_prec_5 <- as.data.frame(sum_a_arthro_sample_data_shannon_prec_5)
sum_a_arthro_sample_data_shannon_prec_5$group <- "Arthropods"
sum_a_arthro_sample_data_shannon_prec_5$diversity <- "Shannon"
sum_a_arthro_sample_data_shannon_prec_5$class <- "High precipitation"
sum_a_arthro_sample_data_shannon_prec_5 <-sum_a_arthro_sample_data_shannon_prec_5[-1,]

sum_a_arthro_sample_data_richness_prec_5 <- capture.output(summary(sub_a_arthro_sample_data_prec_5$Richness),file=NULL,append=FALSE)
sum_a_arthro_sample_data_richness_prec_5 <- as.data.frame(sum_a_arthro_sample_data_richness_prec_5)
sum_a_arthro_sample_data_richness_prec_5$group <- "Arthropods"
sum_a_arthro_sample_data_richness_prec_5$diversity <- "Richness"
sum_a_arthro_sample_data_richness_prec_5$class <- "High precipitation"
sum_a_arthro_sample_data_richness_prec_5 <-sum_a_arthro_sample_data_richness_prec_5[-1,]

sub_a_arthro_sample_data_prec_6 <- a_arthro_sample_data %>% filter(cluster_prec == "Very high")
sum_a_arthro_sample_data_shannon_prec_6 <- capture.output(summary(sub_a_arthro_sample_data_prec_6$Shannon),file=NULL,append=FALSE)
sum_a_arthro_sample_data_shannon_prec_6 <- as.data.frame(sum_a_arthro_sample_data_shannon_prec_6)
sum_a_arthro_sample_data_shannon_prec_6$group <- "Arthropods"
sum_a_arthro_sample_data_shannon_prec_6$diversity <- "Shannon"
sum_a_arthro_sample_data_shannon_prec_6$class <- "Very high precipitation"
sum_a_arthro_sample_data_shannon_prec_6 <-sum_a_arthro_sample_data_shannon_prec_6[-1,]

sum_a_arthro_sample_data_richness_prec_6 <- capture.output(summary(sub_a_arthro_sample_data_prec_6$Richness),file=NULL,append=FALSE)
sum_a_arthro_sample_data_richness_prec_6 <- as.data.frame(sum_a_arthro_sample_data_richness_prec_6)
sum_a_arthro_sample_data_richness_prec_6$group <- "Arthropods"
sum_a_arthro_sample_data_richness_prec_6$diversity <- "Richness"
sum_a_arthro_sample_data_richness_prec_6$class <- "Very high precipitation"
sum_a_arthro_sample_data_richness_prec_6 <-sum_a_arthro_sample_data_richness_prec_6[-1,]




#Annelids

sub_a_anneli_sample_data_season1 <- a_anneli_sample_data %>% filter(Sample_season == "Spring")
sum_a_anneli_richness_season1 <- capture.output(summary(sub_a_anneli_sample_data_season1$Richness),file=NULL,append=FALSE)
sum_a_anneli_richness_season1 <- as.data.frame(sum_a_anneli_richness_season1)
sum_a_anneli_richness_season1$group <- "Annelids"
sum_a_anneli_richness_season1$diversity <- "Richness"
sum_a_anneli_richness_season1$class <- "Season - Spring"
sum_a_anneli_richness_season1 <-sum_a_anneli_richness_season1[-1,]

sum_a_anneli_shannon_season1 <- capture.output(summary(sub_a_anneli_sample_data_season1$Shannon),file=NULL,append=FALSE)
sum_a_anneli_shannon_season1 <- as.data.frame(sum_a_anneli_shannon_season1)
sum_a_anneli_shannon_season1$group <- "Annelids"
sum_a_anneli_shannon_season1$diversity <- "Shannon"
sum_a_anneli_shannon_season1$class <- "Season - Spring"
sum_a_anneli_shannon_season1 <-sum_a_anneli_shannon_season1[-1,]

sub_a_anneli_sample_data_season2 <- a_anneli_sample_data %>% filter(Sample_season == "Summer")
sum_a_anneli_richness_season2 <- capture.output(summary(sub_a_anneli_sample_data_season2$Richness),file=NULL,append=FALSE)
sum_a_anneli_richness_season2 <- as.data.frame(sum_a_anneli_richness_season2)
sum_a_anneli_richness_season2$group <- "Annelids"
sum_a_anneli_richness_season2$diversity <- "Richness"
sum_a_anneli_richness_season2$class <- "Season - Summer"
sum_a_anneli_richness_season2 <-sum_a_anneli_richness_season2[-1,]

sum_a_anneli_shannon_season2 <- capture.output(summary(sub_a_anneli_sample_data_season2$Shannon),file=NULL,append=FALSE)
sum_a_anneli_shannon_season2 <- as.data.frame(sum_a_anneli_shannon_season2)
sum_a_anneli_shannon_season2$group <- "Annelids"
sum_a_anneli_shannon_season2$diversity <- "Shannon"
sum_a_anneli_shannon_season2$class <- "Season - Summer"
sum_a_anneli_shannon_season2 <-sum_a_anneli_shannon_season2[-1,]

sub_a_anneli_sample_data_season3 <- a_anneli_sample_data %>% filter(Sample_season == "Autumn")
sum_a_anneli_richness_season3 <- capture.output(summary(sub_a_anneli_sample_data_season3$Richness),file=NULL,append=FALSE)
sum_a_anneli_richness_season3 <- as.data.frame(sum_a_anneli_richness_season3)
sum_a_anneli_richness_season3$group <- "Annelids"
sum_a_anneli_richness_season3$diversity <- "Richness"
sum_a_anneli_richness_season3$class <- "Season - Autumn"
sum_a_anneli_richness_season3 <-sum_a_anneli_richness_season3[-1,]

sum_a_anneli_shannon_season3 <- capture.output(summary(sub_a_anneli_sample_data_season3$Shannon),file=NULL,append=FALSE)
sum_a_anneli_shannon_season3 <- as.data.frame(sum_a_anneli_shannon_season3)
sum_a_anneli_shannon_season3$group <- "Annelids"
sum_a_anneli_shannon_season3$diversity <- "Shannon"
sum_a_anneli_shannon_season3$class <- "Season - Autumn"
sum_a_anneli_shannon_season3 <-sum_a_anneli_shannon_season3[-1,]

sub_a_anneli_sample_data_depth1 <- a_anneli_sample_data %>% filter(depth_grouped == "Shallow")
sum_a_anneli_richness_depth1 <- capture.output(summary(sub_a_anneli_sample_data_depth1$Richness),file=NULL,append=FALSE)
sum_a_anneli_richness_depth1 <- as.data.frame(sum_a_anneli_richness_depth1)
sum_a_anneli_richness_depth1$group <- "Annelids"
sum_a_anneli_richness_depth1$diversity <- "Richness"
sum_a_anneli_richness_depth1$class <- "Depth <1m"
sum_a_anneli_richness_depth1 <-sum_a_anneli_richness_depth1[-1,]

sum_a_anneli_shannon_depth1 <- capture.output(summary(sub_a_anneli_sample_data_depth1$Shannon),file=NULL,append=FALSE)
sum_a_anneli_shannon_depth1 <- as.data.frame(sum_a_anneli_shannon_depth1)
sum_a_anneli_shannon_depth1$group <- "Annelids"
sum_a_anneli_shannon_depth1$diversity <- "Shannon"
sum_a_anneli_shannon_depth1$class <- "Depth >1m"
sum_a_anneli_shannon_depth1 <-sum_a_anneli_shannon_depth1[-1,]

sub_a_anneli_sample_data_depth2 <- a_anneli_sample_data %>% filter(depth_grouped == "Moderately deep")
sum_a_anneli_richness_depth2 <- capture.output(summary(sub_a_anneli_sample_data_depth2$Richness),file=NULL,append=FALSE)
sum_a_anneli_richness_depth2 <- as.data.frame(sum_a_anneli_richness_depth2)
sum_a_anneli_richness_depth2$group <- "Annelids"
sum_a_anneli_richness_depth2$diversity <- "Richness"
sum_a_anneli_richness_depth2$class <- "Depth 1-2m"
sum_a_anneli_richness_depth2 <-sum_a_anneli_richness_depth2[-1,]

sum_a_anneli_shannon_depth2 <- capture.output(summary(sub_a_anneli_sample_data_depth2$Shannon),file=NULL,append=FALSE)
sum_a_anneli_shannon_depth2 <- as.data.frame(sum_a_anneli_shannon_depth2)
sum_a_anneli_shannon_depth2$group <- "Annelids"
sum_a_anneli_shannon_depth2$diversity <- "Shannon"
sum_a_anneli_shannon_depth2$class <- "Depth 1-2m"
sum_a_anneli_shannon_depth2 <-sum_a_anneli_shannon_depth2[-1,]

sub_a_anneli_sample_data_depth3 <- a_anneli_sample_data %>% filter(depth_grouped == "Deep")
sum_a_anneli_richness_depth3 <- capture.output(summary(sub_a_anneli_sample_data_depth3$Richness),file=NULL,append=FALSE)
sum_a_anneli_richness_depth3 <- as.data.frame(sum_a_anneli_richness_depth3)
sum_a_anneli_richness_depth3$group <- "Annelids"
sum_a_anneli_richness_depth3$diversity <- "Richness"
sum_a_anneli_richness_depth3$class <- "Depth >2m"
sum_a_anneli_richness_depth3 <-sum_a_anneli_richness_depth3[-1,]

sum_a_anneli_shannon_depth3 <- capture.output(summary(sub_a_anneli_sample_data_depth3$Shannon),file=NULL,append=FALSE)
sum_a_anneli_shannon_depth3 <- as.data.frame(sum_a_anneli_shannon_depth3)
sum_a_anneli_shannon_depth3$group <- "Annelids"
sum_a_anneli_shannon_depth3$diversity <- "Shannon"
sum_a_anneli_shannon_depth3$class <- "Depth >2m"
sum_a_anneli_shannon_depth3 <-sum_a_anneli_shannon_depth3[-1,]

sub_a_anneli_sample_data_pH1 <- a_anneli_sample_data %>% filter(pH_grouped == "Acidic")
sum_a_anneli_richness_pH1 <- capture.output(summary(sub_a_anneli_sample_data_pH1$Richness),file=NULL,append=FALSE)
sum_a_anneli_richness_pH1 <- as.data.frame(sum_a_anneli_richness_pH1)
sum_a_anneli_richness_pH1$group <- "Annelids"
sum_a_anneli_richness_pH1$diversity <- "Richness"
sum_a_anneli_richness_pH1$class <- "pH 3.3-6.5"
sum_a_anneli_richness_pH1 <-sum_a_anneli_richness_pH1[-1,]

sum_a_anneli_shannon_pH1 <- capture.output(summary(sub_a_anneli_sample_data_pH1$Shannon),file=NULL,append=FALSE)
sum_a_anneli_shannon_pH1 <- as.data.frame(sum_a_anneli_shannon_pH1)
sum_a_anneli_shannon_pH1$group <- "Annelids"
sum_a_anneli_shannon_pH1$diversity <- "Shannon"
sum_a_anneli_shannon_pH1$class <- "pH 3.3-6.5"
sum_a_anneli_shannon_pH1 <-sum_a_anneli_shannon_pH1[-1,]

sub_a_anneli_sample_data_pH2 <- a_anneli_sample_data %>% filter(pH_grouped == "Neutral")
sum_a_anneli_richness_pH2 <- capture.output(summary(sub_a_anneli_sample_data_pH2$Richness),file=NULL,append=FALSE)
sum_a_anneli_richness_pH2 <- as.data.frame(sum_a_anneli_richness_pH2)
sum_a_anneli_richness_pH2$group <- "Annelids"
sum_a_anneli_richness_pH2$diversity <- "Richness"
sum_a_anneli_richness_pH2$class <- "pH 6.5-7.5"
sum_a_anneli_richness_pH2 <-sum_a_anneli_richness_pH2[-1,]

sum_a_anneli_shannon_pH2 <- capture.output(summary(sub_a_anneli_sample_data_pH2$Shannon),file=NULL,append=FALSE)
sum_a_anneli_shannon_pH2 <- as.data.frame(sum_a_anneli_shannon_pH2)
sum_a_anneli_shannon_pH2$group <- "Annelids"
sum_a_anneli_shannon_pH2$diversity <- "Shannon"
sum_a_anneli_shannon_pH2$class <- "pH 6.5-7.5"
sum_a_anneli_shannon_pH2 <-sum_a_anneli_shannon_pH2[-1,]

sub_a_anneli_sample_data_pH3 <- a_anneli_sample_data %>% filter(pH_grouped == "Alkaline")
sum_a_anneli_richness_pH3 <- capture.output(summary(sub_a_anneli_sample_data_pH3$Richness),file=NULL,append=FALSE)
sum_a_anneli_richness_pH3 <- as.data.frame(sum_a_anneli_richness_pH3)
sum_a_anneli_richness_pH3$group <- "Annelids"
sum_a_anneli_richness_pH3$diversity <- "Richness"
sum_a_anneli_richness_pH3$class <- "pH 7.5+"
sum_a_anneli_richness_pH3 <-sum_a_anneli_richness_pH3[-1,]

sum_a_anneli_shannon_pH3 <- capture.output(summary(sub_a_anneli_sample_data_pH3$Shannon),file=NULL,append=FALSE)
sum_a_anneli_shannon_pH3 <- as.data.frame(sum_a_anneli_shannon_pH3)
sum_a_anneli_shannon_pH3$group <- "Annelids"
sum_a_anneli_shannon_pH3$diversity <- "Shannon"
sum_a_anneli_shannon_pH3$class <- "pH 7.5+"
sum_a_anneli_shannon_pH3 <-sum_a_anneli_shannon_pH3[-1,]

sub_a_anneli_sample_data_erosion1 <- a_anneli_sample_data %>% filter(Erosion_grouped == "Low")
sum_a_anneli_richness_erosion1 <- capture.output(summary(sub_a_anneli_sample_data_erosion1$Richness),file=NULL,append=FALSE)
sum_a_anneli_richness_erosion1 <- as.data.frame(sum_a_anneli_richness_erosion1)
sum_a_anneli_richness_erosion1$group <- "Annelids"
sum_a_anneli_richness_erosion1$diversity <- "Richness"
sum_a_anneli_richness_erosion1$class <- "Erosion 0-1"
sum_a_anneli_richness_erosion1 <-sum_a_anneli_richness_erosion1[-1,]

sum_a_anneli_shannon_erosion1 <- capture.output(summary(sub_a_anneli_sample_data_erosion1$Shannon),file=NULL,append=FALSE)
sum_a_anneli_shannon_erosion1 <- as.data.frame(sum_a_anneli_shannon_erosion1)
sum_a_anneli_shannon_erosion1$group <- "Annelids"
sum_a_anneli_shannon_erosion1$diversity <- "Shannon"
sum_a_anneli_shannon_erosion1$class <- "Erosion 0-1"
sum_a_anneli_shannon_erosion1 <-sum_a_anneli_shannon_erosion1[-1,]

sub_a_anneli_sample_data_erosion2 <- a_anneli_sample_data %>% filter(Erosion_grouped == "Medium high")
sum_a_anneli_richness_erosion2 <- capture.output(summary(sub_a_anneli_sample_data_erosion2$Richness),file=NULL,append=FALSE)
sum_a_anneli_richness_erosion2 <- as.data.frame(sum_a_anneli_richness_erosion2)
sum_a_anneli_richness_erosion2$group <- "Annelids"
sum_a_anneli_richness_erosion2$diversity <- "Richness"
sum_a_anneli_richness_erosion2$class <- "Erosion 1-5"
sum_a_anneli_richness_erosion2 <-sum_a_anneli_richness_erosion2[-1,]

sum_a_anneli_shannon_erosion2 <- capture.output(summary(sub_a_anneli_sample_data_erosion2$Shannon),file=NULL,append=FALSE)
sum_a_anneli_shannon_erosion2 <- as.data.frame(sum_a_anneli_shannon_erosion2)
sum_a_anneli_shannon_erosion2$group <- "Annelids"
sum_a_anneli_shannon_erosion2$diversity <- "Shannon"
sum_a_anneli_shannon_erosion2$class <- "Erosion 1-5"
sum_a_anneli_shannon_erosion2 <-sum_a_anneli_shannon_erosion2[-1,]

sub_a_anneli_sample_data_erosion3 <- a_anneli_sample_data %>% filter(Erosion_grouped == "High")
sum_a_anneli_richness_erosion3 <- capture.output(summary(sub_a_anneli_sample_data_erosion3$Richness),file=NULL,append=FALSE)
sum_a_anneli_richness_erosion3 <- as.data.frame(sum_a_anneli_richness_erosion3)
sum_a_anneli_richness_erosion3$group <- "Annelids"
sum_a_anneli_richness_erosion3$diversity <- "Richness"
sum_a_anneli_richness_erosion3$class <- "Erosion 5-10"
sum_a_anneli_richness_erosion3 <-sum_a_anneli_richness_erosion3[-1,]

sum_a_anneli_shannon_erosion3 <- capture.output(summary(sub_a_anneli_sample_data_erosion3$Shannon),file=NULL,append=FALSE)
sum_a_anneli_shannon_erosion3 <- as.data.frame(sum_a_anneli_shannon_erosion3)
sum_a_anneli_shannon_erosion3$group <- "Annelids"
sum_a_anneli_shannon_erosion3$diversity <- "Shannon"
sum_a_anneli_shannon_erosion3$class <- "Erosion 5-10"
sum_a_anneli_shannon_erosion3 <-sum_a_anneli_shannon_erosion3[-1,]

sub_a_anneli_sample_data_erosion4 <- a_anneli_sample_data %>% filter(Erosion_grouped == "Very high")
sum_a_anneli_richness_erosion4 <- capture.output(summary(sub_a_anneli_sample_data_erosion4$Richness),file=NULL,append=FALSE)
sum_a_anneli_richness_erosion4 <- as.data.frame(sum_a_anneli_richness_erosion4)
sum_a_anneli_richness_erosion4$group <- "Annelids"
sum_a_anneli_richness_erosion4$diversity <- "Richness"
sum_a_anneli_richness_erosion4$class <- "Erosion 10+"
sum_a_anneli_richness_erosion4 <-sum_a_anneli_richness_erosion4[-1,]

sum_a_anneli_shannon_erosion4 <- capture.output(summary(sub_a_anneli_sample_data_erosion4$Shannon),file=NULL,append=FALSE)
sum_a_anneli_shannon_erosion4 <- as.data.frame(sum_a_anneli_shannon_erosion4)
sum_a_anneli_shannon_erosion4$group <- "Annelids"
sum_a_anneli_shannon_erosion4$diversity <- "Shannon"
sum_a_anneli_shannon_erosion4$class <- "Erosion 10+"
sum_a_anneli_shannon_erosion4 <-sum_a_anneli_shannon_erosion4[-1,]


sub_anneli_sample_data_CL1 <- a_anneli_sample_data %>% filter(LU2 == "CL1")
sum_anneli_sample_data_richness_CL1 <- capture.output(summary(sub_anneli_sample_data_CL1$Richness),file=NULL,append=FALSE)
sum_anneli_sample_data_richness_CL1 <- as.data.frame(sum_anneli_sample_data_richness_CL1)
sum_anneli_sample_data_richness_CL1$group <- "Annelids"
sum_anneli_sample_data_richness_CL1$diversity <- "Richness"
sum_anneli_sample_data_richness_CL1$class <- "Land Cover CL1"
sum_anneli_sample_data_richness_CL1 <- sum_anneli_sample_data_richness_CL1[-1,]

sum_anneli_sample_data_shannon_CL1 <- capture.output(summary(sub_anneli_sample_data_CL1$Shannon),file=NULL,append=FALSE)
sum_anneli_sample_data_shannon_CL1 <- as.data.frame(sum_anneli_sample_data_shannon_CL1)
sum_anneli_sample_data_shannon_CL1$group <- "Annelids"
sum_anneli_sample_data_shannon_CL1$diversity <- "Shannon"
sum_anneli_sample_data_shannon_CL1$class <- "Land Cover CL1"
sum_anneli_sample_data_shannon_CL1 <-sum_anneli_sample_data_shannon_CL1[-1,]

sub_anneli_sample_data_CL2 <- a_anneli_sample_data %>% filter(LU2 == "CL2")
sum_anneli_sample_data_richness_CL2 <- capture.output(summary(sub_anneli_sample_data_CL2$Richness),file=NULL,append=FALSE)
sum_anneli_sample_data_richness_CL2 <- as.data.frame(sum_anneli_sample_data_richness_CL2)
sum_anneli_sample_data_richness_CL2$group <- "Annelids"
sum_anneli_sample_data_richness_CL2$diversity <- "Richness"
sum_anneli_sample_data_richness_CL2$class <- "Land Cover CL2"
sum_anneli_sample_data_richness_CL2 <- sum_anneli_sample_data_richness_CL2[-1,]

sum_anneli_sample_data_shannon_CL2 <- capture.output(summary(sub_anneli_sample_data_CL2$Shannon),file=NULL,append=FALSE)
sum_anneli_sample_data_shannon_CL2 <- as.data.frame(sum_anneli_sample_data_shannon_CL2)
sum_anneli_sample_data_shannon_CL2$group <- "Annelids"
sum_anneli_sample_data_shannon_CL2$diversity <- "Shannon"
sum_anneli_sample_data_shannon_CL2$class <- "Land Cover CL2"
sum_anneli_sample_data_shannon_CL2 <-sum_anneli_sample_data_shannon_CL2[-1,]

sub_anneli_sample_data_CL3 <- a_anneli_sample_data %>% filter(LU2 == "CL3")
sum_anneli_sample_data_richness_CL3 <- capture.output(summary(sub_anneli_sample_data_CL3$Richness),file=NULL,append=FALSE)
sum_anneli_sample_data_richness_CL3 <- as.data.frame(sum_anneli_sample_data_richness_CL3)
sum_anneli_sample_data_richness_CL3$group <- "Annelids"
sum_anneli_sample_data_richness_CL3$diversity <- "Richness"
sum_anneli_sample_data_richness_CL3$class <- "Land Cover CL3"
sum_anneli_sample_data_richness_CL3 <- sum_anneli_sample_data_richness_CL3[-1,]

sum_anneli_sample_data_shannon_CL3 <- capture.output(summary(sub_anneli_sample_data_CL3$Shannon),file=NULL,append=FALSE)
sum_anneli_sample_data_shannon_CL3 <- as.data.frame(sum_anneli_sample_data_shannon_CL3)
sum_anneli_sample_data_shannon_CL3$group <- "Annelids"
sum_anneli_sample_data_shannon_CL3$diversity <- "Shannon"
sum_anneli_sample_data_shannon_CL3$class <- "Land Cover CL3"
sum_anneli_sample_data_shannon_CL3 <-sum_anneli_sample_data_shannon_CL3[-1,]

sub_anneli_sample_data_CL4 <- a_anneli_sample_data %>% filter(LU2 == "CL4")
sum_anneli_sample_data_richness_CL4 <- capture.output(summary(sub_anneli_sample_data_CL4$Richness),file=NULL,append=FALSE)
sum_anneli_sample_data_richness_CL4 <- as.data.frame(sum_anneli_sample_data_richness_CL4)
sum_anneli_sample_data_richness_CL4$group <- "Annelids"
sum_anneli_sample_data_richness_CL4$diversity <- "Richness"
sum_anneli_sample_data_richness_CL4$class <- "Land Cover CL4"
sum_anneli_sample_data_richness_CL4 <- sum_anneli_sample_data_richness_CL4[-1,]

sum_anneli_sample_data_shannon_CL4 <- capture.output(summary(sub_anneli_sample_data_CL4$Shannon),file=NULL,append=FALSE)
sum_anneli_sample_data_shannon_CL4 <- as.data.frame(sum_anneli_sample_data_shannon_CL4)
sum_anneli_sample_data_shannon_CL4$group <- "Annelids"
sum_anneli_sample_data_shannon_CL4$diversity <- "Shannon"
sum_anneli_sample_data_shannon_CL4$class <- "Land Cover CL4"
sum_anneli_sample_data_shannon_CL4 <-sum_anneli_sample_data_shannon_CL4[-1,]

sub_anneli_sample_data_GL1 <- a_anneli_sample_data %>% filter(LU2 == "GL1")
sum_anneli_sample_data_richness_GL1 <- capture.output(summary(sub_anneli_sample_data_GL1$Richness),file=NULL,append=FALSE)
sum_anneli_sample_data_richness_GL1 <- as.data.frame(sum_anneli_sample_data_richness_GL1)
sum_anneli_sample_data_richness_GL1$group <- "Annelids"
sum_anneli_sample_data_richness_GL1$diversity <- "Richness"
sum_anneli_sample_data_richness_GL1$class <- "Land Cover GL1"
sum_anneli_sample_data_richness_GL1 <- sum_anneli_sample_data_richness_GL1[-1,]

sum_anneli_sample_data_shannon_GL1 <- capture.output(summary(sub_anneli_sample_data_GL1$Shannon),file=NULL,append=FALSE)
sum_anneli_sample_data_shannon_GL1 <- as.data.frame(sum_anneli_sample_data_shannon_GL1)
sum_anneli_sample_data_shannon_GL1$group <- "Annelids"
sum_anneli_sample_data_shannon_GL1$diversity <- "Shannon"
sum_anneli_sample_data_shannon_GL1$class <- "Land Cover GL1"
sum_anneli_sample_data_shannon_GL1 <-sum_anneli_sample_data_shannon_GL1[-1,]

sub_anneli_sample_data_GL2 <- a_anneli_sample_data %>% filter(LU2 == "GL2")
sum_anneli_sample_data_richness_GL2 <- capture.output(summary(sub_anneli_sample_data_GL2$Richness),file=NULL,append=FALSE)
sum_anneli_sample_data_richness_GL2 <- as.data.frame(sum_anneli_sample_data_richness_GL2)
sum_anneli_sample_data_richness_GL2$group <- "Annelids"
sum_anneli_sample_data_richness_GL2$diversity <- "Richness"
sum_anneli_sample_data_richness_GL2$class <- "Land Cover GL2"
sum_anneli_sample_data_richness_GL2 <- sum_anneli_sample_data_richness_GL2[-1,]

sum_anneli_sample_data_shannon_GL2 <- capture.output(summary(sub_anneli_sample_data_GL2$Shannon),file=NULL,append=FALSE)
sum_anneli_sample_data_shannon_GL2 <- as.data.frame(sum_anneli_sample_data_shannon_GL2)
sum_anneli_sample_data_shannon_GL2$group <- "Annelids"
sum_anneli_sample_data_shannon_GL2$diversity <- "Shannon"
sum_anneli_sample_data_shannon_GL2$class <- "Land Cover GL2"
sum_anneli_sample_data_shannon_GL2 <-sum_anneli_sample_data_shannon_GL2[-1,]

sub_anneli_sample_data_GL3 <- a_anneli_sample_data %>% filter(LU2 == "GL3")
sum_anneli_sample_data_richness_GL3 <- capture.output(summary(sub_anneli_sample_data_GL3$Richness),file=NULL,append=FALSE)
sum_anneli_sample_data_richness_GL3 <- as.data.frame(sum_anneli_sample_data_richness_GL3)
sum_anneli_sample_data_richness_GL3$group <- "Annelids"
sum_anneli_sample_data_richness_GL3$diversity <- "Richness"
sum_anneli_sample_data_richness_GL3$class <- "Land Cover GL3"
sum_anneli_sample_data_richness_GL3 <- sum_anneli_sample_data_richness_GL3[-1,]

sum_anneli_sample_data_shannon_GL3 <- capture.output(summary(sub_anneli_sample_data_GL3$Shannon),file=NULL,append=FALSE)
sum_anneli_sample_data_shannon_GL3 <- as.data.frame(sum_anneli_sample_data_shannon_GL3)
sum_anneli_sample_data_shannon_GL3$group <- "Annelids"
sum_anneli_sample_data_shannon_GL3$diversity <- "Shannon"
sum_anneli_sample_data_shannon_GL3$class <- "Land Cover GL3"
sum_anneli_sample_data_shannon_GL3 <-sum_anneli_sample_data_shannon_GL3[-1,]

sub_anneli_sample_data_WL1 <- a_anneli_sample_data %>% filter(LU2 == "WL1")
sum_anneli_sample_data_richness_WL1 <- capture.output(summary(sub_anneli_sample_data_WL1$Richness),file=NULL,append=FALSE)
sum_anneli_sample_data_richness_WL1 <- as.data.frame(sum_anneli_sample_data_richness_WL1)
sum_anneli_sample_data_richness_WL1$group <- "Annelids"
sum_anneli_sample_data_richness_WL1$diversity <- "Richness"
sum_anneli_sample_data_richness_WL1$class <- "Land Cover WL1"
sum_anneli_sample_data_richness_WL1 <- sum_anneli_sample_data_richness_WL1[-1,]

sum_anneli_sample_data_shannon_WL1 <- capture.output(summary(sub_anneli_sample_data_WL1$Shannon),file=NULL,append=FALSE)
sum_anneli_sample_data_shannon_WL1 <- as.data.frame(sum_anneli_sample_data_shannon_WL1)
sum_anneli_sample_data_shannon_WL1$group <- "Annelids"
sum_anneli_sample_data_shannon_WL1$diversity <- "Shannon"
sum_anneli_sample_data_shannon_WL1$class <- "Land Cover WL1"
sum_anneli_sample_data_shannon_WL1 <-sum_anneli_sample_data_shannon_WL1[-1,]

sub_anneli_sample_data_WL2 <- a_anneli_sample_data %>% filter(LU2 == "WL2")
sum_anneli_sample_data_richness_WL2 <- capture.output(summary(sub_anneli_sample_data_WL2$Richness),file=NULL,append=FALSE)
sum_anneli_sample_data_richness_WL2 <- as.data.frame(sum_anneli_sample_data_richness_WL2)
sum_anneli_sample_data_richness_WL2$group <- "Annelids"
sum_anneli_sample_data_richness_WL2$diversity <- "Richness"
sum_anneli_sample_data_richness_WL2$class <- "Land Cover WL2"
sum_anneli_sample_data_richness_WL2 <- sum_anneli_sample_data_richness_WL2[-1,]

sum_anneli_sample_data_shannon_WL2 <- capture.output(summary(sub_anneli_sample_data_WL2$Shannon),file=NULL,append=FALSE)
sum_anneli_sample_data_shannon_WL2 <- as.data.frame(sum_anneli_sample_data_shannon_WL2)
sum_anneli_sample_data_shannon_WL2$group <- "Annelids"
sum_anneli_sample_data_shannon_WL2$diversity <- "Shannon"
sum_anneli_sample_data_shannon_WL2$class <- "Land Cover WL2"
sum_anneli_sample_data_shannon_WL2 <-sum_anneli_sample_data_shannon_WL2[-1,]

sub_anneli_sample_data_WL3 <- a_anneli_sample_data %>% filter(LU2 == "WL3")
sum_anneli_sample_data_richness_WL3 <- capture.output(summary(sub_anneli_sample_data_WL3$Richness),file=NULL,append=FALSE)
sum_anneli_sample_data_richness_WL3 <- as.data.frame(sum_anneli_sample_data_richness_WL3)
sum_anneli_sample_data_richness_WL3$group <- "Annelids"
sum_anneli_sample_data_richness_WL3$diversity <- "Richness"
sum_anneli_sample_data_richness_WL3$class <- "Land Cover WL3"
sum_anneli_sample_data_richness_WL3 <- sum_anneli_sample_data_richness_WL3[-1,]

sum_anneli_sample_data_shannon_WL3 <- capture.output(summary(sub_anneli_sample_data_WL3$Shannon),file=NULL,append=FALSE)
sum_anneli_sample_data_shannon_WL3 <- as.data.frame(sum_anneli_sample_data_shannon_WL3)
sum_anneli_sample_data_shannon_WL3$group <- "Annelids"
sum_anneli_sample_data_shannon_WL3$diversity <- "Shannon"
sum_anneli_sample_data_shannon_WL3$class <- "Land Cover WL3"
sum_anneli_sample_data_shannon_WL3 <-sum_anneli_sample_data_shannon_WL3[-1,]

sub_anneli_sample_data_WL4 <- a_anneli_sample_data %>% filter(LU2 == "WL4")
sum_anneli_sample_data_richness_WL4 <- capture.output(summary(sub_anneli_sample_data_WL4$Richness),file=NULL,append=FALSE)
sum_anneli_sample_data_richness_WL4 <- as.data.frame(sum_anneli_sample_data_richness_WL4)
sum_anneli_sample_data_richness_WL4$group <- "Annelids"
sum_anneli_sample_data_richness_WL4$diversity <- "Richness"
sum_anneli_sample_data_richness_WL4$class <- "Land Cover WL4"
sum_anneli_sample_data_richness_WL4 <- sum_anneli_sample_data_richness_WL4[-1,]

sum_anneli_sample_data_shannon_WL4 <- capture.output(summary(sub_anneli_sample_data_WL4$Shannon),file=NULL,append=FALSE)
sum_anneli_sample_data_shannon_WL4 <- as.data.frame(sum_anneli_sample_data_shannon_WL4)
sum_anneli_sample_data_shannon_WL4$group <- "Annelids"
sum_anneli_sample_data_shannon_WL4$diversity <- "Shannon"
sum_anneli_sample_data_shannon_WL4$class <- "Land Cover WL4"
sum_anneli_sample_data_shannon_WL4 <-sum_anneli_sample_data_shannon_WL4[-1,]

sub_anneli_sample_data_No_C <- a_anneli_sample_data %>% filter(LU2 == "No_C")
sum_anneli_sample_data_richness_No_C <- capture.output(summary(sub_anneli_sample_data_No_C$Richness),file=NULL,append=FALSE)
sum_anneli_sample_data_richness_No_C <- as.data.frame(sum_anneli_sample_data_richness_No_C)
sum_anneli_sample_data_richness_No_C$group <- "Annelids"
sum_anneli_sample_data_richness_No_C$diversity <- "Richness"
sum_anneli_sample_data_richness_No_C$class <- "Land Cover No_C"
sum_anneli_sample_data_richness_No_C <- sum_anneli_sample_data_richness_No_C[-1,]

sum_anneli_sample_data_shannon_No_C <- capture.output(summary(sub_anneli_sample_data_No_C$Shannon),file=NULL,append=FALSE)
sum_anneli_sample_data_shannon_No_C <- as.data.frame(sum_anneli_sample_data_shannon_No_C)
sum_anneli_sample_data_shannon_No_C$group <- "Annelids"
sum_anneli_sample_data_shannon_No_C$diversity <- "Shannon"
sum_anneli_sample_data_shannon_No_C$class <- "Land Cover No_C"
sum_anneli_sample_data_shannon_No_C <-sum_anneli_sample_data_shannon_No_C[-1,]



sub_a_anneli_sample_data_CL_annual <- a_anneli_sample_data %>% filter(LC1_2018 == "CL_annual")
sum_a_anneli_richness_CL_annual <- capture.output(summary(sub_a_anneli_sample_data_CL_annual$Richness),file=NULL,append=FALSE)
sum_a_anneli_richness_CL_annual <- as.data.frame(sum_a_anneli_richness_CL_annual)
sum_a_anneli_richness_CL_annual$group <- "Annelids"
sum_a_anneli_richness_CL_annual$diversity <- "Richness"
sum_a_anneli_richness_CL_annual$class <- "Land Cover CL_annual"
sum_a_anneli_richness_CL_annual <-sum_a_anneli_richness_CL_annual[-1,]

sum_a_anneli_shannon_CL_annual <- capture.output(summary(sub_a_anneli_sample_data_CL_annual$Shannon),file=NULL,append=FALSE)
sum_a_anneli_shannon_CL_annual <- as.data.frame(sum_a_anneli_shannon_CL_annual)
sum_a_anneli_shannon_CL_annual$group <- "Annelids"
sum_a_anneli_shannon_CL_annual$diversity <- "Shannon"
sum_a_anneli_shannon_CL_annual$class <- "Land Cove CL_annualr"
sum_a_anneli_shannon_CL_annual <-sum_a_anneli_shannon_CL_annual[-1,]

sum_a_anneli_sample_data_CL_permanent <- a_anneli_sample_data %>% filter(LC1_2018 == "CL_permanent")
sum_a_anneli_richness_CL_permanent <- capture.output(summary(sum_a_anneli_sample_data_CL_permanent$Richness),file=NULL,append=FALSE)
sum_a_anneli_richness_CL_permanent <- as.data.frame(sum_a_anneli_richness_CL_permanent)
sum_a_anneli_richness_CL_permanent$group <- "Annelids"
sum_a_anneli_richness_CL_permanent$diversity <- "Richness"
sum_a_anneli_richness_CL_permanent$class <- "Land Cover CL_permanent"
sum_a_anneli_richness_CL_permanent <-sum_a_anneli_richness_CL_permanent[-1,]

sum_a_anneli_shannon_CL_permanent <- capture.output(summary(sum_a_anneli_sample_data_CL_permanent$Shannon),file=NULL,append=FALSE)
sum_a_anneli_shannon_CL_permanent <- as.data.frame(sum_a_anneli_shannon_CL_permanent)
sum_a_anneli_shannon_CL_permanent$group <- "Annelids"
sum_a_anneli_shannon_CL_permanent$diversity <- "Shannon"
sum_a_anneli_shannon_CL_permanent$class <- "Land Cover CL_permanent"
sum_a_anneli_shannon_CL_permanent <-sum_a_anneli_shannon_CL_permanent[-1,]

sub_a_anneli_sample_data_GL_unmanaged <- a_anneli_sample_data %>% filter(LC1_2018 == "GL_unmanaged")
sum_a_anneli_richness_GL_unmanaged <- capture.output(summary(sub_a_anneli_sample_data_GL_unmanaged$Richness),file=NULL,append=FALSE)
sum_a_anneli_richness_GL_unmanaged <- as.data.frame(sum_a_anneli_richness_GL_unmanaged)
sum_a_anneli_richness_GL_unmanaged$group <- "Annelids"
sum_a_anneli_richness_GL_unmanaged$diversity <- "Richness"
sum_a_anneli_richness_GL_unmanaged$class <- "Land Cover GL_unmanaged"
sum_a_anneli_richness_GL_unmanaged <-sum_a_anneli_richness_GL_unmanaged[-1,]

sum_a_anneli_shannon_GL_unmanaged <- capture.output(summary(sub_a_anneli_sample_data_GL_unmanaged$Shannon),file=NULL,append=FALSE)
sum_a_anneli_shannon_GL_unmanaged <- as.data.frame(sum_a_anneli_shannon_GL_unmanaged)
sum_a_anneli_shannon_GL_unmanaged$group <- "Annelids"
sum_a_anneli_shannon_GL_unmanaged$diversity <- "Shannon"
sum_a_anneli_shannon_GL_unmanaged$class <- "Land Cover GL_unmanaged"
sum_a_anneli_shannon_GL_unmanaged <-sum_a_anneli_shannon_GL_unmanaged[-1,]

sub_a_anneli_sample_data_GL_managed <- a_anneli_sample_data %>% filter(LC1_2018 == "GL_managed")
sum_a_anneli_richness_GL_managed <- capture.output(summary(sub_a_anneli_sample_data_GL_managed$Richness),file=NULL,append=FALSE)
sum_a_anneli_richness_GL_managed <- as.data.frame(sum_a_anneli_richness_GL_managed)
sum_a_anneli_richness_GL_managed$group <- "Annelids"
sum_a_anneli_richness_GL_managed$diversity <- "Richness"
sum_a_anneli_richness_GL_managed$class <- "Land Cover GL_managed"
sum_a_anneli_richness_GL_managed <-sum_a_anneli_richness_GL_managed[-1,]

sum_a_anneli_shannon_GL_managed <- capture.output(summary(sub_a_anneli_sample_data_GL_managed$Shannon),file=NULL,append=FALSE)
sum_a_anneli_shannon_GL_managed <- as.data.frame(sum_a_anneli_shannon_GL_managed)
sum_a_anneli_shannon_GL_managed$group <- "Annelids"
sum_a_anneli_shannon_GL_managed$diversity <- "Shannon"
sum_a_anneli_shannon_GL_managed$class <- "Land Cover GL_managed"
sum_a_anneli_shannon_GL_managed <-sum_a_anneli_shannon_GL_managed[-1,]

sub_a_anneli_sample_data_WL_broadleaved <- a_anneli_sample_data %>% filter(LC1_2018 == "WL_broadleaved")
sum_a_anneli_richness_WL_broadleaved <- capture.output(summary(sub_a_anneli_sample_data_WL_broadleaved$Richness),file=NULL,append=FALSE)
sum_a_anneli_richness_WL_broadleaved <- as.data.frame(sum_a_anneli_richness_WL_broadleaved)
sum_a_anneli_richness_WL_broadleaved$group <- "Annelids"
sum_a_anneli_richness_WL_broadleaved$diversity <- "Richness"
sum_a_anneli_richness_WL_broadleaved$class <- "Land Cover WL_broadleaved"
sum_a_anneli_richness_WL_broadleaved <-sum_a_anneli_richness_WL_broadleaved[-1,]

sum_a_anneli_shannon_WL_broadleaved <- capture.output(summary(sub_a_anneli_sample_data_WL_broadleaved$Shannon),file=NULL,append=FALSE)
sum_a_anneli_shannon_WL_broadleaved <- as.data.frame(sum_a_anneli_shannon_WL_broadleaved)
sum_a_anneli_shannon_WL_broadleaved$group <- "Annelids"
sum_a_anneli_shannon_WL_broadleaved$diversity <- "Shannon"
sum_a_anneli_shannon_WL_broadleaved$class <- "Land Cover WL_broadleaved"
sum_a_anneli_shannon_WL_broadleaved <-sum_a_anneli_shannon_WL_broadleaved[-1,]

sub_a_anneli_sample_data_WL_coniferous <- a_anneli_sample_data %>% filter(LC1_2018 == "WL_coniferous")
sum_a_anneli_richness_WL_coniferous <- capture.output(summary(sub_a_anneli_sample_data_WL_coniferous$Richness),file=NULL,append=FALSE)
sum_a_anneli_richness_WL_coniferous <- as.data.frame(sum_a_anneli_richness_WL_coniferous)
sum_a_anneli_richness_WL_coniferous$group <- "Annelids"
sum_a_anneli_richness_WL_coniferous$diversity <- "Richness"
sum_a_anneli_richness_WL_coniferous$class <- "Land Cover WL_coniferous"
sum_a_anneli_richness_WL_coniferous <-sum_a_anneli_richness_WL_coniferous[-1,]

sum_a_anneli_shannon_WL_coniferous <- capture.output(summary(sub_a_anneli_sample_data_WL_coniferous$Shannon),file=NULL,append=FALSE)
sum_a_anneli_shannon_WL_coniferous <- as.data.frame(sum_a_anneli_shannon_WL_coniferous)
sum_a_anneli_shannon_WL_coniferous$group <- "Annelids"
sum_a_anneli_shannon_WL_coniferous$diversity <- "Shannon"
sum_a_anneli_shannon_WL_coniferous$class <- "Land Cover WL_coniferous"
sum_a_anneli_shannon_WL_coniferous <-sum_a_anneli_shannon_WL_coniferous[-1,]

sub_a_anneli_sample_data_temp_1 <- a_anneli_sample_data %>% filter(cluster_temp == "Very low")
sum_a_anneli_sample_data_shannon_temp_1 <- capture.output(summary(sub_a_anneli_sample_data_temp_1$Shannon),file=NULL,append=FALSE)
sum_a_anneli_sample_data_shannon_temp_1 <- as.data.frame(sum_a_anneli_sample_data_shannon_temp_1)
sum_a_anneli_sample_data_shannon_temp_1$group <- "Annelids"
sum_a_anneli_sample_data_shannon_temp_1$diversity <- "Shannon"
sum_a_anneli_sample_data_shannon_temp_1$class <- "Very low temperature"
sum_a_anneli_sample_data_shannon_temp_1 <-sum_a_anneli_sample_data_shannon_temp_1[-1,]

sum_a_anneli_sample_data_richness_temp_1 <- capture.output(summary(sub_a_anneli_sample_data_temp_1$Richness),file=NULL,append=FALSE)
sum_a_anneli_sample_data_richness_temp_1 <- as.data.frame(sum_a_anneli_sample_data_richness_temp_1)
sum_a_anneli_sample_data_richness_temp_1$group <- "Annelids"
sum_a_anneli_sample_data_richness_temp_1$diversity <- "Richness"
sum_a_anneli_sample_data_richness_temp_1$class <- "Very low temperature"
sum_a_anneli_sample_data_richness_temp_1 <-sum_a_anneli_sample_data_richness_temp_1[-1,]

sub_a_anneli_sample_data_temp_2 <- a_anneli_sample_data %>% filter(cluster_temp == "Low")
sum_a_anneli_sample_data_shannon_temp_2 <- capture.output(summary(sub_a_anneli_sample_data_temp_2$Shannon),file=NULL,append=FALSE)
sum_a_anneli_sample_data_shannon_temp_2 <- as.data.frame(sum_a_anneli_sample_data_shannon_temp_2)
sum_a_anneli_sample_data_shannon_temp_2$group <- "Annelids"
sum_a_anneli_sample_data_shannon_temp_2$diversity <- "Shannon"
sum_a_anneli_sample_data_shannon_temp_2$class <- "Low temperature"
sum_a_anneli_sample_data_shannon_temp_2 <-sum_a_anneli_sample_data_shannon_temp_2[-1,]

sum_a_anneli_sample_data_richness_temp_2 <- capture.output(summary(sub_a_anneli_sample_data_temp_2$Richness),file=NULL,append=FALSE)
sum_a_anneli_sample_data_richness_temp_2 <- as.data.frame(sum_a_anneli_sample_data_richness_temp_2)
sum_a_anneli_sample_data_richness_temp_2$group <- "Annelids"
sum_a_anneli_sample_data_richness_temp_2$diversity <- "Richness"
sum_a_anneli_sample_data_richness_temp_2$class <- "Low temperature"
sum_a_anneli_sample_data_richness_temp_2 <-sum_a_anneli_sample_data_richness_temp_2[-1,]

sub_a_anneli_sample_data_temp_3 <- a_anneli_sample_data %>% filter(cluster_temp == "Medium low")
sum_a_anneli_sample_data_shannon_temp_3 <- capture.output(summary(sub_a_anneli_sample_data_temp_3$Shannon),file=NULL,append=FALSE)
sum_a_anneli_sample_data_shannon_temp_3 <- as.data.frame(sum_a_anneli_sample_data_shannon_temp_3)
sum_a_anneli_sample_data_shannon_temp_3$group <- "Annelids"
sum_a_anneli_sample_data_shannon_temp_3$diversity <- "Shannon"
sum_a_anneli_sample_data_shannon_temp_3$class <- "Medium low temperature"
sum_a_anneli_sample_data_shannon_temp_3 <-sum_a_anneli_sample_data_shannon_temp_3[-1,]

sum_a_anneli_sample_data_richness_temp_3 <- capture.output(summary(sub_a_anneli_sample_data_temp_3$Richness),file=NULL,append=FALSE)
sum_a_anneli_sample_data_richness_temp_3 <- as.data.frame(sum_a_anneli_sample_data_richness_temp_3)
sum_a_anneli_sample_data_richness_temp_3$group <- "Annelids"
sum_a_anneli_sample_data_richness_temp_3$diversity <- "Richness"
sum_a_anneli_sample_data_richness_temp_3$class <- "Medium low temperature"
sum_a_anneli_sample_data_richness_temp_3 <-sum_a_anneli_sample_data_richness_temp_3[-1,]

sub_a_anneli_sample_data_temp_4 <- a_anneli_sample_data %>% filter(cluster_temp == "Medium high")
sum_a_anneli_sample_data_shannon_temp_4 <- capture.output(summary(sub_a_anneli_sample_data_temp_4$Shannon),file=NULL,append=FALSE)
sum_a_anneli_sample_data_shannon_temp_4 <- as.data.frame(sum_a_anneli_sample_data_shannon_temp_4)
sum_a_anneli_sample_data_shannon_temp_4$group <- "Annelids"
sum_a_anneli_sample_data_shannon_temp_4$diversity <- "Shannon"
sum_a_anneli_sample_data_shannon_temp_4$class <- "Medium high temperature"
sum_a_anneli_sample_data_shannon_temp_4 <-sum_a_anneli_sample_data_shannon_temp_4[-1,]

sum_a_anneli_sample_data_richness_temp_4 <- capture.output(summary(sub_a_anneli_sample_data_temp_4$Richness),file=NULL,append=FALSE)
sum_a_anneli_sample_data_richness_temp_4 <- as.data.frame(sum_a_anneli_sample_data_richness_temp_4)
sum_a_anneli_sample_data_richness_temp_4$group <- "Annelids"
sum_a_anneli_sample_data_richness_temp_4$diversity <- "Richness"
sum_a_anneli_sample_data_richness_temp_4$class <- "Medium high temperature"
sum_a_anneli_sample_data_richness_temp_4 <-sum_a_anneli_sample_data_richness_temp_4[-1,]

sub_a_anneli_sample_data_temp_5 <- a_anneli_sample_data %>% filter(cluster_temp == "High")
sum_a_anneli_sample_data_shannon_temp_5 <- capture.output(summary(sub_a_anneli_sample_data_temp_5$Shannon),file=NULL,append=FALSE)
sum_a_anneli_sample_data_shannon_temp_5 <- as.data.frame(sum_a_anneli_sample_data_shannon_temp_5)
sum_a_anneli_sample_data_shannon_temp_5$group <- "Annelids"
sum_a_anneli_sample_data_shannon_temp_5$diversity <- "Shannon"
sum_a_anneli_sample_data_shannon_temp_5$class <- "High temperature"
sum_a_anneli_sample_data_shannon_temp_5 <-sum_a_anneli_sample_data_shannon_temp_5[-1,]

sum_a_anneli_sample_data_richness_temp_5 <- capture.output(summary(sub_a_anneli_sample_data_temp_5$Richness),file=NULL,append=FALSE)
sum_a_anneli_sample_data_richness_temp_5 <- as.data.frame(sum_a_anneli_sample_data_richness_temp_5)
sum_a_anneli_sample_data_richness_temp_5$group <- "Annelids"
sum_a_anneli_sample_data_richness_temp_5$diversity <- "Richness"
sum_a_anneli_sample_data_richness_temp_5$class <- "High temperature"
sum_a_anneli_sample_data_richness_temp_5 <-sum_a_anneli_sample_data_richness_temp_5[-1,]

sub_a_anneli_sample_data_temp_6 <- a_anneli_sample_data %>% filter(cluster_temp == "Very high")
sum_a_anneli_sample_data_shannon_temp_6 <- capture.output(summary(sub_a_anneli_sample_data_temp_6$Shannon),file=NULL,append=FALSE)
sum_a_anneli_sample_data_shannon_temp_6 <- as.data.frame(sum_a_anneli_sample_data_shannon_temp_6)
sum_a_anneli_sample_data_shannon_temp_6$group <- "Annelids"
sum_a_anneli_sample_data_shannon_temp_6$diversity <- "Shannon"
sum_a_anneli_sample_data_shannon_temp_6$class <- "Very high temperature"
sum_a_anneli_sample_data_shannon_temp_6 <-sum_a_anneli_sample_data_shannon_temp_6[-1,]

sum_a_anneli_sample_data_richness_temp_6 <- capture.output(summary(sub_a_anneli_sample_data_temp_6$Richness),file=NULL,append=FALSE)
sum_a_anneli_sample_data_richness_temp_6 <- as.data.frame(sum_a_anneli_sample_data_richness_temp_6)
sum_a_anneli_sample_data_richness_temp_6$group <- "Annelids"
sum_a_anneli_sample_data_richness_temp_6$diversity <- "Richness"
sum_a_anneli_sample_data_richness_temp_6$class <- "Very high temperature"
sum_a_anneli_sample_data_richness_temp_6 <-sum_a_anneli_sample_data_richness_temp_6[-1,] 

sub_a_anneli_sample_data_temp_range_1 <- a_anneli_sample_data %>% filter(cluster_temp_range == "Very low")
sum_a_anneli_sample_data_shannon_temp_range_1 <- capture.output(summary(sub_a_anneli_sample_data_temp_range_1$Shannon),file=NULL,append=FALSE)
sum_a_anneli_sample_data_shannon_temp_range_1 <- as.data.frame(sum_a_anneli_sample_data_shannon_temp_range_1)
sum_a_anneli_sample_data_shannon_temp_range_1$group <- "Annelids"
sum_a_anneli_sample_data_shannon_temp_range_1$diversity <- "Shannon"
sum_a_anneli_sample_data_shannon_temp_range_1$class <- "Very low temperature range"
sum_a_anneli_sample_data_shannon_temp_range_1 <-sum_a_anneli_sample_data_shannon_temp_range_1[-1,]

sum_a_anneli_sample_data_richness_temp_range_1 <- capture.output(summary(sub_a_anneli_sample_data_temp_range_1$Richness),file=NULL,append=FALSE)
sum_a_anneli_sample_data_richness_temp_range_1 <- as.data.frame(sum_a_anneli_sample_data_richness_temp_range_1)
sum_a_anneli_sample_data_richness_temp_range_1$group <- "Annelids"
sum_a_anneli_sample_data_richness_temp_range_1$diversity <- "Richness"
sum_a_anneli_sample_data_richness_temp_range_1$class <- "Very low temperature range"
sum_a_anneli_sample_data_richness_temp_range_1 <-sum_a_anneli_sample_data_richness_temp_range_1[-1,]

sub_a_anneli_sample_data_temp_range_2 <- a_anneli_sample_data %>% filter(cluster_temp_range == "Low")
sum_a_anneli_sample_data_shannon_temp_range_2 <- capture.output(summary(sub_a_anneli_sample_data_temp_range_2$Shannon),file=NULL,append=FALSE)
sum_a_anneli_sample_data_shannon_temp_range_2 <- as.data.frame(sum_a_anneli_sample_data_shannon_temp_range_2)
sum_a_anneli_sample_data_shannon_temp_range_2$group <- "Annelids"
sum_a_anneli_sample_data_shannon_temp_range_2$diversity <- "Shannon"
sum_a_anneli_sample_data_shannon_temp_range_2$class <- "Low temperature range"
sum_a_anneli_sample_data_shannon_temp_range_2 <-sum_a_anneli_sample_data_shannon_temp_range_2[-1,]

sum_a_anneli_sample_data_richness_temp_range_2 <- capture.output(summary(sub_a_anneli_sample_data_temp_range_2$Richness),file=NULL,append=FALSE)
sum_a_anneli_sample_data_richness_temp_range_2 <- as.data.frame(sum_a_anneli_sample_data_richness_temp_range_2)
sum_a_anneli_sample_data_richness_temp_range_2$group <- "Annelids"
sum_a_anneli_sample_data_richness_temp_range_2$diversity <- "Richness"
sum_a_anneli_sample_data_richness_temp_range_2$class <- "Low temperature range"
sum_a_anneli_sample_data_richness_temp_range_2 <-sum_a_anneli_sample_data_richness_temp_range_2[-1,]

sub_a_anneli_sample_data_temp_range_3 <- a_anneli_sample_data %>% filter(cluster_temp_range == "Medium low")
sum_a_anneli_sample_data_shannon_temp_range_3 <- capture.output(summary(sub_a_anneli_sample_data_temp_range_3$Shannon),file=NULL,append=FALSE)
sum_a_anneli_sample_data_shannon_temp_range_3 <- as.data.frame(sum_a_anneli_sample_data_shannon_temp_range_3)
sum_a_anneli_sample_data_shannon_temp_range_3$group <- "Annelids"
sum_a_anneli_sample_data_shannon_temp_range_3$diversity <- "Shannon"
sum_a_anneli_sample_data_shannon_temp_range_3$class <- "Medium low temperature range"
sum_a_anneli_sample_data_shannon_temp_range_3 <-sum_a_anneli_sample_data_shannon_temp_range_3[-1,]

sum_a_anneli_sample_data_richness_temp_range_3 <- capture.output(summary(sub_a_anneli_sample_data_temp_range_3$Richness),file=NULL,append=FALSE)
sum_a_anneli_sample_data_richness_temp_range_3 <- as.data.frame(sum_a_anneli_sample_data_richness_temp_range_3)
sum_a_anneli_sample_data_richness_temp_range_3$group <- "Annelids"
sum_a_anneli_sample_data_richness_temp_range_3$diversity <- "Richness"
sum_a_anneli_sample_data_richness_temp_range_3$class <- "Medium low temperature range"
sum_a_anneli_sample_data_richness_temp_range_3 <-sum_a_anneli_sample_data_richness_temp_range_3[-1,]

sub_a_anneli_sample_data_temp_range_4 <- a_anneli_sample_data %>% filter(cluster_temp_range == "Medium high")
sum_a_anneli_sample_data_shannon_temp_range_4 <- capture.output(summary(sub_a_anneli_sample_data_temp_range_4$Shannon),file=NULL,append=FALSE)
sum_a_anneli_sample_data_shannon_temp_range_4 <- as.data.frame(sum_a_anneli_sample_data_shannon_temp_range_4)
sum_a_anneli_sample_data_shannon_temp_range_4$group <- "Annelids"
sum_a_anneli_sample_data_shannon_temp_range_4$diversity <- "Shannon"
sum_a_anneli_sample_data_shannon_temp_range_4$class <- "Medium high temperature range"
sum_a_anneli_sample_data_shannon_temp_range_4 <-sum_a_anneli_sample_data_shannon_temp_range_4[-1,]

sum_a_anneli_sample_data_richness_temp_range_4 <- capture.output(summary(sub_a_anneli_sample_data_temp_range_4$Richness),file=NULL,append=FALSE)
sum_a_anneli_sample_data_richness_temp_range_4 <- as.data.frame(sum_a_anneli_sample_data_richness_temp_range_4)
sum_a_anneli_sample_data_richness_temp_range_4$group <- "Annelids"
sum_a_anneli_sample_data_richness_temp_range_4$diversity <- "Richness"
sum_a_anneli_sample_data_richness_temp_range_4$class <- "Medium high temperature range"
sum_a_anneli_sample_data_richness_temp_range_4 <-sum_a_anneli_sample_data_richness_temp_range_4[-1,]

sub_a_anneli_sample_data_temp_range_5 <- a_anneli_sample_data %>% filter(cluster_temp_range == "High")
sum_a_anneli_sample_data_shannon_temp_range_5 <- capture.output(summary(sub_a_anneli_sample_data_temp_range_5$Shannon),file=NULL,append=FALSE)
sum_a_anneli_sample_data_shannon_temp_range_5 <- as.data.frame(sum_a_anneli_sample_data_shannon_temp_range_5)
sum_a_anneli_sample_data_shannon_temp_range_5$group <- "Annelids"
sum_a_anneli_sample_data_shannon_temp_range_5$diversity <- "Shannon"
sum_a_anneli_sample_data_shannon_temp_range_5$class <- "High temperature range"
sum_a_anneli_sample_data_shannon_temp_range_5 <-sum_a_anneli_sample_data_shannon_temp_range_5[-1,]

sum_a_anneli_sample_data_richness_temp_range_5 <- capture.output(summary(sub_a_anneli_sample_data_temp_range_5$Richness),file=NULL,append=FALSE)
sum_a_anneli_sample_data_richness_temp_range_5 <- as.data.frame(sum_a_anneli_sample_data_richness_temp_range_5)
sum_a_anneli_sample_data_richness_temp_range_5$group <- "Annelids"
sum_a_anneli_sample_data_richness_temp_range_5$diversity <- "Richness"
sum_a_anneli_sample_data_richness_temp_range_5$class <- "High temperature range"
sum_a_anneli_sample_data_richness_temp_range_5 <-sum_a_anneli_sample_data_richness_temp_range_5[-1,]

sub_a_anneli_sample_data_temp_range_6 <- a_anneli_sample_data %>% filter(cluster_temp_range == "Very high")
sum_a_anneli_sample_data_shannon_temp_range_6 <- capture.output(summary(sub_a_anneli_sample_data_temp_range_6$Shannon),file=NULL,append=FALSE)
sum_a_anneli_sample_data_shannon_temp_range_6 <- as.data.frame(sum_a_anneli_sample_data_shannon_temp_range_6)
sum_a_anneli_sample_data_shannon_temp_range_6$group <- "Annelids"
sum_a_anneli_sample_data_shannon_temp_range_6$diversity <- "Shannon"
sum_a_anneli_sample_data_shannon_temp_range_6$class <- "Very high temperature range"
sum_a_anneli_sample_data_shannon_temp_range_6 <-sum_a_anneli_sample_data_shannon_temp_range_6[-1,]

sum_a_anneli_sample_data_richness_temp_range_6 <- capture.output(summary(sub_a_anneli_sample_data_temp_range_6$Richness),file=NULL,append=FALSE)
sum_a_anneli_sample_data_richness_temp_range_6 <- as.data.frame(sum_a_anneli_sample_data_richness_temp_range_6)
sum_a_anneli_sample_data_richness_temp_range_6$group <- "Annelids"
sum_a_anneli_sample_data_richness_temp_range_6$diversity <- "Richness"
sum_a_anneli_sample_data_richness_temp_range_6$class <- "Very high temperature range"
sum_a_anneli_sample_data_richness_temp_range_6 <-sum_a_anneli_sample_data_richness_temp_range_6[-1,]

sub_a_anneli_sample_data_prec_1 <- a_anneli_sample_data %>% filter(cluster_prec == "Very low")
sum_a_anneli_sample_data_shannon_prec_1 <- capture.output(summary(sub_a_anneli_sample_data_prec_1$Shannon),file=NULL,append=FALSE)
sum_a_anneli_sample_data_shannon_prec_1 <- as.data.frame(sum_a_anneli_sample_data_shannon_prec_1)
sum_a_anneli_sample_data_shannon_prec_1$group <- "Annelids"
sum_a_anneli_sample_data_shannon_prec_1$diversity <- "Shannon"
sum_a_anneli_sample_data_shannon_prec_1$class <- "Very low precipitation "
sum_a_anneli_sample_data_shannon_prec_1 <-sum_a_anneli_sample_data_shannon_prec_1[-1,]

sum_a_anneli_sample_data_richness_prec_1 <- capture.output(summary(sub_a_anneli_sample_data_prec_1$Richness),file=NULL,append=FALSE)
sum_a_anneli_sample_data_richness_prec_1 <- as.data.frame(sum_a_anneli_sample_data_richness_prec_1)
sum_a_anneli_sample_data_richness_prec_1$group <- "Annelids"
sum_a_anneli_sample_data_richness_prec_1$diversity <- "Richness"
sum_a_anneli_sample_data_richness_prec_1$class <- "Very low precipitation "
sum_a_anneli_sample_data_richness_prec_1 <-sum_a_anneli_sample_data_richness_prec_1[-1,]

sub_a_anneli_sample_data_prec_2 <- a_anneli_sample_data %>% filter(cluster_prec == "Low")
sum_a_anneli_sample_data_shannon_prec_2 <- capture.output(summary(sub_a_anneli_sample_data_prec_2$Shannon),file=NULL,append=FALSE)
sum_a_anneli_sample_data_shannon_prec_2 <- as.data.frame(sum_a_anneli_sample_data_shannon_prec_2)
sum_a_anneli_sample_data_shannon_prec_2$group <- "Annelids"
sum_a_anneli_sample_data_shannon_prec_2$diversity <- "Shannon"
sum_a_anneli_sample_data_shannon_prec_2$class <- "Low precipitation"
sum_a_anneli_sample_data_shannon_prec_2 <-sum_a_anneli_sample_data_shannon_prec_2[-1,]

sum_a_anneli_sample_data_richness_prec_2 <- capture.output(summary(sub_a_anneli_sample_data_prec_2$Richness),file=NULL,append=FALSE)
sum_a_anneli_sample_data_richness_prec_2 <- as.data.frame(sum_a_anneli_sample_data_richness_prec_2)
sum_a_anneli_sample_data_richness_prec_2$group <- "Annelids"
sum_a_anneli_sample_data_richness_prec_2$diversity <- "Richness"
sum_a_anneli_sample_data_richness_prec_2$class <- "Low precipitation"
sum_a_anneli_sample_data_richness_prec_2 <-sum_a_anneli_sample_data_richness_prec_2[-1,]

sub_a_anneli_sample_data_prec_3 <- a_anneli_sample_data %>% filter(cluster_prec == "Medium low")
sum_a_anneli_sample_data_shannon_prec_3 <- capture.output(summary(sub_a_anneli_sample_data_prec_3$Shannon),file=NULL,append=FALSE)
sum_a_anneli_sample_data_shannon_prec_3 <- as.data.frame(sum_a_anneli_sample_data_shannon_prec_3)
sum_a_anneli_sample_data_shannon_prec_3$group <- "Annelids"
sum_a_anneli_sample_data_shannon_prec_3$diversity <- "Shannon"
sum_a_anneli_sample_data_shannon_prec_3$class <- "Medium low precipitation"
sum_a_anneli_sample_data_shannon_prec_3 <-sum_a_anneli_sample_data_shannon_prec_3[-1,]

sum_a_anneli_sample_data_richness_prec_3 <- capture.output(summary(sub_a_anneli_sample_data_prec_3$Richness),file=NULL,append=FALSE)
sum_a_anneli_sample_data_richness_prec_3 <- as.data.frame(sum_a_anneli_sample_data_richness_prec_3)
sum_a_anneli_sample_data_richness_prec_3$group <- "Annelids"
sum_a_anneli_sample_data_richness_prec_3$diversity <- "Richness"
sum_a_anneli_sample_data_richness_prec_3$class <- "Medium low precipitation"
sum_a_anneli_sample_data_richness_prec_3 <-sum_a_anneli_sample_data_richness_prec_3[-1,]

sub_a_anneli_sample_data_prec_4 <- a_anneli_sample_data %>% filter(cluster_prec == "Medium high")
sum_a_anneli_sample_data_shannon_prec_4 <- capture.output(summary(sub_a_anneli_sample_data_prec_4$Shannon),file=NULL,append=FALSE)
sum_a_anneli_sample_data_shannon_prec_4 <- as.data.frame(sum_a_anneli_sample_data_shannon_prec_4)
sum_a_anneli_sample_data_shannon_prec_4$group <- "Annelids"
sum_a_anneli_sample_data_shannon_prec_4$diversity <- "Shannon"
sum_a_anneli_sample_data_shannon_prec_4$class <- "Medium high precipitation"
sum_a_anneli_sample_data_shannon_prec_4 <-sum_a_anneli_sample_data_shannon_prec_4[-1,]

sum_a_anneli_sample_data_richness_prec_4 <- capture.output(summary(sub_a_anneli_sample_data_prec_4$Richness),file=NULL,append=FALSE)
sum_a_anneli_sample_data_richness_prec_4 <- as.data.frame(sum_a_anneli_sample_data_richness_prec_4)
sum_a_anneli_sample_data_richness_prec_4$group <- "Annelids"
sum_a_anneli_sample_data_richness_prec_4$diversity <- "Richness"
sum_a_anneli_sample_data_richness_prec_4$class <- "Medium high precipitation"
sum_a_anneli_sample_data_richness_prec_4 <-sum_a_anneli_sample_data_richness_prec_4[-1,]

sub_a_anneli_sample_data_prec_5 <- a_anneli_sample_data %>% filter(cluster_prec == "High")
sum_a_anneli_sample_data_shannon_prec_5 <- capture.output(summary(sub_a_anneli_sample_data_prec_5$Shannon),file=NULL,append=FALSE)
sum_a_anneli_sample_data_shannon_prec_5 <- as.data.frame(sum_a_anneli_sample_data_shannon_prec_5)
sum_a_anneli_sample_data_shannon_prec_5$group <- "Annelids"
sum_a_anneli_sample_data_shannon_prec_5$diversity <- "Shannon"
sum_a_anneli_sample_data_shannon_prec_5$class <- "High precipitation"
sum_a_anneli_sample_data_shannon_prec_5 <-sum_a_anneli_sample_data_shannon_prec_5[-1,]

sum_a_anneli_sample_data_richness_prec_5 <- capture.output(summary(sub_a_anneli_sample_data_prec_5$Richness),file=NULL,append=FALSE)
sum_a_anneli_sample_data_richness_prec_5 <- as.data.frame(sum_a_anneli_sample_data_richness_prec_5)
sum_a_anneli_sample_data_richness_prec_5$group <- "Annelids"
sum_a_anneli_sample_data_richness_prec_5$diversity <- "Richness"
sum_a_anneli_sample_data_richness_prec_5$class <- "High precipitation"
sum_a_anneli_sample_data_richness_prec_5 <-sum_a_anneli_sample_data_richness_prec_5[-1,]

sub_a_anneli_sample_data_prec_6 <- a_anneli_sample_data %>% filter(cluster_prec == "Very high")
sum_a_anneli_sample_data_shannon_prec_6 <- capture.output(summary(sub_a_anneli_sample_data_prec_6$Shannon),file=NULL,append=FALSE)
sum_a_anneli_sample_data_shannon_prec_6 <- as.data.frame(sum_a_anneli_sample_data_shannon_prec_6)
sum_a_anneli_sample_data_shannon_prec_6$group <- "Annelids"
sum_a_anneli_sample_data_shannon_prec_6$diversity <- "Shannon"
sum_a_anneli_sample_data_shannon_prec_6$class <- "Very high precipitation"
sum_a_anneli_sample_data_shannon_prec_6 <-sum_a_anneli_sample_data_shannon_prec_6[-1,]

sum_a_anneli_sample_data_richness_prec_6 <- capture.output(summary(sub_a_anneli_sample_data_prec_6$Richness),file=NULL,append=FALSE)
sum_a_anneli_sample_data_richness_prec_6 <- as.data.frame(sum_a_anneli_sample_data_richness_prec_6)
sum_a_anneli_sample_data_richness_prec_6$group <- "Annelids"
sum_a_anneli_sample_data_richness_prec_6$diversity <- "Richness"
sum_a_anneli_sample_data_richness_prec_6$class <- "Very high precipitation"
sum_a_anneli_sample_data_richness_prec_6 <-sum_a_anneli_sample_data_richness_prec_6[-1,]











#rename first column
names(sum_a_nema_sample_data_richness_prec_1)[1] <- "Summary"
names(sum_a_nema_sample_data_richness_prec_2)[1] <- "Summary"
names(sum_a_nema_sample_data_richness_prec_3)[1] <- "Summary"
names(sum_a_nema_sample_data_richness_prec_4)[1] <- "Summary"
names(sum_a_nema_sample_data_richness_prec_5)[1] <- "Summary"
names(sum_a_nema_sample_data_richness_prec_6)[1] <- "Summary"
names(sum_a_nema_sample_data_shannon_prec_1)[1] <- "Summary"
names(sum_a_nema_sample_data_shannon_prec_2)[1] <- "Summary"
names(sum_a_nema_sample_data_shannon_prec_3)[1] <- "Summary"
names(sum_a_nema_sample_data_shannon_prec_4)[1] <- "Summary"
names(sum_a_nema_sample_data_shannon_prec_5)[1] <- "Summary"
names(sum_a_nema_sample_data_shannon_prec_6)[1] <- "Summary"

names(sum_a_roti_sample_data_richness_prec_1)[1] <- "Summary"
names(sum_a_roti_sample_data_richness_prec_2)[1] <- "Summary"
names(sum_a_roti_sample_data_richness_prec_3)[1] <- "Summary"
names(sum_a_roti_sample_data_richness_prec_4)[1] <- "Summary"
names(sum_a_roti_sample_data_richness_prec_5)[1] <- "Summary"
names(sum_a_roti_sample_data_richness_prec_6)[1] <- "Summary"
names(sum_a_roti_sample_data_shannon_prec_1)[1] <- "Summary"
names(sum_a_roti_sample_data_shannon_prec_2)[1] <- "Summary"
names(sum_a_roti_sample_data_shannon_prec_3)[1] <- "Summary"
names(sum_a_roti_sample_data_shannon_prec_4)[1] <- "Summary"
names(sum_a_roti_sample_data_shannon_prec_5)[1] <- "Summary"
names(sum_a_roti_sample_data_shannon_prec_6)[1] <- "Summary"

names(sum_a_tardi_sample_data_shannon_prec_1)[1] <- "Summary"
names(sum_a_tardi_sample_data_shannon_prec_2)[1] <- "Summary"
names(sum_a_tardi_sample_data_shannon_prec_3)[1] <- "Summary"
names(sum_a_tardi_sample_data_shannon_prec_4)[1] <- "Summary"
names(sum_a_tardi_sample_data_shannon_prec_5)[1] <- "Summary"
names(sum_a_tardi_sample_data_shannon_prec_6)[1] <- "Summary"

names(sum_f_sample_data_richness_prec_1)[1] <- "Summary"
names(sum_f_sample_data_richness_prec_2)[1] <- "Summary"
names(sum_f_sample_data_richness_prec_3)[1] <- "Summary"
names(sum_f_sample_data_richness_prec_4)[1] <- "Summary"
names(sum_f_sample_data_richness_prec_5)[1] <- "Summary"
names(sum_f_sample_data_richness_prec_6)[1] <- "Summary"
names(sum_f_sample_data_shannon_prec_1)[1] <- "Summary"
names(sum_f_sample_data_shannon_prec_2)[1] <- "Summary"
names(sum_f_sample_data_shannon_prec_3)[1] <- "Summary"
names(sum_f_sample_data_shannon_prec_4)[1] <- "Summary"
names(sum_f_sample_data_shannon_prec_5)[1] <- "Summary"
names(sum_f_sample_data_shannon_prec_6)[1] <- "Summary"

names(sum_a_tardi_sample_data_richness_prec_1)[1] <- "Summary"
names(sum_a_tardi_sample_data_richness_prec_2)[1] <- "Summary"
names(sum_a_tardi_sample_data_richness_prec_3)[1] <- "Summary"
names(sum_a_tardi_sample_data_richness_prec_4)[1] <- "Summary"
names(sum_a_tardi_sample_data_richness_prec_5)[1] <- "Summary"
names(sum_a_tardi_sample_data_richness_prec_6)[1] <- "Summary"
names(sum_a_tardi_sample_data_shannon_prec_1)[1] <- "Summary"
names(sum_a_tardi_sample_data_shannon_prec_2)[1] <- "Summary"
names(sum_a_tardi_sample_data_shannon_prec_3)[1] <- "Summary"
names(sum_a_tardi_sample_data_shannon_prec_4)[1] <- "Summary"
names(sum_a_tardi_sample_data_shannon_prec_5)[1] <- "Summary"
names(sum_a_tardi_sample_data_shannon_prec_6)[1] <- "Summary"

names(sum_f_sample_data_richness_temp_1)[1] <- "Summary"
names(sum_f_sample_data_richness_temp_2)[1] <- "Summary"
names(sum_f_sample_data_richness_temp_3)[1] <- "Summary"
names(sum_f_sample_data_richness_temp_4)[1] <- "Summary"
names(sum_f_sample_data_richness_temp_5)[1] <- "Summary"
names(sum_f_sample_data_richness_temp_6)[1] <- "Summary"
names(sum_f_sample_data_shannon_temp_1)[1] <- "Summary"
names(sum_f_sample_data_shannon_temp_2)[1] <- "Summary"
names(sum_f_sample_data_shannon_temp_3)[1] <- "Summary"
names(sum_f_sample_data_shannon_temp_4)[1] <- "Summary"
names(sum_f_sample_data_shannon_temp_5)[1] <- "Summary"
names(sum_f_sample_data_shannon_temp_6)[1] <- "Summary"

names(sum_p_sample_data_richness_temp_1)[1] <- "Summary"
names(sum_p_sample_data_richness_temp_2)[1] <- "Summary"
names(sum_p_sample_data_richness_temp_3)[1] <- "Summary"
names(sum_p_sample_data_richness_temp_4)[1] <- "Summary"
names(sum_p_sample_data_richness_temp_5)[1] <- "Summary"
names(sum_p_sample_data_richness_temp_6)[1] <- "Summary"
names(sum_p_sample_data_shannon_temp_1)[1] <- "Summary"
names(sum_p_sample_data_shannon_temp_2)[1] <- "Summary"
names(sum_p_sample_data_shannon_temp_3)[1] <- "Summary"
names(sum_p_sample_data_shannon_temp_4)[1] <- "Summary"
names(sum_p_sample_data_shannon_temp_5)[1] <- "Summary"
names(sum_p_sample_data_shannon_temp_6)[1] <- "Summary"

names(sum_a_nema_sample_data_richness_temp_1)[1] <- "Summary"
names(sum_a_nema_sample_data_richness_temp_2)[1] <- "Summary"
names(sum_a_nema_sample_data_richness_temp_3)[1] <- "Summary"
names(sum_a_nema_sample_data_richness_temp_4)[1] <- "Summary"
names(sum_a_nema_sample_data_richness_temp_5)[1] <- "Summary"
names(sum_a_nema_sample_data_richness_temp_6)[1] <- "Summary"
names(sum_a_nema_sample_data_shannon_temp_1)[1] <- "Summary"
names(sum_a_nema_sample_data_shannon_temp_2)[1] <- "Summary"
names(sum_a_nema_sample_data_shannon_temp_3)[1] <- "Summary"
names(sum_a_nema_sample_data_shannon_temp_4)[1] <- "Summary"
names(sum_a_nema_sample_data_shannon_temp_5)[1] <- "Summary"
names(sum_a_nema_sample_data_shannon_temp_6)[1] <- "Summary"

names(sum_a_roti_sample_data_richness_temp_1)[1] <- "Summary"
names(sum_a_roti_sample_data_richness_temp_2)[1] <- "Summary"
names(sum_a_roti_sample_data_richness_temp_3)[1] <- "Summary"
names(sum_a_roti_sample_data_richness_temp_4)[1] <- "Summary"
names(sum_a_roti_sample_data_richness_temp_5)[1] <- "Summary"
names(sum_a_roti_sample_data_richness_temp_6)[1] <- "Summary"
names(sum_a_roti_sample_data_shannon_temp_1)[1] <- "Summary"
names(sum_a_roti_sample_data_shannon_temp_2)[1] <- "Summary"
names(sum_a_roti_sample_data_shannon_temp_3)[1] <- "Summary"
names(sum_a_roti_sample_data_shannon_temp_4)[1] <- "Summary"
names(sum_a_roti_sample_data_shannon_temp_5)[1] <- "Summary"
names(sum_a_roti_sample_data_shannon_temp_6)[1] <- "Summary"

names(sum_a_tardi_sample_data_richness_temp_1)[1] <- "Summary"
names(sum_a_tardi_sample_data_richness_temp_2)[1] <- "Summary"
names(sum_a_tardi_sample_data_richness_temp_3)[1] <- "Summary"
names(sum_a_tardi_sample_data_richness_temp_4)[1] <- "Summary"
names(sum_a_tardi_sample_data_richness_temp_5)[1] <- "Summary"
names(sum_a_tardi_sample_data_richness_temp_6)[1] <- "Summary"
names(sum_a_tardi_sample_data_shannon_temp_1)[1] <- "Summary"
names(sum_a_tardi_sample_data_shannon_temp_2)[1] <- "Summary"
names(sum_a_tardi_sample_data_shannon_temp_3)[1] <- "Summary"
names(sum_a_tardi_sample_data_shannon_temp_4)[1] <- "Summary"
names(sum_a_tardi_sample_data_shannon_temp_5)[1] <- "Summary"
names(sum_a_tardi_sample_data_shannon_temp_6)[1] <- "Summary"

names(sum_f_sample_data_richness_temp_range_1)[1] <- "Summary"
names(sum_f_sample_data_richness_temp_range_2)[1] <- "Summary"
names(sum_f_sample_data_richness_temp_range_3)[1] <- "Summary"
names(sum_f_sample_data_richness_temp_range_4)[1] <- "Summary"
names(sum_f_sample_data_richness_temp_range_5)[1] <- "Summary"
names(sum_f_sample_data_richness_temp_range_6)[1] <- "Summary"
names(sum_f_sample_data_shannon_temp_range_1)[1] <- "Summary"
names(sum_f_sample_data_shannon_temp_range_2)[1] <- "Summary"
names(sum_f_sample_data_shannon_temp_range_3)[1] <- "Summary"
names(sum_f_sample_data_shannon_temp_range_4)[1] <- "Summary"
names(sum_f_sample_data_shannon_temp_range_5)[1] <- "Summary"
names(sum_f_sample_data_shannon_temp_range_6)[1] <- "Summary"

names(sum_p_sample_data_richness_temp_range_1)[1] <- "Summary"
names(sum_p_sample_data_richness_temp_range_2)[1] <- "Summary"
names(sum_p_sample_data_richness_temp_range_3)[1] <- "Summary"
names(sum_p_sample_data_richness_temp_range_4)[1] <- "Summary"
names(sum_p_sample_data_richness_temp_range_5)[1] <- "Summary"
names(sum_p_sample_data_richness_temp_range_6)[1] <- "Summary"
names(sum_p_sample_data_shannon_temp_range_1)[1] <- "Summary"
names(sum_p_sample_data_shannon_temp_range_2)[1] <- "Summary"
names(sum_p_sample_data_shannon_temp_range_3)[1] <- "Summary"
names(sum_p_sample_data_shannon_temp_range_4)[1] <- "Summary"
names(sum_p_sample_data_shannon_temp_range_5)[1] <- "Summary"
names(sum_p_sample_data_shannon_temp_range_6)[1] <- "Summary"

names(sum_a_roti_sample_data_richness_temp_range_1)[1] <- "Summary"
names(sum_a_roti_sample_data_richness_temp_range_2)[1] <- "Summary"
names(sum_a_roti_sample_data_richness_temp_range_3)[1] <- "Summary"
names(sum_a_roti_sample_data_richness_temp_range_4)[1] <- "Summary"
names(sum_a_roti_sample_data_richness_temp_range_5)[1] <- "Summary"
names(sum_a_roti_sample_data_richness_temp_range_6)[1] <- "Summary"
names(sum_a_roti_sample_data_shannon_temp_range_1)[1] <- "Summary"
names(sum_a_roti_sample_data_shannon_temp_range_2)[1] <- "Summary"
names(sum_a_roti_sample_data_shannon_temp_range_3)[1] <- "Summary"
names(sum_a_roti_sample_data_shannon_temp_range_4)[1] <- "Summary"
names(sum_a_roti_sample_data_shannon_temp_range_5)[1] <- "Summary"
names(sum_a_roti_sample_data_shannon_temp_range_6)[1] <- "Summary"

names(sum_a_tardi_sample_data_richness_temp_range_1)[1] <- "Summary"
names(sum_a_tardi_sample_data_richness_temp_range_2)[1] <- "Summary"
names(sum_a_tardi_sample_data_richness_temp_range_3)[1] <- "Summary"
names(sum_a_tardi_sample_data_richness_temp_range_4)[1] <- "Summary"
names(sum_a_tardi_sample_data_richness_temp_range_5)[1] <- "Summary"
names(sum_a_tardi_sample_data_richness_temp_range_6)[1] <- "Summary"
names(sum_a_tardi_sample_data_shannon_temp_range_1)[1] <- "Summary"
names(sum_a_tardi_sample_data_shannon_temp_range_2)[1] <- "Summary"
names(sum_a_tardi_sample_data_shannon_temp_range_3)[1] <- "Summary"
names(sum_a_tardi_sample_data_shannon_temp_range_4)[1] <- "Summary"
names(sum_a_tardi_sample_data_shannon_temp_range_5)[1] <- "Summary"
names(sum_a_tardi_sample_data_shannon_temp_range_6)[1] <- "Summary"

names(sum_a_nema_sample_data_richness_temp_range_1)[1] <- "Summary"
names(sum_a_nema_sample_data_richness_temp_range_2)[1] <- "Summary"
names(sum_a_nema_sample_data_richness_temp_range_3)[1] <- "Summary"
names(sum_a_nema_sample_data_richness_temp_range_4)[1] <- "Summary"
names(sum_a_nema_sample_data_richness_temp_range_5)[1] <- "Summary"
names(sum_a_nema_sample_data_richness_temp_range_6)[1] <- "Summary"
names(sum_a_nema_sample_data_shannon_temp_range_1)[1] <- "Summary"
names(sum_a_nema_sample_data_shannon_temp_range_2)[1] <- "Summary"
names(sum_a_nema_sample_data_shannon_temp_range_3)[1] <- "Summary"
names(sum_a_nema_sample_data_shannon_temp_range_4)[1] <- "Summary"
names(sum_a_nema_sample_data_shannon_temp_range_5)[1] <- "Summary"
names(sum_a_nema_sample_data_shannon_temp_range_6)[1] <- "Summary"



names(sum_a_nema_richness_season1)[1] <- "Summary"
names(sum_a_nema_richness_season2)[1] <- "Summary"
names(sum_a_nema_richness_season3)[1] <- "Summary"
names(sum_a_nema_shannon_season1)[1] <- "Summary"
names(sum_a_nema_shannon_season2)[1] <- "Summary"
names(sum_a_nema_shannon_season3)[1] <- "Summary"

names(sum_a_roti_richness_season1)[1] <- "Summary"
names(sum_a_roti_richness_season2)[1] <- "Summary"
names(sum_a_roti_richness_season3)[1] <- "Summary"
names(sum_a_roti_shannon_season1)[1] <- "Summary"
names(sum_a_roti_shannon_season2)[1] <- "Summary"
names(sum_a_roti_shannon_season3)[1] <- "Summary"

names(sum_a_tardi_richness_season1)[1] <- "Summary"
names(sum_a_tardi_richness_season2)[1] <- "Summary"
names(sum_a_tardi_richness_season3)[1] <- "Summary"
names(sum_a_tardi_shannon_season1)[1] <- "Summary"
names(sum_a_tardi_shannon_season2)[1] <- "Summary"
names(sum_a_tardi_shannon_season3)[1] <- "Summary"

names(sum_f_richness_season1)[1] <- "Summary"
names(sum_f_richness_season2)[1] <- "Summary"
names(sum_f_richness_season3)[1] <- "Summary"
names(sum_f_shannon_season1)[1] <- "Summary"
names(sum_f_shannon_season2)[1] <- "Summary"
names(sum_f_shannon_season3)[1] <- "Summary"

names(sum_p_richness_season1)[1] <- "Summary"
names(sum_p_richness_season2)[1] <- "Summary"
names(sum_p_richness_season3)[1] <- "Summary"
names(sum_p_shannon_season1)[1] <- "Summary"
names(sum_p_shannon_season2)[1] <- "Summary"
names(sum_p_shannon_season3)[1] <- "Summary"

names(sum_a_nema_richness_erosion1)[1] <- "Summary"
names(sum_a_nema_richness_erosion2)[1] <- "Summary"
names(sum_a_nema_richness_erosion3)[1] <- "Summary"
names(sum_a_nema_richness_erosion4)[1] <- "Summary"
names(sum_a_nema_shannon_erosion1)[1] <- "Summary"
names(sum_a_nema_shannon_erosion2)[1] <- "Summary"
names(sum_a_nema_shannon_erosion3)[1] <- "Summary"
names(sum_a_nema_shannon_erosion4)[1] <- "Summary"

names(sum_a_roti_richness_erosion1)[1] <- "Summary"
names(sum_a_roti_richness_erosion2)[1] <- "Summary"
names(sum_a_roti_richness_erosion3)[1] <- "Summary"
names(sum_a_roti_richness_erosion4)[1] <- "Summary"
names(sum_a_roti_shannon_erosion1)[1] <- "Summary"
names(sum_a_roti_shannon_erosion2)[1] <- "Summary"
names(sum_a_roti_shannon_erosion3)[1] <- "Summary"
names(sum_a_roti_shannon_erosion4)[1] <- "Summary"

names(sum_a_tardi_richness_erosion1)[1] <- "Summary"
names(sum_a_tardi_richness_erosion2)[1] <- "Summary"
names(sum_a_tardi_richness_erosion3)[1] <- "Summary"
names(sum_a_tardi_richness_erosion4)[1] <- "Summary"
names(sum_a_tardi_shannon_erosion1)[1] <- "Summary"
names(sum_a_tardi_shannon_erosion2)[1] <- "Summary"
names(sum_a_tardi_shannon_erosion3)[1] <- "Summary"
names(sum_a_tardi_shannon_erosion4)[1] <- "Summary"

names(sum_f_richness_erosion1)[1] <- "Summary"
names(sum_f_richness_erosion2)[1] <- "Summary"
names(sum_f_richness_erosion3)[1] <- "Summary"
names(sum_f_richness_erosion4)[1] <- "Summary"
names(sum_f_shannon_erosion1)[1] <- "Summary"
names(sum_f_shannon_erosion2)[1] <- "Summary"
names(sum_f_shannon_erosion3)[1] <- "Summary"
names(sum_f_shannon_erosion4)[1] <- "Summary"

names(sum_p_richness_erosion1)[1] <- "Summary"
names(sum_p_richness_erosion2)[1] <- "Summary"
names(sum_p_richness_erosion3)[1] <- "Summary"
names(sum_p_richness_erosion4)[1] <- "Summary"
names(sum_p_shannon_erosion1)[1] <- "Summary"
names(sum_p_shannon_erosion2)[1] <- "Summary"
names(sum_p_shannon_erosion3)[1] <- "Summary"
names(sum_p_shannon_erosion4)[1] <- "Summary"

names(sum_a_nema_richness_depth1)[1] <- "Summary"
names(sum_a_nema_richness_depth2)[1] <- "Summary"
names(sum_a_nema_richness_depth3)[1] <- "Summary"
names(sum_a_nema_shannon_depth1)[1] <- "Summary"
names(sum_a_nema_shannon_depth2)[1] <- "Summary"
names(sum_a_nema_shannon_depth3)[1] <- "Summary"

names(sum_a_roti_richness_depth1)[1] <- "Summary"
names(sum_a_roti_richness_depth2)[1] <- "Summary"
names(sum_a_roti_richness_depth3)[1] <- "Summary"
names(sum_a_roti_shannon_depth1)[1] <- "Summary"
names(sum_a_roti_shannon_depth2)[1] <- "Summary"
names(sum_a_roti_shannon_depth3)[1] <- "Summary"

names(sum_a_tardi_richness_depth1)[1] <- "Summary"
names(sum_a_tardi_richness_depth2)[1] <- "Summary"
names(sum_a_tardi_richness_depth3)[1] <- "Summary"
names(sum_a_tardi_shannon_depth1)[1] <- "Summary"
names(sum_a_tardi_shannon_depth2)[1] <- "Summary"
names(sum_a_tardi_shannon_depth3)[1] <- "Summary"

names(sum_f_richness_depth1)[1] <- "Summary"
names(sum_f_richness_depth2)[1] <- "Summary"
names(sum_f_richness_depth3)[1] <- "Summary"
names(sum_f_shannon_depth1)[1] <- "Summary"
names(sum_f_shannon_depth2)[1] <- "Summary"
names(sum_f_shannon_depth3)[1] <- "Summary"

names(sum_p_richness_depth1)[1] <- "Summary"
names(sum_p_richness_depth2)[1] <- "Summary"
names(sum_p_richness_depth3)[1] <- "Summary"
names(sum_p_shannon_depth1)[1] <- "Summary"
names(sum_p_shannon_depth2)[1] <- "Summary"
names(sum_p_shannon_depth3)[1] <- "Summary"

names(sum_a_nema_richness_pH1)[1] <- "Summary"
names(sum_a_nema_richness_pH2)[1] <- "Summary"
names(sum_a_nema_richness_pH3)[1] <- "Summary"
names(sum_a_nema_shannon_pH1)[1] <- "Summary"
names(sum_a_nema_shannon_pH2)[1] <- "Summary"
names(sum_a_nema_shannon_pH3)[1] <- "Summary"

names(sum_a_roti_richness_pH1)[1] <- "Summary"
names(sum_a_roti_richness_pH2)[1] <- "Summary"
names(sum_a_roti_richness_pH3)[1] <- "Summary"
names(sum_a_roti_shannon_pH1)[1] <- "Summary"
names(sum_a_roti_shannon_pH2)[1] <- "Summary"
names(sum_a_roti_shannon_pH3)[1] <- "Summary"

names(sum_a_tardi_richness_pH1)[1] <- "Summary"
names(sum_a_tardi_richness_pH2)[1] <- "Summary"
names(sum_a_tardi_richness_pH3)[1] <- "Summary"
names(sum_a_tardi_shannon_pH1)[1] <- "Summary"
names(sum_a_tardi_shannon_pH2)[1] <- "Summary"
names(sum_a_tardi_shannon_pH3)[1] <- "Summary"

names(sum_f_richness_pH1)[1] <- "Summary"
names(sum_f_richness_pH2)[1] <- "Summary"
names(sum_f_richness_pH3)[1] <- "Summary"
names(sum_f_shannon_pH1)[1] <- "Summary"
names(sum_f_shannon_pH2)[1] <- "Summary"
names(sum_f_shannon_pH3)[1] <- "Summary"

names(sum_p_richness_pH1)[1] <- "Summary"
names(sum_p_richness_pH2)[1] <- "Summary"
names(sum_p_richness_pH3)[1] <- "Summary"
names(sum_p_shannon_pH1)[1] <- "Summary"
names(sum_p_shannon_pH2)[1] <- "Summary"
names(sum_p_shannon_pH3)[1] <- "Summary"

names(sum_a_roti_richness_CL_annual)[1] <- "Summary"
names(sum_a_roti_richness_CL_permanent)[1] <- "Summary"
names(sum_a_roti_richness_GL_unmanaged)[1] <- "Summary"
names(sum_a_roti_richness_GL_managed)[1] <- "Summary"
names(sum_a_roti_richness_WL_broadleaved)[1] <- "Summary"
names(sum_a_roti_richness_WL_coniferous)[1] <- "Summary"
names(sum_a_roti_shannon_CL_annual)[1] <- "Summary"
names(sum_a_roti_shannon_CL_permanent)[1] <- "Summary"
names(sum_a_roti_shannon_GL_unmanaged)[1] <- "Summary"
names(sum_a_roti_shannon_GL_managed)[1] <- "Summary"
names(sum_a_roti_shannon_WL_broadleaved)[1] <- "Summary"
names(sum_a_roti_shannon_WL_coniferous)[1] <- "Summary"

names(sum_a_nema_richness_CL_annual)[1] <- "Summary"
names(sum_a_nema_richness_CL_permanent)[1] <- "Summary"
names(sum_a_nema_richness_GL_unmanaged)[1] <- "Summary"
names(sum_a_nema_richness_GL_managed)[1] <- "Summary"
names(sum_a_nema_richness_WL_broadleaved)[1] <- "Summary"
names(sum_a_nema_richness_WL_coniferous)[1] <- "Summary"
names(sum_a_nema_shannon_CL_annual)[1] <- "Summary"
names(sum_a_nema_shannon_CL_permanent)[1] <- "Summary"
names(sum_a_nema_shannon_GL_unmanaged)[1] <- "Summary"
names(sum_a_nema_shannon_GL_managed)[1] <- "Summary"
names(sum_a_nema_shannon_WL_broadleaved)[1] <- "Summary"
names(sum_a_nema_shannon_WL_coniferous)[1] <- "Summary"

names(sum_a_tardi_richness_CL_annual)[1] <- "Summary"
names(sum_a_tardi_richness_CL_permanent)[1] <- "Summary"
names(sum_a_tardi_richness_GL_managed)[1] <- "Summary"
names(sum_a_tardi_richness_GL_unmanaged)[1] <- "Summary"
names(sum_a_tardi_richness_WL_broadleaved)[1] <- "Summary"
names(sum_a_tardi_richness_WL_coniferous)[1] <- "Summary"
names(sum_a_tardi_shannon_CL_annual)[1] <- "Summary"
names(sum_a_tardi_shannon_CL_permanent)[1] <- "Summary"
names(sum_a_tardi_shannon_GL_unmanaged)[1] <- "Summary"
names(sum_a_tardi_shannon_GL_managed)[1] <- "Summary"
names(sum_a_tardi_shannon_WL_broadleaved)[1] <- "Summary"
names(sum_a_tardi_shannon_WL_coniferous)[1] <- "Summary"

names(sum_f_sample_data_richness_CL_annual)[1] <- "Summary"
names(sum_f_sample_data_richness_CL_permanent)[1] <- "Summary"
names(sum_f_sample_data_richness_GL_unmanaged)[1] <- "Summary"
names(sum_f_sample_data_richness_GL_managed)[1] <- "Summary"
names(sum_f_sample_data_richness_WL_broadleaved)[1] <- "Summary"
names(sum_f_sample_data_richness_WL_coniferous)[1] <- "Summary"
names(sum_f_sample_data_shannon_CL_annual)[1] <- "Summary"
names(sum_f_sample_data_shannon_CL_permanent)[1] <- "Summary"
names(sum_f_sample_data_shannon_GL_unmanaged)[1] <- "Summary"
names(sum_f_sample_data_shannon_GL_managed)[1] <- "Summary"
names(sum_f_sample_data_shannon_WL_broadleaved)[1] <- "Summary"
names(sum_f_sample_data_shannon_WL_coniferous)[1] <- "Summary"

names(sum_p_sample_data_richness_CL_annual)[1] <- "Summary"
names(sum_p_sample_data_richness_CL_permanent)[1] <- "Summary"
names(sum_p_sample_data_richness_GL_unmanaged)[1] <- "Summary"
names(sum_p_sample_data_richness_GL_managed)[1] <- "Summary"
names(sum_p_sample_data_richness_WL_coniferous)[1] <- "Summary"
names(sum_p_sample_data_richness_WL_broadleaved)[1] <- "Summary"
names(sum_p_sample_data_shannon_CL_annual)[1] <- "Summary"
names(sum_p_sample_data_shannon_CL_permanent)[1] <- "Summary"
names(sum_p_sample_data_shannon_GL_managed)[1] <- "Summary"
names(sum_p_sample_data_shannon_GL_unmanaged)[1] <- "Summary"
names(sum_p_sample_data_shannon_WL_coniferous)[1] <- "Summary"
names(sum_p_sample_data_shannon_WL_broadleaved)[1] <- "Summary"

names(sum_roti_sample_data_richness_CL1)[1] <- "Summary"
names(sum_roti_sample_data_richness_CL2)[1] <- "Summary"
names(sum_roti_sample_data_richness_CL3)[1] <- "Summary"
names(sum_roti_sample_data_richness_CL4)[1] <- "Summary"
names(sum_roti_sample_data_richness_GL1)[1] <- "Summary"
names(sum_roti_sample_data_richness_GL2)[1] <- "Summary"
names(sum_roti_sample_data_richness_GL3)[1] <- "Summary"
names(sum_roti_sample_data_richness_WL1)[1] <- "Summary"
names(sum_roti_sample_data_richness_WL2)[1] <- "Summary"
names(sum_roti_sample_data_richness_WL3)[1] <- "Summary"
names(sum_roti_sample_data_richness_WL4)[1] <- "Summary"
names(sum_roti_sample_data_richness_No_C)[1] <- "Summary"
names(sum_roti_sample_data_shannon_CL1)[1] <- "Summary"
names(sum_roti_sample_data_shannon_CL2)[1] <- "Summary"
names(sum_roti_sample_data_shannon_CL3)[1] <- "Summary"
names(sum_roti_sample_data_shannon_CL4)[1] <- "Summary"
names(sum_roti_sample_data_shannon_GL1)[1] <- "Summary"
names(sum_roti_sample_data_shannon_GL2)[1] <- "Summary"
names(sum_roti_sample_data_shannon_GL3)[1] <- "Summary"
names(sum_roti_sample_data_shannon_WL1)[1] <- "Summary"
names(sum_roti_sample_data_shannon_WL2)[1] <- "Summary"
names(sum_roti_sample_data_shannon_WL3)[1] <- "Summary"
names(sum_roti_sample_data_shannon_WL4)[1] <- "Summary"
names(sum_roti_sample_data_shannon_No_C)[1] <- "Summary"

names(sum_nema_sample_data_richness_CL1)[1] <- "Summary"
names(sum_nema_sample_data_richness_CL2)[1] <- "Summary"
names(sum_nema_sample_data_richness_CL3)[1] <- "Summary"
names(sum_nema_sample_data_richness_CL4)[1] <- "Summary"
names(sum_nema_sample_data_richness_GL1)[1] <- "Summary"
names(sum_nema_sample_data_richness_GL2)[1] <- "Summary"
names(sum_nema_sample_data_richness_GL3)[1] <- "Summary"
names(sum_nema_sample_data_richness_WL1)[1] <- "Summary"
names(sum_nema_sample_data_richness_WL2)[1] <- "Summary"
names(sum_nema_sample_data_richness_WL3)[1] <- "Summary"
names(sum_nema_sample_data_richness_WL4)[1] <- "Summary"
names(sum_nema_sample_data_richness_No_C)[1] <- "Summary"
names(sum_nema_sample_data_shannon_CL1)[1] <- "Summary"
names(sum_nema_sample_data_shannon_CL2)[1] <- "Summary"
names(sum_nema_sample_data_shannon_CL3)[1] <- "Summary"
names(sum_nema_sample_data_shannon_CL4)[1] <- "Summary"
names(sum_nema_sample_data_shannon_GL1)[1] <- "Summary"
names(sum_nema_sample_data_shannon_GL2)[1] <- "Summary"
names(sum_nema_sample_data_shannon_GL3)[1] <- "Summary"
names(sum_nema_sample_data_shannon_WL1)[1] <- "Summary"
names(sum_nema_sample_data_shannon_WL2)[1] <- "Summary"
names(sum_nema_sample_data_shannon_WL3)[1] <- "Summary"
names(sum_nema_sample_data_shannon_WL4)[1] <- "Summary"
names(sum_nema_sample_data_shannon_No_C)[1] <- "Summary"

names(sum_tardi_sample_data_richness_CL1)[1] <- "Summary"
names(sum_tardi_sample_data_richness_CL2)[1] <- "Summary"
names(sum_tardi_sample_data_richness_CL3)[1] <- "Summary"
names(sum_tardi_sample_data_richness_CL4)[1] <- "Summary"
names(sum_tardi_sample_data_richness_GL1)[1] <- "Summary"
names(sum_tardi_sample_data_richness_GL2)[1] <- "Summary"
names(sum_tardi_sample_data_richness_GL3)[1] <- "Summary"
names(sum_tardi_sample_data_richness_WL1)[1] <- "Summary"
names(sum_tardi_sample_data_richness_WL2)[1] <- "Summary"
names(sum_tardi_sample_data_richness_WL3)[1] <- "Summary"
names(sum_tardi_sample_data_richness_WL4)[1] <- "Summary"
names(sum_tardi_sample_data_richness_No_C)[1] <- "Summary"
names(sum_tardi_sample_data_shannon_CL1)[1] <- "Summary"
names(sum_tardi_sample_data_shannon_CL2)[1] <- "Summary"
names(sum_tardi_sample_data_shannon_CL3)[1] <- "Summary"
names(sum_tardi_sample_data_shannon_CL4)[1] <- "Summary"
names(sum_tardi_sample_data_shannon_GL1)[1] <- "Summary"
names(sum_tardi_sample_data_shannon_GL2)[1] <- "Summary"
names(sum_tardi_sample_data_shannon_GL3)[1] <- "Summary"
names(sum_tardi_sample_data_shannon_WL1)[1] <- "Summary"
names(sum_tardi_sample_data_shannon_WL2)[1] <- "Summary"
names(sum_tardi_sample_data_shannon_WL3)[1] <- "Summary"
names(sum_tardi_sample_data_shannon_WL4)[1] <- "Summary"
names(sum_tardi_sample_data_shannon_No_C)[1] <- "Summary"

names(sum_f_sample_data_richness_CL1)[1] <- "Summary"
names(sum_f_sample_data_richness_CL2)[1] <- "Summary"
names(sum_f_sample_data_richness_CL3)[1] <- "Summary"
names(sum_f_sample_data_richness_CL4)[1] <- "Summary"
names(sum_f_sample_data_richness_GL1)[1] <- "Summary"
names(sum_f_sample_data_richness_GL2)[1] <- "Summary"
names(sum_f_sample_data_richness_GL3)[1] <- "Summary"
names(sum_f_sample_data_richness_WL1)[1] <- "Summary"
names(sum_f_sample_data_richness_WL2)[1] <- "Summary"
names(sum_f_sample_data_richness_WL3)[1] <- "Summary"
names(sum_f_sample_data_richness_WL4)[1] <- "Summary"
names(sum_f_sample_data_richness_No_C)[1] <- "Summary"
names(sum_f_sample_data_shannon_CL1)[1] <- "Summary"
names(sum_f_sample_data_shannon_CL2)[1] <- "Summary"
names(sum_f_sample_data_shannon_CL3)[1] <- "Summary"
names(sum_f_sample_data_shannon_CL4)[1] <- "Summary"
names(sum_f_sample_data_shannon_GL1)[1] <- "Summary"
names(sum_f_sample_data_shannon_GL2)[1] <- "Summary"
names(sum_f_sample_data_shannon_GL3)[1] <- "Summary"
names(sum_f_sample_data_shannon_WL1)[1] <- "Summary"
names(sum_f_sample_data_shannon_WL2)[1] <- "Summary"
names(sum_f_sample_data_shannon_WL3)[1] <- "Summary"
names(sum_f_sample_data_shannon_WL4)[1] <- "Summary"
names(sum_f_sample_data_shannon_No_C)[1] <- "Summary"

names(sum_p_sample_data_richness_CL1)[1] <- "Summary"
names(sum_p_sample_data_richness_CL2)[1] <- "Summary"
names(sum_p_sample_data_richness_CL3)[1] <- "Summary"
names(sum_p_sample_data_richness_CL4)[1] <- "Summary"
names(sum_p_sample_data_richness_GL1)[1] <- "Summary"
names(sum_p_sample_data_richness_GL2)[1] <- "Summary"
names(sum_p_sample_data_richness_GL3)[1] <- "Summary"
names(sum_p_sample_data_richness_WL1)[1] <- "Summary"
names(sum_p_sample_data_richness_WL2)[1] <- "Summary"
names(sum_p_sample_data_richness_WL3)[1] <- "Summary"
names(sum_p_sample_data_richness_WL4)[1] <- "Summary"
names(sum_p_sample_data_richness_No_C)[1] <- "Summary"
names(sum_p_sample_data_shannon_CL1)[1] <- "Summary"
names(sum_p_sample_data_shannon_CL2)[1] <- "Summary"
names(sum_p_sample_data_shannon_CL3)[1] <- "Summary"
names(sum_p_sample_data_shannon_CL4)[1] <- "Summary"
names(sum_p_sample_data_shannon_GL1)[1] <- "Summary"
names(sum_p_sample_data_shannon_GL2)[1] <- "Summary"
names(sum_p_sample_data_shannon_GL3)[1] <- "Summary"
names(sum_p_sample_data_shannon_WL1)[1] <- "Summary"
names(sum_p_sample_data_shannon_WL2)[1] <- "Summary"
names(sum_p_sample_data_shannon_WL3)[1] <- "Summary"
names(sum_p_sample_data_shannon_WL4)[1] <- "Summary"
names(sum_p_sample_data_shannon_No_C)[1] <- "Summary"


names(sum_a_arthro_sample_data_richness_prec_1)[1] <- "Summary"
names(sum_a_arthro_sample_data_richness_prec_2)[1] <- "Summary"
names(sum_a_arthro_sample_data_richness_prec_3)[1] <- "Summary"
names(sum_a_arthro_sample_data_richness_prec_4)[1] <- "Summary"
names(sum_a_arthro_sample_data_richness_prec_5)[1] <- "Summary"
names(sum_a_arthro_sample_data_richness_prec_6)[1] <- "Summary"
names(sum_a_arthro_sample_data_shannon_prec_1)[1] <- "Summary"
names(sum_a_arthro_sample_data_shannon_prec_2)[1] <- "Summary"
names(sum_a_arthro_sample_data_shannon_prec_3)[1] <- "Summary"
names(sum_a_arthro_sample_data_shannon_prec_4)[1] <- "Summary"
names(sum_a_arthro_sample_data_shannon_prec_5)[1] <- "Summary"
names(sum_a_arthro_sample_data_shannon_prec_6)[1] <- "Summary"

names(sum_a_arthro_sample_data_richness_temp_1)[1] <- "Summary"
names(sum_a_arthro_sample_data_richness_temp_2)[1] <- "Summary"
names(sum_a_arthro_sample_data_richness_temp_3)[1] <- "Summary"
names(sum_a_arthro_sample_data_richness_temp_4)[1] <- "Summary"
names(sum_a_arthro_sample_data_richness_temp_5)[1] <- "Summary"
names(sum_a_arthro_sample_data_richness_temp_6)[1] <- "Summary"
names(sum_a_arthro_sample_data_shannon_temp_1)[1] <- "Summary"
names(sum_a_arthro_sample_data_shannon_temp_2)[1] <- "Summary"
names(sum_a_arthro_sample_data_shannon_temp_3)[1] <- "Summary"
names(sum_a_arthro_sample_data_shannon_temp_4)[1] <- "Summary"
names(sum_a_arthro_sample_data_shannon_temp_5)[1] <- "Summary"
names(sum_a_arthro_sample_data_shannon_temp_6)[1] <- "Summary"

names(sum_a_arthro_sample_data_richness_temp_range_1)[1] <- "Summary"
names(sum_a_arthro_sample_data_richness_temp_range_2)[1] <- "Summary"
names(sum_a_arthro_sample_data_richness_temp_range_3)[1] <- "Summary"
names(sum_a_arthro_sample_data_richness_temp_range_4)[1] <- "Summary"
names(sum_a_arthro_sample_data_richness_temp_range_5)[1] <- "Summary"
names(sum_a_arthro_sample_data_richness_temp_range_6)[1] <- "Summary"
names(sum_a_arthro_sample_data_shannon_temp_range_1)[1] <- "Summary"
names(sum_a_arthro_sample_data_shannon_temp_range_2)[1] <- "Summary"
names(sum_a_arthro_sample_data_shannon_temp_range_3)[1] <- "Summary"
names(sum_a_arthro_sample_data_shannon_temp_range_4)[1] <- "Summary"
names(sum_a_arthro_sample_data_shannon_temp_range_5)[1] <- "Summary"
names(sum_a_arthro_sample_data_shannon_temp_range_6)[1] <- "Summary"


names(sum_a_arthro_richness_season1)[1] <- "Summary"
names(sum_a_arthro_richness_season2)[1] <- "Summary"
names(sum_a_arthro_richness_season3)[1] <- "Summary"
names(sum_a_arthro_shannon_season1)[1] <- "Summary"
names(sum_a_arthro_shannon_season2)[1] <- "Summary"
names(sum_a_arthro_shannon_season3)[1] <- "Summary"


names(sum_a_arthro_richness_erosion1)[1] <- "Summary"
names(sum_a_arthro_richness_erosion2)[1] <- "Summary"
names(sum_a_arthro_richness_erosion3)[1] <- "Summary"
names(sum_a_arthro_richness_erosion4)[1] <- "Summary"
names(sum_a_arthro_shannon_erosion1)[1] <- "Summary"
names(sum_a_arthro_shannon_erosion2)[1] <- "Summary"
names(sum_a_arthro_shannon_erosion3)[1] <- "Summary"
names(sum_a_arthro_shannon_erosion4)[1] <- "Summary"


names(sum_a_arthro_richness_depth1)[1] <- "Summary"
names(sum_a_arthro_richness_depth2)[1] <- "Summary"
names(sum_a_arthro_richness_depth3)[1] <- "Summary"
names(sum_a_arthro_shannon_depth1)[1] <- "Summary"
names(sum_a_arthro_shannon_depth2)[1] <- "Summary"
names(sum_a_arthro_shannon_depth3)[1] <- "Summary"

names(sum_a_arthro_richness_pH1)[1] <- "Summary"
names(sum_a_arthro_richness_pH2)[1] <- "Summary"
names(sum_a_arthro_richness_pH3)[1] <- "Summary"
names(sum_a_arthro_shannon_pH1)[1] <- "Summary"
names(sum_a_arthro_shannon_pH2)[1] <- "Summary"
names(sum_a_arthro_shannon_pH3)[1] <- "Summary"


names(sum_a_arthro_richness_CL_annual)[1] <- "Summary"
names(sum_a_arthro_richness_CL_permanent)[1] <- "Summary"
names(sum_a_arthro_richness_GL_unmanaged)[1] <- "Summary"
names(sum_a_arthro_richness_GL_managed)[1] <- "Summary"
names(sum_a_arthro_richness_WL_broadleaved)[1] <- "Summary"
names(sum_a_arthro_richness_WL_coniferous)[1] <- "Summary"
names(sum_a_arthro_shannon_CL_annual)[1] <- "Summary"
names(sum_a_arthro_shannon_CL_permanent)[1] <- "Summary"
names(sum_a_arthro_shannon_GL_unmanaged)[1] <- "Summary"
names(sum_a_arthro_shannon_GL_managed)[1] <- "Summary"
names(sum_a_arthro_shannon_WL_broadleaved)[1] <- "Summary"
names(sum_a_arthro_shannon_WL_coniferous)[1] <- "Summary"

names(sum_arthro_sample_data_richness_CL1)[1] <- "Summary"
names(sum_arthro_sample_data_richness_CL2)[1] <- "Summary"
names(sum_arthro_sample_data_richness_CL3)[1] <- "Summary"
names(sum_arthro_sample_data_richness_CL4)[1] <- "Summary"
names(sum_arthro_sample_data_richness_GL1)[1] <- "Summary"
names(sum_arthro_sample_data_richness_GL2)[1] <- "Summary"
names(sum_arthro_sample_data_richness_GL3)[1] <- "Summary"
names(sum_arthro_sample_data_richness_WL1)[1] <- "Summary"
names(sum_arthro_sample_data_richness_WL2)[1] <- "Summary"
names(sum_arthro_sample_data_richness_WL3)[1] <- "Summary"
names(sum_arthro_sample_data_richness_WL4)[1] <- "Summary"
names(sum_arthro_sample_data_richness_No_C)[1] <- "Summary"
names(sum_arthro_sample_data_shannon_CL1)[1] <- "Summary"
names(sum_arthro_sample_data_shannon_CL2)[1] <- "Summary"
names(sum_arthro_sample_data_shannon_CL3)[1] <- "Summary"
names(sum_arthro_sample_data_shannon_CL4)[1] <- "Summary"
names(sum_arthro_sample_data_shannon_GL1)[1] <- "Summary"
names(sum_arthro_sample_data_shannon_GL2)[1] <- "Summary"
names(sum_arthro_sample_data_shannon_GL3)[1] <- "Summary"
names(sum_arthro_sample_data_shannon_WL1)[1] <- "Summary"
names(sum_arthro_sample_data_shannon_WL2)[1] <- "Summary"
names(sum_arthro_sample_data_shannon_WL3)[1] <- "Summary"
names(sum_arthro_sample_data_shannon_WL4)[1] <- "Summary"
names(sum_arthro_sample_data_shannon_No_C)[1] <- "Summary"


names(sum_a_anneli_sample_data_richness_prec_1)[1] <- "Summary"
names(sum_a_anneli_sample_data_richness_prec_2)[1] <- "Summary"
names(sum_a_anneli_sample_data_richness_prec_3)[1] <- "Summary"
names(sum_a_anneli_sample_data_richness_prec_4)[1] <- "Summary"
names(sum_a_anneli_sample_data_richness_prec_5)[1] <- "Summary"
names(sum_a_anneli_sample_data_richness_prec_6)[1] <- "Summary"
names(sum_a_anneli_sample_data_shannon_prec_1)[1] <- "Summary"
names(sum_a_anneli_sample_data_shannon_prec_2)[1] <- "Summary"
names(sum_a_anneli_sample_data_shannon_prec_3)[1] <- "Summary"
names(sum_a_anneli_sample_data_shannon_prec_4)[1] <- "Summary"
names(sum_a_anneli_sample_data_shannon_prec_5)[1] <- "Summary"
names(sum_a_anneli_sample_data_shannon_prec_6)[1] <- "Summary"

names(sum_a_anneli_sample_data_richness_temp_1)[1] <- "Summary"
names(sum_a_anneli_sample_data_richness_temp_2)[1] <- "Summary"
names(sum_a_anneli_sample_data_richness_temp_3)[1] <- "Summary"
names(sum_a_anneli_sample_data_richness_temp_4)[1] <- "Summary"
names(sum_a_anneli_sample_data_richness_temp_5)[1] <- "Summary"
names(sum_a_anneli_sample_data_richness_temp_6)[1] <- "Summary"
names(sum_a_anneli_sample_data_shannon_temp_1)[1] <- "Summary"
names(sum_a_anneli_sample_data_shannon_temp_2)[1] <- "Summary"
names(sum_a_anneli_sample_data_shannon_temp_3)[1] <- "Summary"
names(sum_a_anneli_sample_data_shannon_temp_4)[1] <- "Summary"
names(sum_a_anneli_sample_data_shannon_temp_5)[1] <- "Summary"
names(sum_a_anneli_sample_data_shannon_temp_6)[1] <- "Summary"

names(sum_a_anneli_sample_data_richness_temp_range_1)[1] <- "Summary"
names(sum_a_anneli_sample_data_richness_temp_range_2)[1] <- "Summary"
names(sum_a_anneli_sample_data_richness_temp_range_3)[1] <- "Summary"
names(sum_a_anneli_sample_data_richness_temp_range_4)[1] <- "Summary"
names(sum_a_anneli_sample_data_richness_temp_range_5)[1] <- "Summary"
names(sum_a_anneli_sample_data_richness_temp_range_6)[1] <- "Summary"
names(sum_a_anneli_sample_data_shannon_temp_range_1)[1] <- "Summary"
names(sum_a_anneli_sample_data_shannon_temp_range_2)[1] <- "Summary"
names(sum_a_anneli_sample_data_shannon_temp_range_3)[1] <- "Summary"
names(sum_a_anneli_sample_data_shannon_temp_range_4)[1] <- "Summary"
names(sum_a_anneli_sample_data_shannon_temp_range_5)[1] <- "Summary"
names(sum_a_anneli_sample_data_shannon_temp_range_6)[1] <- "Summary"


names(sum_a_anneli_richness_season1)[1] <- "Summary"
names(sum_a_anneli_richness_season2)[1] <- "Summary"
names(sum_a_anneli_richness_season3)[1] <- "Summary"
names(sum_a_anneli_shannon_season1)[1] <- "Summary"
names(sum_a_anneli_shannon_season2)[1] <- "Summary"
names(sum_a_anneli_shannon_season3)[1] <- "Summary"


names(sum_a_anneli_richness_erosion1)[1] <- "Summary"
names(sum_a_anneli_richness_erosion2)[1] <- "Summary"
names(sum_a_anneli_richness_erosion3)[1] <- "Summary"
names(sum_a_anneli_richness_erosion4)[1] <- "Summary"
names(sum_a_anneli_shannon_erosion1)[1] <- "Summary"
names(sum_a_anneli_shannon_erosion2)[1] <- "Summary"
names(sum_a_anneli_shannon_erosion3)[1] <- "Summary"
names(sum_a_anneli_shannon_erosion4)[1] <- "Summary"


names(sum_a_anneli_richness_depth1)[1] <- "Summary"
names(sum_a_anneli_richness_depth2)[1] <- "Summary"
names(sum_a_anneli_richness_depth3)[1] <- "Summary"
names(sum_a_anneli_shannon_depth1)[1] <- "Summary"
names(sum_a_anneli_shannon_depth2)[1] <- "Summary"
names(sum_a_anneli_shannon_depth3)[1] <- "Summary"

names(sum_a_anneli_richness_pH1)[1] <- "Summary"
names(sum_a_anneli_richness_pH2)[1] <- "Summary"
names(sum_a_anneli_richness_pH3)[1] <- "Summary"
names(sum_a_anneli_shannon_pH1)[1] <- "Summary"
names(sum_a_anneli_shannon_pH2)[1] <- "Summary"
names(sum_a_anneli_shannon_pH3)[1] <- "Summary"


names(sum_a_anneli_richness_CL_annual)[1] <- "Summary"
names(sum_a_anneli_richness_CL_permanent)[1] <- "Summary"
names(sum_a_anneli_richness_GL_unmanaged)[1] <- "Summary"
names(sum_a_anneli_richness_GL_managed)[1] <- "Summary"
names(sum_a_anneli_richness_WL_broadleaved)[1] <- "Summary"
names(sum_a_anneli_richness_WL_coniferous)[1] <- "Summary"
names(sum_a_anneli_shannon_CL_annual)[1] <- "Summary"
names(sum_a_anneli_shannon_CL_permanent)[1] <- "Summary"
names(sum_a_anneli_shannon_GL_unmanaged)[1] <- "Summary"
names(sum_a_anneli_shannon_GL_managed)[1] <- "Summary"
names(sum_a_anneli_shannon_WL_broadleaved)[1] <- "Summary"
names(sum_a_anneli_shannon_WL_coniferous)[1] <- "Summary"

names(sum_anneli_sample_data_richness_CL1)[1] <- "Summary"
names(sum_anneli_sample_data_richness_CL2)[1] <- "Summary"
names(sum_anneli_sample_data_richness_CL3)[1] <- "Summary"
names(sum_anneli_sample_data_richness_CL4)[1] <- "Summary"
names(sum_anneli_sample_data_richness_GL1)[1] <- "Summary"
names(sum_anneli_sample_data_richness_GL2)[1] <- "Summary"
names(sum_anneli_sample_data_richness_GL3)[1] <- "Summary"
names(sum_anneli_sample_data_richness_WL1)[1] <- "Summary"
names(sum_anneli_sample_data_richness_WL2)[1] <- "Summary"
names(sum_anneli_sample_data_richness_WL3)[1] <- "Summary"
names(sum_anneli_sample_data_richness_WL4)[1] <- "Summary"
names(sum_anneli_sample_data_richness_No_C)[1] <- "Summary"
names(sum_anneli_sample_data_shannon_CL1)[1] <- "Summary"
names(sum_anneli_sample_data_shannon_CL2)[1] <- "Summary"
names(sum_anneli_sample_data_shannon_CL3)[1] <- "Summary"
names(sum_anneli_sample_data_shannon_CL4)[1] <- "Summary"
names(sum_anneli_sample_data_shannon_GL1)[1] <- "Summary"
names(sum_anneli_sample_data_shannon_GL2)[1] <- "Summary"
names(sum_anneli_sample_data_shannon_GL3)[1] <- "Summary"
names(sum_anneli_sample_data_shannon_WL1)[1] <- "Summary"
names(sum_anneli_sample_data_shannon_WL2)[1] <- "Summary"
names(sum_anneli_sample_data_shannon_WL3)[1] <- "Summary"
names(sum_anneli_sample_data_shannon_WL4)[1] <- "Summary"
names(sum_anneli_sample_data_shannon_No_C)[1] <- "Summary"







LC_summmary_table <- rbind(sum_p_sample_data_richness_CL_annual,sum_p_sample_data_richness_CL_permanent,sum_p_sample_data_richness_GL_unmanaged,sum_p_sample_data_richness_GL_managed,sum_p_sample_data_richness_WL_coniferous,sum_p_sample_data_richness_WL_broadleaved,sum_p_sample_data_shannon_CL_annual,sum_p_sample_data_shannon_CL_permanent,sum_p_sample_data_shannon_GL_unmanaged,sum_p_sample_data_richness_GL_managed,sum_p_sample_data_shannon_WL_coniferous,sum_p_sample_data_shannon_WL_broadleaved,
                           sum_f_sample_data_richness_CL_annual,sum_f_sample_data_richness_CL_permanent,sum_f_sample_data_richness_GL_unmanaged,sum_f_sample_data_richness_GL_managed,sum_f_sample_data_richness_WL_coniferous,sum_f_sample_data_richness_WL_broadleaved,sum_f_sample_data_shannon_CL_annual,sum_f_sample_data_shannon_CL_permanent,sum_f_sample_data_shannon_GL_unmanaged,sum_f_sample_data_shannon_GL_managed,sum_f_sample_data_shannon_WL_coniferous,sum_f_sample_data_shannon_WL_broadleaved,
                           sum_a_roti_richness_CL_annual,sum_a_roti_richness_CL_permanent,sum_a_roti_richness_GL_unmanaged,sum_a_roti_richness_GL_managed,sum_a_roti_richness_WL_coniferous,sum_a_roti_richness_WL_broadleaved,sum_a_roti_shannon_CL_annual,sum_a_roti_shannon_CL_permanent,sum_a_roti_shannon_GL_unmanaged,sum_a_roti_shannon_GL_managed,sum_a_roti_shannon_WL_coniferous,sum_a_roti_shannon_WL_broadleaved,
                           sum_a_tardi_richness_CL_annual,sum_a_tardi_richness_CL_permanent,sum_a_tardi_richness_GL_unmanaged,sum_a_tardi_richness_GL_managed,sum_a_tardi_richness_WL_coniferous,sum_a_tardi_richness_WL_broadleaved,sum_a_tardi_shannon_CL_annual,sum_a_tardi_shannon_CL_permanent,sum_a_tardi_shannon_GL_unmanaged,sum_a_tardi_shannon_GL_managed,sum_a_tardi_shannon_WL_coniferous,sum_a_tardi_shannon_WL_broadleaved,
                           sum_a_nema_richness_CL_annual,sum_a_nema_richness_CL_permanent,sum_a_nema_richness_GL_unmanaged,sum_a_nema_richness_GL_managed,sum_a_nema_richness_WL_coniferous,sum_a_nema_richness_WL_broadleaved,sum_a_nema_shannon_CL_annual,sum_a_nema_shannon_CL_permanent,sum_a_nema_shannon_GL_unmanaged,sum_a_nema_shannon_GL_managed,sum_a_nema_shannon_WL_coniferous,sum_a_nema_shannon_WL_broadleaved,
                           sum_a_arthro_richness_CL_annual,sum_a_arthro_richness_CL_permanent,sum_a_arthro_richness_GL_unmanaged,sum_a_arthro_richness_GL_managed,sum_a_arthro_richness_WL_coniferous,sum_a_arthro_richness_WL_broadleaved,sum_a_arthro_shannon_CL_annual,sum_a_arthro_shannon_CL_permanent,sum_a_arthro_shannon_GL_unmanaged,sum_a_arthro_shannon_GL_managed,sum_a_arthro_shannon_WL_coniferous,sum_a_arthro_shannon_WL_broadleaved,
                           sum_a_anneli_richness_CL_annual,sum_a_anneli_richness_CL_permanent,sum_a_anneli_richness_GL_unmanaged,sum_a_anneli_richness_GL_managed,sum_a_anneli_richness_WL_coniferous,sum_a_anneli_richness_WL_broadleaved,sum_a_anneli_shannon_CL_annual,sum_a_anneli_shannon_CL_permanent,sum_a_anneli_shannon_GL_unmanaged,sum_a_anneli_shannon_GL_managed,sum_a_anneli_shannon_WL_coniferous,sum_a_anneli_shannon_WL_broadleaved)

LC_summmary_table <- LC_summmary_table[,c(2,3,4,1)]

complete_summmary_table <- bind_rows(sum_p_sample_data_richness_CL_annual,sum_p_sample_data_richness_CL_permanent,sum_p_sample_data_richness_GL_managed,sum_p_sample_data_richness_GL_unmanaged,sum_p_sample_data_richness_WL_coniferous,sum_p_sample_data_richness_WL_broadleaved,sum_p_sample_data_shannon_CL_annual,sum_p_sample_data_shannon_CL_permanent,sum_p_sample_data_shannon_GL_managed,sum_p_sample_data_shannon_GL_unmanaged,sum_p_sample_data_shannon_WL_coniferous,sum_p_sample_data_shannon_WL_broadleaved,
                                     sum_p_sample_data_richness_CL1,sum_p_sample_data_richness_CL2,sum_p_sample_data_richness_CL3,sum_p_sample_data_richness_CL4,sum_p_sample_data_richness_GL1,sum_p_sample_data_richness_GL2,sum_p_sample_data_richness_GL3,sum_p_sample_data_richness_WL1,sum_p_sample_data_richness_WL2,sum_p_sample_data_richness_WL3,sum_p_sample_data_richness_WL4,sum_p_sample_data_richness_No_C,sum_p_sample_data_shannon_CL1,sum_p_sample_data_shannon_CL2,sum_p_sample_data_shannon_CL3,sum_p_sample_data_shannon_CL4,sum_p_sample_data_shannon_GL1,sum_p_sample_data_shannon_GL2,sum_p_sample_data_shannon_GL3,sum_p_sample_data_shannon_WL1,sum_p_sample_data_shannon_WL2,sum_p_sample_data_shannon_WL3,sum_p_sample_data_shannon_WL4,sum_p_sample_data_shannon_No_C,
                                     sum_p_sample_data_richness_temp_1,sum_p_sample_data_richness_temp_2,sum_p_sample_data_richness_temp_3,sum_p_sample_data_richness_temp_4,sum_p_sample_data_richness_temp_5,sum_p_sample_data_richness_temp_6,sum_p_sample_data_richness_temp_range_1,sum_p_sample_data_richness_temp_range_2,sum_p_sample_data_richness_temp_range_3,sum_p_sample_data_richness_temp_range_4,sum_p_sample_data_richness_temp_range_5,sum_p_sample_data_richness_temp_range_6,sum_p_sample_data_shannon_temp_1,sum_p_sample_data_shannon_temp_2,sum_p_sample_data_shannon_temp_3,sum_p_sample_data_shannon_temp_4,sum_p_sample_data_shannon_temp_5,sum_p_sample_data_shannon_temp_6,sum_p_sample_data_shannon_temp_range_1,sum_p_sample_data_shannon_temp_range_2,sum_p_sample_data_shannon_temp_range_3,sum_p_sample_data_shannon_temp_range_4,sum_p_sample_data_shannon_temp_range_5,sum_p_sample_data_shannon_temp_range_6,
                                     sum_a_roti_richness_CL_annual,sum_a_roti_richness_CL_permanent,sum_a_roti_richness_GL_managed,sum_a_roti_richness_GL_unmanaged,sum_a_roti_richness_WL_coniferous,sum_a_roti_richness_WL_broadleaved,sum_a_roti_shannon_CL_annual,sum_a_roti_shannon_CL_permanent,sum_a_roti_shannon_GL_managed,sum_a_roti_shannon_GL_unmanaged,sum_a_roti_shannon_WL_coniferous,sum_a_roti_shannon_WL_broadleaved,
                                     sum_a_anneli_richness_CL_annual,sum_a_anneli_richness_CL_permanent,sum_a_anneli_richness_GL_managed,sum_a_anneli_richness_GL_unmanaged,sum_a_anneli_richness_WL_coniferous,sum_a_anneli_richness_WL_broadleaved,sum_a_anneli_shannon_CL_annual,sum_a_anneli_shannon_CL_permanent,sum_a_anneli_shannon_GL_managed,sum_a_anneli_shannon_GL_unmanaged,sum_a_anneli_shannon_WL_coniferous,sum_a_anneli_shannon_WL_broadleaved,
                                     sum_a_arthro_richness_CL_annual,sum_a_arthro_richness_CL_permanent,sum_a_arthro_richness_GL_managed,sum_a_arthro_richness_GL_unmanaged,sum_a_arthro_richness_WL_coniferous,sum_a_arthro_richness_WL_broadleaved,sum_a_arthro_shannon_CL_annual,sum_a_arthro_shannon_CL_permanent,sum_a_arthro_shannon_GL_managed,sum_a_arthro_shannon_GL_unmanaged,sum_a_arthro_shannon_WL_coniferous,sum_a_arthro_shannon_WL_broadleaved,
                                     sum_f_sample_data_richness_CL_annual,sum_f_sample_data_richness_CL_permanent,sum_f_sample_data_richness_GL_managed,sum_f_sample_data_richness_GL_unmanaged,sum_f_sample_data_richness_WL_coniferous,sum_f_sample_data_richness_WL_broadleaved,sum_f_sample_data_shannon_CL_annual,sum_f_sample_data_shannon_CL_permanent,sum_f_sample_data_shannon_GL_managed,sum_f_sample_data_shannon_GL_unmanaged,sum_f_sample_data_shannon_WL_coniferous,sum_f_sample_data_shannon_WL_broadleaved,
                                     sum_f_sample_data_richness_CL1,sum_f_sample_data_richness_CL2,sum_f_sample_data_richness_CL3,sum_f_sample_data_richness_CL4,sum_f_sample_data_richness_GL1,sum_f_sample_data_richness_GL2,sum_f_sample_data_richness_GL3,sum_f_sample_data_richness_WL1,sum_f_sample_data_richness_WL2,sum_f_sample_data_richness_WL3,sum_f_sample_data_richness_WL4,sum_f_sample_data_richness_No_C,sum_f_sample_data_shannon_CL1,sum_f_sample_data_shannon_CL2,sum_f_sample_data_shannon_CL3,sum_f_sample_data_shannon_CL4,sum_f_sample_data_shannon_GL1,sum_f_sample_data_shannon_GL2,sum_f_sample_data_shannon_GL3,sum_f_sample_data_shannon_WL1,sum_f_sample_data_shannon_WL2,sum_f_sample_data_shannon_WL3,sum_f_sample_data_shannon_WL4,sum_f_sample_data_shannon_No_C,
                                     sum_f_sample_data_richness_prec_1,sum_f_sample_data_richness_prec_2,sum_f_sample_data_richness_prec_3,sum_f_sample_data_richness_prec_4,sum_f_sample_data_richness_prec_5,sum_f_sample_data_richness_prec_6,sum_f_sample_data_richness_temp_1,sum_f_sample_data_richness_temp_2,sum_f_sample_data_richness_temp_3,sum_f_sample_data_richness_temp_4,sum_f_sample_data_richness_temp_5,sum_f_sample_data_richness_temp_6,sum_f_sample_data_richness_temp_range_1,sum_f_sample_data_richness_temp_range_2,sum_f_sample_data_richness_temp_range_3,sum_f_sample_data_richness_temp_range_4,sum_f_sample_data_richness_temp_range_5,sum_f_sample_data_richness_temp_range_6,sum_f_sample_data_shannon_prec_1,sum_f_sample_data_shannon_prec_2,sum_f_sample_data_shannon_prec_3,sum_f_sample_data_shannon_prec_4,sum_f_sample_data_shannon_prec_5,sum_f_sample_data_shannon_prec_6,sum_f_sample_data_shannon_temp_1,sum_f_sample_data_shannon_temp_2,sum_f_sample_data_shannon_temp_3,sum_f_sample_data_shannon_temp_4,sum_f_sample_data_shannon_temp_5,sum_f_sample_data_shannon_temp_6,sum_f_sample_data_shannon_temp_range_1,sum_f_sample_data_shannon_temp_range_2,sum_f_sample_data_shannon_temp_range_3,sum_f_sample_data_shannon_temp_range_4,sum_f_sample_data_shannon_temp_range_5,sum_f_sample_data_shannon_temp_range_6,
                                     sum_a_roti_richness_CL_annual,sum_a_roti_richness_CL_permanent,sum_a_roti_richness_GL_managed,sum_a_roti_richness_GL_unmanaged,sum_a_roti_richness_WL_coniferous,sum_a_roti_richness_WL_broadleaved,sum_a_roti_shannon_CL_annual,sum_a_roti_shannon_CL_permanent,sum_a_roti_shannon_GL_managed,sum_a_roti_shannon_GL_unmanaged,sum_a_roti_shannon_WL_coniferous,sum_a_roti_shannon_WL_broadleaved,
                                     sum_roti_sample_data_richness_CL1,sum_roti_sample_data_richness_CL2,sum_roti_sample_data_richness_CL3,sum_roti_sample_data_richness_CL4,sum_roti_sample_data_richness_GL1,sum_roti_sample_data_richness_GL2,sum_roti_sample_data_richness_GL3,sum_roti_sample_data_richness_WL1,sum_roti_sample_data_richness_WL2,sum_roti_sample_data_richness_WL3,sum_roti_sample_data_richness_WL4,sum_roti_sample_data_richness_No_C,sum_roti_sample_data_shannon_CL1,sum_roti_sample_data_shannon_CL2,sum_roti_sample_data_shannon_CL3,sum_roti_sample_data_shannon_CL4,sum_roti_sample_data_shannon_GL1,sum_roti_sample_data_shannon_GL2,sum_roti_sample_data_shannon_GL3,sum_roti_sample_data_shannon_WL1,sum_roti_sample_data_shannon_WL2,sum_roti_sample_data_shannon_WL3,sum_roti_sample_data_shannon_WL4,sum_roti_sample_data_shannon_No_C,
                                     sum_a_roti_sample_data_richness_temp_1,sum_a_roti_sample_data_richness_temp_2,sum_a_roti_sample_data_richness_temp_3,sum_a_roti_sample_data_richness_temp_4,sum_a_roti_sample_data_richness_temp_5,sum_a_roti_sample_data_richness_temp_6,sum_a_roti_sample_data_richness_temp_range_1,sum_a_roti_sample_data_richness_temp_range_2,sum_a_roti_sample_data_richness_temp_range_3,sum_a_roti_sample_data_richness_temp_range_4,sum_a_roti_sample_data_richness_temp_range_5,sum_a_roti_sample_data_richness_temp_range_6,sum_a_roti_sample_data_shannon_temp_1,sum_a_roti_sample_data_shannon_temp_2,sum_a_roti_sample_data_shannon_temp_3,sum_a_roti_sample_data_shannon_temp_4,sum_a_roti_sample_data_shannon_temp_5,sum_a_roti_sample_data_shannon_temp_6,sum_a_roti_sample_data_shannon_temp_range_1,sum_a_roti_sample_data_shannon_temp_range_2,sum_a_roti_sample_data_shannon_temp_range_3,sum_a_roti_sample_data_shannon_temp_range_4,sum_a_roti_sample_data_shannon_temp_range_5,sum_a_roti_sample_data_shannon_temp_range_6,
                                     sum_a_roti_richness_CL_annual,sum_a_roti_richness_CL_permanent,sum_a_roti_richness_GL_managed,sum_a_roti_richness_GL_unmanaged,sum_a_roti_richness_WL_coniferous,sum_a_roti_richness_WL_broadleaved,sum_a_roti_shannon_CL_annual,sum_a_roti_shannon_CL_permanent,sum_a_roti_shannon_GL_managed,sum_a_roti_shannon_GL_unmanaged,sum_a_roti_shannon_WL_coniferous,sum_a_roti_shannon_WL_broadleaved,
                                     sum_a_arthro_richness_CL_annual,sum_a_arthro_richness_CL_permanent,sum_a_arthro_richness_GL_managed,sum_a_arthro_richness_GL_unmanaged,sum_a_arthro_richness_WL_coniferous,sum_a_arthro_richness_WL_broadleaved,sum_a_arthro_shannon_CL_annual,sum_a_arthro_shannon_CL_permanent,sum_a_arthro_shannon_GL_managed,sum_a_arthro_shannon_GL_unmanaged,sum_a_arthro_shannon_WL_coniferous,sum_a_arthro_shannon_WL_broadleaved,
                                     sum_arthro_sample_data_richness_CL1,sum_arthro_sample_data_richness_CL2,sum_arthro_sample_data_richness_CL3,sum_arthro_sample_data_richness_CL4,sum_arthro_sample_data_richness_GL1,sum_arthro_sample_data_richness_GL2,sum_arthro_sample_data_richness_GL3,sum_arthro_sample_data_richness_WL1,sum_arthro_sample_data_richness_WL2,sum_arthro_sample_data_richness_WL3,sum_arthro_sample_data_richness_WL4,sum_arthro_sample_data_richness_No_C,sum_arthro_sample_data_shannon_CL1,sum_arthro_sample_data_shannon_CL2,sum_arthro_sample_data_shannon_CL3,sum_arthro_sample_data_shannon_CL4,sum_arthro_sample_data_shannon_GL1,sum_arthro_sample_data_shannon_GL2,sum_arthro_sample_data_shannon_GL3,sum_arthro_sample_data_shannon_WL1,sum_arthro_sample_data_shannon_WL2,sum_arthro_sample_data_shannon_WL3,sum_arthro_sample_data_shannon_WL4,sum_arthro_sample_data_shannon_No_C,
                                     sum_a_arthro_sample_data_richness_temp_1,sum_a_arthro_sample_data_richness_temp_2,sum_a_arthro_sample_data_richness_temp_3,sum_a_arthro_sample_data_richness_temp_4,sum_a_arthro_sample_data_richness_temp_5,sum_a_arthro_sample_data_richness_temp_6,sum_a_arthro_sample_data_richness_temp_range_1,sum_a_arthro_sample_data_richness_temp_range_2,sum_a_arthro_sample_data_richness_temp_range_3,sum_a_arthro_sample_data_richness_temp_range_4,sum_a_arthro_sample_data_richness_temp_range_5,sum_a_arthro_sample_data_richness_temp_range_6,sum_a_arthro_sample_data_shannon_temp_1,sum_a_arthro_sample_data_shannon_temp_2,sum_a_arthro_sample_data_shannon_temp_3,sum_a_arthro_sample_data_shannon_temp_4,sum_a_arthro_sample_data_shannon_temp_5,sum_a_arthro_sample_data_shannon_temp_6,sum_a_arthro_sample_data_shannon_temp_range_1,sum_a_arthro_sample_data_shannon_temp_range_2,sum_a_arthro_sample_data_shannon_temp_range_3,sum_a_arthro_sample_data_shannon_temp_range_4,sum_a_arthro_sample_data_shannon_temp_range_5,sum_a_arthro_sample_data_shannon_temp_range_6,
                                     sum_a_arthro_richness_CL_annual,sum_a_arthro_richness_CL_permanent,sum_a_arthro_richness_GL_managed,sum_a_arthro_richness_GL_unmanaged,sum_a_arthro_richness_WL_coniferous,sum_a_arthro_richness_WL_broadleaved,sum_a_arthro_shannon_CL_annual,sum_a_arthro_shannon_CL_permanent,sum_a_arthro_shannon_GL_managed,sum_a_arthro_shannon_GL_unmanaged,sum_a_arthro_shannon_WL_coniferous,sum_a_arthro_shannon_WL_broadleaved,
                                     sum_a_anneli_richness_CL_annual,sum_a_anneli_richness_CL_permanent,sum_a_anneli_richness_GL_managed,sum_a_anneli_richness_GL_unmanaged,sum_a_anneli_richness_WL_coniferous,sum_a_anneli_richness_WL_broadleaved,sum_a_anneli_shannon_CL_annual,sum_a_anneli_shannon_CL_permanent,sum_a_anneli_shannon_GL_managed,sum_a_anneli_shannon_GL_unmanaged,sum_a_anneli_shannon_WL_coniferous,sum_a_anneli_shannon_WL_broadleaved,
                                     sum_anneli_sample_data_richness_CL1,sum_anneli_sample_data_richness_CL2,sum_anneli_sample_data_richness_CL3,sum_anneli_sample_data_richness_CL4,sum_anneli_sample_data_richness_GL1,sum_anneli_sample_data_richness_GL2,sum_anneli_sample_data_richness_GL3,sum_anneli_sample_data_richness_WL1,sum_anneli_sample_data_richness_WL2,sum_anneli_sample_data_richness_WL3,sum_anneli_sample_data_richness_WL4,sum_anneli_sample_data_richness_No_C,sum_anneli_sample_data_shannon_CL1,sum_anneli_sample_data_shannon_CL2,sum_anneli_sample_data_shannon_CL3,sum_anneli_sample_data_shannon_CL4,sum_anneli_sample_data_shannon_GL1,sum_anneli_sample_data_shannon_GL2,sum_anneli_sample_data_shannon_GL3,sum_anneli_sample_data_shannon_WL1,sum_anneli_sample_data_shannon_WL2,sum_anneli_sample_data_shannon_WL3,sum_anneli_sample_data_shannon_WL4,sum_anneli_sample_data_shannon_No_C,
                                     sum_a_anneli_sample_data_richness_temp_1,sum_a_anneli_sample_data_richness_temp_2,sum_a_anneli_sample_data_richness_temp_3,sum_a_anneli_sample_data_richness_temp_4,sum_a_anneli_sample_data_richness_temp_5,sum_a_anneli_sample_data_richness_temp_6,sum_a_anneli_sample_data_richness_temp_range_1,sum_a_anneli_sample_data_richness_temp_range_2,sum_a_anneli_sample_data_richness_temp_range_3,sum_a_anneli_sample_data_richness_temp_range_4,sum_a_anneli_sample_data_richness_temp_range_5,sum_a_anneli_sample_data_richness_temp_range_6,sum_a_anneli_sample_data_shannon_temp_1,sum_a_anneli_sample_data_shannon_temp_2,sum_a_anneli_sample_data_shannon_temp_3,sum_a_anneli_sample_data_shannon_temp_4,sum_a_anneli_sample_data_shannon_temp_5,sum_a_anneli_sample_data_shannon_temp_6,sum_a_anneli_sample_data_shannon_temp_range_1,sum_a_anneli_sample_data_shannon_temp_range_2,sum_a_anneli_sample_data_shannon_temp_range_3,sum_a_anneli_sample_data_shannon_temp_range_4,sum_a_anneli_sample_data_shannon_temp_range_5,sum_a_anneli_sample_data_shannon_temp_range_6,
                                     sum_a_anneli_richness_CL_annual,sum_a_anneli_richness_CL_permanent,sum_a_anneli_richness_GL_managed,sum_a_anneli_richness_GL_unmanaged,sum_a_anneli_richness_WL_coniferous,sum_a_anneli_richness_WL_broadleaved,sum_a_anneli_shannon_CL_annual,sum_a_anneli_shannon_CL_permanent,sum_a_anneli_shannon_GL_managed,sum_a_anneli_shannon_GL_unmanaged,sum_a_anneli_shannon_WL_coniferous,sum_a_anneli_shannon_WL_broadleaved,
                                     sum_a_tardi_richness_CL_annual,sum_a_tardi_richness_CL_permanent,sum_a_tardi_richness_GL_managed,sum_a_tardi_richness_GL_unmanaged,sum_a_tardi_richness_WL_coniferous,sum_a_tardi_richness_WL_broadleaved,sum_a_tardi_shannon_CL_annual,sum_a_tardi_shannon_CL_permanent,sum_a_tardi_shannon_GL_managed,sum_a_tardi_shannon_GL_unmanaged,sum_a_tardi_shannon_WL_coniferous,sum_a_tardi_shannon_WL_broadleaved,
                                     sum_tardi_sample_data_richness_CL1,sum_tardi_sample_data_richness_CL2,sum_tardi_sample_data_richness_CL3,sum_tardi_sample_data_richness_CL4,sum_tardi_sample_data_richness_GL1,sum_tardi_sample_data_richness_GL2,sum_tardi_sample_data_richness_GL3,sum_tardi_sample_data_richness_WL1,sum_tardi_sample_data_richness_WL2,sum_tardi_sample_data_richness_WL3,sum_tardi_sample_data_richness_WL4,sum_tardi_sample_data_richness_No_C,sum_tardi_sample_data_shannon_CL1,sum_tardi_sample_data_shannon_CL2,sum_tardi_sample_data_shannon_CL3,sum_tardi_sample_data_shannon_CL4,sum_tardi_sample_data_shannon_GL1,sum_tardi_sample_data_shannon_GL2,sum_tardi_sample_data_shannon_GL3,sum_tardi_sample_data_shannon_WL1,sum_tardi_sample_data_shannon_WL2,sum_tardi_sample_data_shannon_WL3,sum_tardi_sample_data_shannon_WL4,sum_tardi_sample_data_shannon_No_C,
                                     sum_a_tardi_sample_data_richness_temp_1,sum_a_tardi_sample_data_richness_temp_2,sum_a_tardi_sample_data_richness_temp_3,sum_a_tardi_sample_data_richness_temp_4,sum_a_tardi_sample_data_richness_temp_5,sum_a_tardi_sample_data_richness_temp_6,sum_a_tardi_sample_data_richness_temp_range_1,sum_a_tardi_sample_data_richness_temp_range_2,sum_a_tardi_sample_data_richness_temp_range_3,sum_a_tardi_sample_data_richness_temp_range_4,sum_a_tardi_sample_data_richness_temp_range_5,sum_a_tardi_sample_data_richness_temp_range_6,sum_a_tardi_sample_data_shannon_prec_1,sum_a_tardi_sample_data_shannon_prec_2,sum_a_tardi_sample_data_shannon_prec_3,sum_a_tardi_sample_data_shannon_prec_4,sum_a_tardi_sample_data_shannon_prec_5,sum_a_tardi_sample_data_shannon_prec_6,sum_a_tardi_sample_data_shannon_temp_1,sum_a_tardi_sample_data_shannon_temp_2,sum_a_tardi_sample_data_shannon_temp_3,sum_a_tardi_sample_data_shannon_temp_4,sum_a_tardi_sample_data_shannon_temp_5,sum_a_tardi_sample_data_shannon_temp_6,sum_a_tardi_sample_data_shannon_temp_range_1,sum_a_tardi_sample_data_shannon_temp_range_2,sum_a_tardi_sample_data_shannon_temp_range_3,sum_a_tardi_sample_data_shannon_temp_range_4,sum_a_tardi_sample_data_shannon_temp_range_5,sum_a_tardi_sample_data_shannon_temp_range_6,
                                     sum_a_roti_richness_CL_annual,sum_a_roti_richness_CL_permanent,sum_a_roti_richness_GL_managed,sum_a_roti_richness_GL_unmanaged,sum_a_roti_richness_WL_coniferous,sum_a_roti_richness_WL_broadleaved,sum_a_roti_shannon_CL_annual,sum_a_roti_shannon_CL_permanent,sum_a_roti_shannon_GL_managed,sum_a_roti_shannon_GL_unmanaged,sum_a_roti_shannon_WL_coniferous,sum_a_roti_shannon_WL_broadleaved,
                                     sum_a_anneli_richness_CL_annual,sum_a_anneli_richness_CL_permanent,sum_a_anneli_richness_GL_managed,sum_a_anneli_richness_GL_unmanaged,sum_a_anneli_richness_WL_coniferous,sum_a_anneli_richness_WL_broadleaved,sum_a_anneli_shannon_CL_annual,sum_a_anneli_shannon_CL_permanent,sum_a_anneli_shannon_GL_managed,sum_a_anneli_shannon_GL_unmanaged,sum_a_anneli_shannon_WL_coniferous,sum_a_anneli_shannon_WL_broadleaved,
                                     sum_a_arthro_richness_CL_annual,sum_a_arthro_richness_CL_permanent,sum_a_arthro_richness_GL_managed,sum_a_arthro_richness_GL_unmanaged,sum_a_arthro_richness_WL_coniferous,sum_a_arthro_richness_WL_broadleaved,sum_a_arthro_shannon_CL_annual,sum_a_arthro_shannon_CL_permanent,sum_a_arthro_shannon_GL_managed,sum_a_arthro_shannon_GL_unmanaged,sum_a_arthro_shannon_WL_coniferous,sum_a_arthro_shannon_WL_broadleaved,
                                     sum_a_nema_richness_CL_annual,sum_a_nema_richness_CL_permanent,sum_a_nema_richness_GL_managed,sum_a_nema_richness_GL_unmanaged,sum_a_nema_richness_WL_coniferous,sum_a_nema_richness_WL_broadleaved,sum_a_nema_shannon_CL_annual,sum_a_nema_shannon_CL_permanent,sum_a_nema_shannon_GL_managed,sum_a_nema_shannon_GL_unmanaged,sum_a_nema_shannon_WL_coniferous,sum_a_nema_shannon_WL_broadleaved,
                                     sum_nema_sample_data_richness_CL1,sum_nema_sample_data_richness_CL2,sum_nema_sample_data_richness_CL3,sum_nema_sample_data_richness_CL4,sum_nema_sample_data_richness_GL1,sum_nema_sample_data_richness_GL2,sum_nema_sample_data_richness_GL3,sum_nema_sample_data_richness_WL1,sum_nema_sample_data_richness_WL2,sum_nema_sample_data_richness_WL3,sum_nema_sample_data_richness_WL4,sum_nema_sample_data_richness_No_C,sum_nema_sample_data_shannon_CL1,sum_nema_sample_data_shannon_CL2,sum_nema_sample_data_shannon_CL3,sum_nema_sample_data_shannon_CL4,sum_nema_sample_data_shannon_GL1,sum_nema_sample_data_shannon_GL2,sum_nema_sample_data_shannon_GL3,sum_nema_sample_data_shannon_WL1,sum_nema_sample_data_shannon_WL2,sum_nema_sample_data_shannon_WL3,sum_nema_sample_data_shannon_WL4,sum_nema_sample_data_shannon_No_C,
                                     sum_a_nema_sample_data_richness_temp_1,sum_a_nema_sample_data_richness_temp_2,sum_a_nema_sample_data_richness_temp_3,sum_a_nema_sample_data_richness_temp_4,sum_a_nema_sample_data_richness_temp_5,sum_a_nema_sample_data_richness_temp_6,sum_a_nema_sample_data_richness_temp_range_1,sum_a_nema_sample_data_richness_temp_range_2,sum_a_nema_sample_data_richness_temp_range_3,sum_a_nema_sample_data_richness_temp_range_4,sum_a_nema_sample_data_richness_temp_range_5,sum_a_nema_sample_data_richness_temp_range_6,sum_a_nema_sample_data_shannon_temp_1,sum_a_nema_sample_data_shannon_temp_2,sum_a_nema_sample_data_shannon_temp_3,sum_a_nema_sample_data_shannon_temp_4,sum_a_nema_sample_data_shannon_temp_5,sum_a_nema_sample_data_shannon_temp_6,sum_a_nema_sample_data_shannon_temp_range_1,sum_a_nema_sample_data_shannon_temp_range_2,sum_a_nema_sample_data_shannon_temp_range_3,sum_a_nema_sample_data_shannon_temp_range_4,sum_a_nema_sample_data_shannon_temp_range_5,sum_a_nema_sample_data_shannon_temp_range_6,
                                     sum_a_roti_richness_CL_annual,sum_a_roti_richness_CL_permanent,sum_a_roti_richness_GL_managed,sum_a_roti_richness_GL_unmanaged,sum_a_roti_richness_WL_coniferous,sum_a_roti_richness_WL_broadleaved,sum_a_roti_shannon_CL_annual,sum_a_roti_shannon_CL_permanent,sum_a_roti_shannon_GL_managed,sum_a_roti_shannon_GL_unmanaged,sum_a_roti_shannon_WL_coniferous,sum_a_roti_shannon_WL_broadleaved,
                                     sum_a_anneli_richness_CL_annual,sum_a_anneli_richness_CL_permanent,sum_a_anneli_richness_GL_managed,sum_a_anneli_richness_GL_unmanaged,sum_a_anneli_richness_WL_coniferous,sum_a_anneli_richness_WL_broadleaved,sum_a_anneli_shannon_CL_annual,sum_a_anneli_shannon_CL_permanent,sum_a_anneli_shannon_GL_managed,sum_a_anneli_shannon_GL_unmanaged,sum_a_anneli_shannon_WL_coniferous,sum_a_anneli_shannon_WL_broadleaved,
                                     sum_a_arthro_richness_CL_annual,sum_a_arthro_richness_CL_permanent,sum_a_arthro_richness_GL_managed,sum_a_arthro_richness_GL_unmanaged,sum_a_arthro_richness_WL_coniferous,sum_a_arthro_richness_WL_broadleaved,sum_a_arthro_shannon_CL_annual,sum_a_arthro_shannon_CL_permanent,sum_a_arthro_shannon_GL_managed,sum_a_arthro_shannon_GL_unmanaged,sum_a_arthro_shannon_WL_coniferous,sum_a_arthro_shannon_WL_broadleaved,
                                     sum_f_richness_depth1,sum_f_richness_depth2,sum_f_richness_depth3,sum_f_shannon_depth1,sum_f_shannon_depth2,sum_f_shannon_depth3,
                                     sum_p_richness_depth1,sum_p_richness_depth2,sum_p_richness_depth3,sum_p_shannon_depth1,sum_p_shannon_depth2,sum_p_shannon_depth3,
                                     sum_a_roti_richness_depth1,sum_a_roti_richness_depth2,sum_a_roti_richness_depth3,sum_a_roti_shannon_depth1,sum_a_roti_shannon_depth2,sum_a_roti_shannon_depth3,
                                     sum_a_anneli_richness_depth1,sum_a_anneli_richness_depth2,sum_a_anneli_richness_depth3,sum_a_anneli_shannon_depth1,sum_a_anneli_shannon_depth2,sum_a_anneli_shannon_depth3,
                                     sum_a_arthro_richness_depth1,sum_a_arthro_richness_depth2,sum_a_arthro_richness_depth3,sum_a_arthro_shannon_depth1,sum_a_arthro_shannon_depth2,sum_a_arthro_shannon_depth3,
                                     sum_a_tardi_richness_depth1,sum_a_tardi_richness_depth2,sum_a_tardi_richness_depth3,sum_a_tardi_shannon_depth1,sum_a_tardi_shannon_depth2,sum_a_tardi_shannon_depth3,
                                     sum_a_nema_richness_depth1,sum_a_nema_richness_depth2,sum_a_nema_richness_depth3,sum_a_nema_shannon_depth1,sum_a_nema_shannon_depth2,sum_a_nema_shannon_depth3,
                                     sum_f_richness_erosion1,sum_f_richness_erosion2,sum_f_richness_erosion3,sum_f_richness_erosion4,sum_f_shannon_erosion1,sum_f_shannon_erosion2,sum_f_shannon_erosion3,sum_f_shannon_erosion4,
                                     sum_p_richness_erosion1,sum_p_richness_erosion2,sum_p_richness_erosion3,sum_p_richness_erosion4,sum_p_shannon_erosion1,sum_p_shannon_erosion2,sum_p_shannon_erosion3,sum_p_shannon_erosion4,
                                     sum_a_roti_richness_erosion1,sum_a_roti_richness_erosion2,sum_a_roti_richness_erosion3,sum_a_roti_richness_erosion4,sum_a_roti_shannon_erosion1,sum_a_roti_shannon_erosion2,sum_a_roti_shannon_erosion3,sum_a_roti_shannon_erosion4,
                                     sum_a_anneli_richness_erosion1,sum_a_anneli_richness_erosion2,sum_a_anneli_richness_erosion3,sum_a_anneli_richness_erosion4,sum_a_anneli_shannon_erosion1,sum_a_anneli_shannon_erosion2,sum_a_anneli_shannon_erosion3,sum_a_anneli_shannon_erosion4,
                                     sum_a_arthro_richness_erosion1,sum_a_arthro_richness_erosion2,sum_a_arthro_richness_erosion3,sum_a_arthro_richness_erosion4,sum_a_arthro_shannon_erosion1,sum_a_arthro_shannon_erosion2,sum_a_arthro_shannon_erosion3,sum_a_arthro_shannon_erosion4,
                                     sum_a_tardi_richness_erosion1,sum_a_tardi_richness_erosion2,sum_a_tardi_richness_erosion3,sum_a_tardi_richness_erosion4,sum_a_tardi_shannon_erosion1,sum_a_tardi_shannon_erosion2,sum_a_tardi_shannon_erosion3,sum_a_tardi_shannon_erosion4,
                                     sum_a_nema_richness_erosion1,sum_a_nema_richness_erosion2,sum_a_nema_richness_erosion3,sum_a_nema_richness_erosion4,sum_a_nema_shannon_erosion1,sum_a_nema_shannon_erosion2,sum_a_nema_shannon_erosion3,sum_a_nema_shannon_erosion4,
                                     sum_f_richness_season1,sum_f_richness_season2,sum_f_richness_season3,sum_f_shannon_season1,sum_f_shannon_season2,sum_f_shannon_season3,
                                     sum_p_richness_season1,sum_p_richness_season2,sum_p_richness_season3,sum_p_shannon_season1,sum_p_shannon_season2,sum_p_shannon_season3,
                                     sum_a_roti_richness_season1,sum_a_roti_richness_season2,sum_a_roti_richness_season3,sum_a_roti_shannon_season1,sum_a_roti_shannon_season2,sum_a_roti_shannon_season3,
                                     sum_a_anneli_richness_season1,sum_a_anneli_richness_season2,sum_a_anneli_richness_season3,sum_a_anneli_shannon_season1,sum_a_anneli_shannon_season2,sum_a_anneli_shannon_season3,
                                     sum_a_arthro_richness_season1,sum_a_arthro_richness_season2,sum_a_arthro_richness_season3,sum_a_arthro_shannon_season1,sum_a_arthro_shannon_season2,sum_a_arthro_shannon_season3,
                                     sum_a_tardi_richness_season1,sum_a_tardi_richness_season2,sum_a_tardi_richness_season3,sum_a_tardi_shannon_season1,sum_a_tardi_shannon_season2,sum_a_tardi_shannon_season3,
                                     sum_a_nema_richness_season1,sum_a_nema_richness_season2,sum_a_nema_richness_season3,sum_a_nema_shannon_season1,sum_a_nema_shannon_season2,sum_a_nema_shannon_season3,
                                     sum_f_richness_pH1,sum_f_richness_pH2,sum_f_richness_pH3,sum_f_shannon_pH1,sum_f_shannon_pH2,sum_f_shannon_pH3,sum_p_richness_pH1,sum_p_richness_pH2,sum_p_richness_pH3,sum_p_shannon_pH1,sum_p_shannon_pH2,sum_p_shannon_pH3,sum_a_roti_richness_pH1,
                                     sum_a_roti_richness_pH2,sum_a_roti_richness_pH3,sum_a_roti_shannon_pH1,sum_a_roti_shannon_pH2,sum_a_roti_shannon_pH3,
                                     sum_a_anneli_richness_pH2,sum_a_anneli_richness_pH3,sum_a_anneli_shannon_pH1,sum_a_anneli_shannon_pH2,sum_a_anneli_shannon_pH3,
                                     sum_a_arthro_richness_pH2,sum_a_arthro_richness_pH3,sum_a_arthro_shannon_pH1,sum_a_arthro_shannon_pH2,sum_a_arthro_shannon_pH3,
                                     sum_a_nema_richness_pH1,sum_a_nema_richness_pH2,sum_a_nema_richness_pH3,sum_a_nema_shannon_pH1,sum_a_nema_shannon_pH2,sum_a_nema_shannon_pH3,sum_a_tardi_richness_pH1,sum_a_tardi_richness_pH2,sum_a_tardi_richness_pH3,sum_a_tardi_shannon_pH1,sum_a_tardi_shannon_pH2,sum_a_tardi_shannon_pH3)

write.csv(complete_summmary_table, file = "complete_summmary_table1805.csv", row.names = TRUE)





