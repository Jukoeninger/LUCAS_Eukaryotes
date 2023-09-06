#############################################
#Beta diversity (Figure 5, Figure S10, S11)
#Set graphical themes
theme_jk <- function(){
  theme_bw() +
    theme(text = element_text(size=12, family = "sans"),
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 14, vjust = -3),
          plot.margin = unit(c(1, 1, 1, 1), units =, "cm"),
          plot.title = element_text(size = 14, vjust = 1, hjust = 0),
          legend.text = element_text(size = 16), 
          legend.key.size = unit(0,"cm"),
          legend.key.width= unit(0, "cm"),
          legend.title = element_text(size=16),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 0.5, linetype = "blank"))
}

theme_jk5 <- function(){
  theme_bw() +
    theme(text = element_text(size=12, family = "sans"),
          axis.text = element_text(size = 8), 
          axis.title = element_text(size = 12, vjust = -3),
          plot.margin = unit(c(1, 1, 1, 1), units =, "cm"),
          plot.title = element_text(size = 12, vjust = 1, hjust = 0),
          legend.text = element_text(size = 8), 
          legend.key.size = unit(1,"cm"),
          legend.key.width= unit(1, "cm"),
          legend.title = element_text(size=10),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 0.5, linetype = "blank"))
}

theme_jk4 <- function(){
  theme_bw() +
    theme(text = element_text(size=12, family = "sans"),
          axis.text = element_text(size = 8), 
          axis.title = element_text(size = 12, vjust = -3),
          plot.margin = unit(c(1, 1, 1, 1), units =, "cm"),
          plot.title = element_text(size = 12, vjust = 1, hjust = 0),
          legend.text = element_text(size = 8), 
          legend.key.size = unit(2,"cm"),
          legend.key.width= unit(2, "cm"),
          legend.title = element_text(size=10),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 0.5, linetype = "blank"))
}

######################
#dbRDA BETADIVERESITY MULTIVARIATE ANALYSIS
######################
#########################################################
########################################################
#Fungi
otutable_f_sqrt <- sqrt(t(asv_f_SRS)) #square root transform ASV table
otutable_f_sqrt <- as.data.frame(otutable_f_sqrt) #transform to dataframe
sampledf_f <- data.frame(f_sample_data) #load sample data
hell_f <- vegan::decostand(otutable_f_sqrt,"hellinger") #transform ASV table to Hellinger distances
BCdist_f <- vegdist(hell_f, "bray") #Create Bray-Curtis distances

#choose relevant variables based on VIF (<11) and AIC
sampledf_f1 <- data.frame(scale(sampledf_f[,c("Aridity_in_sample_month","AvPrec2000_2018","Bulk_density","Carbonates","Coarse_fragments","Clay","C.N","Cmic","Electrical_conductivity","Mean_annual_temperature","K","lon","lat","P","pH","Precipitation_seasonality","Precipitation_of_coldest_quarter","Precipitation_of_warmest_quarter","Sand","Soil_depth","SSM","Temperature_in_sample_month")]))
sampledf_f1$LC1_2018 <- sampledf_f$LC1_2018 #adding factor variables

#rename variables for better read
sampledf_f1$LC1_2018[sampledf_f1$LC1_2018 =="CL_annual"] <- "CL_1"
sampledf_f1$LC1_2018[sampledf_f1$LC1_2018 =="CL_permanent"] <- "CL_2"
sampledf_f1$LC1_2018[sampledf_f1$LC1_2018 =="GL_managed"] <- "GL_1"
sampledf_f1$LC1_2018[sampledf_f1$LC1_2018 =="GL_unmanaged"] <- "GL_2"
sampledf_f1$LC1_2018[sampledf_f1$LC1_2018 =="WL_broadleaved"] <- "WL_1"
sampledf_f1$LC1_2018[sampledf_f1$LC1_2018 =="WL_coniferous"] <- "WL_2"
db1_rda.f <-  dbrda((BCdist_f) ~ ., data = sampledf_f1, add = FALSE) #create dbRDA

vif.cca(db1_rda.f) #check for variance inflation factor, filter values <11
extractAIC(db1_rda.f)[2] #check for the Akaike information criterion

#forward selection
db0_rda.f <- dbrda((BCdist_f) ~1, data = sampledf_f1, add = FALSE) #create dbRDA for forward selection
adjR2.dbrda_f <- RsquareAdj(db0_rda.f)$adj.r.squared #define adjusted r2 value
forward_f <- ordiR2step(db0_rda.f,scope=formula(db1_rda.f),direction="forward",pstep=1000) #forward selection of variables
forward_sel_f <- forward_f$anova #create vector with selected variables

#rename selected variables
names(sampledf_f1)[names(sampledf_f1) =="AvPrec2000_2018"] <- "AvPrec2000_2018"
names(sampledf_f1)[names(sampledf_f1) =="Temperature_in_sample_month"] <- "Temp_sample_month"
names(sampledf_f1)[names(sampledf_f1) =="P"] <- "Phosphorus"
names(sampledf_f1)[names(sampledf_f1) =="Mean_annual_temperature"] <- "AvTemp1970_2000"
names(sampledf_f1)[names(sampledf_f1) =="K"] <- "Potassium"
names(sampledf_f1)[names(sampledf_f1) =="C.N"] <- "C_N_ratio"
names(sampledf_f1)[names(sampledf_f1) =="SSM"] <- "Surface_soil_moisture"
names(sampledf_f1)[names(sampledf_f1) =="lon"] <- "Longitude"
names(sampledf_f1)[names(sampledf_f1) =="lat"] <- "Latitude"
names(sampledf_f1)[names(sampledf_f1) =="Precipitation_of_warmest_quarter"] <- "Prec_w_1970_2000"
names(sampledf_f1)[names(sampledf_f1) =="Precipitation_of_coldest_quarter"] <- "Prec_c_1970_2000"
names(sampledf_f1)[names(sampledf_f1) =="Precipitation_seasonality"] <- "Prec_season_1970_2000"
names(sampledf_f1)[names(sampledf_f1) =="LC1_2018"] <- "ES_"

#create dbRDA with selected variables
dbRDA_f <- dbrda(BCdist_f ~ ES_ + Temp_sample_month + AvPrec2000_2018 + Prec_c_1970_2000 + Prec_w_1970_2000 + Prec_season_1970_2000 + AvTemp1970_2000 + Coarse_fragments + Soil_depth + Cmic + Latitude + Longitude + pH + C_N_ratio + Phosphorus + Carbonates + Sand + Surface_soil_moisture + Electrical_conductivity + Clay + Potassium, data=as.data.frame(sampledf_f1), add = FALSE)

#For plotting, rename ecosystem types legend text
sampledf_f1$ES_[sampledf_f1$ES_ =="CL_1"] <- "CL annual"
sampledf_f1$ES_[sampledf_f1$ES_ =="CL_2"] <- "CL permanent"
sampledf_f1$ES_[sampledf_f1$ES_ =="GL_1"] <- "GL managed"
sampledf_f1$ES_[sampledf_f1$ES_ =="GL_2"] <- "GL unmanaged"
sampledf_f1$ES_[sampledf_f1$ES_ =="WL_1"] <-"WL broadleaved"
sampledf_f1$ES_[sampledf_f1$ES_ =="WL_2"] <- "WL coniferous"

#plot dbRDA
plot_ES1_f <- ggord(dbRDA_f, sampledf_f1$ES_, repel=T, hull = F,max.overlaps = 25, grp_title = "Ecosystem type",cols = c("#FEFB01","orange","green","#40BAA4","#FD79AD","#FE036A"),ellipse = F, poly= F,  ext=0.99,polylntyp="solid", addsize=-2, arrow = 0.2, txt = 3.8,veccol ="#4C4E52",veclsz = 0.7, vec_ext = 0.3)+theme_jk()
plot_ES1_f

# #Variation Partitioning including dbMEM
sampledf_f1$LU <- sampledf_f$LU2 # add further factor variable
ES_f <- sampledf_f1[,c(23,24)] #group ecosystem variables
# 
sampledf_soil_f <- sampledf_f1[,c(3:9,11,14,15,19,20,21)] #group soil variables
sampledf_climate_f <- sampledf_f1[,c(1,2,10,16:18,22)] #group climate variables
# 
# longlats_f <- data.frame(long = sampledf_f1$Longitude, lat = sampledf_f1$Latitude,row.names=rownames(sampledf_f1)) #extract coordinates for dbMEM

#conduct dbMEM
f.dbmem.tmp <- dbmem(longlats_f, silent = FALSE)
f.dbmem <- as.data.frame(f.dbmem.tmp)
f.h.det <- resid(lm(as.matrix(hell_f) ~ ., data = longlats_f))
(f.dbmem.fwd <- forward.sel(f.h.det, as.matrix(f.dbmem)))
(nb.sig.dbmem <- nrow(f.dbmem.fwd)) # Number of signif. dbMEM
# Identity of the significant dbMEM in increasing order
(dbmem.sign.f <- sort(f.dbmem.fwd[ ,2]))
# Write the significant dbMEM to a new object
dbmem.red.f <- f.dbmem[ ,c(dbmem.sign.f)]
spatial_var_f <- cbind(longlats_f, dbmem.red.f) #define spatial variable group
saveRDS(spatial_var_f, "/home/koenjul/spatial_var_f.rds")
spatial_var_f <- readRDS("/home/koenjul/spatial_var_f.rds")
f.varpart2<-varpart(BCdist_f, spatial_var_f, sampledf_climate_f,ES_f, sampledf_soil_f) #create varpart
plot(f.varpart2,digits=2,main="Fungi variation partitioning",cex.names=2.5,Xnames = c("D","C","ES","Soil"),bg = c('green', "yellow",'brown',"blue")) #plot varpart


##########################################################
#########################################################
#PROTISTS
otutable_p_sqrt <- sqrt(t(asv_p_SRS)) #square root transform ASV table
otutable_p_sqrt <- as.data.frame(otutable_p_sqrt) #transform to dataframe
sampledf_p <- data.frame(p_sample_data) #load sample data
hell_p <- vegan::decostand(otutable_p_sqrt,"hellinger") #transform ASV table to Hellinger distances
BCdist_p <- vegdist(hell_p, "bray") #Create Bray-Curtis distances

#select variables based on VIF, AIC to prevent correlation
sampledf_p1 <- data.frame(scale(sampledf_p[,c("Aridity_in_sample_month","AvPrec2000_2018","Bulk_density","Carbonates","Coarse_fragments","Clay","C.N","Cmic","Electrical_conductivity","Mean_annual_temperature","K","lon","lat","P","pH","Precipitation_seasonality","Precipitation_of_coldest_quarter","Precipitation_of_warmest_quarter","Sand","Soil_depth","SSM","Temperature_in_sample_month")]))
sampledf_p1$LC1_2018 <- sampledf_p$LC1_2018 #adding factor variables

#rename ecosystem variables for better readability
sampledf_p1$LC1_2018[sampledf_p1$LC1_2018 =="CL_annual"] <- "CL_1"
sampledf_p1$LC1_2018[sampledf_p1$LC1_2018 =="CL_permanent"] <- "CL_2"
sampledf_p1$LC1_2018[sampledf_p1$LC1_2018 =="GL_managed"] <- "GL_1"
sampledf_p1$LC1_2018[sampledf_p1$LC1_2018 =="GL_unmanaged"] <- "GL_2"
sampledf_p1$LC1_2018[sampledf_p1$LC1_2018 =="WL_broadleaved"] <- "WL_1"
sampledf_p1$LC1_2018[sampledf_p1$LC1_2018 =="WL_coniferous"] <- "WL_2"

db1_rda.p <-  dbrda((BCdist_p) ~ ., data = sampledf_p1, add = FALSE) #create dbRDA 
vif.cca(db1_rda.p) #check for variance inflation factor, filter values <11
extractAIC(db1_rda.p)[2] #check for the Akaike information criterion

#forward selection
db0_rda.p <- dbrda((BCdist_p) ~1, data = sampledf_p1, add = FALSE) #create dbRDA for forward selection
adjR2.dbrda_p <- RsquareAdj(db0_rda.p)$adj.r.squared #define adjusted r2 value
forward_p <- ordiR2step(db0_rda.p,scope=formula(db1_rda.p),direction="forward",pstep=1000) #forward selection of variables
forward_sel_p <- forward_p$anova #create vector with selected variables 

#rename selected variables
names(sampledf_p1)[names(sampledf_p1) =="AvPrec2000_2018"] <- "AvPrec2000_2018"
names(sampledf_p1)[names(sampledf_p1) =="Temperature_in_sample_month"] <- "Temp_sample_month"
names(sampledf_p1)[names(sampledf_p1) =="P"] <- "Phosphorus"
names(sampledf_p1)[names(sampledf_p1) =="Mean_annual_temperature"] <- "AvTemp1970_2000"
names(sampledf_p1)[names(sampledf_p1) =="K"] <- "Potassium"
names(sampledf_p1)[names(sampledf_p1) =="C.N"] <- "C_N_ratio"
names(sampledf_p1)[names(sampledf_p1) =="SSM"] <- "Surface_soil_moisture"
names(sampledf_p1)[names(sampledf_p1) =="lon"] <- "Longitude"
names(sampledf_p1)[names(sampledf_p1) =="lat"] <- "Latitude"
names(sampledf_p1)[names(sampledf_p1) =="Precipitation_of_warmest_quarter"] <- "Prec_w_1970_2000"
names(sampledf_p1)[names(sampledf_p1) =="Precipitation_of_coldest_quarter"] <- "Prec_c_1970_2000"
names(sampledf_p1)[names(sampledf_p1) =="Precipitation_seasonality"] <- "Prec_season_1970_2000"
names(sampledf_p1)[names(sampledf_p1) =="LC1_2018"] <- "ES_"

#create dbRDA with selected variables
dbRDA_p <- dbrda(BCdist_p ~ ES_ + Soil_depth + Temp_sample_month + Longitude + Latitude + AvTemp1970_2000 + Coarse_fragments + Carbonates + Phosphorus + pH + Electrical_conductivity + C_N_ratio + Clay + Sand + Potassium + Surface_soil_moisture + Cmic + Prec_w_1970_2000 + Prec_c_1970_2000 + Prec_season_1970_2000 + AvPrec2000_2018, data=as.data.frame(sampledf_p1), add = FALSE)

#For plotting, rename ecosystem types legend text
sampledf_p1$ES_[sampledf_p1$ES_ =="CL_1"] <- "CL annual"
sampledf_p1$ES_[sampledf_p1$ES_ =="CL_2"] <- "CL permanent"
sampledf_p1$ES_[sampledf_p1$ES_ =="GL_1"] <- "GL managed"
sampledf_p1$ES_[sampledf_p1$ES_ =="GL_2"] <- "GL unmanaged"
sampledf_p1$ES_[sampledf_p1$ES_ =="WL_1"] <-"WL broadleaved" 
sampledf_p1$ES_[sampledf_p1$ES_ =="WL_2"] <- "WL coniferous"

plot_ES1_p <- ggord(dbRDA_p, sampledf_p1$ES_, repel=T, hull = F,max.overlaps = 25, grp_title = "Ecosystem type",cols = c("#FEFB01","orange","green","#40BAA4","#FD79AD","#FE036A"),ellipse = F, poly= F,  ext=0.99,polylntyp="solid", addsize=-2, arrow = 0.2, txt = 3.8,veccol ="#4C4E52",veclsz = 0.7, vec_ext = 0.3)+theme_jk() 
plot_ES1_p

#Variation Partitioning including dbMEM
sampledf_p1$LU2 <- sampledf_p$LU2 #adding further factor variable
ES_p <- sampledf_p1[,c(23,24)] #group ecosystem variables

sampledf_soil_p <- sampledf_p1[,c(3:9,11,14,15,19:21)] #group soil variables
sampledf_climate_p <- sampledf_p1[,c(1,2,10,16,17,18,22)] #group climate variables

longlats_p <- data.frame(long = sampledf_p1$Longitude, lat = sampledf_p1$Latitude,row.names=rownames(sampledf_p1)) #extract coordinates for dbMEM
#conduct dbMEM
p.dbmem.tmp <- dbmem(longlats_p, silent = FALSE)
p.dbmem <- as.data.frame(p.dbmem.tmp)
p.h.det <- resid(lm(as.matrix(hell_p) ~ ., data = longlats_p))
(p.dbmem.fwd <- forward.sel(p.h.det, as.matrix(p.dbmem)))
(nb.sig.dbmem <- nrow(p.dbmem.fwd)) # Number of signif. dbMEM
# Identity of the significant dbMEM in increasing order
(dbmem.sign.p <- sort(p.dbmem.fwd[ ,2]))
# Write the significant dbMEM to a new object
dbmem.red.p <- p.dbmem[ ,c(dbmem.sign.p)]
spatial_var_p <- cbind(longlats_p, dbmem.red.p) #define spatial variable group
load RDS spatial group
saveRDS(spatial_var_p, "/home/koenjul/spatial_var_p.rds")
spatial_var_p <- readRDS("/home/koenjul/spatial_var_p.rds")

#plot VarPart Venn diagram (Fig. S11)
p.varpart2 <-
  varpart(BCdist_p, spatial_var_p, sampledf_climate_p,ES_p, sampledf_soil_p) #create varpart
plot(p.varpart2,digits=2,Xnames = c("D","C","ES","Soil"),bg = c('green', "yellow",'brown',"blue")) #plot varpart


########################
#Rotifers
otutable_roti_sqrt <- sqrt(t(asv_roti_SRS)) #square root transform ASV table
otutable_roti_sqrt <- as.data.frame(otutable_roti_sqrt) #transform to dataframe
sampledf_roti <- data.frame(a_roti_sample_data) #load sample data
hell_roti <- vegan::decostand(otutable_roti_sqrt,"hellinger") #transform ASV table to Hellinger distances
BCdist_roti <- vegdist(hell_roti, "bray") #Create Bray-Curtis distances

#select variables based on VIF, AIC to prevent correlation
sampledf_roti1 <- data.frame(scale(sampledf_roti[,c("Aridity_in_sample_month","AvPrec2000_2018","Bulk_density","Carbonates","Coarse_fragments","Clay","C.N","Cmic","Electrical_conductivity","Mean_annual_temperature","K","lon","lat","P","pH","Precipitation_seasonality","Precipitation_of_coldest_quarter","Precipitation_of_warmest_quarter","Sand","Soil_depth","SSM","Temperature_in_sample_month")])) 
sampledf_roti1$LC1_2018 <- sampledf_roti$LC1_2018 #adding factor variable

#change names for better readability
sampledf_roti1$LC1_2018[sampledf_roti1$LC1_2018 =="CL_annual"] <- "CL_1"
sampledf_roti1$LC1_2018[sampledf_roti1$LC1_2018 =="CL_permanent"] <- "CL_2"
sampledf_roti1$LC1_2018[sampledf_roti1$LC1_2018 =="GL_managed"] <- "GL_1"
sampledf_roti1$LC1_2018[sampledf_roti1$LC1_2018 =="GL_unmanaged"] <- "GL_2"
sampledf_roti1$LC1_2018[sampledf_roti1$LC1_2018 =="WL_broadleaved"] <- "WL_1"
sampledf_roti1$LC1_2018[sampledf_roti1$LC1_2018 =="WL_coniferous"] <- "WL_2"

#sampledf_roti1$LC1 <- sampledf_roti$LC1 #adding factor variables
db1_rda.roti <-  dbrda(BCdist_roti ~ ., data = sampledf_roti1, add = FALSE) #create dbRDA with all variables
vif.cca(db1_rda.roti) #check for variance inflation factor, filter values <11
extractAIC(db1_rda.roti)[2] #check for the Akaike information criterion

#forward selection
db0_rda.roti <- dbrda((BCdist_roti) ~1, data = sampledf_roti1, add = FALSE) #create dbRDA for forward selection
adjR2.dbrda_roti <- RsquareAdj(db0_rda.roti)$adj.r.squared #define adjusted r2 value
forward_roti <- ordiR2step(db0_rda.roti,scope=formula(db1_rda.roti),direction="forward",pstep=1000) #forward selection of variables
forward_sel_roti <- forward_roti$anova #create vector with selected variables 

#rename selected variables
names(sampledf_roti1)[names(sampledf_roti1) =="Temperature_in_sample_month"] <- "Temp_sample_month"
names(sampledf_roti1)[names(sampledf_roti1) =="P"] <- "Phosphorus"
names(sampledf_roti1)[names(sampledf_roti1) =="Mean_annual_temperature"] <- "AvTemp1970_2000"
names(sampledf_roti1)[names(sampledf_roti1) =="K"] <- "Potassium"
names(sampledf_roti1)[names(sampledf_roti1) =="C.N"] <- "C_N_ratio"
names(sampledf_roti1)[names(sampledf_roti1) =="SSM"] <- "Surface_soil_moisture"
names(sampledf_roti1)[names(sampledf_roti1) =="lon"] <- "Longitude"
names(sampledf_roti1)[names(sampledf_roti1) =="lat"] <- "Latitude"
names(sampledf_roti1)[names(sampledf_roti1) =="Precipitation_of_warmest_quarter"] <- "Prec_w_1970_2000"
names(sampledf_roti1)[names(sampledf_roti1) =="Precipitation_of_coldest_quarter"] <- "Prec_c_1970_2000"
names(sampledf_roti1)[names(sampledf_roti1) =="Precipitation_seasonality"] <- "Prec_season_1970_2000"
names(sampledf_roti1)[names(sampledf_roti1) =="LC1_2018"] <- "ES_"

#create dbRDA with selected variables
dbRDA_roti <- dbrda(BCdist_roti ~ ES_ + Surface_soil_moisture + Bulk_density + Soil_depth + Prec_w_1970_2000 + Prec_season_1970_2000 + Prec_c_1970_2000 + AvTemp1970_2000 + Temp_sample_month + Longitude + Cmic + Latitude + Potassium + Electrical_conductivity + pH + C_N_ratio + Clay + Sand + Soil_depth + Carbonates + Phosphorus, data=as.data.frame(sampledf_roti1),add = FALSE)

#For plotting, rename ecosystem types legend text
sampledf_roti1$ES_[sampledf_roti1$ES_ =="CL_1"] <- "CL annual"
sampledf_roti1$ES_[sampledf_roti1$ES_ =="CL_2"] <- "CL permanent"
sampledf_roti1$ES_[sampledf_roti1$ES_ =="GL_1"] <- "GL managed"
sampledf_roti1$ES_[sampledf_roti1$ES_ =="GL_2"] <- "GL unmanaged"
sampledf_roti1$ES_[sampledf_roti1$ES_ =="WL_1"] <-"WL broadleaved" 
sampledf_roti1$ES_[sampledf_roti1$ES_ =="WL_2"] <- "WL coniferous"

#plot dnRDA
plot_ES1_roti <- ggord(dbRDA_roti, sampledf_roti1$ES_, repel=T, hull = F,max.overlaps = 35, grp_title = "Ecosystem type",cols = c("#FEFB01","orange","green","#40BAA4","#FD79AD","#FE036A"),ellipse = F, poly= F,  ext=0.99,polylntyp="solid", addsize=-2, arrow = 0.2, txt = 3.8,veccol ="#4C4E52",veclsz = 0.7, vec_ext = 0.3)+theme_jk()
plot_ES1_roti 

#Varpart including dbMEM
sampledf_roti1$LU2 <- sampledf_roti$LU2 #adding furhter factor variables
ES_roti<- sampledf_roti1[,c(23,24)] #define ecosystem type variables
sampledf_soil_roti <- sampledf_roti1[,c(3:9,11,14,15,19,20,21)] #group soil variables
sampledf_climate_roti <- sampledf_roti1[,c(1,2,10,16,17,18,22)] #group climate variables

longlats_roti <- data.frame(long = sampledf_roti1$Longitude, lat = sampledf_roti1$Latitude,row.names=rownames(sampledf_roti1)) #extract coordinates for dbMEM
#conduct dbMEM
roti.dbmem.tmp <- dbmem(longlats_roti, silent = FALSE)
roti.dbmem <- as.data.frame(roti.dbmem.tmp)
roti.h.det <- resid(lm(as.matrix(hell_roti) ~ ., data = longlats_roti))
(roti.dbmem.fwd <- forward.sel(roti.h.det, as.matrix(roti.dbmem)))
(nb.sig.dbmem <- nrow(roti.dbmem.fwd)) # Number of signif. dbMEM
# Identity of the significant dbMEM in increasing order
(dbmem.sign.roti <- sort(roti.dbmem.fwd[ ,2]))
# Write the significant dbMEM to a new object
dbmem.red.roti <- roti.dbmem[ ,c(dbmem.sign.roti)]
spatial_roti <- cbind(dbmem.red.roti,longlats_roti) #define spatial variable group
saveRDS(spatial_roti, "/home/koenjul/spatial_roti.rds")
spatial_roti <- readRDS("/home/koenjul/spatial_roti.rds")

#plot VarPart Venn diagram (Fig. S11)
roti.varpart2 <-
varpart(BCdist_roti, spatial_roti, sampledf_climate_roti,ES_roti,sampledf_soil_roti) #create varpart
plot(roti.varpart2,digits=2,Xnames = c("D","C","ES","Soil"),bg = c('green', "yellow",'brown',"blue")) #plot varpart; D referring to "geographical distance"

##########################################################
#########################################################
#tardigrades
otutable_tardi_sqrt <- sqrt(t(asv_tardi_SRS)) #square root transform ASV table
sampledf_tardi <- data.frame(a_tardi_sample_data) #load sample data
hell_tardi <- vegan::decostand(otutable_tardi_sqrt,"hellinger") #transform ASV table to Hellinger distances
BCdist_tardi <- vegdist(hell_tardi, "bray") #Create Bray-Curtis distances

#select variables based on VIF, AIC to prevent correlation
sampledf_tardi1 <- data.frame(scale(sampledf_tardi[,c("Aridity_in_sample_month","AvPrec2000_2018","Bulk_density","Carbonates","Coarse_fragments","Clay","C.N","Cmic","Electrical_conductivity","Mean_annual_temperature","K","lon","lat","P","pH","Precipitation_seasonality","Precipitation_of_coldest_quarter","Precipitation_of_warmest_quarter","Sand","Soil_depth","SSM","Temperature_in_sample_month")]))
sampledf_tardi1$LC1_2018 <- sampledf_tardi$LC1_2018 #adding factor variables

#rename ecosystem variables for better readability
sampledf_tardi1$LC1_2018[sampledf_tardi1$LC1_2018 =="CL_annual"] <- "CL_1"
sampledf_tardi1$LC1_2018[sampledf_tardi1$LC1_2018 =="CL_permanent"] <- "CL_2"
sampledf_tardi1$LC1_2018[sampledf_tardi1$LC1_2018 =="GL_managed"] <- "GL_1"
sampledf_tardi1$LC1_2018[sampledf_tardi1$LC1_2018 =="GL_unmanaged"] <- "GL_2"
sampledf_tardi1$LC1_2018[sampledf_tardi1$LC1_2018 =="WL_broadleaved"] <- "WL_1"
sampledf_tardi1$LC1_2018[sampledf_tardi1$LC1_2018 =="WL_coniferous"] <- "WL_2"

db1_rda.tardi <-  dbrda((BCdist_tardi) ~., data = sampledf_tardi1, add = FALSE) #create dbRDA
vif.cca(db1_rda.tardi)#check for variance inflation factor, filter values <11
extractAIC(db1_rda.tardi)[2] #check for the Akaike information criterion

#forward selection
db0_rda.tardi <- dbrda((BCdist_tardi) ~1, data = sampledf_tardi1, add = FALSE) #create dbRDA for forward selection
adjR2.dbrda_tardi <- RsquareAdj(db0_rda.tardi)$adj.r.squared #define adjusted r2 value
forward_tardi <- ordiR2step(db0_rda.tardi,scope=formula(db1_rda.tardi),direction="forward",pstep=1000) #forward selection of variables
forward_sel_tardi <- forward_tardi$anova #create vector with selected variables 

#rename selected variables
names(sampledf_tardi1)[names(sampledf_tardi1) =="AvPrec2000_2018"] <- "AvPrec2000_2018"
names(sampledf_tardi1)[names(sampledf_tardi1) =="Temperature_in_sample_month"] <- "Temp_sample_month"
names(sampledf_tardi1)[names(sampledf_tardi1) =="P"] <- "Phosphorus"
names(sampledf_tardi1)[names(sampledf_tardi1) =="Mean_annual_temperature"] <- "AvTemp1970_2000"
names(sampledf_tardi1)[names(sampledf_tardi1) =="K"] <- "Potassium"
names(sampledf_tardi1)[names(sampledf_tardi1) =="C.N"] <- "C_N_ratio"
names(sampledf_tardi1)[names(sampledf_tardi1) =="SSM"] <- "Surface_soil_moisture"
names(sampledf_tardi1)[names(sampledf_tardi1) =="lon"] <- "Longitude"
names(sampledf_tardi1)[names(sampledf_tardi1) =="lat"] <- "Latitude"
names(sampledf_tardi1)[names(sampledf_tardi1) =="Precipitation_of_warmest_quarter"] <- "Prec_w_1970_2000"
names(sampledf_tardi1)[names(sampledf_tardi1) =="Precipitation_of_coldest_quarter"] <- "Prec_c_1970_2000"
names(sampledf_tardi1)[names(sampledf_tardi1) =="Precipitation_seasonality"] <- "Prec_season_1970_2000"
names(sampledf_tardi1)[names(sampledf_tardi1) =="LC1_2018"] <- "ES_"

#create dbRDA with selected variables
dbRDA_tardi <- dbrda(BCdist_tardi ~ ES_ + AvPrec2000_2018 + AvTemp1970_2000 + Prec_w_1970_2000 + Prec_season_1970_2000 + Surface_soil_moisture + Longitude + Latitude + Coarse_fragments + Cmic + pH + Electrical_conductivity + C_N_ratio + Clay + Carbonates + Phosphorus + Sand + Soil_depth, data=as.data.frame(sampledf_tardi1), add = FALSE)

#For plotting, rename ecosystem types legend text
sampledf_tardi1$ES_[sampledf_tardi1$ES_ =="CL_1"] <- "CL annual"
sampledf_tardi1$ES_[sampledf_tardi1$ES_ =="CL_2"] <- "CL permanent"
sampledf_tardi1$ES_[sampledf_tardi1$ES_ =="GL_1"] <- "GL managed"
sampledf_tardi1$ES_[sampledf_tardi1$ES_ =="GL_2"] <- "GL unmanaged"
sampledf_tardi1$ES_[sampledf_tardi1$ES_ =="WL_1"] <-"WL broadleaved" 
sampledf_tardi1$ES_[sampledf_tardi1$ES_ =="WL_2"] <- "WL coniferous"
plot_ES1_tardi <- ggord(dbRDA_tardi, sampledf_tardi1$ES_, repel=T, hull = F,max.overlaps = 35, grp_title = "Ecosystem type",cols = c("#FEFB01","orange","green","#40BAA4","#FD79AD","#FE036A"),ellipse = F, poly= F,  ext=0.99,polylntyp="solid", addsize=-2, arrow = 0.2, txt = 3.8,veccol ="#4C4E52",veclsz = 0.7, vec_ext = 0.3)+theme_jk() 
plot_ES1_tardi

#VarPart including dbMEM
sampledf_tardi1$LU2 <- sampledf_tardi$LU2 #adding furthre factor variable
ES_tardi <- sampledf_tardi1[,c(23,24)] #group ecosystem variables

sampledf_soil_tardi <- sampledf_tardi1[,c(3:9,11,14,15,19,20,21)] #group soil variables
sampledf_climate_tardi <- sampledf_tardi1[,c(1,2,10,16,17,18,22)] #group climate variables

longlats_tardi <- data.frame(long = sampledf_tardi1$Longitude, lat = sampledf_tardi1$Latitude,row.names=rownames(sampledf_tardi1)) #extract coordinates for dbMEM
#conduct dbMEM
tardi.dbmem.tmp <- dbmem(longlats_tardi, silent = FALSE)
tardi.dbmem <- as.data.frame(tardi.dbmem.tmp)
tardi.h.det <- resid(lm(as.matrix(hell_tardi) ~ ., data = longlats_tardi))
(tardi.dbmem.fwd <- forward.sel(tardi.h.det, as.matrix(tardi.dbmem)))
(nb.sig.dbmem <- nrow(tardi.dbmem.fwd)) # Number of signif. dbMEM
# Identity of the significant dbMEM in increasing order
(dbmem.sign.tardi <- sort(tardi.dbmem.fwd[ ,2]))
# Write the significant dbMEM to a new object
dbmem.red.tardi <- tardi.dbmem[ ,c(dbmem.sign.tardi)]
spatial_vars_tardi <- cbind(dbmem.red.tardi, longlats_tardi) #define spatial variable group
saveRDS(spatial_vars_tardi, "/home/koenjul/spatial_vars_tardi.rds")
spatial_vars_tardi <- readRDS("/home/koenjul/spatial_vars_tardi.rds")

#plot VarPart Venn diagram (Fig. S11)
tardi.varpart2 <-
  varpart(BCdist_tardi, spatial_vars_tardi,sampledf_climate_tardi,ES_tardi,sampledf_soil_tardi) #create varpart
plot(tardi.varpart2,digits=2,Xnames = c("D","C","ES","Soil"),bg = c('green', "yellow",'brown',"blue")) #plot varpart


##########################################################
#########################################################
#nema
otutable_nema_sqrt <- sqrt(t(asv_nema_SRS)) #square root transform ASV table
otutable_nema_sqrt <- as.data.frame(otutable_nema_sqrt) #transform to dataframe
sampledf_nema <- data.frame(a_nema_sample_data) #load sample data
hell_nema <- vegan::decostand(otutable_nema_sqrt,"hellinger") #transform ASV table to Hellinger distances
BCdist_nema <- vegdist(hell_nema, "bray") #Create Bray-Curtis distances

#select variables based on VIF, AIC to prevent correlation
sampledf_nema1 <- data.frame(scale(sampledf_nema[,c("AvPrec2000_2018","Bulk_density","Carbonates","Coarse_fragments","Clay","C.N","Cmic","Electrical_conductivity","Mean_annual_temperature","K","lon","lat","P","pH","Precipitation_seasonality","Sand","Soil_depth","SSM","Temperature_in_sample_month")]))
sampledf_nema1$LC1_2018 <- sampledf_nema$LC1_2018 #adding factor variables

#rename for better readability in plot
sampledf_nema1$LC1_2018[sampledf_nema1$LC1_2018 =="CL_annual"] <- "CL_1"
sampledf_nema1$LC1_2018[sampledf_nema1$LC1_2018 =="CL_permanent"] <- "CL_2"
sampledf_nema1$LC1_2018[sampledf_nema1$LC1_2018 =="GL_managed"] <- "GL_1"
sampledf_nema1$LC1_2018[sampledf_nema1$LC1_2018 =="GL_unmanaged"] <- "GL_2"
sampledf_nema1$LC1_2018[sampledf_nema1$LC1_2018 =="WL_broadleaved"] <- "WL_1"
sampledf_nema1$LC1_2018[sampledf_nema1$LC1_2018 =="WL_coniferous"] <- "WL_2"

db1_rda.nema <-  dbrda((BCdist_nema) ~., data = sampledf_nema1, add = FALSE) #create dbRDA 
vif.cca(db1_rda.nema) #check for variance inflation factor, filter values <11
extractAIC(db1_rda.nema)[2] #check for the Akaike information criterion

#forward selection
db0_rda.nema <-  dbrda((BCdist_nema) ~1, data = sampledf_nema1, add = FALSE) #create dbRDA for forward selection
adjR2.dbrda_nema <- RsquareAdj(db0_rda.nema)$adj.r.squared #define adjusted r2 value
forward_nema <- ordiR2step(db0_rda.nema,scope=formula(db1_rda.nema),direction="forward",pstep=1000) #forward selection of variables
forward_sel_nema <- forward_nema$anova #create vector with selected variables 

#rename selected variables
names(sampledf_nema1)[names(sampledf_nema1) =="P"] <- "Phosphorus"
names(sampledf_nema1)[names(sampledf_nema1) =="Mean_annual_temperature"] <- "AvTemp1970_2000"
names(sampledf_nema1)[names(sampledf_nema1) =="C.N"] <- "C_N_ratio"
names(sampledf_nema1)[names(sampledf_nema1) =="SSM"] <- "Surface_soil_moisture"
names(sampledf_nema1)[names(sampledf_nema1) =="lon"] <- "Longitude"
names(sampledf_nema1)[names(sampledf_nema1) =="lat"] <- "Latitude"
names(sampledf_nema1)[names(sampledf_nema1) =="Precipitation_seasonality"] <- "Prec_season_1970_2000"
names(sampledf_nema1)[names(sampledf_nema1) =="LC1_2018"] <- "ES_"

#create dbRDA with selected variables
dbRDA_nema <- dbrda(BCdist_nema ~  AvTemp1970_2000 + Prec_season_1970_2000 + ES_ + Latitude + Longitude + Surface_soil_moisture + Electrical_conductivity + pH + Carbonates + C_N_ratio + Phosphorus, data=as.data.frame(sampledf_nema1), add = FALSE)

#For plotting, rename ecosystem types legend text
sampledf_nema1$ES_[sampledf_nema1$ES_ =="CL_1"] <- "CL annual"
sampledf_nema1$ES_[sampledf_nema1$ES_ =="CL_2"] <- "CL permanent"
sampledf_nema1$ES_[sampledf_nema1$ES_ =="GL_1"] <- "GL managed"
sampledf_nema1$ES_[sampledf_nema1$ES_ =="GL_2"] <- "GL unmanaged"
sampledf_nema1$ES_[sampledf_nema1$ES_ =="WL_1"] <-"WL broadleaved" 
sampledf_nema1$ES_[sampledf_nema1$ES_ =="WL_2"] <- "WL coniferous"

#plot dbRDA
plot_ES1_nema <- ggord(dbRDA_nema, sampledf_nema1$ES_, repel=T, hull = F,max.overlaps = 25, grp_title = "Ecosystem type",cols = c("#FEFB01","orange","green","#40BAA4","#FD79AD","#FE036A"),ellipse = F, poly= F,  ext=0.99,polylntyp="solid", addsize=-2, arrow = 0.2, txt = 3.8,veccol ="#4C4E52",veclsz = 0.7, vec_ext = 0.3)+theme_jk() 
plot_ES1_nema

#VarPart including dbMEM
sampledf_nema1$LU2 <- sampledf_nema$LU2 #adding further factor variables
ES_nema <- sampledf_nema1[,c(20,21)] #group ecosystem variables
sampledf_soil_nema <- sampledf_nema1[,c(2:8,10,13,14,16:18)] #group soil variables
sampledf_climate_nema <- sampledf_nema1[,c(1,9,15,19)] #group climate variables

longlats_nema <- data.frame(long = sampledf_nema1$Longitude, lat = sampledf_nema1$Latitude,row.names=rownames(sampledf_nema1)) #extract coordinates for dbMEM
#conduct dbMEM
nema.dbmem.tmp <- dbmem(longlats_nema, silent = FALSE)
nema.dbmem <- as.data.frame(nema.dbmem.tmp)
nema.h.det <- resid(lm(as.matrix(hell_nema) ~ ., data = longlats_nema))
(nema.dbmem.fwd <- forward.sel(nema.h.det, as.matrix(nema.dbmem)))
(nb.sig.dbmem <- nrow(nema.dbmem.fwd)) # Number of signif. dbMEM
# Identity of the significant dbMEM in increasing order
(dbmem.sign.nema <- sort(nema.dbmem.fwd[ ,2]))
# Write the significant dbMEM to a new object
dbmem.red.nema <- nema.dbmem[ ,c(dbmem.sign.nema)]

spatial_var_nema <- cbind(longlats_nema, dbmem.red.nema) #define spatial variable group
saveRDS(spatial_var_nema, "/home/koenjul/spatial_var_nema.rds")
spatial_var_nema <- readRDS("/home/koenjul/spatial_var_nema.rds")

#plot VarPart Venn diagram (Fig. S11)
nema.varpart2 <-
  varpart(BCdist_nema, spatial_var_nema, sampledf_climate_nema,ES_nema, sampledf_soil_nema) #create varpart
plot(nema.varpart2,digits=2,Xnames = c("D","C","ES","Soil"),bg = c('green', "yellow",'brown',"blue")) #plot varpart



##########################################################
#########################################################
#arthropods
otutable_arthro_sqrt <- sqrt(t(asv_arthro_SRS)) #square root transform ASV table
sampledf_arthro <- data.frame(a_arthro_sample_data) #load sample data
hell_arthro <- vegan::decostand(otutable_arthro_sqrt,"hellinger") #transform ASV table to Hellinger distances
BCdist_arthro <- vegdist(hell_arthro, "bray") #Create Bray-Curtis distances

#select variables based on VIF, AIC to prevent correlation
sampledf_arthro1 <- data.frame(scale(sampledf_arthro[,c("Aridity_in_sample_month","AvPrec2000_2018","Bulk_density","Carbonates","Coarse_fragments","Clay","C.N","Cmic","Electrical_conductivity","Mean_annual_temperature","K","lon","lat","P","pH","Precipitation_seasonality","Precipitation_of_coldest_quarter","Precipitation_of_warmest_quarter","Sand","Soil_depth","SSM","Temperature_in_sample_month")]))
sampledf_arthro1$LC1_2018 <- sampledf_arthro$LC1_2018 #adding factor variables

#rename ecosystem variables for better readability
sampledf_arthro1$LC1_2018[sampledf_arthro1$LC1_2018 =="CL_annual"] <- "CL_1"
sampledf_arthro1$LC1_2018[sampledf_arthro1$LC1_2018 =="CL_permanent"] <- "CL_2"
sampledf_arthro1$LC1_2018[sampledf_arthro1$LC1_2018 =="GL_managed"] <- "GL_1"
sampledf_arthro1$LC1_2018[sampledf_arthro1$LC1_2018 =="GL_unmanaged"] <- "GL_2"
sampledf_arthro1$LC1_2018[sampledf_arthro1$LC1_2018 =="WL_broadleaved"] <- "WL_1"
sampledf_arthro1$LC1_2018[sampledf_arthro1$LC1_2018 =="WL_coniferous"] <- "WL_2"

db1_rda.arthro <-  dbrda((BCdist_arthro) ~., data = sampledf_arthro1, add = FALSE) #create dbRDA
vif.cca(db1_rda.arthro)#check for variance inflation factor, filter values <11
extractAIC(db1_rda.arthro)[2] #check for the Akaike information criterion

#forward selection
db0_rda.arthro <- dbrda((BCdist_arthro) ~1, data = sampledf_arthro1, add = FALSE) #create dbRDA for forward selection
adjR2.dbrda_arthro <- RsquareAdj(db0_rda.arthro)$adj.r.squared #define adjusted r2 value
forward_arthro <- ordiR2step(db0_rda.arthro,scope=formula(db1_rda.arthro),direction="forward",pstep=1000) #forward selection of variables
forward_sel_arthro <- forward_arthro$anova #create vector with selected variables 

#rename selected variables
names(sampledf_arthro1)[names(sampledf_arthro1) =="AvPrec2000_2018"] <- "AvPrec2000_2018"
names(sampledf_arthro1)[names(sampledf_arthro1) =="Temperature_in_sample_month"] <- "Temp_sample_month"
names(sampledf_arthro1)[names(sampledf_arthro1) =="P"] <- "Phosphorus"
names(sampledf_arthro1)[names(sampledf_arthro1) =="Mean_annual_temperature"] <- "AvTemp1970_2000"
names(sampledf_arthro1)[names(sampledf_arthro1) =="K"] <- "Potassium"
names(sampledf_arthro1)[names(sampledf_arthro1) =="C.N"] <- "C_N_ratio"
names(sampledf_arthro1)[names(sampledf_arthro1) =="SSM"] <- "Surface_soil_moisture"
names(sampledf_arthro1)[names(sampledf_arthro1) =="lon"] <- "Longitude"
names(sampledf_arthro1)[names(sampledf_arthro1) =="lat"] <- "Latitude"
names(sampledf_arthro1)[names(sampledf_arthro1) =="Precipitation_of_warmest_quarter"] <- "Prec_w_1970_2000"
names(sampledf_arthro1)[names(sampledf_arthro1) =="Precipitation_of_coldest_quarter"] <- "Prec_c_1970_2000"
names(sampledf_arthro1)[names(sampledf_arthro1) =="Precipitation_seasonality"] <- "Prec_season_1970_2000"
names(sampledf_arthro1)[names(sampledf_arthro1) =="LC1_2018"] <- "ES_"

#create dbRDA with selected variables
dbRDA_arthro <- dbrda(BCdist_arthro ~ ES_ + AvPrec2000_2018 + AvTemp1970_2000 + Prec_w_1970_2000 + Prec_season_1970_2000 + Surface_soil_moisture + Longitude + Latitude + Coarse_fragments + Cmic + pH + Electrical_conductivity + C_N_ratio + Clay + Carbonates + Phosphorus + Sand + Soil_depth, data=as.data.frame(sampledf_arthro1), add = FALSE)

#For plotting, rename ecosystem types legend text
sampledf_arthro1$ES_[sampledf_arthro1$ES_ =="CL_1"] <- "CL annual"
sampledf_arthro1$ES_[sampledf_arthro1$ES_ =="CL_2"] <- "CL permanent"
sampledf_arthro1$ES_[sampledf_arthro1$ES_ =="GL_1"] <- "GL managed"
sampledf_arthro1$ES_[sampledf_arthro1$ES_ =="GL_2"] <- "GL unmanaged"
sampledf_arthro1$ES_[sampledf_arthro1$ES_ =="WL_1"] <-"WL broadleaved" 
sampledf_arthro1$ES_[sampledf_arthro1$ES_ =="WL_2"] <- "WL coniferous"
plot_ES1_arthro <- ggord(dbRDA_arthro, sampledf_arthro1$ES_, repel=T, hull = F,max.overlaps = 30, grp_title = "Ecosystem type",cols = c("#FEFB01","orange","green","#40BAA4","#FD79AD","#FE036A"),ellipse = F, poly= F,  ext=0.99,polylntyp="solid", addsize=-2, arrow = 0.2, txt = 3.8,veccol ="#4C4E52",veclsz = 0.7, vec_ext = 0.3)+theme_jk() 
plot_ES1_arthro

#VarPart including dbMEM
sampledf_arthro1$LU2 <- sampledf_arthro$LU2 #adding furthre factor variable
ES_arthro <- sampledf_arthro1[,c(23,24)] #group ecosystem variables

sampledf_soil_arthro <- sampledf_arthro1[,c(3:9,11,14,15,19,20,21)] #group soil variables
sampledf_climate_arthro <- sampledf_arthro1[,c(1,2,10,16,17,18,22)] #group climate variables

# longlats_arthro <- data.frame(long = sampledf_arthro1$Longitude, lat = sampledf_arthro1$Latitude,row.names=rownames(sampledf_arthro1)) #extract coordinates for dbMEM
# #conduct dbMEM
# arthro.dbmem.tmp <- dbmem(longlats_arthro, silent = FALSE)
# arthro.dbmem <- as.data.frame(arthro.dbmem.tmp)
# arthro.h.det <- resid(lm(as.matrix(hell_arthro) ~ ., data = longlats_arthro))
# (arthro.dbmem.fwd <- forward.sel(arthro.h.det, as.matrix(arthro.dbmem)))
# (nb.sig.dbmem <- nrow(arthro.dbmem.fwd)) # Number of signif. dbMEM
# # Identity of the significant dbMEM in increasing order
# (dbmem.sign.arthro <- sort(arthro.dbmem.fwd[ ,2]))
# # Write the significant dbMEM to a new object
# dbmem.red.arthro <- arthro.dbmem[ ,c(dbmem.sign.arthro)]
# spatial_vars_arthro <- cbind(dbmem.red.arthro, longlats_arthro) #define spatial variable group
# saveRDS(spatial_vars_arthro,"/home/koenjul/spatial_vars_arthro.rds")
spatial_vars_arthro <- readRDS("/home/koenjul/spatial_var_arthro.rds")

#plot VarPart Venn diagram (Fig. S11)
arthro.varpart2 <-
  varpart(BCdist_arthro, spatial_vars_arthro,sampledf_climate_arthro,ES_arthro,sampledf_soil_arthro) #create varpart
plot(arthro.varpart2,digits=2,Xnames = c("D","C","ES","Soil"),bg = c('green', "yellow",'brown',"blue")) #plot varpart




##########################################################
#########################################################
#annelids
otutable_anneli_sqrt <- sqrt(t(asv_anneli_SRS)) #square root transform ASV table
sampledf_anneli <- data.frame(a_anneli_sample_data) #load sample data
hell_anneli <- vegan::decostand(otutable_anneli_sqrt,"hellinger") #transform ASV table to Hellinger distances
BCdist_anneli <- vegdist(hell_anneli, "bray") #Create Bray-Curtis distances

#select variables based on VIF, AIC to prevent correlation
sampledf_anneli1 <- data.frame(scale(sampledf_anneli[,c("Aridity_in_sample_month","AvPrec2000_2018","Bulk_density","Carbonates","Coarse_fragments","Clay","C.N","Cmic","Electrical_conductivity","Mean_annual_temperature","K","lon","lat","P","pH","Precipitation_seasonality","Precipitation_of_coldest_quarter","Precipitation_of_warmest_quarter","Sand","Soil_depth","SSM","Temperature_in_sample_month")]))
sampledf_anneli1$LC1_2018 <- sampledf_anneli$LC1_2018 #adding factor variables

#rename ecosystem variables for better readability
sampledf_anneli1$LC1_2018[sampledf_anneli1$LC1_2018 =="CL_annual"] <- "CL_1"
sampledf_anneli1$LC1_2018[sampledf_anneli1$LC1_2018 =="CL_permanent"] <- "CL_2"
sampledf_anneli1$LC1_2018[sampledf_anneli1$LC1_2018 =="GL_managed"] <- "GL_1"
sampledf_anneli1$LC1_2018[sampledf_anneli1$LC1_2018 =="GL_unmanaged"] <- "GL_2"
sampledf_anneli1$LC1_2018[sampledf_anneli1$LC1_2018 =="WL_broadleaved"] <- "WL_1"
sampledf_anneli1$LC1_2018[sampledf_anneli1$LC1_2018 =="WL_coniferous"] <- "WL_2"

db1_rda.anneli <-  dbrda((BCdist_anneli) ~., data = sampledf_anneli1, add = FALSE) #create dbRDA
vif.cca(db1_rda.anneli)#check for variance inflation factor, filter values <11
extractAIC(db1_rda.anneli)[2] #check for the Akaike information criterion

#forward selection
db0_rda.anneli <- dbrda((BCdist_anneli) ~1, data = sampledf_anneli1, add = FALSE) #create dbRDA for forward selection
adjR2.dbrda_anneli <- RsquareAdj(db0_rda.anneli)$adj.r.squared #define adjusted r2 value
forward_anneli <- ordiR2step(db0_rda.anneli,scope=formula(db1_rda.anneli),direction="forward",pstep=1000) #forward selection of variables
forward_sel_anneli <- forward_anneli$anova #create vector with selected variables 

#rename selected variables
names(sampledf_anneli1)[names(sampledf_anneli1) =="AvPrec2000_2018"] <- "AvPrec2000_2018"
names(sampledf_anneli1)[names(sampledf_anneli1) =="Temperature_in_sample_month"] <- "Temp_sample_month"
names(sampledf_anneli1)[names(sampledf_anneli1) =="P"] <- "Phosphorus"
names(sampledf_anneli1)[names(sampledf_anneli1) =="Mean_annual_temperature"] <- "AvTemp1970_2000"
names(sampledf_anneli1)[names(sampledf_anneli1) =="K"] <- "Potassium"
names(sampledf_anneli1)[names(sampledf_anneli1) =="C.N"] <- "C_N_ratio"
names(sampledf_anneli1)[names(sampledf_anneli1) =="SSM"] <- "Surface_soil_moisture"
names(sampledf_anneli1)[names(sampledf_anneli1) =="lon"] <- "Longitude"
names(sampledf_anneli1)[names(sampledf_anneli1) =="lat"] <- "Latitude"
names(sampledf_anneli1)[names(sampledf_anneli1) =="Precipitation_of_warmest_quarter"] <- "Prec_w_1970_2000"
names(sampledf_anneli1)[names(sampledf_anneli1) =="Precipitation_of_coldest_quarter"] <- "Prec_c_1970_2000"
names(sampledf_anneli1)[names(sampledf_anneli1) =="Precipitation_seasonality"] <- "Prec_season_1970_2000"
names(sampledf_anneli1)[names(sampledf_anneli1) =="LC1_2018"] <- "ES_"

#create dbRDA with selected variables
dbRDA_anneli <- dbrda(BCdist_anneli ~ ES_ + AvPrec2000_2018 + AvTemp1970_2000 + Prec_w_1970_2000 + Prec_season_1970_2000 + Surface_soil_moisture + Longitude + Latitude + Coarse_fragments + Cmic + pH + Electrical_conductivity + C_N_ratio + Clay + Carbonates + Phosphorus + Sand + Soil_depth, data=as.data.frame(sampledf_anneli1), add = FALSE)

#For plotting, rename ecosystem types legend text
sampledf_anneli1$ES_[sampledf_anneli1$ES_ =="CL_1"] <- "CL annual"
sampledf_anneli1$ES_[sampledf_anneli1$ES_ =="CL_2"] <- "CL permanent"
sampledf_anneli1$ES_[sampledf_anneli1$ES_ =="GL_1"] <- "GL managed"
sampledf_anneli1$ES_[sampledf_anneli1$ES_ =="GL_2"] <- "GL unmanaged"
sampledf_anneli1$ES_[sampledf_anneli1$ES_ =="WL_1"] <-"WL broadleaved" 
sampledf_anneli1$ES_[sampledf_anneli1$ES_ =="WL_2"] <- "WL coniferous"
plot_ES1_anneli <- ggord(dbRDA_anneli, sampledf_anneli1$ES_, repel=T, hull = F,max.overlaps = 25, grp_title = "Ecosystem type",cols = c("#FEFB01","orange","green","#40BAA4","#FD79AD","#FE036A"),ellipse = F, poly= F,  ext=0.99,polylntyp="solid", addsize=-2, arrow = 0.2, txt = 3.8,veccol ="#4C4E52",veclsz = 0.7, vec_ext = 0.3)+theme_jk() 
plot_ES1_anneli

#VarPart including dbMEM
sampledf_anneli1$LU2 <- sampledf_anneli$LU2 #adding furthre factor variable
ES_anneli <- sampledf_anneli1[,c(23,24)] #group ecosystem variables

sampledf_soil_anneli <- sampledf_anneli1[,c(3:9,11,14,15,19,20,21)] #group soil variables
sampledf_climate_anneli <- sampledf_anneli1[,c(1,2,10,16,17,18,22)] #group climate variables

longlats_anneli <- data.frame(long = sampledf_anneli1$Longitude, lat = sampledf_anneli1$Latitude,row.names=rownames(sampledf_anneli1)) #extract coordinates for dbMEM
#conduct dbMEM
anneli.dbmem.tmp <- dbmem(longlats_anneli, silent = FALSE)
anneli.dbmem <- as.data.frame(anneli.dbmem.tmp)
anneli.h.det <- resid(lm(as.matrix(hell_anneli) ~ ., data = longlats_anneli))
(anneli.dbmem.fwd <- forward.sel(anneli.h.det, as.matrix(anneli.dbmem)))
(nb.sig.dbmem <- nrow(anneli.dbmem.fwd)) # Number of signif. dbMEM
# Identity of the significant dbMEM in increasing order
(dbmem.sign.anneli <- sort(anneli.dbmem.fwd[ ,2]))
# Write the significant dbMEM to a new object
dbmem.red.anneli <- anneli.dbmem[ ,c(dbmem.sign.anneli)]
spatial_vars_anneli <- cbind(dbmem.red.anneli, longlats_anneli) #define spatial variable group
saveRDS(spatial_vars_anneli,"/home/koenjul/spatial_vars_anneli.rds")
spatial_vars_anneli <- readRDS("/home/koenjul/spatial_var_anneli.rds")

#plot VarPart Venn diagram (Fig. S11)
anneli.varpart2 <-
  varpart(BCdist_anneli, spatial_vars_anneli,sampledf_climate_anneli,ES_anneli,sampledf_soil_anneli) #create varpart
plot(anneli.varpart2,digits=2,Xnames = c("D","C","ES","Soil"),bg = c('green', "yellow",'brown',"blue")) #plot varpart



#Plot Figure S11
svg("2305_varpart_anneli.svg",width=8,height=7)
#plot
plot(anneli.varpart2,digits=2,Xnames = c("D","C","ES","Soil"),bg = c('green', "yellow",'brown',"blue")) #plot varpart
dev.off()

svg("2305_varpart_arthro.svg",width=8,height=7)
#plot
plot(arthro.varpart2,digits=2,Xnames = c("D","C","ES","Soil"),bg = c('green', "yellow",'brown',"blue")) #plot varpart
dev.off()

svg("2305_varpart_nema.svg",width=8,height=7)
#plot
plot(nema.varpart2,digits=2,Xnames = c("D","C","ES","Soil"),bg = c('green', "yellow",'brown',"blue")) #plot varpart
dev.off()

svg("2305_varpart_tardi.svg",width=8,height=7)
#plot
plot(tardi.varpart2,digits=2,Xnames = c("D","C","ES","Soil"),bg = c('green', "yellow",'brown',"blue")) #plot varpart
dev.off()

svg("2305_varpart_roti.svg",width=8,height=7)
#plot
plot(roti.varpart2,digits=2,Xnames = c("D","C","ES","Soil"),bg = c('green', "yellow",'brown',"blue")) #plot varpart
dev.off()

svg("2305_varpart_p.svg",width=8,height=7)
#plot
plot(p.varpart2,digits=2,Xnames = c("D","C","ES","Soil"),bg = c('green', "yellow",'brown',"blue")) #plot varpart
dev.off()

svg("2305_varpart_f.svg",width=8,height=7)
#plot
plot(f.varpart2,digits=2,Xnames = c("D","C","ES","Soil"),bg = c('green', "yellow",'brown',"blue")) #plot varpart
dev.off()


#combined dbRDA plot for all eukaryotic groups
#Plot figure 5
svg("2305_Fig5ABCDEF.svg",width=15,height=30)
#plot
ggarrange( plot_ES1_f,plot_ES1_p,plot_ES1_roti,plot_ES1_tardi, plot_ES1_nema, plot_ES1_arthro, plot_ES1_anneli,
           labels = c("(a)", "(b)","(c)","(d)","(e)","(f)","(g)"),align="v",common.legend = TRUE,
           ncol = 2, nrow = 4)
dev.off() # Close the graphics device



#Figure 5F
#VarPart stacked barplot
first_column <- c("Geographical distance", "Climatic variables", "Ecosystem type", "Soil properties","Shared","Unexplained")
second_column <- c(1.4,0.3,1.3,1.5,6.4,89.1)
third_column <- "Fungi"
df_f <- data.frame(first_column, second_column,third_column)

first_column <- c("Geographical distance", "Climatic variables", "Ecosystem type", "Soil properties","Shared","Unexplained")
second_column <- c(1,0.2,0.7,1,3.2,93.9)
third_column <- "Protists"
df_p <- data.frame(first_column, second_column,third_column)

first_column <- c("Geographical distance", "Climatic variables", "Ecosystem type", "Soil properties","Shared","Unexplained")
second_column <- c(1.1,0.3,0.7,1.2,2.2,94.7)
third_column <- "Rotifers"
df_roti <- data.frame(first_column, second_column,third_column)

first_column <- c("Geographical distance", "Climatic variables", "Ecosystem type", "Soil properties","Shared","Unexplained")
second_column <- c(0.9,0.4,0.9,1.8,3.4,92.6)
third_column <- "Tardigrades"
df_tardi <- data.frame(first_column, second_column,third_column)

first_column <- c("Geographical distance", "Climatic variables", "Ecosystem type", "Soil properties","Shared","Unexplained")
second_column <- c(1.4,0.3,1,1.2,2.5,92.9)
third_column <- "Nematodes"
df_nema <- data.frame(first_column, second_column,third_column)

first_column <- c("Geographical distance", "Climatic variables", "Ecosystem type", "Soil properties","Shared","Unexplained")
second_column <- c(1.2,0.2,0.8,1.6,3.8,91.9)
third_column <- "Arthropods"
df_arthro <- data.frame(first_column, second_column,third_column)

first_column <- c("Geographical distance", "Climatic variables", "Ecosystem type", "Soil properties","Shared","Unexplained")
second_column <- c(2.3,0.1,4.1,3.2,20,70.3)
third_column <- "Annelids"
df_anneli <- data.frame(first_column, second_column,third_column)

df <- rbind(df_f,df_roti,df_tardi,df_nema,df_p,df_arthro,df_anneli)
colnames(df)[1] <- "Env_driver"
colnames(df)[2] <- "VarPart"
colnames(df)[3] <- "group"

df$Env_driver <- factor(df$Env_driver, levels=c("Geographical distance","Climatic variables","Ecosystem type","Soil properties","Shared","Unexplained"))
df$group <- factor(df$group, levels=c("Fungi","Protists","Rotifers", "Tardigrades", "Nematodes","Arthropods","Annelids"))
df_long <- df %>% tidyr::gather(group, value, -Env_driver, -VarPart,-group)
df_long <- df_long[,-1]
plot_df <- ggplot(df,aes(x=group,y=VarPart)) 
svg("2305_Fig5F.svg",width=14,height=8)
#plot
plot_df + geom_bar(aes(fill=Env_driver),stat="identity") +
  ylab("\nVariation explained") +
  scale_fill_manual(values=c("orange","#A7C7E7","#77dd77","#836953","#cb99c9","grey")) + 
  theme_jk5()+ theme(text = element_text(size=12, family = "sans"),
                                                  axis.title.x=element_text(size=16,angle=90, vjust = 0.5, hjust=1),
                                                   panel.grid.major = element_blank(),
                                                   panel.grid.minor = element_blank(),
                                                   panel.background = element_blank(),
                                                   panel.grid  = element_blank(),
                                                   axis.text = element_text(size = 22),
                                                   axis.title = element_text(size = 22, hjust = 0.5),
                                                   plot.margin = unit(c(1, 1, 1, 1), units =, "cm"),
                                                   legend.text = element_text(size = 20),
                                                   legend.key.size = unit(1,"cm"),
                                                   legend.key.width= unit(1, "cm"),
                                                   legend.title = element_blank(),
                                                    legend.position="bottom",
                                                   legend.background = element_rect(color = "black",
                                                                                    fill = "transparent",
                                                                                    size = 0.5, linetype = "blank"))

dev.off() # Close the graphics device



######################################################################
#Testing homogeneity of group dispersion
# Table for betadisper results
#Step 1: multivariate homogeneity of group dispersions analysis - Betadisper
betadisper_roti_LC1 <- betadisper(BCdist_roti, a_roti_sample_data$LC1_2018)
betadisper_roti_LC1_aov <- anova(betadisper_roti_LC1)
table_roti_LC1_aov <- as.data.frame(betadisper_roti_LC1_aov)
table_roti_LC1_aov$Group <- "Rotifers"
table_roti_LC1_aov$group <- "Ecosystem type"
# #Permutest
permutest(betadisper_roti_LC1, permutations = 999, pairwise = TRUE)
#Extract distances to centorid betadisper
betadisper_roti_LC1_groups <- as.data.frame.table(betadisper_roti_LC1$group.distances)
betadisper_roti_LC1_groups$Group <- "Rotifers"
betadisper_roti_LC1_groups$group <- "Ecosystem type"

betadisper_roti_LU2 <- betadisper(BCdist_roti, a_roti_sample_data$LU2)
betadisper_roti_LU2_aov <- anova(betadisper_roti_LU2)
table_roti_LU2_aov <- as.data.frame(betadisper_roti_LU2_aov)
table_roti_LU2_aov$Group <- "Rotifers"
table_roti_LU2_aov$group <- "Plant cover"
# #Permutest
permutest(betadisper_roti_LU2, permutations = 999, pairwise = TRUE)
#Extract distances to centorid betadisper
betadisper_roti_LU2_groups <- as.data.frame.table(betadisper_roti_LU2$group.distances)
betadisper_roti_LU2_groups$Group <- "Rotifers"
betadisper_roti_LU2_groups$group <- "Plant cover"

betadisper_roti_pH <- betadisper(BCdist_roti, a_roti_sample_data$pH_grouped)
betadisper_roti_pH_aov <-anova(betadisper_roti_pH)
table_roti_pH_aov <- as.data.frame(betadisper_roti_pH_aov)
table_roti_pH_aov$Group <- "Rotifers"
table_roti_pH_aov$group <- "pH"
#Permutest
permutest(betadisper_roti_pH, permutations = 999, pairwise = TRUE)
betadisper_roti_pH_groups <- as.data.frame.table(betadisper_roti_pH$group.distances)
betadisper_roti_pH_groups$Group <- "Rotifers"
betadisper_roti_pH_groups$group <- "pH"

betadisper_roti_Season <- betadisper(BCdist_roti, a_roti_sample_data$Sample_season)
betadisper_roti_Season_aov <-anova(betadisper_roti_Season)
table_roti_Season_aov <- as.data.frame(betadisper_roti_Season_aov)
table_roti_Season_aov$Group <- "Rotifers"
table_roti_Season_aov$group <- "Sampling season"
betadisper_roti_Season_groups <- as.data.frame.table(betadisper_roti_Season$group.distances)
betadisper_roti_Season_groups$Group <- "Rotifers"
betadisper_roti_Season_groups$group <- "Sampling season"
#Permutest
permutest(betadisper_roti_Season, permutations = 999, pairwise = TRUE)

betadisper_roti_depth <- betadisper(BCdist_roti, a_roti_sample_data$depth_grouped)
betadisper_roti_depth_aov <- anova(betadisper_roti_depth)
table_roti_depth_aov <- as.data.frame(betadisper_roti_depth_aov)
table_roti_depth_aov$Group <- "Rotifers"
table_roti_depth_aov$group <- "depth"
betadisper_roti_depth_groups <- as.data.frame.table(betadisper_roti_depth$group.distances)
betadisper_roti_depth_groups$Group <- "Rotifers"
betadisper_roti_depth_groups$group <- "Soil depth"
#Permutest
permutest(betadisper_roti_depth, permutations = 999, pairwise = TRUE)

betadisper_roti_erosion <- betadisper(BCdist_roti, a_roti_sample_data$Erosion_grouped)
betadisper_roti_erosion_aov <- anova(betadisper_roti_erosion)
table_roti_erosion_aov <- as.data.frame(betadisper_roti_erosion_aov)
table_roti_erosion_aov$Group <- "Rotifers"
table_roti_erosion_aov$group <- "erosion"
betadisper_roti_erosion_groups <- as.data.frame.table(betadisper_roti_erosion$group.distances)
betadisper_roti_erosion_groups$Group <- "Rotifers"
betadisper_roti_erosion_groups$group <- "Erosion risk"
#Permutest
permutest(betadisper_roti_erosion, permutations = 999, pairwise = TRUE)

betadisper_roti_temp <- betadisper(BCdist_roti, a_roti_sample_data$cluster_temp)
betadisper_roti_temp_aov <-anova(betadisper_roti_temp)
table_roti_temp_aov <- as.data.frame(betadisper_roti_temp_aov)
table_roti_temp_aov$Group <- "Rotifers"
table_roti_temp_aov$group <- "Temperature"
betadisper_roti_temp_groups <- as.data.frame.table(betadisper_roti_temp$group.distances)
betadisper_roti_temp_groups$Group <- "Rotifers"
betadisper_roti_temp_groups$group <- "Temperature"
#Permutest
permutest(betadisper_roti_temp, permutations = 999, pairwise = TRUE)

betadisper_roti_temp_range <- betadisper(BCdist_roti, a_roti_sample_data$cluster_temp_range)
betadisper_roti_temp_range_aov <-anova(betadisper_roti_temp_range)
table_roti_temp_range_aov <- as.data.frame(betadisper_roti_temp_range_aov)
table_roti_temp_range_aov$Group <- "Rotifers"
table_roti_temp_range_aov$group <- "Temperature range"
betadisper_roti_temp_range_groups <- as.data.frame.table(betadisper_roti_temp_range$group.distances)
betadisper_roti_temp_range_groups$Group <- "Rotifers"
betadisper_roti_temp_range_groups$group <- "Temperature Range"
#Permutest
permutest(betadisper_roti_temp_range, permutations = 999, pairwise = TRUE)

betadisper_roti_prec <- betadisper(BCdist_roti, a_roti_sample_data$cluster_prec)
betadisper_roti_prec_aov <- anova(betadisper_roti_prec)
table_roti_prec_aov <- as.data.frame(betadisper_roti_prec_aov)
table_roti_prec_aov$Group <- "Rotifers"
table_roti_prec_aov$group <- "Precipitation"
betadisper_roti_prec_groups <- as.data.frame.table(betadisper_roti_prec$group.distances)
betadisper_roti_prec_groups$Group <- "Rotifers"
betadisper_roti_prec_groups$group <- "Precipitation"
#Permutest
permutest(betadisper_roti_prec, permutations = 999, pairwise = TRUE)

betadisper_roti_Cmic <- betadisper(BCdist_roti, a_roti_sample_data$cluster_Cmic)
betadisper_roti_Cmic_aov <-anova(betadisper_roti_Cmic)
table_roti_Cmic_aov <- as.data.frame(betadisper_roti_Cmic_aov)
table_roti_Cmic_aov$Group <- "Rotifers"
table_roti_Cmic_aov$group <- "Cmic"
betadisper_roti_Cmic_groups <- as.data.frame.table(betadisper_roti_Cmic$group.distances)
betadisper_roti_Cmic_groups$Group <- "Rotifers"
betadisper_roti_Cmic_groups$group <- "Cmic"
#Permutest
permutest(betadisper_roti_Cmic, permutations = 999, pairwise = TRUE)

#Tardigrades
betadisper_tardi_LC1 <- betadisper(BCdist_tardi, a_tardi_sample_data$LC1_2018)
betadisper_tardi_LC1_aov <- anova(betadisper_tardi_LC1)
table_tardi_LC1_aov <- as.data.frame(betadisper_tardi_LC1_aov)
table_tardi_LC1_aov$Group <- "Tardigrades"
table_tardi_LC1_aov$group <- "Ecosystem type"
#Extract distances to centorid betadisper
betadisper_tardi_LC1_groups <- as.data.frame.table(betadisper_tardi_LC1$group.distances)
betadisper_tardi_LC1_groups$Group <- "Tardigrades"
betadisper_tardi_LC1_groups$group <- "Ecosystem type"
#Permutest 
permutest(betadisper_tardi_LC1, permutations = 999, pairwise = TRUE)

betadisper_tardi_LU2 <- betadisper(BCdist_tardi, a_tardi_sample_data$LU2)
betadisper_tardi_LU2_aov <- anova(betadisper_tardi_LU2)
table_tardi_LU2_aov <- as.data.frame(betadisper_tardi_LU2_aov)
table_tardi_LU2_aov$Group <- "Tardigrades"
table_tardi_LU2_aov$group <- "Plant cover"
# #Permutest
permutest(betadisper_tardi_LU2, permutations = 999, pairwise = TRUE)
#Extract distances to centorid betadisper
betadisper_tardi_LU2_groups <- as.data.frame.table(betadisper_tardi_LU2$group.distances)
betadisper_tardi_LU2_groups$Group <- "Tardigrades"
betadisper_tardi_LU2_groups$group <- "Plant cover"

betadisper_tardi_pH <- betadisper(BCdist_tardi, a_tardi_sample_data$pH_grouped)
betadisper_tardi_pH_aov <- anova(betadisper_tardi_pH)
table_tardi_pH_aov <- as.data.frame(betadisper_tardi_pH_aov)
table_tardi_pH_aov$Group <- "Tardigrades"
table_tardi_pH_aov$group <- "pH"
betadisper_tardi_pH_groups <- as.data.frame.table(betadisper_tardi_pH$group.distances)
betadisper_tardi_pH_groups$Group <- "Tardigrades"
betadisper_tardi_pH_groups$group <- "pH"
#Permutest
permutest(betadisper_tardi_pH, permutations = 999, pairwise = TRUE)

betadisper_tardi_Season <- betadisper(BCdist_tardi, a_tardi_sample_data$Sample_season)
betadisper_tardi_Season_aov <-anova(betadisper_tardi_Season)
table_tardi_Season_aov <- as.data.frame(betadisper_tardi_Season_aov)
table_tardi_Season_aov$Group <- "Tardigrades"
table_tardi_Season_aov$group <- "Sampling season"
betadisper_tardi_Season_groups <- as.data.frame.table(betadisper_tardi_Season$group.distances)
betadisper_tardi_Season_groups$Group <- "Tardigrades"
betadisper_tardi_Season_groups$group <- "Sampling season"

#Permutest
permutest(betadisper_tardi_Season, permutations = 999, pairwise = TRUE)

betadisper_tardi_depth <-betadisper(BCdist_tardi, a_tardi_sample_data$depth_grouped)
betadisper_tardi_depth_aov <-anova(betadisper_tardi_depth)
table_tardi_depth_aov <- as.data.frame(betadisper_tardi_depth_aov)
table_tardi_depth_aov$Group <- "Tardigrades"
table_tardi_depth_aov$group <- "Depth"
betadisper_tardi_depth_groups <- as.data.frame.table(betadisper_tardi_depth$group.distances)
betadisper_tardi_depth_groups$Group <- "Tardigrades"
betadisper_tardi_depth_groups$group <- "Soil depth"
#Permutest
permutest(betadisper_tardi_depth, permutations = 999, pairwise = TRUE)

betadisper_tardi_erosion <- betadisper(BCdist_tardi, a_tardi_sample_data$Erosion_grouped)
betadisper_tardi_erosion_aov <- anova(betadisper_tardi_erosion)
table_tardi_erosion_aov <- as.data.frame(betadisper_tardi_erosion_aov)
table_tardi_erosion_aov$Group <- "Tardigrades"
table_tardi_erosion_aov$group <- "Erosion"
betadisper_tardi_erosion_groups <- as.data.frame.table(betadisper_tardi_erosion$group.distances)
betadisper_tardi_erosion_groups$Group <- "Tardigrades"
betadisper_tardi_erosion_groups$group <- "Erosion risk"
#Permutest
permutest(betadisper_tardi_erosion, permutations = 999, pairwise = TRUE)

betadisper_tardi_temp <- betadisper(BCdist_tardi, a_tardi_sample_data$cluster_temp)
betadisper_tardi_temp_aov <- anova(betadisper_tardi_temp)
table_tardi_temp_aov <- as.data.frame(betadisper_tardi_temp_aov)
table_tardi_temp_aov$Group <- "Tardigrades"
table_tardi_temp_aov$group <- "Temperature"
betadisper_tardi_temp_groups <- as.data.frame.table(betadisper_tardi_temp$group.distances)
betadisper_tardi_temp_groups$Group <- "Tardigrades"
betadisper_tardi_temp_groups$group <- "Temperature"
#Permutest
permutest(betadisper_tardi_temp, permutations = 999, pairwise = TRUE)

betadisper_tardi_temp_range <- betadisper(BCdist_tardi, a_tardi_sample_data$cluster_temp_range)
betadisper_tardi_temp_range_aov <-anova(betadisper_tardi_temp_range)
table_tardi_temp_range_aov <- as.data.frame(betadisper_tardi_temp_range_aov)
table_tardi_temp_range_aov$Group <- "Tardigrades"
table_tardi_temp_range_aov$group <- "Temperature range"
betadisper_tardi_temp_range_groups <- as.data.frame.table(betadisper_tardi_temp_range$group.distances)
betadisper_tardi_temp_range_groups$Group <- "Tardigrades"
betadisper_tardi_temp_range_groups$group <- "Temperature Range"
#Permutest
permutest(betadisper_tardi_temp_range, permutations = 999, pairwise = TRUE)

betadisper_tardi_prec <- betadisper(BCdist_tardi, a_tardi_sample_data$cluster_prec)
betadisper_tardi_prec_aov <- anova(betadisper_tardi_prec)
table_tardi_prec_aov <- as.data.frame(betadisper_tardi_prec_aov)
table_tardi_prec_aov$Group <- "Tardigrades"
table_tardi_prec_aov$group <- "Precipitation"
betadisper_tardi_prec_groups <- as.data.frame.table(betadisper_tardi_prec$group.distances)
betadisper_tardi_prec_groups$Group <- "Tardigrades"
betadisper_tardi_prec_groups$group <- "Precipitation"
#Permutest
permutest(betadisper_tardi_prec, permutations = 999, pairwise = TRUE)

betadisper_tardi_Cmic <- betadisper(BCdist_tardi, a_tardi_sample_data$cluster_Cmic)
betadisper_tardi_Cmic_aov <- anova(betadisper_tardi_Cmic)
table_tardi_Cmic_aov <- as.data.frame(betadisper_tardi_Cmic_aov)
table_tardi_Cmic_aov$Group <- "Tardigrades"
table_tardi_Cmic_aov$group <- "Cmic"
betadisper_tardi_Cmic_groups <- as.data.frame.table(betadisper_tardi_Cmic$group.distances)
betadisper_tardi_Cmic_groups$Group <- "Tardigrades"
betadisper_tardi_Cmic_groups$group <- "Cmic"
#Permutest
permutest(betadisper_tardi_Cmic, permutations = 999, pairwise = TRUE)

#nematodes
betadisper_nema_LC1 <- betadisper(BCdist_nema, a_nema_sample_data$LC1_2018)
betadisper_nema_LC1_aov <-anova(betadisper_nema_LC1)
table_nema_LC1_aov <- as.data.frame(betadisper_nema_LC1_aov)
table_nema_LC1_aov$Group <- "Nematodes"
table_nema_LC1_aov$group <- "Ecosystem type"
#Permutest
permutest(betadisper_nema_LC1, permutations = 999, pairwise = TRUE)
#Extract distances to centorid betadisper
betadisper_nema_LC1_groups <- as.data.frame.table(betadisper_nema_LC1$group.distances)
betadisper_nema_LC1_groups$Group <- "Nematodes"
betadisper_nema_LC1_groups$group <- "Ecosystem type"

betadisper_nema_LU2 <- betadisper(BCdist_nema, a_nema_sample_data$LU2)
betadisper_nema_LU2_aov <- anova(betadisper_nema_LU2)
table_nema_LU2_aov <- as.data.frame(betadisper_nema_LU2_aov)
table_nema_LU2_aov$Group <- "Nematodes"
table_nema_LU2_aov$group <- "Plant cover"
# #Permutest
permutest(betadisper_nema_LU2, permutations = 999, pairwise = TRUE)
#Extract distances to centorid betadisper
betadisper_nema_LU2_groups <- as.data.frame.table(betadisper_nema_LU2$group.distances)
betadisper_nema_LU2_groups$Group <- "Nematodes"
betadisper_nema_LU2_groups$group <- "Plant cover"

betadisper_nema_pH <- betadisper(BCdist_nema, a_nema_sample_data$pH_grouped)
betadisper_nema_pH_aov <-anova(betadisper_nema_pH)
table_nema_pH_aov <- as.data.frame(betadisper_nema_pH_aov)
table_nema_pH_aov$Group <- "Nematodes"
table_nema_pH_aov$group <- "pH"
anova(betadisper_nema_pH)
betadisper_nema_pH_groups <- as.data.frame.table(betadisper_nema_pH$group.distances)
betadisper_nema_pH_groups$Group <- "Nematodes"
betadisper_nema_pH_groups$group <- "pH"
#Permutest
permutest(betadisper_nema_pH, permutations = 999, pairwise = TRUE)

betadisper_nema_Season <- betadisper(BCdist_nema, a_nema_sample_data$Sample_season)
betadisper_nema_Season_aov <-anova(betadisper_nema_Season)
table_nema_Season_aov <- as.data.frame(betadisper_nema_Season_aov)
table_nema_Season_aov$Group <- "Nematodes"
table_nema_Season_aov$group <- "Season"
betadisper_nema_Season_groups <- as.data.frame.table(betadisper_nema_Season$group.distances)
betadisper_nema_Season_groups$Group <- "Nematodes"
betadisper_nema_Season_groups$group <- "Sampling season"
#Permutest
permutest(betadisper_nema_Season, permutations = 999, pairwise = TRUE)

betadisper_nema_depth <-betadisper(BCdist_nema, a_nema_sample_data$depth_grouped)
betadisper_nema_depth_aov <- anova(betadisper_nema_depth)
table_nema_depth_aov <- as.data.frame(betadisper_nema_depth_aov)
table_nema_depth_aov$Group <- "Nematodes"
table_nema_depth_aov$group <- "Depth"
betadisper_nema_depth_groups <- as.data.frame.table(betadisper_nema_depth$group.distances)
betadisper_nema_depth_groups$Group <- "Nematodes"
betadisper_nema_depth_groups$group <- "Soil depth"
#Permutest
permutest(betadisper_nema_depth, permutations = 999, pairwise = TRUE)

betadisper_nema_erosion <- betadisper(BCdist_nema, a_nema_sample_data$Erosion_grouped)
betadisper_nema_erosion_aov <- anova(betadisper_nema_erosion)
table_nema_erosion_aov <- as.data.frame(betadisper_nema_erosion_aov)
table_nema_erosion_aov$Group <- "Nematodes"
table_nema_erosion_aov$group <- "Erosion"
permutest(betadisper_nema_erosion, permutations = 999, pairwise = TRUE)
betadisper_nema_erosion_groups <- as.data.frame.table(betadisper_nema_erosion$group.distances)
betadisper_nema_erosion_groups$Group <- "Nematodes"
betadisper_nema_erosion_groups$group <- "Erosion risk"

betadisper_nema_temp <- betadisper(BCdist_nema, a_nema_sample_data$cluster_temp)
betadisper_nema_temp_aov <- anova(betadisper_nema_temp)
table_nema_temp_aov <- as.data.frame(betadisper_nema_temp_aov)
table_nema_temp_aov$Group <- "Nematodes"
table_nema_temp_aov$group <- "Temperature"
betadisper_nema_temp_groups <- as.data.frame.table(betadisper_nema_temp$group.distances)
betadisper_nema_temp_groups$Group <- "Nematodes"
betadisper_nema_temp_groups$group <- "Temperature"
#Permutest
permutest(betadisper_nema_temp, permutations = 999, pairwise = TRUE)

betadisper_nema_temp_range <- betadisper(BCdist_nema, a_nema_sample_data$cluster_temp_range)
betadisper_nema_temp_range_aov <- anova(betadisper_nema_temp_range)
table_nema_temp_range_aov <- as.data.frame(betadisper_nema_temp_range_aov)
table_nema_temp_range_aov$Group <- "Nematodes"
table_nema_temp_range_aov$group <- "Temperature range"
betadisper_nema_temp_range_groups <- as.data.frame.table(betadisper_nema_temp_range$group.distances)
betadisper_nema_temp_range_groups$Group <- "Nematodes"
betadisper_nema_temp_range_groups$group <- "Temperature Range"
#Permutest
permutest(betadisper_nema_temp_range, permutations = 999, pairwise = TRUE)

betadisper_nema_prec <- betadisper(BCdist_nema, a_nema_sample_data$cluster_prec)
betadisper_nema_prec_aov <-anova(betadisper_nema_prec)
table_nema_prec_aov <- as.data.frame(betadisper_nema_prec_aov)
table_nema_prec_aov$Group <- "Nematodes"
table_nema_prec_aov$group <- "Precipitation"
betadisper_nema_prec_groups <- as.data.frame.table(betadisper_nema_prec$group.distances)
betadisper_nema_prec_groups$Group <- "Nematodes"
betadisper_nema_prec_groups$group <- "Precipitation"
#Permutest
permutest(betadisper_nema_prec, permutations = 999, pairwise = TRUE)

betadisper_nema_Cmic <- betadisper(BCdist_nema, a_nema_sample_data$cluster_Cmic)
betadisper_nema_Cmic_aov <-anova(betadisper_nema_Cmic)
table_nema_Cmic_aov <- as.data.frame(betadisper_nema_Cmic_aov)
table_nema_Cmic_aov$Group <- "Nematodes"
table_nema_Cmic_aov$group <- "Cmic"
betadisper_nema_Cmic_groups <- as.data.frame.table(betadisper_nema_Cmic$group.distances)
betadisper_nema_Cmic_groups$Group <- "Nematodes"
betadisper_nema_Cmic_groups$group <- "Cmic"
#Permutest
permutest(betadisper_nema_Cmic, permutations = 999, pairwise = TRUE)


#Arthropods
betadisper_arthro_LC1 <- betadisper(BCdist_arthro, a_arthro_sample_data$LC1_2018)
betadisper_arthro_LC1_aov <-anova(betadisper_arthro_LC1)
table_arthro_LC1_aov <- as.data.frame(betadisper_arthro_LC1_aov)
table_arthro_LC1_aov$Group <- "Arthropods"
table_arthro_LC1_aov$group <- "Ecosystem type"
#Permutest
permutest(betadisper_arthro_LC1, permutations = 999, pairwise = TRUE)
#Extract distances to centorid betadisper
betadisper_arthro_LC1_groups <- as.data.frame.table(betadisper_arthro_LC1$group.distances)
betadisper_arthro_LC1_groups$Group <- "Arthropods"
betadisper_arthro_LC1_groups$group <- "Ecosystem type"

betadisper_arthro_LU2 <- betadisper(BCdist_arthro, a_arthro_sample_data$LU2)
betadisper_arthro_LU2_aov <- anova(betadisper_arthro_LU2)
table_arthro_LU2_aov <- as.data.frame(betadisper_arthro_LU2_aov)
table_arthro_LU2_aov$Group <- "Arthropods"
table_arthro_LU2_aov$group <- "Plant cover"
# #Permutest
permutest(betadisper_arthro_LU2, permutations = 999, pairwise = TRUE)
#Extract distances to centorid betadisper
betadisper_arthro_LU2_groups <- as.data.frame.table(betadisper_arthro_LU2$group.distances)
betadisper_arthro_LU2_groups$Group <- "Arthropods"
betadisper_arthro_LU2_groups$group <- "Plant cover"

betadisper_arthro_pH <- betadisper(BCdist_arthro, a_arthro_sample_data$pH_grouped)
betadisper_arthro_pH_aov <-anova(betadisper_arthro_pH)
table_arthro_pH_aov <- as.data.frame(betadisper_arthro_pH_aov)
table_arthro_pH_aov$Group <- "Arthropods"
table_arthro_pH_aov$group <- "pH"
anova(betadisper_arthro_pH)
betadisper_arthro_pH_groups <- as.data.frame.table(betadisper_arthro_pH$group.distances)
betadisper_arthro_pH_groups$Group <- "Arthropods"
betadisper_arthro_pH_groups$group <- "pH"
#Permutest
permutest(betadisper_arthro_pH, permutations = 999, pairwise = TRUE)

betadisper_arthro_Season <- betadisper(BCdist_arthro, a_arthro_sample_data$Sample_season)
betadisper_arthro_Season_aov <-anova(betadisper_arthro_Season)
table_arthro_Season_aov <- as.data.frame(betadisper_arthro_Season_aov)
table_arthro_Season_aov$Group <- "Arthropods"
table_arthro_Season_aov$group <- "Season"
betadisper_arthro_Season_groups <- as.data.frame.table(betadisper_arthro_Season$group.distances)
betadisper_arthro_Season_groups$Group <- "Arthropods"
betadisper_arthro_Season_groups$group <- "Sampling season"
#Permutest
permutest(betadisper_arthro_Season, permutations = 999, pairwise = TRUE)

betadisper_arthro_depth <-betadisper(BCdist_arthro, a_arthro_sample_data$depth_grouped)
betadisper_arthro_depth_aov <- anova(betadisper_arthro_depth)
table_arthro_depth_aov <- as.data.frame(betadisper_arthro_depth_aov)
table_arthro_depth_aov$Group <- "Arthropods"
table_arthro_depth_aov$group <- "Depth"
betadisper_arthro_depth_groups <- as.data.frame.table(betadisper_arthro_depth$group.distances)
betadisper_arthro_depth_groups$Group <- "Arthropods"
betadisper_arthro_depth_groups$group <- "Soil depth"
#Permutest
permutest(betadisper_arthro_depth, permutations = 999, pairwise = TRUE)

betadisper_arthro_erosion <- betadisper(BCdist_arthro, a_arthro_sample_data$Erosion_grouped)
betadisper_arthro_erosion_aov <- anova(betadisper_arthro_erosion)
table_arthro_erosion_aov <- as.data.frame(betadisper_arthro_erosion_aov)
table_arthro_erosion_aov$Group <- "Arthropods"
table_arthro_erosion_aov$group <- "Erosion"
permutest(betadisper_arthro_erosion, permutations = 999, pairwise = TRUE)
betadisper_arthro_erosion_groups <- as.data.frame.table(betadisper_arthro_erosion$group.distances)
betadisper_arthro_erosion_groups$Group <- "Arthropods"
betadisper_arthro_erosion_groups$group <- "Erosion risk"

betadisper_arthro_temp <- betadisper(BCdist_arthro, a_arthro_sample_data$cluster_temp)
betadisper_arthro_temp_aov <- anova(betadisper_arthro_temp)
table_arthro_temp_aov <- as.data.frame(betadisper_arthro_temp_aov)
table_arthro_temp_aov$Group <- "Arthropods"
table_arthro_temp_aov$group <- "Temperature"
betadisper_arthro_temp_groups <- as.data.frame.table(betadisper_arthro_temp$group.distances)
betadisper_arthro_temp_groups$Group <- "Arthropods"
betadisper_arthro_temp_groups$group <- "Temperature"
#Permutest
permutest(betadisper_arthro_temp, permutations = 999, pairwise = TRUE)

betadisper_arthro_temp_range <- betadisper(BCdist_arthro, a_arthro_sample_data$cluster_temp_range)
betadisper_arthro_temp_range_aov <- anova(betadisper_arthro_temp_range)
table_arthro_temp_range_aov <- as.data.frame(betadisper_arthro_temp_range_aov)
table_arthro_temp_range_aov$Group <- "Arthropods"
table_arthro_temp_range_aov$group <- "Temperature range"
betadisper_arthro_temp_range_groups <- as.data.frame.table(betadisper_arthro_temp_range$group.distances)
betadisper_arthro_temp_range_groups$Group <- "Arthropods"
betadisper_arthro_temp_range_groups$group <- "Temperature Range"
#Permutest
permutest(betadisper_arthro_temp_range, permutations = 999, pairwise = TRUE)

betadisper_arthro_prec <- betadisper(BCdist_arthro, a_arthro_sample_data$cluster_prec)
betadisper_arthro_prec_aov <-anova(betadisper_arthro_prec)
table_arthro_prec_aov <- as.data.frame(betadisper_arthro_prec_aov)
table_arthro_prec_aov$Group <- "Arthropods"
table_arthro_prec_aov$group <- "Precipitation"
betadisper_arthro_prec_groups <- as.data.frame.table(betadisper_arthro_prec$group.distances)
betadisper_arthro_prec_groups$Group <- "Arthropods"
betadisper_arthro_prec_groups$group <- "Precipitation"
#Permutest
permutest(betadisper_arthro_prec, permutations = 999, pairwise = TRUE)

betadisper_arthro_Cmic <- betadisper(BCdist_arthro, a_arthro_sample_data$cluster_Cmic)
betadisper_arthro_Cmic_aov <-anova(betadisper_arthro_Cmic)
table_arthro_Cmic_aov <- as.data.frame(betadisper_arthro_Cmic_aov)
table_arthro_Cmic_aov$Group <- "Arthropods"
table_arthro_Cmic_aov$group <- "Cmic"
betadisper_arthro_Cmic_groups <- as.data.frame.table(betadisper_arthro_Cmic$group.distances)
betadisper_arthro_Cmic_groups$Group <- "Arthropods"
betadisper_arthro_Cmic_groups$group <- "Cmic"
#Permutest
permutest(betadisper_arthro_Cmic, permutations = 999, pairwise = TRUE)


#Annelids
betadisper_anneli_LC1 <- betadisper(BCdist_anneli, a_anneli_sample_data$LC1_2018)
betadisper_anneli_LC1_aov <-anova(betadisper_anneli_LC1)
table_anneli_LC1_aov <- as.data.frame(betadisper_anneli_LC1_aov)
table_anneli_LC1_aov$Group <- "Annelids"
table_anneli_LC1_aov$group <- "Ecosystem type"
#Permutest
permutest(betadisper_anneli_LC1, permutations = 999, pairwise = TRUE)
#Extract distances to centorid betadisper
betadisper_anneli_LC1_groups <- as.data.frame.table(betadisper_anneli_LC1$group.distances)
betadisper_anneli_LC1_groups$Group <- "Annelids"
betadisper_anneli_LC1_groups$group <- "Ecosystem type"

betadisper_anneli_LU2 <- betadisper(BCdist_anneli, a_anneli_sample_data$LU2)
betadisper_anneli_LU2_aov <- anova(betadisper_anneli_LU2)
table_anneli_LU2_aov <- as.data.frame(betadisper_anneli_LU2_aov)
table_anneli_LU2_aov$Group <- "Annelids"
table_anneli_LU2_aov$group <- "Plant cover"
# #Permutest
permutest(betadisper_anneli_LU2, permutations = 999, pairwise = TRUE)
#Extract distances to centorid betadisper
betadisper_anneli_LU2_groups <- as.data.frame.table(betadisper_anneli_LU2$group.distances)
betadisper_anneli_LU2_groups$Group <- "Annelids"
betadisper_anneli_LU2_groups$group <- "Plant cover"

betadisper_anneli_pH <- betadisper(BCdist_anneli, a_anneli_sample_data$pH_grouped)
betadisper_anneli_pH_aov <-anova(betadisper_anneli_pH)
table_anneli_pH_aov <- as.data.frame(betadisper_anneli_pH_aov)
table_anneli_pH_aov$Group <- "Annelids"
table_anneli_pH_aov$group <- "pH"
anova(betadisper_anneli_pH)
betadisper_anneli_pH_groups <- as.data.frame.table(betadisper_anneli_pH$group.distances)
betadisper_anneli_pH_groups$Group <- "Annelids"
betadisper_anneli_pH_groups$group <- "pH"
#Permutest
permutest(betadisper_anneli_pH, permutations = 999, pairwise = TRUE)

betadisper_anneli_Season <- betadisper(BCdist_anneli, a_anneli_sample_data$Sample_season)
betadisper_anneli_Season_aov <-anova(betadisper_anneli_Season)
table_anneli_Season_aov <- as.data.frame(betadisper_anneli_Season_aov)
table_anneli_Season_aov$Group <- "Annelids"
table_anneli_Season_aov$group <- "Season"
betadisper_anneli_Season_groups <- as.data.frame.table(betadisper_anneli_Season$group.distances)
betadisper_anneli_Season_groups$Group <- "Annelids"
betadisper_anneli_Season_groups$group <- "Sampling season"
#Permutest
permutest(betadisper_anneli_Season, permutations = 999, pairwise = TRUE)

betadisper_anneli_depth <-betadisper(BCdist_anneli, a_anneli_sample_data$depth_grouped)
betadisper_anneli_depth_aov <- anova(betadisper_anneli_depth)
table_anneli_depth_aov <- as.data.frame(betadisper_anneli_depth_aov)
table_anneli_depth_aov$Group <- "Annelids"
table_anneli_depth_aov$group <- "Depth"
betadisper_anneli_depth_groups <- as.data.frame.table(betadisper_anneli_depth$group.distances)
betadisper_anneli_depth_groups$Group <- "Annelids"
betadisper_anneli_depth_groups$group <- "Soil depth"
#Permutest
permutest(betadisper_anneli_depth, permutations = 999, pairwise = TRUE)

betadisper_anneli_erosion <- betadisper(BCdist_anneli, a_anneli_sample_data$Erosion_grouped)
betadisper_anneli_erosion_aov <- anova(betadisper_anneli_erosion)
table_anneli_erosion_aov <- as.data.frame(betadisper_anneli_erosion_aov)
table_anneli_erosion_aov$Group <- "Annelids"
table_anneli_erosion_aov$group <- "Erosion"
permutest(betadisper_anneli_erosion, permutations = 999, pairwise = TRUE)
betadisper_anneli_erosion_groups <- as.data.frame.table(betadisper_anneli_erosion$group.distances)
betadisper_anneli_erosion_groups$Group <- "Annelids"
betadisper_anneli_erosion_groups$group <- "Erosion risk"

betadisper_anneli_temp <- betadisper(BCdist_anneli, a_anneli_sample_data$cluster_temp)
betadisper_anneli_temp_aov <- anova(betadisper_anneli_temp)
table_anneli_temp_aov <- as.data.frame(betadisper_anneli_temp_aov)
table_anneli_temp_aov$Group <- "Annelids"
table_anneli_temp_aov$group <- "Temperature"
betadisper_anneli_temp_groups <- as.data.frame.table(betadisper_anneli_temp$group.distances)
betadisper_anneli_temp_groups$Group <- "Annelids"
betadisper_anneli_temp_groups$group <- "Temperature"
#Permutest
permutest(betadisper_anneli_temp, permutations = 999, pairwise = TRUE)

betadisper_anneli_temp_range <- betadisper(BCdist_anneli, a_anneli_sample_data$cluster_temp_range)
betadisper_anneli_temp_range_aov <- anova(betadisper_anneli_temp_range)
table_anneli_temp_range_aov <- as.data.frame(betadisper_anneli_temp_range_aov)
table_anneli_temp_range_aov$Group <- "Annelids"
table_anneli_temp_range_aov$group <- "Temperature range"
betadisper_anneli_temp_range_groups <- as.data.frame.table(betadisper_anneli_temp_range$group.distances)
betadisper_anneli_temp_range_groups$Group <- "Annelids"
betadisper_anneli_temp_range_groups$group <- "Temperature Range"
#Permutest
permutest(betadisper_anneli_temp_range, permutations = 999, pairwise = TRUE)

betadisper_anneli_prec <- betadisper(BCdist_anneli, a_anneli_sample_data$cluster_prec)
betadisper_anneli_prec_aov <-anova(betadisper_anneli_prec)
table_anneli_prec_aov <- as.data.frame(betadisper_anneli_prec_aov)
table_anneli_prec_aov$Group <- "Annelids"
table_anneli_prec_aov$group <- "Precipitation"
betadisper_anneli_prec_groups <- as.data.frame.table(betadisper_anneli_prec$group.distances)
betadisper_anneli_prec_groups$Group <- "Annelids"
betadisper_anneli_prec_groups$group <- "Precipitation"
#Permutest
permutest(betadisper_anneli_prec, permutations = 999, pairwise = TRUE)

betadisper_anneli_Cmic <- betadisper(BCdist_anneli, a_anneli_sample_data$cluster_Cmic)
betadisper_anneli_Cmic_aov <-anova(betadisper_anneli_Cmic)
table_anneli_Cmic_aov <- as.data.frame(betadisper_anneli_Cmic_aov)
table_anneli_Cmic_aov$Group <- "Annelids"
table_anneli_Cmic_aov$group <- "Cmic"
betadisper_anneli_Cmic_groups <- as.data.frame.table(betadisper_anneli_Cmic$group.distances)
betadisper_anneli_Cmic_groups$Group <- "Annelids"
betadisper_anneli_Cmic_groups$group <- "Cmic"
#Permutest
permutest(betadisper_anneli_Cmic, permutations = 999, pairwise = TRUE)



#Protists
betadisper_p_LC1 <- betadisper(BCdist_p, p_sample_data$LC1_2018)
betadisper_p_LC1_aov <- anova(betadisper_p_LC1)
table_p_LC1_aov <- as.data.frame(betadisper_p_LC1_aov)
table_p_LC1_aov$Group <- "Protists"
table_p_LC1_aov$group <- "Ecosystem type"
#Extract distances to centorid betadisper
betadisper_p_LC1_groups <- as.data.frame.table(betadisper_p_LC1$group.distances)
betadisper_p_LC1_groups$Group <- "Protists"
betadisper_p_LC1_groups$group <- "Ecosystem type"
#Permutest
permutest(betadisper_p_LC1, permutations = 999, pairwise = TRUE)

betadisper_p_LU2 <- betadisper(BCdist_p, p_sample_data$LU2)
betadisper_p_LU2_aov <- anova(betadisper_p_LU2)
table_p_LU2_aov <- as.data.frame(betadisper_p_LU2_aov)
table_p_LU2_aov$Group <- "Protists"
table_p_LU2_aov$group <- "Plant cover"
# #Permutest
permutest(betadisper_p_LU2, permutations = 999, pairwise = TRUE)
#Extract distances to centorid betadisper
betadisper_p_LU2_groups <- as.data.frame.table(betadisper_p_LU2$group.distances)
betadisper_p_LU2_groups$Group <- "Protists"
betadisper_p_LU2_groups$group <- "Plant cover"

betadisper_p_pH <- betadisper(BCdist_p, p_sample_data$pH_grouped)
betadisper_p_pH_aov <-anova(betadisper_p_pH)
table_p_pH_aov <- as.data.frame(betadisper_p_pH_aov)
table_p_pH_aov$Group <- "Protists"
table_p_pH_aov$group <- "pH"
betadisper_p_pH_groups <- as.data.frame.table(betadisper_p_pH$group.distances)
betadisper_p_pH_groups$Group <- "Protists"
betadisper_p_pH_groups$group <- "pH"
#Permutest
permutest(betadisper_p_pH, permutations = 999, pairwise = TRUE)

betadisper_p_Season <- betadisper(BCdist_p, p_sample_data$Sample_season)
betadisper_p_Season_aov <-anova(betadisper_p_Season)
table_p_Season_aov <- as.data.frame(betadisper_p_Season_aov)
table_p_Season_aov$Group <- "Protists"
table_p_Season_aov$group <- "Season"
betadisper_p_Season_groups <- as.data.frame.table(betadisper_p_Season$group.distances)
betadisper_p_Season_groups$Group <- "Protists"
betadisper_p_Season_groups$group <- "Sampling season"

#Permutest
permutest(betadisper_p_Season, permutations = 999, pairwise = TRUE)

betadisper_p_depth <-betadisper(BCdist_p, p_sample_data$depth_grouped)
betadisper_p_depth_aov <-anova(betadisper_p_depth)
table_p_depth_aov <- as.data.frame(betadisper_p_depth_aov)
table_p_depth_aov$Group <- "Protists"
table_p_depth_aov$group <- "Depth"
betadisper_p_depth_groups <- as.data.frame.table(betadisper_p_depth$group.distances)
betadisper_p_depth_groups$Group <- "Protists"
betadisper_p_depth_groups$group <- "Soil depth"
#Permutest
permutest(betadisper_p_depth, permutations = 999, pairwise = TRUE)

betadisper_p_erosion <- betadisper(BCdist_p, p_sample_data$Erosion_grouped)
betadisper_p_erosion_aov <-anova(betadisper_p_erosion)
table_p_erosion_aov <- as.data.frame(betadisper_p_erosion_aov)
table_p_erosion_aov$Group <- "Protists"
table_p_erosion_aov$group <- "Erosion"
betadisper_p_erosion_groups <- as.data.frame.table(betadisper_p_erosion$group.distances)
betadisper_p_erosion_groups$Group <- "Protists"
betadisper_p_erosion_groups$group <- "Erosion risk"
#Permutest
permutest(betadisper_p_erosion, permutations = 999, pairwise = TRUE)

betadisper_p_temp <- betadisper(BCdist_p, p_sample_data$cluster_temp)
betadisper_p_temp_aov <-anova(betadisper_p_temp)
table_p_temp_aov <- as.data.frame(betadisper_p_temp_aov)
table_p_temp_aov$Group <- "Protists"
table_p_temp_aov$group <- "Temperature"
betadisper_p_temp_groups <- as.data.frame.table(betadisper_p_temp$group.distances)
betadisper_p_temp_groups$Group <- "Protists"
betadisper_p_temp_groups$group <- "Temperature"
#Permutest
permutest(betadisper_p_temp, permutations = 999, pairwise = TRUE)

betadisper_p_temp_range <- betadisper(BCdist_p, p_sample_data$cluster_temp_range)
betadisper_p_temp_range_aov <-anova(betadisper_p_temp_range)
table_p_temp_range_aov <- as.data.frame(betadisper_p_temp_range_aov)
table_p_temp_range_aov$Group <- "Protists"
table_p_temp_range_aov$group <- "Temperature range"
betadisper_p_temp_range_groups <- as.data.frame.table(betadisper_p_temp_range$group.distances)
betadisper_p_temp_range_groups$Group <- "Protists"
betadisper_p_temp_range_groups$group <- "Temperature Range"
#Permutest
permutest(betadisper_p_temp_range, permutations = 999, pairwise = TRUE)

betadisper_p_prec <- betadisper(BCdist_p, p_sample_data$cluster_prec)
betadisper_p_prec_aov <-anova(betadisper_p_prec)
table_p_prec_aov <- as.data.frame(betadisper_p_prec_aov)
table_p_prec_aov$Group <- "Protists"
table_p_prec_aov$group <- "Precipitation"
betadisper_p_prec_groups <- as.data.frame.table(betadisper_p_prec$group.distances)
betadisper_p_prec_groups$Group <- "Protists"
betadisper_p_prec_groups$group <- "Precipitation"
#Permutest
permutest(betadisper_p_prec, permutations = 999, pairwise = TRUE)

betadisper_p_Cmic <- betadisper(BCdist_p, p_sample_data$cluster_Cmic)
betadisper_p_Cmic_aov <-anova(betadisper_p_Cmic)
table_p_Cmic_aov <- as.data.frame(betadisper_p_Cmic_aov)
table_p_Cmic_aov$Group <- "Protists"
table_p_Cmic_aov$group <- "Cmic"
betadisper_p_Cmic_groups <- as.data.frame.table(betadisper_p_Cmic$group.distances)
betadisper_p_Cmic_groups$Group <- "Protists"
betadisper_p_Cmic_groups$group <- "Cmic"
#Permutest
permutest(betadisper_p_Cmic, permutations = 999, pairwise = TRUE)

#Fungi
betadisper_f_LC1 <- betadisper(BCdist_f, f_sample_data$LC1_2018)
betadisper_f_LC1_aov <- anova(betadisper_f_LC1)
table_f_LC1_aov <- as.data.frame(betadisper_f_LC1_aov)
table_f_LC1_aov$Group <- "Fungi"
table_f_LC1_aov$group <- "Ecosystem type"
#Extract distances to centorid betadisper
betadisper_f_LC1_groups <- as.data.frame.table(betadisper_f_LC1$group.distances)
betadisper_f_LC1_groups$Group <- "Fungi"
betadisper_f_LC1_groups$group <- "Ecosystem type"
#Permutest 
permutest_f_LC1 <- permutest(betadisper_f_LC1, permutations = 999, pairwise = TRUE)
permutest_f_LC1
betadisper_f_pH <- betadisper(BCdist_f, f_sample_data$pH_grouped)
betadisper_f_pH_aov <-anova(betadisper_f_pH)
table_f_pH_aov <- as.data.frame(betadisper_f_pH_aov)
table_f_pH_aov$Group <- "Fungi"
table_f_pH_aov$group <- "pH"
betadisper_f_pH_groups <- as.data.frame.table(betadisper_f_pH$group.distances)
betadisper_f_pH_groups$Group <- "Fungi"
betadisper_f_pH_groups$group <- "pH"
#Permutest
permutest(betadisper_f_pH, permutations = 999, pairwise = TRUE)

betadisper_f_LU2 <- betadisper(BCdist_f, f_sample_data$LU2)
betadisper_f_LU2_aov <- anova(betadisper_f_LU2)
table_f_LU2_aov <- as.data.frame(betadisper_f_LU2_aov)
table_f_LU2_aov$Group <- "Fungi"
table_f_LU2_aov$group <- "Plant cover"
# #Permutest
permutest(betadisper_f_LU2, permutations = 999, pairwise = TRUE)
#Extract distances to centorid betadisper
betadisper_f_LU2_groups <- as.data.frame.table(betadisper_f_LU2$group.distances)
betadisper_f_LU2_groups$Group <- "Fungi"
betadisper_f_LU2_groups$group <- "Plant cover"

betadisper_f_Season <- betadisper(BCdist_f, f_sample_data$Sample_season)
betadisper_f_Season_aov <-anova(betadisper_f_Season)
table_f_Season_aov <- as.data.frame(betadisper_f_Season_aov)
table_f_Season_aov$Group <- "Fungi"
table_f_Season_aov$group <- "Season"
betadisper_f_Season_groups <- as.data.frame.table(betadisper_f_Season$group.distances)
betadisper_f_Season_groups$Group <- "Fungi"
betadisper_f_Season_groups$group <- "Sampling season"
#Permutest
permutest(betadisper_f_Season, permutations = 999, pairwise = TRUE)

betadisper_f_depth <- betadisper(BCdist_f, f_sample_data$depth_grouped)
betadisper_f_depth_aov <-anova(betadisper_f_depth)
table_f_depth_aov <- as.data.frame(betadisper_f_depth_aov)
table_f_depth_aov$Group <- "Fungi"
table_f_depth_aov$group <- "Depth"
betadisper_f_depth_groups <- as.data.frame.table(betadisper_f_depth$group.distances)
betadisper_f_depth_groups$Group <- "Fungi"
betadisper_f_depth_groups$group <- "Soil depth"
#Permutest
permutest(betadisper_f_depth, permutations = 999, pairwise = TRUE)

betadisper_f_erosion <- betadisper(BCdist_f, f_sample_data$Erosion_grouped)
betadisper_f_erosion_aov <-anova(betadisper_f_erosion)
table_f_erosion_aov <- as.data.frame(betadisper_f_erosion_aov)
table_f_erosion_aov$Group <- "Fungi"
table_f_erosion_aov$group <- "Erosion"
betadisper_f_erosion_groups <- as.data.frame.table(betadisper_f_erosion$group.distances)
betadisper_f_erosion_groups$Group <- "Fungi"
betadisper_f_erosion_groups$group <- "Erosion risk"
#Permutest
permutest(betadisper_f_erosion, permutations = 999, pairwise = TRUE)

betadisper_f_temp <- betadisper(BCdist_f, f_sample_data$cluster_temp)
betadisper_f_temp_aov <-anova(betadisper_f_temp)
table_f_temp_aov <- as.data.frame(betadisper_f_temp_aov)
table_f_temp_aov$Group <- "Fungi"
table_f_temp_aov$group <- "Temperature"
betadisper_f_temp_groups <- as.data.frame.table(betadisper_f_temp$group.distances)
betadisper_f_temp_groups$Group <- "Fungi"
betadisper_f_temp_groups$group <- "Temperature"
#Permutest
permutest(betadisper_f_temp, permutations = 999, pairwise = TRUE)

betadisper_f_temp_range <- betadisper(BCdist_f, f_sample_data$cluster_temp_range)
betadisper_f_temp_range_aov <-anova(betadisper_f_temp_range)
table_f_temp_range_aov <- as.data.frame(betadisper_f_temp_range_aov)
table_f_temp_range_aov$Group <- "Fungi"
table_f_temp_range_aov$group <- "Temperature range"
betadisper_f_temp_range_groups <- as.data.frame.table(betadisper_f_temp_range$group.distances)
betadisper_f_temp_range_groups$Group <- "Fungi"
betadisper_f_temp_range_groups$group <- "Temperature Range"
#Permutest
permutest(betadisper_f_temp_range, permutations = 999, pairwise = TRUE)

betadisper_f_prec <- betadisper(BCdist_f, f_sample_data$cluster_prec)
betadisper_f_prec_aov <-anova(betadisper_f_prec)
table_f_prec_aov <- as.data.frame(betadisper_f_prec_aov)
table_f_prec_aov$Group <- "Fungi"
table_f_prec_aov$group <- "Precipitation"
betadisper_f_prec_groups <- as.data.frame.table(betadisper_f_prec$group.distances)
betadisper_f_prec_groups$Group <- "Fungi"
betadisper_f_prec_groups$group <- "Precipitation"
#Permutest
permutest(betadisper_f_prec, permutations = 999, pairwise = TRUE)

betadisper_f_Cmic <- betadisper(BCdist_f, f_sample_data$cluster_Cmic)
betadisper_f_Cmic_aov <-anova(betadisper_f_Cmic)
table_f_Cmic_aov <- as.data.frame(betadisper_f_Cmic_aov)
table_f_Cmic_aov$Group <- "Fungi"
table_f_Cmic_aov$group <- "Cmic"
betadisper_f_Cmic_groups <- as.data.frame.table(betadisper_f_Cmic$group.distances)
betadisper_f_Cmic_groups$Group <- "Fungi"
betadisper_f_Cmic_groups$group <- "Cmic"
#Permutest
permutest(betadisper_f_Cmic, permutations = 999, pairwise = TRUE)


#betadisper distance
betadisper_distances <- rbind(betadisper_f_LC1_groups,betadisper_f_LU2_groups,betadisper_f_pH_groups,betadisper_f_Season_groups,betadisper_f_erosion_groups,betadisper_f_depth_groups,betadisper_f_prec_groups,betadisper_f_temp_range_groups,betadisper_f_temp_groups,betadisper_f_Cmic_groups,
                              betadisper_p_LC1_groups,betadisper_p_LU2_groups,betadisper_p_pH_groups,betadisper_p_Season_groups,betadisper_p_erosion_groups,betadisper_p_depth_groups,betadisper_p_prec_groups,betadisper_p_temp_range_groups,betadisper_p_temp_groups,betadisper_p_Cmic_groups,
                              betadisper_roti_LC1_groups,betadisper_roti_LU2_groups,betadisper_roti_pH_groups,betadisper_roti_Season_groups,betadisper_roti_erosion_groups,betadisper_roti_depth_groups,betadisper_roti_prec_groups,betadisper_roti_temp_range_groups,betadisper_roti_temp_groups,betadisper_roti_Cmic_groups,
                              betadisper_nema_LC1_groups,betadisper_nema_LU2_groups,betadisper_nema_pH_groups,betadisper_nema_Season_groups,betadisper_nema_erosion_groups,betadisper_nema_depth_groups,betadisper_nema_prec_groups,betadisper_nema_temp_range_groups,betadisper_nema_temp_groups,betadisper_nema_Cmic_groups,
                              betadisper_tardi_LC1_groups,betadisper_tardi_LU2_groups,betadisper_tardi_pH_groups,betadisper_tardi_Season_groups,betadisper_tardi_erosion_groups,betadisper_tardi_depth_groups,betadisper_tardi_prec_groups,betadisper_tardi_temp_range_groups,betadisper_tardi_temp_groups,betadisper_tardi_Cmic_groups,
                              betadisper_arthro_LC1_groups,betadisper_arthro_LU2_groups,betadisper_arthro_pH_groups,betadisper_arthro_Season_groups,betadisper_arthro_erosion_groups,betadisper_arthro_depth_groups,betadisper_arthro_prec_groups,betadisper_arthro_temp_range_groups,betadisper_arthro_temp_groups,betadisper_arthro_Cmic_groups,
                              betadisper_anneli_LC1_groups,betadisper_anneli_LU2_groups,betadisper_anneli_pH_groups,betadisper_anneli_Season_groups,betadisper_anneli_erosion_groups,betadisper_anneli_depth_groups,betadisper_anneli_prec_groups,betadisper_anneli_temp_range_groups,betadisper_anneli_temp_groups,betadisper_anneli_Cmic_groups)
write.csv(betadisper_distances, file = "2405_betadisper_distances.csv", row.names = F)

anova_betadisper_table <-rbind(table_f_LC1_aov,table_f_LU2_aov,table_f_pH_aov,table_f_erosion_aov,table_f_depth_aov,table_f_prec_aov,table_f_temp_aov,table_f_temp_range_aov,table_f_Cmic_aov,
                        table_p_LC1_aov,table_p_LU2_aov,table_p_pH_aov,table_p_erosion_aov,table_p_depth_aov,table_p_prec_aov,table_p_temp_aov,table_p_temp_range_aov,table_p_Cmic_aov,
                       table_roti_LC1_aov,table_roti_LU2_aov,table_roti_pH_aov,table_roti_erosion_aov,table_roti_depth_aov,table_roti_prec_aov,table_roti_temp_aov,table_roti_temp_range_aov,table_roti_Cmic_aov,
                        table_nema_LC1_aov,table_nema_LU2_aov,table_nema_pH_aov,table_nema_erosion_aov,table_nema_depth_aov,table_nema_prec_aov,table_nema_temp_aov,table_nema_temp_range_aov,table_nema_Cmic_aov,
                        table_tardi_LC1_aov,table_tardi_LU2_aov,table_tardi_pH_aov,table_tardi_erosion_aov,table_tardi_depth_aov,table_tardi_prec_aov,table_tardi_temp_aov,table_tardi_temp_range_aov,table_tardi_Cmic_aov,
                       table_arthro_LC1_aov,table_arthro_LU2_aov,table_arthro_pH_aov,table_arthro_erosion_aov,table_arthro_depth_aov,table_arthro_prec_aov,table_arthro_temp_aov,table_arthro_temp_range_aov,table_arthro_Cmic_aov,
                       table_anneli_LC1_aov,table_anneli_LU2_aov,table_anneli_pH_aov,table_anneli_erosion_aov,table_anneli_depth_aov,table_anneli_prec_aov,table_anneli_temp_aov,table_anneli_temp_range_aov,table_anneli_Cmic_aov)
write.csv(anova_betadisper_table, file = "2405_anova_betadisper_table.csv", row.names = F)


betadisper_distances$Group <- factor(betadisper_distances$Group, levels=c("Fungi","Protists","Rotifers", "Tardigrades", "Nematodes","Arthropods","Annelids" ))
Cmic_cols <- c("Very low"="grey","Low"="brown","Medium low"="#C71585","Medium high"="#FF69B4","High"="blue","Very high"="green")
Depth_cols <- c("Shallow"="green","Moderately deep"="#FEFB01","Deep"="orange")
Erosion_cols <- c("Low"="green", "Medium low"="#FEFB01","High"="orange","Very high"="#FF99CC")
LC_cols <- c("CL_annual"="#FEFB01", "CL_permanent"="orange","GL_managed"="green","GL_unmanaged"="#40BAA4","WL_broadleaved"="#FD79AD","WL_coniferous"="#FE036A")
pH_cols <- c("Acidic"="#FF99CC", "Neutral"="green", "Alkaline"="#59bfff")
Prec_cols <- c("Very low"="red","Low"="orange","Medium low"="yellow","Medium high"="#71C671","High"="#388E8E","Very high"="#7D9EC0")
Season_cols <- c("Spring"="green","Autumn"="orange","Summer"="#FEFB01")
Temp_cols <- c("Very low"="#7D9EC0","Low"="#388E8E","Medium low"="#71C671","Medium high"="yellow","High temperature"="orange","Very high"="red")
Temp_range_cols <- c("Very low"="#7D9EC0","Low"="#388E8E","Medium low"="#71C671","Medium high"="yellow","High"="orange","Very high"="red")
LU2_cols <- c("CL1"="#FEFB01","CL2"="#ECC979","CL3"="orange","CL4"="#ECA579","GL1"="green","GL2"="#40BAA4","GL3"="#1F8276","WL1"="#eec2c9","WL2"="#FD79AD","WL3"="#EC79D7","WL4"="#FE036A","No_C"="grey")

betadisper_distances$Group <- factor(betadisper_distances$Group, levels=c("Fungi","Protists","Rotifers", "Tardigrades", "Nematodes","Arthropods","Annelids"))
Cmic_betadisper_distances <- betadisper_distances %>% filter(group == "Cmic")
pH_betadisper_distances <- betadisper_distances %>% filter(group == "pH")
Depth_betadisper_distances <- betadisper_distances %>% filter(group == "Soil depth")
Erosion_betadisper_distances <- betadisper_distances %>% filter(group == "Erosion risk")
Season_betadisper_distances <- betadisper_distances %>% filter(group == "Sampling season")
LC1_betadisper_distances <- betadisper_distances %>% filter(group == "Ecosystem type")
Prec_betadisper_distances <- betadisper_distances %>% filter(group == "Precipitation")
Temp_range_betadisper_distances <- betadisper_distances %>% filter(group == "Temperature Range")
Temp_betadisper_distances <- betadisper_distances %>% filter(group == "Temperature")
LU2_betadisper_distances <- betadisper_distances %>% filter(group == "Plant cover")

Prec_betadisper_distances$Var1 <- factor(Prec_betadisper_distances$Var1, levels=c("Very low", "Low", "Medium low","Medium high", "High", "Very high"))
Temp_betadisper_distances$Var1 <- factor(Temp_betadisper_distances$Var1, levels=c("Very low", "Low", "Medium low","Medium high", "High temperature", "Very high"))
Temp_range_betadisper_distances$Var1 <- factor(Temp_range_betadisper_distances$Var1, levels=c("Very low", "Low", "Medium low","Medium high", "High", "Very high"))
Cmic_betadisper_distances$Var1 <- factor(Cmic_betadisper_distances$Var1, levels=c("Very low", "Low", "Medium low","Medium high", "Highn", "Very high"))
Erosion_betadisper_distances$Var1 <- factor(Erosion_betadisper_distances$Var1, levels=c("Low", "Medium low", "High","Very high"))
Season_betadisper_distances$Var1 <- factor(Season_betadisper_distances$Var1, levels=c("Spring", "Summer", "Autumn"))
pH_betadisper_distances$Var1 <- factor(pH_betadisper_distances$Var1, levels=c("Acidic", "Neutral", "Alkaline"))
Depth_betadisper_distances$Var1 <- factor(Depth_betadisper_distances$Var1, levels=c("Shallow", "Moderately deep", "Deep"))


#Cmic
Cmic_betadisper_plot <- ggplot(Cmic_betadisper_distances, aes(x = Group, y = Freq)) + 
  geom_boxplot(aes(colour = Var1),width = 0.6, outlier.shape = NA, alpha = 0.5) +
  scale_colour_manual("Variable classes", values = Cmic_cols)+
  xlab("Var1") +
  ylab("Distance")  +
  ggtitle("Cmic")+
  theme_jk()+
  theme(axis.text.x =element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x=element_blank())

pH_betadisper_plot <- ggplot(pH_betadisper_distances, aes(x = Group, y = Freq)) + 
  geom_boxplot(aes(colour = Var1),width = 0.6, outlier.shape = NA, alpha = 0.5) +
  scale_colour_manual("Variable classes", values = pH_cols)+
  xlab("Var1") +
  ylab("Distance")  +
  ggtitle("pH")+
  theme_jk()+
  theme(axis.text.x =element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x=element_blank())

Depth_betadisper_plot <- ggplot(Depth_betadisper_distances, aes(x = Group, y = Freq)) + 
  geom_boxplot(aes(colour = Var1),width = 0.6, outlier.shape = NA, alpha = 0.5) +
  scale_colour_manual("Variable classes", values = Depth_cols)+
  xlab("Var1") +
  ylab("Distance")  +
  ggtitle("Soil depth")+
  theme_jk()+
  theme(axis.text.x =element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x=element_blank())

Erosion_betadisper_plot <- ggplot(Erosion_betadisper_distances, aes(x = Group, y = Freq)) + 
  geom_boxplot(aes(colour = Var1),width = 0.6, outlier.shape = NA, alpha = 0.5) +
  scale_colour_manual("Variable classes", values = Erosion_cols)+
  xlab("Var1") +
  ylab("Distance")  +
  ggtitle("Erosion")+
  theme_jk()+
  theme(axis.text.x =element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x=element_blank())

LC1_betadisper_plot <- ggplot(LC1_betadisper_distances, aes(x = Group, y = Freq)) + 
  geom_boxplot(aes(colour = Var1),width = 0.6, outlier.shape = NA, alpha = 0.5) +
  scale_colour_manual("Variable classes", values = LC_cols)+
  xlab("Var1") +
  ylab("Distance")  +
  ggtitle("Ecosystem type")+
  theme_jk()+
  theme(axis.text.x =element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x=element_blank())

LU2_betadisper_plot <- ggplot(LU2_betadisper_distances, aes(x = Group, y = Freq)) + 
  geom_boxplot(aes(colour = Var1),width = 0.6, outlier.shape = NA, alpha = 0.5) +
  scale_colour_manual("Variable classes", values = LU2_cols)+
  xlab("Var1") +
  ylab("Distance")  +
  ggtitle("Plant cover")+
  theme_jk()+
  theme(axis.text.x =element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x=element_blank())

Season_betadisper_plot <- ggplot(Season_betadisper_distances, aes(x = Group, y = Freq)) + 
  geom_boxplot(aes(colour = Var1),width = 0.6, outlier.shape = NA, alpha = 0.5) +
  scale_colour_manual("Variable classes", values = Season_cols)+
  xlab("Var1") +
  ylab("Distance")  +
  ggtitle("Sampling season")+
  theme_jk()+
  theme(axis.text.x =element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x=element_blank())

Temp_betadisper_plot <- ggplot(Temp_betadisper_distances, aes(x = Group, y = Freq)) + 
  geom_boxplot(aes(colour = Var1),width = 0.6, outlier.shape = NA, alpha = 0.5) +
  scale_colour_manual("Variable classes", values = Temp_cols)+
  xlab("Var1") +
  ylab("Distance")  +  
  ggtitle("Temperature")+
  theme_jk()+
  theme(axis.text.x =element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x=element_blank())

Temp_range_betadisper_plot <- ggplot(Temp_range_betadisper_distances, aes(x = Group, y = Freq)) + 
  geom_boxplot(aes(colour = Var1),width = 0.6, outlier.shape = NA, alpha = 0.5) +
  scale_colour_manual("Variable classes", values =Temp_range_cols)+
  xlab("Var1") +
  ylab("Distance")  +
  ggtitle("Temperature range")+
  theme_jk()+
  theme(axis.text.x =element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x=element_blank())

Prec_betadisper_plot <- ggplot(Prec_betadisper_distances, aes(x = Group, y = Freq)) + 
  geom_boxplot(aes(colour = Var1),width = 0.6, outlier.shape = NA, alpha = 0.5) +
  scale_colour_manual("Variable classes", values = Prec_cols)+
  xlab("Var1") +
  ylab("Distance")  +
  ggtitle("Precipitation")+
  theme_jk()+
  theme(axis.text.x =element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x=element_blank())

#plot Figure S10
svg("2405_FigS10.svg",width=19,height=12)
ggarrange(LC1_betadisper_plot, LU2_betadisper_plot, pH_betadisper_plot, Season_betadisper_plot, Depth_betadisper_plot, Erosion_betadisper_plot,Temp_betadisper_plot,Temp_range_betadisper_plot,Prec_betadisper_plot,
          labels = c("(a)", "(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)"),align="v",hjust=c(-1.5),
          ncol = 3, nrow = 3)
dev.off() # Close the graphics device

#save png plot
png("2405_FigS10.png",width=1600,height=1300)
ggarrange(LC1_betadisper_plot, LU2_betadisper_plot, pH_betadisper_plot, Season_betadisper_plot, Depth_betadisper_plot, Erosion_betadisper_plot,Temp_betadisper_plot,Temp_range_betadisper_plot,Prec_betadisper_plot,
          labels = c("(a)", "(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)"),align="v",hjust=c(-1.5),
          ncol = 3, nrow = 3)
dev.off() # Close the graphics device

svg("2405_FigS10_small.svg",width=16,height=13)
ggarrange(LC1_betadisper_plot, LU2_betadisper_plot, pH_betadisper_plot, Season_betadisper_plot, Depth_betadisper_plot, Erosion_betadisper_plot,Temp_betadisper_plot,Temp_range_betadisper_plot,Prec_betadisper_plot,
          labels = c("(a)", "(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)"),align="v",hjust=c(-1.5),
          ncol = 3, nrow = 3)
dev.off() # Close the graphics device


#--> PERMUTEST DISPERSION
#Step 2: Adonis (TESTING IF GROUP COMPOSITION ARE SIMILAR OR NOT)
adonis_roti_LC1 <- adonis2(BCdist_roti ~ a_roti_sample_data$LC1_2018)
adonis_roti_LC1$variable <-"Ecosystem type"
adonis_roti_LC1$group <- "Rotifers"
adonis_roti_LU2 <- adonis2(BCdist_roti ~ a_roti_sample_data$LU2)
adonis_roti_LU2$variable <-"Plant cover"
adonis_roti_LU2$group <- "Rotifers"
adonis_roti_pH <- adonis2(BCdist_roti ~ a_roti_sample_data$pH_grouped)
adonis_roti_pH$variable <-"pH"
adonis_roti_pH$group <- "Rotifers"
adonis_roti_depth <- adonis2(BCdist_roti ~ a_roti_sample_data$depth_grouped)
adonis_roti_depth$variable <-"Soil depth"
adonis_roti_depth$group <- "Rotifers"
adonis_roti_erosion <- adonis2(BCdist_roti ~ a_roti_sample_data$Erosion_grouped)
adonis_roti_erosion$variable <-"Erosion risk"
adonis_roti_erosion$group <- "Rotifers"
adonis_roti_prec <- adonis2(BCdist_roti ~ a_roti_sample_data$cluster_prec)
adonis_roti_prec$variable <-"Precipitation"
adonis_roti_prec$group <- "Rotifers"
adonis_roti_temp <- adonis2(BCdist_roti ~ a_roti_sample_data$cluster_temp)
adonis_roti_temp$variable <-"Temperature"
adonis_roti_temp$group <- "Rotifers"
adonis_roti_temp_range <- adonis2(BCdist_roti ~ a_roti_sample_data$cluster_temp_range)
adonis_roti_temp_range$variable <-"Temperature_range"
adonis_roti_temp_range$group <- "Rotifers"
adonis_roti_Cmic <- adonis2(BCdist_roti ~ a_roti_sample_data$cluster_Cmic)
adonis_roti_Cmic$variable <-"Cmic"
adonis_roti_Cmic$group <- "Rotifers"

adonis_nema_LC1 <- adonis2(BCdist_nema ~ a_nema_sample_data$LC1_2018)
adonis_nema_LC1$variable <-"Ecosystem type"
adonis_nema_LC1$group <- "Nematodes"
adonis_nema_LU2 <- adonis2(BCdist_nema ~ a_nema_sample_data$LU2)
adonis_nema_LU2$variable <-"Plant cover"
adonis_nema_LU2$group <- "Nematodes"
adonis_nema_pH <- adonis2(BCdist_nema ~ a_nema_sample_data$pH_grouped)
adonis_nema_pH$variable <-"pH"
adonis_nema_pH$group <- "Nematodes"
adonis_nema_depth <- adonis2(BCdist_nema ~ a_nema_sample_data$depth_grouped)
adonis_nema_depth$variable <-"Soil depth"
adonis_nema_depth$group <- "Nematodes"
adonis_nema_erosion <- adonis2(BCdist_nema ~ a_nema_sample_data$Erosion_grouped)
adonis_nema_erosion$variable <-"Erosion risk"
adonis_nema_erosion$group <- "Nematodes"
adonis_nema_prec <- adonis2(BCdist_nema ~ a_nema_sample_data$cluster_prec)
adonis_nema_prec$variable <-"Precipitation"
adonis_nema_prec$group <- "Nematodes"
adonis_nema_temp <- adonis2(BCdist_nema ~ a_nema_sample_data$cluster_temp)
adonis_nema_temp$variable <-"Temperature"
adonis_nema_temp$group <- "Nematodes"
adonis_nema_temp_range <- adonis2(BCdist_nema ~ a_nema_sample_data$cluster_temp_range)
adonis_nema_temp_range$variable <-"Temperature_range"
adonis_nema_temp_range$group <- "Nematodes"
adonis_nema_Cmic <- adonis2(BCdist_nema ~ a_nema_sample_data$cluster_Cmic)
adonis_nema_Cmic$variable <-"Cmic"
adonis_nema_Cmic$group <- "Nematodes"

adonis_tardi_LC1 <- adonis2(BCdist_tardi ~ a_tardi_sample_data$LC1_2018)
adonis_tardi_LC1$variable <-"Ecosystem type"
adonis_tardi_LC1$group <- "Tardigrades"
adonis_tardi_LU2 <- adonis2(BCdist_tardi ~ a_tardi_sample_data$LU2)
adonis_tardi_LU2$variable <-"Plant cover"
adonis_tardi_LU2$group <- "Tardigrades"
adonis_tardi_pH <-adonis2(BCdist_tardi ~ a_tardi_sample_data$pH_grouped)
adonis_tardi_pH$variable <-"pH"
adonis_tardi_pH$group <- "Tardigrades"
adonis_tardi_depth <- adonis2(BCdist_tardi ~ a_tardi_sample_data$depth_grouped)
adonis_tardi_depth$variable <-"Soil depth"
adonis_tardi_depth$group <- "Tardigrades"
adonis_tardi_erosion <- adonis2(BCdist_tardi ~ a_tardi_sample_data$Erosion_grouped)
adonis_tardi_erosion$variable <-"Erosion risk"
adonis_tardi_erosion$group <- "Tardigrades"
adonis_tardi_prec <- adonis2(BCdist_tardi ~ a_tardi_sample_data$cluster_prec)
adonis_tardi_prec$variable <-"Precipitation"
adonis_tardi_prec$group <- "Tardigrades"
adonis_tardi_temp <- adonis2(BCdist_tardi ~ a_tardi_sample_data$cluster_temp)
adonis_tardi_temp$variable <-"Temperature"
adonis_tardi_temp$group <- "Tardigrades"
adonis_tardi_temp_range <- adonis2(BCdist_tardi ~ a_tardi_sample_data$cluster_temp_range)
adonis_tardi_temp_range$variable <-"Temperature range"
adonis_tardi_temp_range$group <- "Tardigrades"
adonis_tardi_Cmic <- adonis2(BCdist_tardi ~ a_tardi_sample_data$cluster_Cmic)
adonis_tardi_Cmic$variable <-"Cmic"
adonis_tardi_Cmic$group <- "Tardigrades"

adonis_arthro_LC1 <- adonis2(BCdist_arthro ~ a_arthro_sample_data$LC1_2018)
adonis_arthro_LC1$variable <-"Ecosystem type"
adonis_arthro_LC1$group <- "Arthropods"
adonis_arthro_LU2 <- adonis2(BCdist_arthro ~ a_arthro_sample_data$LU2)
adonis_arthro_LU2$variable <-"Plant cover"
adonis_arthro_LU2$group <- "Arthropods"
adonis_arthro_pH <-adonis2(BCdist_arthro ~ a_arthro_sample_data$pH_grouped)
adonis_arthro_pH$variable <-"pH"
adonis_arthro_pH$group <- "Arthropods"
adonis_arthro_depth <- adonis2(BCdist_arthro ~ a_arthro_sample_data$depth_grouped)
adonis_arthro_depth$variable <-"Soil depth"
adonis_arthro_depth$group <- "Arthropods"
adonis_arthro_erosion <- adonis2(BCdist_arthro ~ a_arthro_sample_data$Erosion_grouped)
adonis_arthro_erosion$variable <-"Erosion risk"
adonis_arthro_erosion$group <- "Arthropods"
adonis_arthro_prec <- adonis2(BCdist_arthro ~ a_arthro_sample_data$cluster_prec)
adonis_arthro_prec$variable <-"Precipitation"
adonis_arthro_prec$group <- "Arthropods"
adonis_arthro_temp <- adonis2(BCdist_arthro ~ a_arthro_sample_data$cluster_temp)
adonis_arthro_temp$variable <-"Temperature"
adonis_arthro_temp$group <- "Arthropods"
adonis_arthro_temp_range <- adonis2(BCdist_arthro ~ a_arthro_sample_data$cluster_temp_range)
adonis_arthro_temp_range$variable <-"Temperature range"
adonis_arthro_temp_range$group <- "Arthropods"
adonis_arthro_Cmic <- adonis2(BCdist_arthro ~ a_arthro_sample_data$cluster_Cmic)
adonis_arthro_Cmic$variable <-"Cmic"
adonis_arthro_Cmic$group <- "Arthropods"

adonis_anneli_LC1 <- adonis2(BCdist_anneli ~ a_anneli_sample_data$LC1_2018)
adonis_anneli_LC1$variable <-"Ecosystem type"
adonis_anneli_LC1$group <- "Annelids"
adonis_anneli_LU2 <- adonis2(BCdist_anneli ~ a_anneli_sample_data$LU2)
adonis_anneli_LU2$variable <-"Plant cover"
adonis_anneli_LU2$group <- "Annelids"
adonis_anneli_pH <-adonis2(BCdist_anneli ~ a_anneli_sample_data$pH_grouped)
adonis_anneli_pH$variable <-"pH"
adonis_anneli_pH$group <- "Annelids"
adonis_anneli_depth <- adonis2(BCdist_anneli ~ a_anneli_sample_data$depth_grouped)
adonis_anneli_depth$variable <-"Soil depth"
adonis_anneli_depth$group <- "Annelids"
adonis_anneli_erosion <- adonis2(BCdist_anneli ~ a_anneli_sample_data$Erosion_grouped)
adonis_anneli_erosion$variable <-"Erosion risk"
adonis_anneli_erosion$group <- "Annelids"
adonis_anneli_prec <- adonis2(BCdist_anneli ~ a_anneli_sample_data$cluster_prec)
adonis_anneli_prec$variable <-"Precipitation"
adonis_anneli_prec$group <- "Annelids"
adonis_anneli_temp <- adonis2(BCdist_anneli ~ a_anneli_sample_data$cluster_temp)
adonis_anneli_temp$variable <-"Temperature"
adonis_anneli_temp$group <- "Annelids"
adonis_anneli_temp_range <- adonis2(BCdist_anneli ~ a_anneli_sample_data$cluster_temp_range)
adonis_anneli_temp_range$variable <-"Temperature range"
adonis_anneli_temp_range$group <- "Annelids"
adonis_anneli_Cmic <- adonis2(BCdist_anneli ~ a_anneli_sample_data$cluster_Cmic)
adonis_anneli_Cmic$variable <-"Cmic"
adonis_anneli_Cmic$group <- "Annelids"


adonis_p_LC1 <- adonis2(BCdist_p ~ p_sample_data$LC1_2018)
adonis_p_LC1$variable <-"Ecosystem type"
adonis_p_LC1$group <- "Protists"
adonis_p_LU2 <- adonis2(BCdist_p ~ p_sample_data$LU2)
adonis_p_LU2$variable <-"Plant cover"
adonis_p_LU2$group <- "Protists"
adonis_p_pH <- adonis2(BCdist_p ~ p_sample_data$pH_grouped)
adonis_p_pH$variable <-"pH"
adonis_p_pH$group <- "Protists"
adonis_p_depth <- adonis2(BCdist_p ~ p_sample_data$depth_grouped)
adonis_p_depth$variable <-"Soil depth"
adonis_p_depth$group <- "Protists"
adonis_p_erosion <- adonis2(BCdist_p ~ p_sample_data$Erosion_grouped)
adonis_p_erosion$variable <-"Erosion risk"
adonis_p_erosion$group <- "Protists"
adonis_p_prec <- adonis2(BCdist_p ~ p_sample_data$cluster_prec)
adonis_p_prec$variable <-"Precipitation"
adonis_p_prec$group <- "Protists"
adonis_p_temp <- adonis2(BCdist_p ~ p_sample_data$cluster_temp)
adonis_p_temp$variable <-"Temperatrure"
adonis_p_temp$group <- "Protists"
adonis_p_temp_range <- adonis2(BCdist_p ~ p_sample_data$cluster_temp_range)
adonis_p_temp_range$variable <-"Temperature range"
adonis_p_temp_range$group <- "Protists"
adonis_p_Cmic <- adonis2(BCdist_p ~ p_sample_data$cluster_Cmic)
adonis_p_Cmic$variable <-"Cmic"
adonis_p_Cmic$group <- "Protists"

adonis_f_LC1 <- adonis2(BCdist_f ~ f_sample_data$LC1_2018)
adonis_f_LC1$variable <-"Ecosystem type"
adonis_f_LC1$group <- "Fungi"
adonis_f_LU2 <- adonis2(BCdist_f ~ f_sample_data$LU2)
adonis_f_LU2$variable <-"Plant cover"
adonis_f_LU2$group <- "Fungi"
adonis_f_pH <- adonis2(BCdist_f ~ f_sample_data$pH_grouped)
adonis_f_pH$variable <-"pH"
adonis_f_pH$group <- "Fungi"
adonis_f_depth <- adonis2(BCdist_f ~ f_sample_data$depth_grouped)
adonis_f_depth$variable <-"Soil depth"
adonis_f_depth$group <- "Fungi"
adonis_f_erosion <- adonis2(BCdist_f ~ f_sample_data$Erosion_grouped)
adonis_f_erosion$variable <-"Erosion risk"
adonis_f_erosion$group <- "Fungi"
adonis_f_prec <- adonis2(BCdist_f ~ f_sample_data$cluster_prec)
adonis_f_prec$variable <-"Precipitation"
adonis_f_prec$group <- "Fungi"
adonis_f_temp <- adonis2(BCdist_f ~ f_sample_data$cluster_temp)
adonis_f_temp$variable <-"Temperature"
adonis_f_temp$group <- "Fungi"
adonis_f_temp_range <- adonis2(BCdist_f ~ f_sample_data$cluster_temp_range)
adonis_f_temp_range$variable <-"Temperature range"
adonis_f_temp_range$group <- "Fungi"
adonis_f_Cmic <- adonis2(BCdist_f ~ f_sample_data$cluster_Cmic)
adonis_f_Cmic$variable <-"Cmic"
adonis_f_Cmic$group <- "Fungi"

#create adonis dataframe 
adonis_ALL <- rbind(adonis_f_LC1,adonis_f_LU2,adonis_f_pH,adonis_f_depth,adonis_f_erosion,adonis_f_prec,adonis_f_temp,adonis_f_temp_range,adonis_f_Cmic,
                    adonis_p_LC1,adonis_p_LU2,adonis_p_pH,adonis_p_depth,adonis_p_erosion,adonis_p_prec,adonis_p_temp,adonis_p_temp_range,adonis_p_Cmic,
                    adonis_roti_LC1,adonis_roti_LU2,adonis_roti_pH,adonis_roti_depth,adonis_roti_erosion,adonis_roti_prec,adonis_roti_temp,adonis_roti_temp_range,adonis_roti_Cmic,
                    adonis_tardi_LC1,adonis_tardi_LU2,adonis_tardi_pH,adonis_tardi_depth,adonis_tardi_erosion,adonis_tardi_prec,adonis_tardi_temp,adonis_tardi_temp_range,adonis_tardi_Cmic,
                    adonis_nema_LC1,adonis_nema_LU2, adonis_nema_pH,adonis_nema_depth,adonis_nema_erosion,adonis_nema_prec,adonis_nema_temp,adonis_nema_temp_range,adonis_nema_Cmic,
                    adonis_arthro_LC1,adonis_arthro_LU2,adonis_arthro_pH,adonis_arthro_depth,adonis_arthro_erosion,adonis_arthro_prec,adonis_arthro_temp,adonis_arthro_temp_range,adonis_arthro_Cmic,
                    adonis_anneli_LC1,adonis_anneli_LU2,adonis_anneli_pH,adonis_anneli_depth,adonis_anneli_erosion,adonis_anneli_prec,adonis_anneli_temp,adonis_anneli_temp_range,adonis_anneli_Cmic)
                  
                    
                      
write.csv(adonis_ALL, file = "adonis_1605.csv", row.names = TRUE)





