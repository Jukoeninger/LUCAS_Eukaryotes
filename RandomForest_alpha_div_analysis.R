#############################################
#Random forest observed richness, Figure 5, Figure S8
#Set graphical themes
theme_jk5 <- function(){
  theme_bw() +
    theme(text = element_text(size=12, family = "sans"),
          axis.text = element_text(size = 12,face="bold"),
          axis.text.x = element_text(angle=45,vjust=1,hjust=1),
          axis.title = element_text(size = 14, vjust = -3,face="bold"),
          plot.margin = unit(c(1, 1, 1, 1), units =, "cm"),
          plot.title = element_text(size = 12, vjust = 1, hjust = 0),
          legend.text = element_text(size = 12),
          legend.key.size = unit(1,"cm"),
          legend.key.width= unit(1, "cm"),
          legend.title = element_text(size=10),
          legend.background = element_rect(color = "black",
                                           fill = "transparent",
                                           size = 0.5, linetype = "blank"))
}

#import ecosystem history
ecosystem_history <- read.csv("~/Euk_18S_data/ecosystem_history.csv", header=1, row.names=1)
ecosystem_history <- ecosystem_history[,c(1,4:6)]

set.seed(123)

#rotifers
otutable_roti <-ranacapa::vegan_otu(lc_a_roti)
sampledf_roti <- data.frame(sample_data(lc_a_roti))
sampledf_roti <- merge(sampledf_roti,ecosystem_history,by="BARCODE_ID",all=F)
rownames(sampledf_roti) <-sampledf_roti[,"BARCODE_ID"]

#feature selection
sampledf_roti2_obs <- sampledf_roti[,c("Richness","Elevation","Aridity_in_sample_month","AvTemp2000_2018","AvPrec2000_2018","Erosion_risk","Slope", "Aspect","Sample_month","Temperature_in_sample_month","Aridity","bas","Cmic","qO2","Water_content","LC1_2018","Coarse_fragments","Clay","Sand","pH","Electrical_conductivity","P","K","LU2","C.N","Soil_depth","Bulk_density","Mean_annual_temperature","Mean_diurnal_range","Isothermality","Temperature_seasonality","Max_temperature_of_warmest.month","Min_temperature_of_coldest.month","Temperature_annual_range","Mean_temperature_of_the_wettest_quarter","Mean_temperature_of_driest_quarter","Mean_temperature_of_warmest_quarter","Mean_temperature_of_coldest_quarter","Total_annual_precipitation","Precipitation_of_wettest_month","Precipitation_of_driest_month","Precipitation_seasonality","Precipitation_of_wettest_quarter","Precipitation_of_driest_quarter","Precipitation_of_warmest_quarter","Precipitation_of_coldest_quarter","SSM","SSMA","LC1_LUCAS15","LC1_LUCAS09","LC1")]
colnames(sampledf_roti2_obs)[colnames(sampledf_roti2_obs) == "LC1_2018"] <- "Ecosystem_type"
colnames(sampledf_roti2_obs)[colnames(sampledf_roti2_obs) == "LU2"] <- "Plant_cover"
colnames(sampledf_roti2_obs)[colnames(sampledf_roti2_obs) == "SSM"] <- "Surface_moisture"
colnames(sampledf_roti2_obs)[colnames(sampledf_roti2_obs) == "LC1_LUCAS15"] <- "Ecosystem_type_2015"
colnames(sampledf_roti2_obs)[colnames(sampledf_roti2_obs) == "LC1_LUCAS09"] <- "Ecosystem_type_2009"
sampledf_roti5_obs <-na.exclude(sampledf_roti2_obs)#feature selection

sampledf_roti2_obs <- scale(sampledf_roti5_obs[,-c(49:51,24,16)]) #scale numerical variables
sampledf_roti2_obs <- as.data.frame(sampledf_roti2_obs)
sampledf_roti5_obs <- cbind(sampledf_roti2_obs,sampledf_roti5_obs[,c(49:51,24,16)])
train.index_roti = sample(seq_len(nrow(sampledf_roti5_obs)), size = floor(0.8 * nrow(sampledf_roti5_obs)))
train.obs_roti <- sampledf_roti5_obs[train.index_roti,]
test.obs_roti <- sampledf_roti5_obs[-train.index_roti,]

#rfe
ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   number=10,
                   verbose = FALSE)

rfe_roti <- rfe(x=train.obs_roti[,-(1)],y=train.obs_roti[,1],
                sizes = c(2:51),
                rfeControl = ctrl)

rfe_roti$optVariables
train.obs_roti <- train.obs_roti[,c("Richness",rfe_roti$optVariables)]


#tuning model using caret
myControl <- trainControl(
  method = "cv",
  number = 10,
  savePredictions=TRUE,
  classProbs = TRUE, # IMPORTANT!
  verboseIter = FALSE
)
x <- train.obs_roti[,c(!colnames(train.obs_roti) %in% "Richness")]
mtry <- sqrt(ncol(x))
tunegrid <- expand.grid(.mtry=mtry)
set.seed(123)

test_caret_roti <- train(Richness~., data=train.obs_roti, method="rf", metric="rmse",
                         tuneGrid=tunegrid, ntree = 1000, trControl=myControl)
test_caret_roti$bestTune <- as.numeric(test_caret_roti$bestTune)

#ranger
rf_roti_obs <- ranger(Richness ~., data=train.obs_roti, mtry=test_caret_roti$bestTune,num.trees=900, min.node.size=5,
                      quantreg = TRUE,keep.inbag = TRUE,importance='impurity',
                      oob.error=TRUE,verbose=TRUE) #16.9%
rf_roti_obs
saveRDS(rf_roti_obs, "/home/koenjul/rf_roti_obs.rds")
rf_roti_obs <-readRDS("/home/koenjul/rf_roti_obs.rds")
v<-as.vector(rf_roti_obs$variable.importance)
w<-as.vector(names(rf_roti_obs$variable.importance))
vw<-as.data.frame(cbind(v,w));names(vw);vw$v<-as.numeric(vw$v)
vw2_roti <- vw[with(vw,order(-v)),]
vw3_roti <- vw2_roti[1:50,]


#tardigrades
otutable_tardi <-ranacapa::vegan_otu(lc_a_tardi)
sampledf_tardi <- data.frame(sample_data(lc_a_tardi))
#merge with intensity data
sampledf_tardi <- merge(sampledf_tardi,ecosystem_history,by="BARCODE_ID",all=F)
rownames(sampledf_tardi) <-sampledf_tardi[,"BARCODE_ID"]

#feature selection
sampledf_tardi2_obs <- sampledf_tardi[,c("Richness","Elevation","Aridity_in_sample_month","AvTemp2000_2018","AvPrec2000_2018","Erosion_risk","Slope", "Aspect","Sample_month","Temperature_in_sample_month","Aridity","bas","Cmic","qO2","Water_content","LC1_2018","Coarse_fragments","Clay","Sand","pH","Electrical_conductivity","P","K","LU2","C.N","Soil_depth","Bulk_density","Mean_annual_temperature","Mean_diurnal_range","Isothermality","Temperature_seasonality","Max_temperature_of_warmest.month","Min_temperature_of_coldest.month","Temperature_annual_range","Mean_temperature_of_the_wettest_quarter","Mean_temperature_of_driest_quarter","Mean_temperature_of_warmest_quarter","Mean_temperature_of_coldest_quarter","Total_annual_precipitation","Precipitation_of_wettest_month","Precipitation_of_driest_month","Precipitation_seasonality","Precipitation_of_wettest_quarter","Precipitation_of_driest_quarter","Precipitation_of_warmest_quarter","Precipitation_of_coldest_quarter","SSM","SSMA","LC1_LUCAS15","LC1_LUCAS09","LC1")]
colnames(sampledf_tardi2_obs)[colnames(sampledf_tardi2_obs) == "LC1_2018"] <- "Ecosystem_type"
colnames(sampledf_tardi2_obs)[colnames(sampledf_tardi2_obs) == "LU2"] <- "Plant_cover"
colnames(sampledf_tardi2_obs)[colnames(sampledf_tardi2_obs) == "SSM"] <- "Surface_moisture"
colnames(sampledf_tardi2_obs)[colnames(sampledf_tardi2_obs) == "LC1_LUCAS15"] <- "Ecosystem_type_2015"
colnames(sampledf_tardi2_obs)[colnames(sampledf_tardi2_obs) == "LC1_LUCAS09"] <- "Ecosystem_type_2009"
sampledf_tardi5_obs <-na.exclude(sampledf_tardi2_obs)#feature selection

sampledf_tardi3_obs <- scale(sampledf_tardi5_obs[,-c(49:51,24,16)]) #scale numerical variables
sampledf_tardi3_obs <- as.data.frame(sampledf_tardi3_obs)
sampledf_tardi5_ob <- cbind(sampledf_tardi3_obs,sampledf_tardi5_obs[,c(49:51,24,16)])
train.index_tardi = sample(seq_len(nrow(sampledf_tardi5_obs)), size = floor(0.8 * nrow(sampledf_tardi5_obs)))
train.obs_tardi <- sampledf_tardi5_obs[train.index_tardi,]
test.obs_tardi <- sampledf_tardi5_obs[-train.index_tardi,]

#rfe
ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   number= 10,
                   verbose = FALSE)

rfe_tardi <- rfe(x=train.obs_tardi[,-(1)],y=train.obs_tardi[,1],
             sizes = c(2:51),
             rfeControl = ctrl)

rfe_tardi$optVariables #show selected variables - all variables selected
train.obs_tardi <- train.obs_tardi[,c("Richness",rfe_tardi$optVariables)]

#tuning model using caret
myControl <- trainControl(
  method = "cv",
  number = 10,
  savePredictions=TRUE,
  classProbs = TRUE, # IMPORTANT!
  verboseIter = FALSE
)
x <- train.obs_tardi[,c(!colnames(train.obs_tardi) %in% "Richness")]
mtry <- sqrt(ncol(x))
tunegrid <- expand.grid(.mtry=mtry)
set.seed(123)

test_caret_tardi <- train(Richness~., data=train.obs_tardi, method="rf", metric="rmse",
                          tuneGrid=tunegrid, ntree = 1000, trControl=myControl,na.action=na.pass)
test_caret_tardi$bestTune <- as.numeric(test_caret_tardi$bestTune)

#ranger
rf_tardi_obs <- ranger(Richness ~., data=train.obs_tardi, mtry=test_caret_tardi$bestTune,num.trees=900, min.node.size=5,
                       quantreg = TRUE,keep.inbag = TRUE,importance='impurity',
                       oob.error=TRUE,verbose=TRUE) #13.1%
rf_tardi_obs
saveRDS(rf_tardi_obs, "/home/koenjul/rf_tardi_obs.rds")
rf_tardi_obs <-readRDS("/home/koenjul/rf_tardi_obs.rds")
v<-as.vector(rf_tardi_obs$variable.importance)
w<-as.vector(names(rf_tardi_obs$variable.importance))
vw<-as.data.frame(cbind(v,w));names(vw);vw$v<-as.numeric(vw$v)
vw2_tardi <- vw[with(vw,order(-v)),]
vw3_tardi <- vw2_tardi[1:50,]



#nematodes
otutable_nema <-ranacapa::vegan_otu(lc_a_nema)
sampledf_nema <- data.frame(sample_data(lc_a_nema))
#merge with intensity data
sampledf_nema <- merge(sampledf_nema,ecosystem_history,by="BARCODE_ID",all=F)
rownames(sampledf_nema) <-sampledf_nema[,"BARCODE_ID"]

#feature selection
sampledf_nema2_obs <- sampledf_nema[,c("Richness","Elevation","Aridity_in_sample_month","AvTemp2000_2018","AvPrec2000_2018","Erosion_risk","Slope", "Aspect","Sample_month","Temperature_in_sample_month","Aridity","bas","Cmic","qO2","Water_content","LC1_2018","Coarse_fragments","Clay","Sand","pH","Electrical_conductivity","P","K","LU2","C.N","Soil_depth","Bulk_density","Mean_annual_temperature","Mean_diurnal_range","Isothermality","Temperature_seasonality","Max_temperature_of_warmest.month","Min_temperature_of_coldest.month","Temperature_annual_range","Mean_temperature_of_the_wettest_quarter","Mean_temperature_of_driest_quarter","Mean_temperature_of_warmest_quarter","Mean_temperature_of_coldest_quarter","Total_annual_precipitation","Precipitation_of_wettest_month","Precipitation_of_driest_month","Precipitation_seasonality","Precipitation_of_wettest_quarter","Precipitation_of_driest_quarter","Precipitation_of_warmest_quarter","Precipitation_of_coldest_quarter","SSM","SSMA","LC1_LUCAS15","LC1_LUCAS09","LC1")]
colnames(sampledf_nema2_obs)[colnames(sampledf_nema2_obs) == "LC1_2018"] <- "Ecosystem_type"
colnames(sampledf_nema2_obs)[colnames(sampledf_nema2_obs) == "LU2"] <- "Plant_cover"
colnames(sampledf_nema2_obs)[colnames(sampledf_nema2_obs) == "SSM"] <- "Surface_moisture"
colnames(sampledf_nema2_obs)[colnames(sampledf_nema2_obs) == "LC1_LUCAS15"] <- "Ecosystem_type_2015"
colnames(sampledf_nema2_obs)[colnames(sampledf_nema2_obs) == "LC1_LUCAS09"] <- "Ecosystem_type_2009"
sampledf_nema5_obs <-na.exclude(sampledf_nema2_obs)#feature selection

sampledf_nema3_obs <- scale(sampledf_nema5_obs[,-c(49:51,24,16)]) #scale numerical variables
sampledf_nema3_obs <- as.data.frame(sampledf_nema3_obs)
sampledf_nema5_ob <- cbind(sampledf_nema3_obs,sampledf_nema5_obs[,c(49:51,24,16)])
train.index_nema = sample(seq_len(nrow(sampledf_nema5_obs)), size = floor(0.8 * nrow(sampledf_nema5_obs)))
train.obs_nema <- sampledf_nema5_obs[train.index_nema,]
test.obs_nema <- sampledf_nema5_obs[-train.index_nema,]

#rfe
ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   number= 10,
                   verbose = FALSE)

rfe_nema <- rfe(x=train.obs_nema[,-(1)],y=train.obs_nema[,1],
                 sizes = c(2:51),
                 rfeControl = ctrl)

rfe_nema$optVariables #show selected variables - all variables selected
train.obs_nema <- train.obs_nema[,c("Richness",rfe_nema$optVariables)]

#tuning model using caret
myControl <- trainControl(
  method = "cv",
  number = 10,
  savePredictions=TRUE,
  classProbs = TRUE, # IMPORTANT!
  verboseIter = FALSE
)
x <- train.obs_nema[,c(!colnames(train.obs_nema) %in% "Richness")]
mtry <- sqrt(ncol(x))
tunegrid <- expand.grid(.mtry=mtry)
set.seed(123)

test_caret_nema <- train(Richness~., data=train.obs_nema, method="rf", metric="rmse",
                          tuneGrid=tunegrid, ntree = 1000, trControl=myControl,na.action=na.pass)
test_caret_nema$bestTune <- as.numeric(test_caret_nema$bestTune)

#ranger
rf_nema_obs <- ranger(Richness ~., data=train.obs_nema, mtry=test_caret_nema$bestTune,num.trees=900, min.node.size=5,
                       quantreg = TRUE,keep.inbag = TRUE,importance='impurity',
                       oob.error=TRUE,verbose=TRUE) #13.1%
rf_nema_obs
saveRDS(rf_nema_obs, "/home/koenjul/rf_nema_obs.rds")
rf_nema_obs <-readRDS("/home/koenjul/rf_nema_obs.rds")
v<-as.vector(rf_nema_obs$variable.importance)
w<-as.vector(names(rf_nema_obs$variable.importance))
vw<-as.data.frame(cbind(v,w));names(vw);vw$v<-as.numeric(vw$v)
vw2_nema <- vw[with(vw,order(-v)),]
vw3_nema <- vw2_nema[2:50,]




#arthropods
otutable_arthro <-ranacapa::vegan_otu(lc_a_arthro)
sampledf_arthro <- data.frame(sample_data(lc_a_arthro))
#merge with intensity data
sampledf_arthro <- merge(sampledf_arthro,ecosystem_history,by="BARCODE_ID",all=F)
rownames(sampledf_arthro) <-sampledf_arthro[,"BARCODE_ID"]

#feature selection
sampledf_arthro2_obs <- sampledf_arthro[,c("Richness","Elevation","Aridity_in_sample_month","AvTemp2000_2018","AvPrec2000_2018","Erosion_risk","Slope", "Aspect","Sample_month","Temperature_in_sample_month","Aridity","bas","Cmic","qO2","Water_content","LC1_2018","Coarse_fragments","Clay","Sand","pH","Electrical_conductivity","P","K","LU2","C.N","Soil_depth","Bulk_density","Mean_annual_temperature","Mean_diurnal_range","Isothermality","Temperature_seasonality","Max_temperature_of_warmest.month","Min_temperature_of_coldest.month","Temperature_annual_range","Mean_temperature_of_the_wettest_quarter","Mean_temperature_of_driest_quarter","Mean_temperature_of_warmest_quarter","Mean_temperature_of_coldest_quarter","Total_annual_precipitation","Precipitation_of_wettest_month","Precipitation_of_driest_month","Precipitation_seasonality","Precipitation_of_wettest_quarter","Precipitation_of_driest_quarter","Precipitation_of_warmest_quarter","Precipitation_of_coldest_quarter","SSM","SSMA","LC1_LUCAS15","LC1_LUCAS09","LC1")]
colnames(sampledf_arthro2_obs)[colnames(sampledf_arthro2_obs) == "LC1_2018"] <- "Ecosystem_type"
colnames(sampledf_arthro2_obs)[colnames(sampledf_arthro2_obs) == "LU2"] <- "Plant_cover"
colnames(sampledf_arthro2_obs)[colnames(sampledf_arthro2_obs) == "SSM"] <- "Surface_moisture"
colnames(sampledf_arthro2_obs)[colnames(sampledf_arthro2_obs) == "LC1_LUCAS15"] <- "Ecosystem_type_2015"
colnames(sampledf_arthro2_obs)[colnames(sampledf_arthro2_obs) == "LC1_LUCAS09"] <- "Ecosystem_type_2009"
sampledf_arthro5_obs <-na.exclude(sampledf_arthro2_obs)#feature selection

sampledf_arthro3_obs <- scale(sampledf_arthro5_obs[,-c(49:51,24,16)]) #scale numerical variables
sampledf_arthro3_obs <- as.data.frame(sampledf_arthro3_obs)
sampledf_arthro5_ob <- cbind(sampledf_arthro3_obs,sampledf_arthro5_obs[,c(49:51,24,16)])
train.index_arthro = sample(seq_len(nrow(sampledf_arthro5_obs)), size = floor(0.8 * nrow(sampledf_arthro5_obs)))
train.obs_arthro <- sampledf_arthro5_obs[train.index_arthro,]
test.obs_arthro <- sampledf_arthro5_obs[-train.index_arthro,]

#rfe
ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   number= 10,
                   verbose = FALSE)

rfe_arthro <- rfe(x=train.obs_arthro[,-(1)],y=train.obs_arthro[,1],
                 sizes = c(2:51),
                 rfeControl = ctrl)

rfe_arthro$optVariables #show selected variables - all variables selected
train.obs_arthro <- train.obs_arthro[,c("Richness",rfe_arthro$optVariables)]

#tuning model using caret
myControl <- trainControl(
  method = "cv",
  number = 10,
  savePredictions=TRUE,
  classProbs = TRUE, # IMPORTANT!
  verboseIter = FALSE
)
x <- train.obs_arthro[,c(!colnames(train.obs_arthro) %in% "Richness")]
mtry <- sqrt(ncol(x))
tunegrid <- expand.grid(.mtry=mtry)
set.seed(123)

test_caret_arthro <- train(Richness~., data=train.obs_arthro, method="rf", metric="rmse",
                          tuneGrid=tunegrid, ntree = 1000, trControl=myControl,na.action=na.pass)
test_caret_arthro$bestTune <- as.numeric(test_caret_arthro$bestTune)

#ranger
rf_arthro_obs <- ranger(Richness ~., data=train.obs_arthro, mtry=test_caret_arthro$bestTune,num.trees=900, min.node.size=5,
                       quantreg = TRUE,keep.inbag = TRUE,importance='impurity',
                       oob.error=TRUE,verbose=TRUE) #13.1%
rf_arthro_obs
saveRDS(rf_arthro_obs, "/home/koenjul/rf_arthro_obs.rds")
rf_arthro_obs <-readRDS("/home/koenjul/rf_arthro_obs.rds")
v<-as.vector(rf_arthro_obs$variable.importance)
w<-as.vector(names(rf_arthro_obs$variable.importance))
vw<-as.data.frame(cbind(v,w));names(vw);vw$v<-as.numeric(vw$v)
vw2_arthro <- vw[with(vw,order(-v)),]
vw3_arthro <- vw2_arthro[1:5,]



#annelids
otutable_anneli <-ranacapa::vegan_otu(lc_a_anneli)
sampledf_anneli <- data.frame(sample_data(lc_a_anneli))
sampledf_anneli <- merge(sampledf_anneli,ecosystem_history,by="BARCODE_ID",all=F)
rownames(sampledf_anneli) <-sampledf_anneli[,"BARCODE_ID"]
#Richness
sampledf_anneli2_obs <- sampledf_anneli[,c("Richness","Elevation","AvTemp2000_2018","AvPrec2000_2018","Erosion_risk","Slope", "Aspect","Sample_month","Precipitation_in_sample_month","Temperature_in_sample_month","Aridity","bas","Cmic","qO2","Water_content","LC1_2018","Coarse_fragments","Clay","Sand","pH","Electrical_conductivity","P","K","LU2","C.N","Soil_depth","Bulk_density","Mean_annual_temperature","Mean_diurnal_range","Isothermality","Temperature_seasonality","Max_temperature_of_warmest.month","Min_temperature_of_coldest.month","Temperature_annual_range","Mean_temperature_of_the_wettest_quarter","Mean_temperature_of_driest_quarter","Mean_temperature_of_warmest_quarter","Mean_temperature_of_coldest_quarter","Total_annual_precipitation","Precipitation_of_wettest_month","Precipitation_of_driest_month","Precipitation_seasonality","Precipitation_of_wettest_quarter","Precipitation_of_driest_quarter","Precipitation_of_warmest_quarter","Precipitation_of_coldest_quarter","SSM","SSMA","LC1_LUCAS15","LC1_LUCAS09","LC1")]
colnames(sampledf_anneli2_obs)[colnames(sampledf_anneli2_obs) == "LC1_2018"] <- "Ecosystem_type"
colnames(sampledf_anneli2_obs)[colnames(sampledf_anneli2_obs) == "LU2"] <- "Plant_cover"
colnames(sampledf_anneli2_obs)[colnames(sampledf_anneli2_obs) == "SSM"] <- "Surface_moisture"
colnames(sampledf_anneli2_obs)[colnames(sampledf_anneli2_obs) == "LC1_LUCAS15"] <- "Ecosystem_type_2015"
colnames(sampledf_anneli2_obs)[colnames(sampledf_anneli2_obs) == "LC1_LUCAS09"] <- "Ecosystem_type_2009"
sampledf_anneli5_obs <-na.exclude(sampledf_anneli2_obs)#feature selection
sampledf_anneli2_obs <- scale(sampledf_anneli5_obs[,-c(49:51,24,16)]) #scale numerical variables
sampledf_anneli2_obs <- as.data.frame(sampledf_anneli2_obs)
sampledf_anneli5_obs <-cbind(sampledf_anneli2_obs,sampledf_anneli5_obs[,c(49:51,24,16)])

train.index_anneli = sample(seq_len(nrow(sampledf_anneli5_obs)), size = floor(0.8 * nrow(sampledf_anneli5_obs))) 
train.obs_anneli <- sampledf_anneli5_obs[train.index_anneli,]
test.obs_anneli <- sampledf_anneli5_obs[-train.index_anneli,]

#rfe
ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   number=10,
                   verbose = FALSE)

rfe_anneli <- rfe(x=train.obs_anneli[,-(1)],y=train.obs_anneli[,1],
             sizes = c(2:51),
             rfeControl = ctrl)

rfe_anneli$optVariables #show selected variables
train.obs_anneli <- train.obs_anneli[,c("Richness",rfe_anneli$optVariables)]

#tuning model using caret
myControl <- trainControl(
  method = "cv",
  number = 10,
  savePredictions=TRUE,
  classProbs = TRUE, # IMPORTANT!
  verboseIter = FALSE
)
x <- train.obs_anneli[,c(!colnames(train.obs_anneli) %in% "Richness")]
mtry <- sqrt(ncol(x))
tunegrid <- expand.grid(.mtry=mtry)
set.seed(123)

test_caret_anneli <- train(Richness~., data=train.obs_anneli, method="rf", metric="rmse",
                          tuneGrid=tunegrid, ntree = 1000, trControl=myControl)
test_caret_anneli$bestTune <- as.numeric(test_caret_anneli$bestTune)

#ranger
rf_anneli_obs <- ranger(Richness ~., data=train.obs_anneli, mtry=test_caret_anneli$bestTune,num.trees=450, min.node.size=5,
                       quantreg = TRUE,keep.inbag = TRUE,importance='impurity',
                       oob.error=TRUE,verbose=TRUE) #4.4%
rf_anneli_obs
saveRDS(rf_anneli_obs, "/home/koenjul/rf_anneli_obs.rds")
rf_anneli_obs <-readRDS("/home/koenjul/rf_anneli_obs.rds")
v<-as.vector(rf_anneli_obs$variable.importance)
w<-as.vector(names(rf_anneli_obs$variable.importance))
vw<-as.data.frame(cbind(v,w));names(vw);vw$v<-as.numeric(vw$v)
vw2_anneli <- vw[with(vw,order(-v)),]
vw3_anneli <- vw2_anneli[1:37,]



#Protists
otutable_p<-ranacapa::vegan_otu(lc_p)
sampledf_p <- data.frame(sample_data(lc_p))
sampledf_p <- merge(sampledf_p,ecosystem_history,by="BARCODE_ID",all=F)
rownames(sampledf_p) <-sampledf_p[,"BARCODE_ID"]

colnames(sampledf_p)
sampledf_p2_obs <- sampledf_p[,c("Richness","Elevation","AvTemp2000_2018","AvPrec2000_2018","Erosion_risk","Slope", "Aspect","Sample_month","Precipitation_in_sample_month","Temperature_in_sample_month","Aridity","bas","Cmic","qO2","Water_content","LC1_2018","Coarse_fragments","Clay","Sand","pH","Electrical_conductivity","P","K","LU2","C.N","Soil_depth","Bulk_density","Mean_annual_temperature","Mean_diurnal_range","Isothermality","Temperature_seasonality","Max_temperature_of_warmest.month","Min_temperature_of_coldest.month","Temperature_annual_range","Mean_temperature_of_the_wettest_quarter","Mean_temperature_of_driest_quarter","Mean_temperature_of_warmest_quarter","Mean_temperature_of_coldest_quarter","Total_annual_precipitation","Precipitation_of_wettest_month","Precipitation_of_driest_month","Precipitation_seasonality","Precipitation_of_wettest_quarter","Precipitation_of_driest_quarter","Precipitation_of_warmest_quarter","Precipitation_of_coldest_quarter","SSM","SSMA","LC1_LUCAS15","LC1_LUCAS09","LC1")]
colnames(sampledf_p2_obs)[colnames(sampledf_p2_obs) == "LC1_2018"] <- "Ecosystem_type"
colnames(sampledf_p2_obs)[colnames(sampledf_p2_obs) == "LU2"] <- "Plant_cover"
colnames(sampledf_p2_obs)[colnames(sampledf_p2_obs) == "SSM"] <- "Surface_moisture"
colnames(sampledf_p2_obs)[colnames(sampledf_p2_obs) == "LC1_LUCAS15"] <- "Ecosystem_type_2015"
colnames(sampledf_p2_obs)[colnames(sampledf_p2_obs) == "LC1_LUCAS09"] <- "Ecosystem_type_2009"
sampledf_p5_obs <-na.exclude(sampledf_p2_obs)#feature selection
sampledf_p2_obs <- scale(sampledf_p5_obs[,-c(49:51,24,16)]) #scale numerical variables
sampledf_p2_obs <- as.data.frame(sampledf_p2_obs)
sampledf_p5_obs <-cbind(sampledf_p2_obs,sampledf_p5_obs[,c(49:51,24,16)])

train.index_p = sample(seq_len(nrow(sampledf_p5_obs)), size = floor(0.8 * nrow(sampledf_p5_obs))) 
train.obs_p <- sampledf_p5_obs[train.index_p,]
test.obs_p <- sampledf_p5_obs[-train.index_p,]

#rfe
ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   number=10,
                   verbose = FALSE)

rfe_p <- rfe(x=train.obs_p[,-(1)],y=train.obs_p[,1],
                 sizes = c(2:51),
                 rfeControl = ctrl)

rfe_p$optVariables #show selected variables - all variables selected
train.obs_p <- train.obs_p[,c("Richness",rfe_p$optVariables)]

#tuning model using caret
myControl <- trainControl(
  method = "cv",
  number = 10,
  savePredictions=TRUE,
  classProbs = TRUE, # IMPORTANT!
  verboseIter = FALSE
)
x <- train.obs_p[,c(!colnames(train.obs_p) %in% "Richness")]
mtry <- sqrt(ncol(x))
tunegrid <- expand.grid(.mtry=mtry)
set.seed(123)

test_caret_p <- train(Richness~., data=train.obs_p, method="rf", metric="rmse",
                      tuneGrid=tunegrid, ntree = 1000, trControl=myControl)
test_caret_p$bestTune <- as.numeric(test_caret_p$bestTune)

#ranger
rf_p_obs <- ranger(Richness ~., data=train.obs_p, mtry=test_caret_p$bestTune,num.trees=900, min.node.size=5,
                   quantreg = TRUE,keep.inbag = TRUE,importance='impurity',
                   oob.error=TRUE,verbose=TRUE) #16.2%
rf_p_obs
saveRDS(rf_p_obs, "/home/koenjul/rf_p_obs.rds")
rf_p_obs <-readRDS("/home/koenjul/rf_p_obs.rds")
v<-as.vector(rf_p_obs$variable.importance)
w<-as.vector(names(rf_p_obs$variable.importance))
vw<-as.data.frame(cbind(v,w));names(vw);vw$v<-as.numeric(vw$v)
vw2_p <- vw[with(vw,order(-v)),]
vw3_p <- vw2_p[1:50,]



#Fungi
otutable_f <-ranacapa::vegan_otu(lc_f)
sampledf_f <- data.frame(sample_data(lc_f))
sampledf_f <- merge(sampledf_f,ecosystem_history,by="BARCODE_ID",all=F)
rownames(sampledf_f) <-sampledf_f[,"BARCODE_ID"]
sampledf_f2_obs <- sampledf_f[,c("Richness","Elevation","AvTemp2000_2018","AvPrec2000_2018","Erosion_risk","Slope", "Aspect","Sample_month","Precipitation_in_sample_month","Temperature_in_sample_month","Aridity","bas","Cmic","qO2","Water_content","LC1_2018","Coarse_fragments","Clay","Sand","pH","Electrical_conductivity","P","K","LU2","C.N","Soil_depth","Bulk_density","Mean_annual_temperature","Mean_diurnal_range","Isothermality","Temperature_seasonality","Max_temperature_of_warmest.month","Min_temperature_of_coldest.month","Temperature_annual_range","Mean_temperature_of_the_wettest_quarter","Mean_temperature_of_driest_quarter","Mean_temperature_of_warmest_quarter","Mean_temperature_of_coldest_quarter","Total_annual_precipitation","Precipitation_of_wettest_month","Precipitation_of_driest_month","Precipitation_seasonality","Precipitation_of_wettest_quarter","Precipitation_of_driest_quarter","Precipitation_of_warmest_quarter","Precipitation_of_coldest_quarter","SSM","SSMA","LC1_LUCAS15","LC1_LUCAS09","LC1")]
colnames(sampledf_f2_obs)[colnames(sampledf_f2_obs) == "LC1_2018"] <- "Ecosystem_type"
colnames(sampledf_f2_obs)[colnames(sampledf_f2_obs) == "LU2"] <- "Plant_cover"
colnames(sampledf_f2_obs)[colnames(sampledf_f2_obs) == "SSM"] <- "Surface_moisture"
colnames(sampledf_f2_obs)[colnames(sampledf_f2_obs) == "LC1_LUCAS15"] <- "Ecosystem_type_2015"
colnames(sampledf_f2_obs)[colnames(sampledf_f2_obs) == "LC1_LUCAS09"] <- "Ecosystem_type_2009"
sampledf_f5_obs <-na.exclude(sampledf_f2_obs)#feature selection
sampledf_f2_obs <- scale(sampledf_f5_obs[,-c(49:51,24,16)]) #scale numerical variables
sampledf_f2_obs <- as.data.frame(sampledf_f2_obs)
sampledf_f5_obs <-cbind(sampledf_f2_obs,sampledf_f5_obs[,c(49:51,24,16)])

train.index_f = sample(seq_len(nrow(sampledf_f5_obs)), size = floor(0.8 * nrow(sampledf_f5_obs))) 
train.obs_f <- sampledf_f5_obs[train.index_f,]
test.obs_f <- sampledf_f5_obs[-train.index_f,]

#rfe
ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   number=10,
                   verbose = FALSE)

rfe_f <- rfe(x=train.obs_f[,-(1)],y=train.obs_f[,1],
             sizes = c(2:51),
             rfeControl = ctrl)

rfe_f$optVariables #show selected variables
train.obs_f <- train.obs_f[,c("Richness",rfe_f$optVariables)]

#tuning model using caret
myControl <- trainControl(
  method = "cv",
  number = 10,
  savePredictions=TRUE,
  classProbs = TRUE, # IMPORTANT!
  verboseIter = FALSE
)
x <- train.obs_f[,c(!colnames(train.obs_f) %in% "Richness")]
mtry <- sqrt(ncol(x))
tunegrid <- expand.grid(.mtry=mtry)
set.seed(123)

test_caret_f <- train(Richness~., data=train.obs_f, method="rf", metric="rmse",
                      tuneGrid=tunegrid, ntree = 1000, trControl=myControl)
test_caret_f$bestTune <- as.numeric(test_caret_f$bestTune)

#ranger
set.seed(123)
rf_f_obs <- ranger(Richness ~., data=train.obs_f, mtry=test_caret_f$bestTune,num.trees=900, min.node.size=115,
                   quantreg = TRUE,keep.inbag = TRUE,importance='impurity',
                   oob.error=TRUE,verbose=TRUE) #20.3%
rf_f_obs
saveRDS(rf_f_obs, "/home/koenjul/rf_f_obs.rds")
rf_f_obs <-readRDS("/home/koenjul/rf_f_obs.rds")
v<-as.vector(rf_f_obs$variable.importance)
w<-as.vector(names(rf_f_obs$variable.importance))
vw<-as.data.frame(cbind(v,w));names(vw);vw$v<-as.numeric(vw$v)
vw2_f <- vw[with(vw,order(-v)),]
vw3_f <- vw2_f[1:13,]



## Diagram showing var part explained alpha diversity
#Richness
first_column <- c("Variance")
second_column <- c(22.5)
third_column <- "Fungi"
df_f <- data.frame(first_column, second_column,third_column)

first_column <- c("Variance")
second_column <- c(15.4)
third_column <- "Protists"
df_p <- data.frame(first_column, second_column,third_column)

first_column <- c("Variance")
second_column <- c(2.9)
third_column <- "Rotifers"
df_roti <- data.frame(first_column, second_column,third_column)

first_column <- c("Variance")
second_column <- c(0.7)
third_column <- "Tardigrades"
df_tardi <- data.frame(first_column, second_column,third_column)

first_column <- c("Variance")
second_column <- c(12.2)
third_column <- "Nematodes"
df_nema <- data.frame(first_column, second_column,third_column)

first_column <- c("Variance")
second_column <- c(12.7)
third_column <- "Arthropods"
df_arthro <- data.frame(first_column, second_column,third_column)

first_column <- c("Variance")
second_column <- c(5.4)
third_column <- "Annelids"
df_anneli <- data.frame(first_column, second_column,third_column)




df2 <- rbind(df_f,df_p,df_roti,df_tardi,df_nema,df_arthro,df_anneli)
colnames(df2)[1] <- "Env_driver"
colnames(df2)[2] <- "Var_exp"
colnames(df2)[3] <- "group"

df2_long <- df2 %>% tidyr::gather(group, value, -Env_driver, -Var_exp,-group)
df2_long$group <- factor(df2_long$group, levels=c("Fungi","Protists","Rotifers","Tardigrades","Nematodes","Arthropods","Annelids"))

#PLOT
plot_df_obs <- ggplot(df2_long,aes(x=group,y=Var_exp)) 
plot_df_obs <- plot_df_obs + geom_bar(aes(fill=group),stat="identity") +
  ylab("\nVariation explained in %") +
  ggtitle("Observed richness")+
  ylim(0,28)+
  scale_fill_manual(values=c("#708238","#E47250","#5A4A6F","#9D5A6C","#A7C7E7","#FF99CC","orange")) + 
  theme_jk5()+theme(axis.text.x = element_text(size=13,angle=90, vjust = 0.5, hjust=1),
                    axis.title.x = element_blank(),
                    axis.ticks = element_blank(),
                    axis.text.y = element_text(size = 14),
                    legend.text = element_text(size = 14),
                    legend.title = element_blank(),
                    plot.title = element_text(size = 14, vjust = 1, hjust = 0))


##############################################################################
###########################################################
#Plot variable importance 

#add column with eukaryotic group
vw3_roti$group <- "Rotifers"
vw3_anneli$group <- "Annelids"
vw3_tardi$group <- "Tardigrades"
vw3_p$group <- "Protists"
vw3_f$group <- "Fungi"
vw3_arthro$group <- "Arthropods"
vw3_nema$group <- "Nematodes"

#create dataframe with all eukaryotic organisms' richness
obs_imp_var_long <- rbind(vw3_f,vw3_p,vw3_roti,vw3_tardi,vw3_nema,vw3_arthro,vw3_anneli)
obs_imp_var_long <- as.data.frame(obs_imp_var_long)
#rename variables
obs_imp_var_long$w[obs_imp_var_long$w =="AvTemp2000_2018"] <- "Mean temperature 2000-2018"
obs_imp_var_long$w[obs_imp_var_long$w =="AvPrec2000_2018"] <- "Mean precipitation 2000-2018"
obs_imp_var_long$w[obs_imp_var_long$w =="Sample_month"] <- "Sample month"
obs_imp_var_long$w[obs_imp_var_long$w =="Precipitation_in_sample_month"] <- "Precipitation in sample month"
obs_imp_var_long$w[obs_imp_var_long$w =="Temperature_in_sample_month"] <- "Temperature in sample month"
obs_imp_var_long$w[obs_imp_var_long$w =="Aridity_in_sample_month"] <- "Aridity in sample month"
obs_imp_var_long$w[obs_imp_var_long$w =="bas"] <- "Basal respiration"
obs_imp_var_long$w[obs_imp_var_long$w =="Cmic"] <- "tardibial biomass carbon"
obs_imp_var_long$w[obs_imp_var_long$w =="qO2"] <- "Respiration quotient"
obs_imp_var_long$w[obs_imp_var_long$w =="Water_content"] <- "Soil water content"
obs_imp_var_long$w[obs_imp_var_long$w =="Coarse_fragments"] <- "Coarse fragments"
obs_imp_var_long$w[obs_imp_var_long$w =="Electrical_conductivity"] <- "Electrical conductivity"
obs_imp_var_long$w[obs_imp_var_long$w =="P"] <- "Phosphorus"
obs_imp_var_long$w[obs_imp_var_long$w =="K"] <- "Potassium"
obs_imp_var_long$w[obs_imp_var_long$w =="C.N"] <- "C:N ratio"
obs_imp_var_long$w[obs_imp_var_long$w =="Erosion_risk"] <- "Erosion risk"
obs_imp_var_long$w[obs_imp_var_long$w =="Soil_depth"] <- "Soil depth"
obs_imp_var_long$w[obs_imp_var_long$w =="Bulk_density"] <- "Bulk density"
obs_imp_var_long$w[obs_imp_var_long$w =="Isothermality"] <- "Isothermality 1970-2000"
obs_imp_var_long$w[obs_imp_var_long$w =="Mean_annual_temperature"] <- "Mean annual temperature 1970-2000"
obs_imp_var_long$w[obs_imp_var_long$w =="Mean_diurnal_range"] <- "Mean diurnal range 1970-2000"
obs_imp_var_long$w[obs_imp_var_long$w =="Temperature_seasonality"] <- "Temperature seasonality 1970-2000"
obs_imp_var_long$w[obs_imp_var_long$w =="Max_temperature_of_warmest.month"] <- "Max temperature of warmest month 1970-2000"
obs_imp_var_long$w[obs_imp_var_long$w =="Min_temperature_of_coldest.month"] <- "Min temperature of coldest month 1970-2000"
obs_imp_var_long$w[obs_imp_var_long$w =="Temperature_annual_range"] <- "Temperature annual range 1970-2000"
obs_imp_var_long$w[obs_imp_var_long$w =="Mean_temperature_of_the_wettest_quarter"] <- "Mean temperature of wettest quarter 1970-2000"
obs_imp_var_long$w[obs_imp_var_long$w =="Mean_temperature_of_driest_quarter"] <- "Mean temperature of driest quarter 1970-2000"
obs_imp_var_long$w[obs_imp_var_long$w =="Mean_temperature_of_warmest_quarter"] <- "Mean temperature of warmest quarter 1970-2000"
obs_imp_var_long$w[obs_imp_var_long$w =="Mean_temperature_of_coldest_quarter"] <- "Mean temperature of coldest quarter 1970-2000"
obs_imp_var_long$w[obs_imp_var_long$w =="Total_annual_precipitation"] <- "Total annual precipitation 1970-2000"
obs_imp_var_long$w[obs_imp_var_long$w =="Precipitation_of_wettest_month"] <- "Precipitation of wettest month 1970-2000"
obs_imp_var_long$w[obs_imp_var_long$w =="Precipitation_of_driest_month"] <- "Precipitation of driest month 1970-2000"
obs_imp_var_long$w[obs_imp_var_long$w =="Precipitation_seasonality"] <- "Precipitation seasonality 1970-2000"
obs_imp_var_long$w[obs_imp_var_long$w =="Precipitation_of_wettest_quarter"] <- "Precipitation of wettest quarter 1970-2000"
obs_imp_var_long$w[obs_imp_var_long$w =="Precipitation_of_driest_quarter"] <- "Precipitation of driest quarter 1970-2000"
obs_imp_var_long$w[obs_imp_var_long$w =="Precipitation_of_warmest_quarter"] <- "Precipitation of warmest quarter 1970-2000"
obs_imp_var_long$w[obs_imp_var_long$w =="Precipitation_of_coldest_quarter"] <- "Precipitation of coldest quarter 1970-2000"
obs_imp_var_long$w[obs_imp_var_long$w =="SSMA"] <- "Surface moisture abnormality"
obs_imp_var_long$w[obs_imp_var_long$w =="Aridity"] <- "Aridity 1970-2000"
obs_imp_var_long$w[obs_imp_var_long$w =="Ecosystem_type"] <- "Ecosystem type 2018"
obs_imp_var_long$w[obs_imp_var_long$w =="Plant_cover"] <- "Plant cover"
obs_imp_var_long$w[obs_imp_var_long$w =="Surface_moisture"] <- "Surface moisture"
obs_imp_var_long$w[obs_imp_var_long$w =="Ecosystem_type_2015"] <- "Ecosystem type 2015"
obs_imp_var_long$w[obs_imp_var_long$w =="Ecosystem_type_2009"] <- "Ecosystem type 2009"
obs_imp_var_long$w[obs_imp_var_long$w =="LC1"] <- "Intensity gradient 2009-2018"

#set order of variables
obs_imp_var_long$group <- factor(obs_imp_var_long$group, levels=c("Fungi","Protists","Rotifers","Tardigrades","Nematodes","Arthropods","Annelids"))
obs_imp_var_long$w <- factor(obs_imp_var_long$w, levels=c("Ecosystem type 2018","Ecosystem type 2015","Ecosystem type 2009","Intensity gradient 2009-2018","Plant cover","Aspect","Elevation","Slope","Basal respiration","tardibial biomass carbon","Respiration quotient","Sample month","Bulk density","Carbonates","Clay", "C:N ratio", "Coarse fragments","Electrical conductivity","Erosion risk", "Potassium","Phosphorus", "pH","Sand","Soil depth","Surface moisture","Surface moisture abnormality","Soil water content","Aridity in sample month","Precipitation in sample month","Temperature in sample month","Mean temperature 2000-2018","Mean precipitation 2000-2018","Aridity 1970-2000","Isothermality 1970-2000","Max temperature of warmest month 1970-2000","Mean annual temperature 1970-2000","Mean diurnal range 1970-2000","Mean temperature of coldest quarter 1970-2000","Mean temperature of driest quarter 1970-2000","Mean temperature of wettest quarter 1970-2000","Mean temperature of warmest quarter 1970-2000","Min temperature of coldest month 1970-2000","Precipitation of coldest quarter 1970-2000","Precipitation of driest month 1970-2000","Precipitation of driest quarter 1970-2000","Precipitation of warmest quarter 1970-2000","Precipitation of wettest month 1970-2000","Precipitation of wettest quarter 1970-2000","Precipitation seasonality 1970-2000","Temperature annual range 1970-2000","Temperature seasonality 1970-2000","Total annual precipitation 1970-2000"))
obs_imp_var_long <- obs_imp_var_long %>% tidyr::gather(group, value, -v, -w,-group)

#setting colors
colH <- c("#387338","#387338","#387338","#387338","#387338",
          "#a0a7b0","#a0a7b0","#a0a7b0",
          "#b352a9","#b352a9","#b352a9",
          "#e5e83a",
          "#b38952","#b38952","#b38952","#b38952","#b38952","#b38952","#b38952","#b38952","#b38952","#b38952","#b38952","#b38952","#b38952","#b38952",
          "#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe") 

#PLOT FIGURE 4
plot_obs_imp_var_long <- ggplot(obs_imp_var_long,aes(x=v,y=w))
svg("1605_Fig4.svg",width=10,height=8)
plot_obs_imp_var_long + geom_bar(aes(fill=w),stat="identity") + 
  facet_wrap(~group,nrow=1, scales = "free_x")+
  ylab("\nVariation explained") +
  ggtitle("Observed ASV richness")+
  xlab("\nVariable importance") +
  scale_fill_manual(values=c(colH)) +  
  theme_jk5()+theme(legend.position="none",
                    axis.text.y = element_text(size = 10),
                    axis.text.x = element_text(size=10,angle=90, vjust = 0.5, hjust=1),
                    strip.text = element_text(face = "bold",size=9))
dev.off() # Close the graphics device











#############################################
#Random forest Shannon index, Figure S8, S9
#Set graphical themes
theme_jk5 <- function(){
  theme_bw() +
    theme(text = element_text(size=12, family = "sans"),
          axis.text = element_text(size = 12,face="bold"),
          axis.text.x = element_text(angle=45,vjust=1,hjust=1),
          axis.title = element_text(size = 14, vjust = -3,face="bold"),
          plot.margin = unit(c(1, 1, 1, 1), units =, "cm"),
          plot.title = element_text(size = 12, vjust = 1, hjust = 0),
          legend.text = element_text(size = 12),
          legend.key.size = unit(1,"cm"),
          legend.key.width= unit(1, "cm"),
          legend.title = element_text(size=10),
          legend.background = element_rect(color = "black",
                                           fill = "transparent",
                                           size = 0.5, linetype = "blank"))
}

#import ecosystem history
ecosystem_history <- read.csv("~/Downloads/ecosystem_history1605.csv", header=1, row.names=1)
ecosystem_history <- ecosystem_history[,c(1,4,5,6)]

#Random Forest 
set.seed(123)

#Rotifers
otutable_roti <- ranacapa::vegan_otu(lc_a_roti)
sampledf_roti <- data.frame(sample_data(lc_a_roti))


#merge with ecosystem history
sampledf_roti <- merge(sampledf_roti,ecosystem_history,by="BARCODE_ID",all=F)
rownames(sampledf_roti) <-sampledf_roti[,"BARCODE_ID"]

#Shannon index
#feature selection
sampledf_roti2_H <- sampledf_roti[,c("Shannon","Elevation","Aridity","Aridity_in_sample_month","AvTemp2000_2018","AvPrec2000_2018","Erosion_risk","Slope", "Aspect","Sample_month","Temperature_in_sample_month","bas","Cmic","qO2","Water_content","LC1_2018","Coarse_fragments","Clay","Sand","pH","Electrical_conductivity","P","K","LU2","C.N","Soil_depth","Bulk_density","Mean_annual_temperature","Mean_diurnal_range","Isothermality","Temperature_seasonality","Max_temperature_of_warmest.month","Min_temperature_of_coldest.month","Temperature_annual_range","Mean_temperature_of_the_wettest_quarter","Mean_temperature_of_driest_quarter","Mean_temperature_of_warmest_quarter","Mean_temperature_of_coldest_quarter","Total_annual_precipitation","Precipitation_of_wettest_month","Precipitation_of_driest_month","Precipitation_seasonality","Precipitation_of_wettest_quarter","Precipitation_of_driest_quarter","Precipitation_of_warmest_quarter","Precipitation_of_coldest_quarter","SSM","SSMA","LC1_LUCAS15","LC1_LUCAS09","LC1")]
colnames(sampledf_roti2_H)[colnames(sampledf_roti2_H) == "LC1_2018"] <- "Ecosystem_type"
colnames(sampledf_roti2_H)[colnames(sampledf_roti2_H) == "LU2"] <- "Plant_cover"
colnames(sampledf_roti2_H)[colnames(sampledf_roti2_H) == "SSM"] <- "Surface_moisture"
colnames(sampledf_roti2_H)[colnames(sampledf_roti2_H) == "LC1_LUCAS15"] <- "Ecosystem_type_2015"
colnames(sampledf_roti2_H)[colnames(sampledf_roti2_H) == "LC1_LUCAS09"] <- "Ecosystem_type_2009"
sampledf_roti5_H <-na.exclude(sampledf_roti2_H)#feature selection
sampledf_roti3_H <- scale(sampledf_roti5_H[,-c(49:51,24,16)]) #scale numerical variables
sampledf_roti3_H <- as.data.frame(sampledf_roti3_H)
sampledf_roti5_H <- cbind(sampledf_roti3_H,sampledf_roti5_H[,c(49:51,24,16)])
train.index_roti = sample(seq_len(nrow(sampledf_roti5_H)), size = floor(0.8 * nrow(sampledf_roti5_H)))
train.H_roti <- sampledf_roti5_H[train.index_roti,]
test.H_roti <- sampledf_roti5_H[-train.index_roti,]

#rfe
ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   number= 10,
                   verbose = FALSE)
rfe_roti_H <- rfe(x=train.H_roti[,-(1)],y=train.H_roti[,1],
                  sizes = c(2:51),
                  rfeControl = ctrl)
rfe_roti_H$optVariables #show selected variables
train.H_roti <- train.H_roti[,c("Shannon",rfe_roti_H$optVariables)]


#tuning model using caret
myControl <- trainControl(
  method = "cv",
  number = 10,
  savePredictions=TRUE,
  classProbs = TRUE, # IMPORTANT!
  verboseIter = FALSE
)

x <- train.H_roti[,c(!colnames(train.H_roti) %in% "Shannon")]
mtry <- sqrt(ncol(x))
tunegrid <- expand.grid(.mtry=mtry)
set.seed(123)

test_caret_roti <- train(Shannon~., data=train.H_roti, method="rf", metric="rmse",
                         tuneGrid=tunegrid, ntree = 1000, trControl=myControl,na.action=na.pass)
test_caret_roti$bestTune <- as.numeric(test_caret_roti$bestTune)

#ranger
rf_roti_H <- ranger(Shannon ~., data=train.H_roti, mtry=test_caret_roti$bestTune,num.trees=900, min.node.size=250,
                    quantreg = TRUE,keep.inbag = TRUE,importance='impurity',
                    oob.error=TRUE,verbose=TRUE) #2.2%
rf_roti_H
saveRDS(rf_roti_H, "/home/koenjul/rf_roti_H.rds")
rf_roti_H <-readRDS("/home/koenjul/rf_roti_H.rds")
v<-as.vector(rf_roti_H$variable.importance)
w<-as.vector(names(rf_roti_H$variable.importance))
vw<-as.data.frame(cbind(v,w));names(vw);vw$v<-as.numeric(vw$v)
vw2_roti <- vw[with(vw,order(-v)),]
vw3_roti_H <- vw2_roti[1:43,]

#Tardigrades
otutable_tardi <-ranacapa::vegan_otu(lc_a_tardi)
sampledf_tardi <- data.frame(sample_data(lc_a_tardi))
sampledf_tardi <- merge(sampledf_tardi,ecosystem_history,by="BARCODE_ID",all=F)
rownames(sampledf_tardi) <-sampledf_tardi[,"BARCODE_ID"]

#Shannon index
#feature selection
sampledf_tardi2_H <- sampledf_tardi[,c("Shannon","Aridity","Elevation","Aridity_in_sample_month","AvTemp2000_2018","AvPrec2000_2018","Erosion_risk","Slope", "Aspect","Sample_month","Temperature_in_sample_month","bas","Cmic","qO2","Water_content","LC1_2018","Coarse_fragments","Clay","Sand","pH","Electrical_conductivity","P","K","LU2","C.N","Soil_depth","Bulk_density","Mean_annual_temperature","Mean_diurnal_range","Isothermality","Temperature_seasonality","Max_temperature_of_warmest.month","Min_temperature_of_coldest.month","Temperature_annual_range","Mean_temperature_of_the_wettest_quarter","Mean_temperature_of_driest_quarter","Mean_temperature_of_warmest_quarter","Mean_temperature_of_coldest_quarter","Total_annual_precipitation","Precipitation_of_wettest_month","Precipitation_of_driest_month","Precipitation_seasonality","Precipitation_of_wettest_quarter","Precipitation_of_driest_quarter","Precipitation_of_warmest_quarter","Precipitation_of_coldest_quarter","SSM","SSMA","LC1_LUCAS15","LC1_LUCAS09","LC1")]
colnames(sampledf_tardi2_H)[colnames(sampledf_tardi2_H) == "LC1_2018"] <- "Ecosystem_type"
colnames(sampledf_tardi2_H)[colnames(sampledf_tardi2_H) == "LU2"] <- "Plant_cover"
colnames(sampledf_tardi2_H)[colnames(sampledf_tardi2_H) == "SSM"] <- "Surface_moisture"
colnames(sampledf_tardi2_H)[colnames(sampledf_tardi2_H) == "LC1_LUCAS15"] <- "Ecosystem_type_2015"
colnames(sampledf_tardi2_H)[colnames(sampledf_tardi2_H) == "LC1_LUCAS09"] <- "Ecosystem_type_2009"
sampledf_tardi5_H <-na.exclude(sampledf_tardi2_H)#feature selection

sampledf_tardi2_H <- scale(sampledf_tardi5_H[,-c(49:51,24,16)]) #scale numerical variables
sampledf_tardi2_H <- as.data.frame(sampledf_tardi2_H)
sampledf_tardi5_H <- cbind(sampledf_tardi2_H,sampledf_tardi5_H[,c(49:51,24,16)])
train.index_tardi = sample(seq_len(nrow(sampledf_tardi5_H)), size = floor(0.8 * nrow(sampledf_tardi5_H))) 
train.H_tardi <- sampledf_tardi5_H[train.index_tardi,]
test.H_tardi <- sampledf_tardi5_H[-train.index_tardi,]

#rfe
ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   number=10,
                   verbose = FALSE)

rfe_tardi_H <- rfe(x=train.H_tardi[,-(1)],y=train.H_tardi[,1],
                   sizes = c(2:51),
                   rfeControl = ctrl)

rfe_tardi_H$optVariables
train.H_tardi <- train.H_tardi[,c("Shannon",rfe_tardi_H$optVariables)]

#tuning model using caret
myControl <- trainControl(
  method = "cv",
  number = 10,
  savePredictions=TRUE,
  classProbs = TRUE, # IMPORTANT!
  verboseIter = FALSE
)

x <- train.H_tardi[,c(!colnames(train.H_tardi) %in% "Shannon")]
mtry <- sqrt(ncol(x))
tunegrid <- expand.grid(.mtry=mtry)
set.seed(123)

test_caret_tardi <- train(Shannon~., data=train.H_tardi, method="rf", metric="rmse",
                          tuneGrid=tunegrid, ntree = 1000, trControl=myControl)
test_caret_tardi$bestTune <- as.numeric(test_caret_tardi$bestTune)

#ranger
rf_tardi_H <- ranger(Shannon ~., data=train.H_tardi, mtry=test_caret_tardi$bestTune,num.trees=900, min.node.size=5,
                     quantreg = TRUE,keep.inbag = TRUE,importance='impurity',
                     oob.error=TRUE,verbose=TRUE) #7.9%
rf_tardi_H
saveRDS(rf_tardi_H, "/home/koenjul/rf_roti_H.rds")
rf_tardi_H <-readRDS("/home/koenjul/rf_roti_H.rds")
v<-as.vector(rf_tardi_H$variable.importance)
w<-as.vector(names(rf_tardi_H$variable.importance))
vw<-as.data.frame(cbind(v,w));names(vw);vw$v<-as.numeric(vw$v)
vw2_tardi <- vw[with(vw,order(-v)),]
vw3_tardi_H <- vw2_tardi[1:50,]



#Nematodes
otutable_nema <-ranacapa::vegan_otu(lc_a_nema)
sampledf_nema <- data.frame(sample_data(lc_a_nema))
sampledf_nema <- merge(sampledf_nema,ecosystem_history,by="BARCODE_ID",all=F)
rownames(sampledf_nema) <-sampledf_nema[,"BARCODE_ID"]

sampledf_nema2_H <- sampledf_nema[,c("Shannon","Aridity","Elevation","Aridity_in_sample_month","AvTemp2000_2018","AvPrec2000_2018","Erosion_risk","Slope", "Aspect","Sample_month","Temperature_in_sample_month","bas","Cmic","qO2","Water_content","LC1_2018","Coarse_fragments","Clay","Sand","pH","Electrical_conductivity","P","K","LU2","C.N","Soil_depth","Bulk_density","Mean_annual_temperature","Mean_diurnal_range","Isothermality","Temperature_seasonality","Max_temperature_of_warmest.month","Min_temperature_of_coldest.month","Temperature_annual_range","Mean_temperature_of_the_wettest_quarter","Mean_temperature_of_driest_quarter","Mean_temperature_of_warmest_quarter","Mean_temperature_of_coldest_quarter","Total_annual_precipitation","Precipitation_of_wettest_month","Precipitation_of_driest_month","Precipitation_seasonality","Precipitation_of_wettest_quarter","Precipitation_of_driest_quarter","Precipitation_of_warmest_quarter","Precipitation_of_coldest_quarter","SSM","SSMA","LC1_LUCAS15","LC1_LUCAS09","LC1")]
colnames(sampledf_nema2_H)[colnames(sampledf_nema2_H) == "LC1_2018"] <- "Ecosystem_type"
colnames(sampledf_nema2_H)[colnames(sampledf_nema2_H) == "LU2"] <- "Plant_cover"
colnames(sampledf_nema2_H)[colnames(sampledf_nema2_H) == "SSM"] <- "Surface_moisture"
colnames(sampledf_nema2_H)[colnames(sampledf_nema2_H) == "LC1_LUCAS15"] <- "Ecosystem_type_2015"
colnames(sampledf_nema2_H)[colnames(sampledf_nema2_H) == "LC1_LUCAS09"] <- "Ecosystem_type_2009"
sampledf_nema5_H <-na.exclude(sampledf_nema2_H)#feature selection
sampledf_nema2_H <- scale(sampledf_nema5_H[,-c(49:51,24,16)]) #scale numerical variables
sampledf_nema2_H <- as.data.frame(sampledf_nema2_H)
sampledf_nema5_H <-cbind(sampledf_nema2_H,sampledf_nema5_H[,c(49:51,24,16)])

train.index_nema = sample(seq_len(nrow(sampledf_nema5_H)), size = floor(0.8 * nrow(sampledf_nema5_H)))
train.H_nema <- sampledf_nema5_H[train.index_nema,]
test.H_nema <- sampledf_nema5_H[-train.index_nema,]

#rfe
ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   number=10,
                   verbose = FALSE)

rfe_nema_H <- rfe(x=train.H_nema[,-(1)],y=train.H_nema[,1],
                  sizes = c(2:50),
                  rfeControl = ctrl)

rfe_nema_H$optVariables
train.H_nema <- train.H_nema[,c("Shannon",rfe_nema_H$optVariables)]

#tuning model using caret
myControl <- trainControl(
  method = "cv",
  number = 10,
  savePredictions=TRUE,
  classProbs = TRUE, # IMPORTANT!
  verboseIter = FALSE
)
x <- train.H_nema[,c(!colnames(train.H_nema) %in% "Shannon")]
mtry <- sqrt(ncol(x))
tunegrid <- expand.grid(.mtry=mtry)
set.seed(123)

test_caret_nema <- train(Shannon~., data=train.H_nema, method="rf", metric="rmse",
                         tuneGrid=tunegrid, ntree = 1000, trControl=myControl)
test_caret_nema$bestTune <- as.numeric(test_caret_nema$bestTune)

#ranger
rf_nema_H <- ranger(Shannon ~., data=train.H_nema, mtry=test_caret_nema$bestTune,num.trees=450, min.node.size=155,
                    quantreg = TRUE,keep.inbag = TRUE,importance='impurity',
                    oob.error=TRUE,verbose=TRUE) #2.2%
rf_nema_H
saveRDS(rf_nema_H, "/home/koenjul/rf_nema_H.rds")
rf_nema_H <-readRDS("/home/koenjul/rf_nema_H.rds")
v<-as.vector(rf_nema_H$variable.importance)
w<-as.vector(names(rf_nema_H$variable.importance))
vw<-as.data.frame(cbind(v,w));names(vw);vw$v<-as.numeric(vw$v)
vw2_nema <- vw[with(vw,order(-v)),]
vw3_nema_H <- vw2_nema[1:48,]


#arthropods
otutable_arthro <-ranacapa::vegan_otu(lc_a_arthro)
sampledf_arthro <- data.frame(sample_data(lc_a_arthro))
sampledf_arthro <- merge(sampledf_arthro,ecosystem_history,by="BARCODE_ID",all=F)
rownames(sampledf_arthro) <-sampledf_arthro[,"BARCODE_ID"]

sampledf_arthro2_H <- sampledf_arthro[,c("Shannon","Aridity","Elevation","Aridity_in_sample_month","AvTemp2000_2018","AvPrec2000_2018","Erosion_risk","Slope", "Aspect","Sample_month","Temperature_in_sample_month","bas","Cmic","qO2","Water_content","LC1_2018","Coarse_fragments","Clay","Sand","pH","Electrical_conductivity","P","K","LU2","C.N","Soil_depth","Bulk_density","Mean_annual_temperature","Mean_diurnal_range","Isothermality","Temperature_seasonality","Max_temperature_of_warmest.month","Min_temperature_of_coldest.month","Temperature_annual_range","Mean_temperature_of_the_wettest_quarter","Mean_temperature_of_driest_quarter","Mean_temperature_of_warmest_quarter","Mean_temperature_of_coldest_quarter","Total_annual_precipitation","Precipitation_of_wettest_month","Precipitation_of_driest_month","Precipitation_seasonality","Precipitation_of_wettest_quarter","Precipitation_of_driest_quarter","Precipitation_of_warmest_quarter","Precipitation_of_coldest_quarter","SSM","SSMA","LC1_LUCAS15","LC1_LUCAS09","LC1")]
colnames(sampledf_arthro2_H)[colnames(sampledf_arthro2_H) == "LC1_2018"] <- "Ecosystem_type"
colnames(sampledf_arthro2_H)[colnames(sampledf_arthro2_H) == "LU2"] <- "Plant_cover"
colnames(sampledf_arthro2_H)[colnames(sampledf_arthro2_H) == "SSM"] <- "Surface_moisture"
colnames(sampledf_arthro2_H)[colnames(sampledf_arthro2_H) == "LC1_LUCAS15"] <- "Ecosystem_type_2015"
colnames(sampledf_arthro2_H)[colnames(sampledf_arthro2_H) == "LC1_LUCAS09"] <- "Ecosystem_type_2009"
sampledf_arthro5_H <-na.exclude(sampledf_arthro2_H)#feature selection
sampledf_arthro2_H <- scale(sampledf_arthro5_H[,-c(49:51,24,16)]) #scale numerical variables
sampledf_arthro2_H <- as.data.frame(sampledf_arthro2_H)
sampledf_arthro5_H <-cbind(sampledf_arthro2_H,sampledf_arthro5_H[,c(49:51,24,16)])

train.index_arthro = sample(seq_len(nrow(sampledf_arthro5_H)), size = floor(0.8 * nrow(sampledf_arthro5_H)))
train.H_arthro <- sampledf_arthro5_H[train.index_arthro,]
test.H_arthro <- sampledf_arthro5_H[-train.index_arthro,]

#rfe
ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   number=10,
                   verbose = FALSE)

rfe_arthro_H <- rfe(x=train.H_arthro[,-(1)],y=train.H_arthro[,1],
                    sizes = c(2:50),
                    rfeControl = ctrl)

rfe_arthro_H$optVariables
train.H_arthro <- train.H_arthro[,c("Shannon",rfe_arthro_H$optVariables)]

#tuning model using caret
myControl <- trainControl(
  method = "cv",
  number = 10,
  savePredictions=TRUE,
  classProbs = TRUE, # IMPORTANT!
  verboseIter = FALSE
)
x <- train.H_arthro[,c(!colnames(train.H_arthro) %in% "Shannon")]
mtry <- sqrt(ncol(x))
tunegrid <- expand.grid(.mtry=mtry)
set.seed(123)

test_caret_arthro <- train(Shannon~., data=train.H_arthro, method="rf", metric="rmse",
                           tuneGrid=tunegrid, ntree = 1000, trControl=myControl)
test_caret_arthro$bestTune <- as.numeric(test_caret_arthro$bestTune)

#ranger
rf_arthro_H <- ranger(Shannon ~., data=train.H_arthro, mtry=test_caret_arthro$bestTune,num.trees=450, min.node.size=155,
                      quantreg = TRUE,keep.inbag = TRUE,importance='impurity',
                      oob.error=TRUE,verbose=TRUE) #2.2%
rf_arthro_H
saveRDS(rf_arthro_H, "/home/koenjul/rf_arthro_H.rds")
rf_arthro_H <-readRDS("/home/koenjul/rf_arthro_H.rds")
v<-as.vector(rf_arthro_H$variable.importance)
w<-as.vector(names(rf_arthro_H$variable.importance))
vw<-as.data.frame(cbind(v,w));names(vw);vw$v<-as.numeric(vw$v)
vw2_arthro <- vw[with(vw,order(-v)),]
vw3_arthro_H <- vw2_arthro[1:47,]



#Annelids
otutable_anneli <-ranacapa::vegan_otu(lc_a_anneli)
sampledf_anneli <- data.frame(sample_data(lc_a_anneli))
sampledf_anneli <- merge(sampledf_anneli,ecosystem_history,by="BARCODE_ID",all=F)
rownames(sampledf_anneli) <-sampledf_anneli[,"BARCODE_ID"]

sampledf_anneli2_H <- sampledf_anneli[,c("Shannon","Aridity","Elevation","Aridity_in_sample_month","AvTemp2000_2018","AvPrec2000_2018","Erosion_risk","Slope", "Aspect","Sample_month","Temperature_in_sample_month","bas","Cmic","qO2","Water_content","LC1_2018","Coarse_fragments","Clay","Sand","pH","Electrical_conductivity","P","K","LU2","C.N","Soil_depth","Bulk_density","Mean_annual_temperature","Mean_diurnal_range","Isothermality","Temperature_seasonality","Max_temperature_of_warmest.month","Min_temperature_of_coldest.month","Temperature_annual_range","Mean_temperature_of_the_wettest_quarter","Mean_temperature_of_driest_quarter","Mean_temperature_of_warmest_quarter","Mean_temperature_of_coldest_quarter","Total_annual_precipitation","Precipitation_of_wettest_month","Precipitation_of_driest_month","Precipitation_seasonality","Precipitation_of_wettest_quarter","Precipitation_of_driest_quarter","Precipitation_of_warmest_quarter","Precipitation_of_coldest_quarter","SSM","SSMA","LC1_LUCAS15","LC1_LUCAS09","LC1")]
colnames(sampledf_anneli2_H)[colnames(sampledf_anneli2_H) == "LC1_2018"] <- "Ecosystem_type"
colnames(sampledf_anneli2_H)[colnames(sampledf_anneli2_H) == "LU2"] <- "Plant_cover"
colnames(sampledf_anneli2_H)[colnames(sampledf_anneli2_H) == "SSM"] <- "Surface_moisture"
colnames(sampledf_anneli2_H)[colnames(sampledf_anneli2_H) == "LC1_LUCAS15"] <- "Ecosystem_type_2015"
colnames(sampledf_anneli2_H)[colnames(sampledf_anneli2_H) == "LC1_LUCAS09"] <- "Ecosystem_type_2009"
sampledf_anneli5_H <-na.exclude(sampledf_anneli2_H)#feature selection
sampledf_anneli2_H <- scale(sampledf_anneli5_H[,-c(49:51,24,16)]) #scale numerical variables
sampledf_anneli2_H <- as.data.frame(sampledf_anneli2_H)
sampledf_anneli5_H <-cbind(sampledf_anneli2_H,sampledf_anneli5_H[,c(49:51,24,16)])

train.index_anneli = sample(seq_len(nrow(sampledf_anneli5_H)), size = floor(0.8 * nrow(sampledf_anneli5_H)))
train.H_anneli <- sampledf_anneli5_H[train.index_anneli,]
test.H_anneli <- sampledf_anneli5_H[-train.index_anneli,]

#rfe
ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   number=10,
                   verbose = FALSE)

rfe_anneli_H <- rfe(x=train.H_anneli[,-(1)],y=train.H_anneli[,1],
                    sizes = c(2:50),
                    rfeControl = ctrl)

rfe_anneli_H$optVariables
train.H_anneli <- train.H_anneli[,c("Shannon",rfe_anneli_H$optVariables)]

#tuning model using caret
myControl <- trainControl(
  method = "cv",
  number = 10,
  savePredictions=TRUE,
  classProbs = TRUE, # IMPORTANT!
  verboseIter = FALSE
)
x <- train.H_anneli[,c(!colnames(train.H_anneli) %in% "Shannon")]
mtry <- sqrt(ncol(x))
tunegrid <- expand.grid(.mtry=mtry)
set.seed(123)

test_caret_anneli <- train(Shannon~., data=train.H_anneli, method="rf", metric="rmse",
                           tuneGrid=tunegrid, ntree = 1000, trControl=myControl)
test_caret_anneli$bestTune <- as.numeric(test_caret_anneli$bestTune)

#ranger
rf_anneli_H <- ranger(Shannon ~., data=train.H_anneli, mtry=test_caret_anneli$bestTune,num.trees=450, min.node.size=155,
                      quantreg = TRUE,keep.inbag = TRUE,importance='impurity',
                      oob.error=TRUE,verbose=TRUE) #2.2%
rf_anneli_H
saveRDS(rf_anneli_H, "/home/koenjul/rf_anneli_H.rds")
rf_anneli_H <-readRDS("/home/koenjul/rf_anneli_H.rds")
v<-as.vector(rf_anneli_H$variable.importance)
w<-as.vector(names(rf_anneli_H$variable.importance))
vw<-as.data.frame(cbind(v,w));names(vw);vw$v<-as.numeric(vw$v)
vw2_anneli <- vw[with(vw,order(-v)),]
vw3_anneli_H <- vw2_anneli[1:50,]






#Protists
otutable_p<-ranacapa::vegan_otu(lc_p)
sampledf_p <- data.frame(sample_data(lc_p))
sampledf_p <- merge(sampledf_p,ecosystem_history,by="BARCODE_ID",all=F)
rownames(sampledf_p) <-sampledf_p[,"BARCODE_ID"]

colnames(sampledf_p)
sampledf_p2_H <- sampledf_p[,c("Shannon","Elevation","AvTemp2000_2018","AvPrec2000_2018","Erosion_risk","Slope", "Aspect","Sample_month","Precipitation_in_sample_month","Temperature_in_sample_month","Aridity","bas","Cmic","qO2","Water_content","LC1_2018","Coarse_fragments","Clay","Sand","pH","Electrical_conductivity","P","K","LU2","C.N","Soil_depth","Bulk_density","Mean_annual_temperature","Mean_diurnal_range","Isothermality","Temperature_seasonality","Max_temperature_of_warmest.month","Min_temperature_of_coldest.month","Temperature_annual_range","Mean_temperature_of_the_wettest_quarter","Mean_temperature_of_driest_quarter","Mean_temperature_of_warmest_quarter","Mean_temperature_of_coldest_quarter","Total_annual_precipitation","Precipitation_of_wettest_month","Precipitation_of_driest_month","Precipitation_seasonality","Precipitation_of_wettest_quarter","Precipitation_of_driest_quarter","Precipitation_of_warmest_quarter","Precipitation_of_coldest_quarter","SSM","SSMA","LC1_LUCAS15","LC1_LUCAS09","LC1")]
colnames(sampledf_p2_H)[colnames(sampledf_p2_H) == "LC1_2018"] <- "Ecosystem_type"
colnames(sampledf_p2_H)[colnames(sampledf_p2_H) == "LU2"] <- "Plant_cover"
colnames(sampledf_p2_H)[colnames(sampledf_p2_H) == "SSM"] <- "Surface_moisture"
colnames(sampledf_p2_H)[colnames(sampledf_p2_H) == "LC1_LUCAS15"] <- "Ecosystem_type_2015"
colnames(sampledf_p2_H)[colnames(sampledf_p2_H) == "LC1_LUCAS09"] <- "Ecosystem_type_2009"
sampledf_p5_H <-na.exclude(sampledf_p2_H)#feature selection
sampledf_p2_H <- scale(sampledf_p5_H[,-c(49:51,24,16)]) #scale numerical variables
sampledf_p2_H <- as.data.frame(sampledf_p2_H)
sampledf_p5_H <-cbind(sampledf_p2_H,sampledf_p5_H[,c(49:51,24,16)])

train.index_p = sample(seq_len(nrow(sampledf_p5_H)), size = floor(0.8 * nrow(sampledf_p5_H))) 
train.H_p <- sampledf_p5_H[train.index_p,]
test.H_p <- sampledf_p5_H[-train.index_p,]

#rfe
ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   number=10,
                   verbose = FALSE)

rfe_p_H <- rfe(x=train.H_p[,-(1)],y=train.H_p[,1],
               sizes = c(2:51),
               rfeControl = ctrl)

rfe_p_H$optVariables
train.H_p <- train.H_p[,c("Shannon",rfe_p_H$optVariables)]

#tuning model using caret
myControl <- trainControl(
  method = "cv",
  number = 10,
  savePredictions=TRUE,
  classProbs = TRUE, # IMPORTANT!
  verboseIter = FALSE
)

x <- train.H_p[,c(!colnames(train.H_p) %in% "Shannon")]

mtry <- sqrt(ncol(x))
tunegrid <- expand.grid(.mtry=mtry)
set.seed(123)

test_caret_p <- train(Shannon~., data=train.H_p, method="rf", metric="rmse",
                      tuneGrid=tunegrid, ntree = 1000, trControl=myControl)
test_caret_p$bestTune <- as.numeric(test_caret_p$bestTune)

#ranger
rf_p_H <- ranger(Shannon ~., data=train.H_p, mtry=test_caret_p$bestTune,num.trees=900, min.node.size=5,
                 quantreg = TRUE,keep.inbag = TRUE,importance='impurity',
                 oob.error=TRUE,verbose=TRUE) #8.4%
rf_p_H
saveRDS(rf_p_H, "/home/koenjul/rf_p_H.rds")
rf_p_H <-readRDS("/home/koenjul/rf_p_H.rds")
v<-as.vector(rf_p_H$variable.importance)
w<-as.vector(names(rf_p_H$variable.importance))
vw<-as.data.frame(cbind(v,w));names(vw);vw$v<-as.numeric(vw$v)
vw2_p <- vw[with(vw,order(-v)),]
vw3_p_H <- vw2_p[1:50,]



#Fungi
otutable_f <-ranacapa::vegan_otu(lc_f)
sampledf_f <- data.frame(sample_data(lc_f))
sampledf_f <- merge(sampledf_f,ecosystem_history,by="BARCODE_ID",all=F)
rownames(sampledf_f) <-sampledf_f[,"BARCODE_ID"]

sampledf_f2_H <- sampledf_f[,c("Shannon","Elevation","AvTemp2000_2018","AvPrec2000_2018","Erosion_risk","Slope", "Aspect","Sample_month","Precipitation_in_sample_month","Temperature_in_sample_month","Aridity","bas","Cmic","qO2","Water_content","LC1_2018","Coarse_fragments","Clay","Sand","pH","Electrical_conductivity","P","K","LU2","C.N","Soil_depth","Bulk_density","Mean_annual_temperature","Mean_diurnal_range","Isothermality","Temperature_seasonality","Max_temperature_of_warmest.month","Min_temperature_of_coldest.month","Temperature_annual_range","Mean_temperature_of_the_wettest_quarter","Mean_temperature_of_driest_quarter","Mean_temperature_of_warmest_quarter","Mean_temperature_of_coldest_quarter","Total_annual_precipitation","Precipitation_of_wettest_month","Precipitation_of_driest_month","Precipitation_seasonality","Precipitation_of_wettest_quarter","Precipitation_of_driest_quarter","Precipitation_of_warmest_quarter","Precipitation_of_coldest_quarter","SSM","SSMA","LC1_LUCAS15","LC1_LUCAS09","LC1")]
colnames(sampledf_f2_H)[colnames(sampledf_f2_H) == "LC1_2018"] <- "Ecosystem_type"
colnames(sampledf_f2_H)[colnames(sampledf_f2_H) == "LU2"] <- "Plant_cover"
colnames(sampledf_f2_H)[colnames(sampledf_f2_H) == "SSM"] <- "Surface_moisture"
colnames(sampledf_f2_H)[colnames(sampledf_f2_H) == "LC1_LUCAS15"] <- "Ecosystem_type_2015"
colnames(sampledf_f2_H)[colnames(sampledf_f2_H) == "LC1_LUCAS09"] <- "Ecosystem_type_2009"
sampledf_f5_H <-na.exclude(sampledf_f2_H)#feature selection
sampledf_f2_H <- scale(sampledf_f5_H[,-c(49:51,24,16)]) #scale numerical variables
sampledf_f2_H <- as.data.frame(sampledf_f2_H)
sampledf_f5_H <-cbind(sampledf_f2_H,sampledf_f5_H[,c(49:51,24,16)])

train.index_f = sample(seq_len(nrow(sampledf_f5_H)), size = floor(0.8 * nrow(sampledf_f5_H))) 
train.H_f <- sampledf_f5_H[train.index_f,]
test.H_f <- sampledf_f5_H[-train.index_f,]

#rfe
ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   number=10,
                   verbose = FALSE)

rfe_f_H <- rfe(x=train.H_f[,-(1)],y=train.H_f[,1],
               sizes = c(2:51),
               rfeControl = ctrl)
rfe_f_H$optVariables #show selected variables
train.H_f <- train.H_f[,c("Shannon",rfe_f_H$optVariables)]

#tuning model using caret
myControl <- trainControl(
  method = "cv",
  number = 10,
  savePredictions=TRUE,
  classProbs = TRUE, # IMPORTANT!
  verboseIter = FALSE
)
x <- train.H_f[,c(!colnames(train.H_f) %in% "Shannon")]
mtry <- sqrt(ncol(x))
tunegrid <- expand.grid(.mtry=mtry)
set.seed(123)

test_caret_f <- train(Shannon~., data=train.H_f, method="rf", metric="rmse",
                      tuneGrid=tunegrid, ntree = 1000, trControl=myControl)
test_caret_f$bestTune <- as.numeric(test_caret_f$bestTune)

#ranger
set.seed(123)
rf_f_H <- ranger(Shannon ~., data=train.H_f, mtry=test_caret_f$bestTune,num.trees=700, min.node.size=400,
                 quantreg = TRUE,keep.inbag = TRUE,importance='impurity',
                 oob.error=TRUE,verbose=TRUE) #7.6%
rf_f_H
saveRDS(rf_f_H, "/home/koenjul/rf_f_H.rds")
rf_f_H <-readRDS("/home/koenjul/rf_f_H.rds")
v<-as.vector(rf_f_H$variable.importance)
w<-as.vector(names(rf_f_H$variable.importance))
vw<-as.data.frame(cbind(v,w));names(vw);vw$v<-as.numeric(vw$v)
vw2_f <- vw[with(vw,order(-v)),]
vw3_f_H <- vw2_f[1:49,]



## Diagram showing var part explained Shannon index
first_column <- c("Variance")
second_column <- c(6.4)
third_column <- "Fungi"
df_f <- data.frame(first_column, second_column,third_column)

first_column <- c("Variance")
second_column <- c(11.3)
third_column <- "Protists"
df_p <- data.frame(first_column, second_column,third_column)

first_column <- c("Variance")
second_column <- c(0)
third_column <- "Rotifers"
df_roti <- data.frame(first_column, second_column,third_column)

first_column <- c("Variance")
second_column <- c(0.6)
third_column <- "Tardigrades"
df_tardi <- data.frame(first_column, second_column,third_column)

first_column <- c("Variance")
second_column <- c(2.7)
third_column <- "Nematodes"
df_nema <- data.frame(first_column, second_column,third_column)

first_column <- c("Variance")
second_column <- c(4.9)
third_column <- "Arthropods"
df_arthro <- data.frame(first_column, second_column,third_column)

first_column <- c("Variance")
second_column <- c(3.8)
third_column <- "Annelids"
df_anneli <- data.frame(first_column, second_column,third_column)

df2 <- rbind(df_roti,df_tardi,df_nema,df_p,df_f,df_anneli,df_arthro)
colnames(df2)[1] <- "Env_driver"
colnames(df2)[2] <- "Var_exp"
colnames(df2)[3] <- "group"

df2_long <- df2 %>% tidyr::gather(group, value, -Env_driver, -Var_exp,-group)
df2_long$group <- factor(df2_long$group, levels=c("Fungi","Protists","Rotifers", "Tardigrades","Nematodes","Arthropods","Annelids"))

plot_df_H <- ggplot(df2_long,aes(x=group,y=Var_exp)) 
plot_df_H <- plot_df_H + geom_bar(aes(fill=group),stat="identity") +
  ylab("\nVariation explained in %") +
  ggtitle("Shannon index")+
  ylim(0,15)+
  scale_fill_manual(values=c("#708238","#E47250","#5A4A6F","#9D5A6C","#A7C7E7","#FF99CC","orange")) +
  theme_jk5()+theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1,size=13),
                    axis.title.x = element_blank(),
                    axis.ticks = element_blank(),
                    axis.text.y = element_text(size = 14),
                    legend.text = element_text(size = 14),
                    legend.title = element_blank(),
                    plot.title = element_text(size = 14, vjust = 1, hjust = 0))

#plot Alpha diversity (richness + Shannon index), Figure S8
svg("1605_FigS8.svg",width=12,height=9)
ggarrange(plot_df_obs, plot_df_H,common.legend =T,ncol=2,legend="bottom",labels = c("A","B"),nrow=2)
dev.off() # Close the graphics device


##############################################################################
###########################################################
#Variable importance 

#add column with eukaryotic group
vw3_tardi_H$group <- "Tardigrades"
vw3_nema_H$group <- "Nematodes"
vw3_roti_H$group <- "Rotifers"
vw3_p_H$group <- "Protists"
vw3_f_H$group <- "Fungi"
vw3_anneli_H$group <- "Annelids"
vw3_arthro_H$group <- "Arthropods"

#create dataframe combining all eukaryotic organisms' Shannon index
H_imp_var_long <- rbind(vw3_f_H,vw3_p_H,vw3_tardi_H,vw3_nema_H,vw3_arthro_H,vw3_anneli_H)
H_imp_var_long <- as.data.frame(H_imp_var_long)

#rename variables
H_imp_var_long$w[H_imp_var_long$w =="AvTemp2000_2018"] <- "Mean temperature 2000-2018"
H_imp_var_long$w[H_imp_var_long$w =="AvPrec2000_2018"] <- "Mean precipitation 2000-2018"
H_imp_var_long$w[H_imp_var_long$w =="Sample_month"] <- "Sample month"
H_imp_var_long$w[H_imp_var_long$w =="Precipitation_in_sample_month"] <- "Precipitation in sample month"
H_imp_var_long$w[H_imp_var_long$w =="Temperature_in_sample_month"] <- "Temperature in sample month"
H_imp_var_long$w[H_imp_var_long$w =="Aridity_in_sample_month"] <- "Aridity in sample month"
H_imp_var_long$w[H_imp_var_long$w =="bas"] <- "Basal respiration"
H_imp_var_long$w[H_imp_var_long$w =="Cmic"] <- "rotibial biomass carbon"
H_imp_var_long$w[H_imp_var_long$w =="qO2"] <- "Respiration quotient"
H_imp_var_long$w[H_imp_var_long$w =="Water_content"] <- "Soil water content"
H_imp_var_long$w[H_imp_var_long$w =="Coarse_fragments"] <- "Coarse fragments"
H_imp_var_long$w[H_imp_var_long$w =="Electrical_conductivity"] <- "Electrical conductivity"
H_imp_var_long$w[H_imp_var_long$w =="P"] <- "Phosphorus"
H_imp_var_long$w[H_imp_var_long$w =="K"] <- "Potassium"
H_imp_var_long$w[H_imp_var_long$w =="C.N"] <- "C:N ratio"
H_imp_var_long$w[H_imp_var_long$w =="Erosion_risk"] <- "Erosion risk"
H_imp_var_long$w[H_imp_var_long$w =="Soil_depth"] <- "Soil depth"
H_imp_var_long$w[H_imp_var_long$w =="Bulk_density"] <- "Bulk density"
H_imp_var_long$w[H_imp_var_long$w =="Isothermality"] <- "Isothermality 1970-2000"
H_imp_var_long$w[H_imp_var_long$w =="Mean_annual_temperature"] <- "Mean annual temperature 1970-2000"
H_imp_var_long$w[H_imp_var_long$w =="Mean_diurnal_range"] <- "Mean diurnal range 1970-2000"
H_imp_var_long$w[H_imp_var_long$w =="Temperature_seasonality"] <- "Temperature seasonality 1970-2000"
H_imp_var_long$w[H_imp_var_long$w =="Max_temperature_of_warmest.month"] <- "Max temperature of warmest month 1970-2000"
H_imp_var_long$w[H_imp_var_long$w =="Min_temperature_of_coldest.month"] <- "Min temperature of coldest month 1970-2000"
H_imp_var_long$w[H_imp_var_long$w =="Temperature_annual_range"] <- "Temperature annual range 1970-2000"
H_imp_var_long$w[H_imp_var_long$w =="Mean_temperature_of_the_wettest_quarter"] <- "Mean temperature of wettest quarter 1970-2000"
H_imp_var_long$w[H_imp_var_long$w =="Mean_temperature_of_driest_quarter"] <- "Mean temperature of driest quarter 1970-2000"
H_imp_var_long$w[H_imp_var_long$w =="Mean_temperature_of_warmest_quarter"] <- "Mean temperature of warmest quarter 1970-2000"
H_imp_var_long$w[H_imp_var_long$w =="Mean_temperature_of_coldest_quarter"] <- "Mean temperature of coldest quarter 1970-2000"
H_imp_var_long$w[H_imp_var_long$w =="Total_annual_precipitation"] <- "Total annual precipitation 1970-2000"
H_imp_var_long$w[H_imp_var_long$w =="Precipitation_of_wettest_month"] <- "Precipitation of wettest month 1970-2000"
H_imp_var_long$w[H_imp_var_long$w =="Precipitation_of_driest_month"] <- "Precipitation of driest month 1970-2000"
H_imp_var_long$w[H_imp_var_long$w =="Precipitation_seasonality"] <- "Precipitation seasonality 1970-2000"
H_imp_var_long$w[H_imp_var_long$w =="Precipitation_of_wettest_quarter"] <- "Precipitation of wettest quarter 1970-2000"
H_imp_var_long$w[H_imp_var_long$w =="Precipitation_of_driest_quarter"] <- "Precipitation of driest quarter 1970-2000"
H_imp_var_long$w[H_imp_var_long$w =="Precipitation_of_warmest_quarter"] <- "Precipitation of warmest quarter 1970-2000"
H_imp_var_long$w[H_imp_var_long$w =="Precipitation_of_coldest_quarter"] <- "Precipitation of coldest quarter 1970-2000"
H_imp_var_long$w[H_imp_var_long$w =="SSMA"] <- "Surface moisture abnormality"
H_imp_var_long$w[H_imp_var_long$w =="Aridity"] <- "Aridity 1970-2000"
H_imp_var_long$w[H_imp_var_long$w =="Ecosystem_type"] <- "Ecosystem type 2018"
H_imp_var_long$w[H_imp_var_long$w =="Plant_cover"] <- "Plant cover"
H_imp_var_long$w[H_imp_var_long$w =="Surface_moisture"] <- "Surface moisture"
H_imp_var_long$w[H_imp_var_long$w =="Ecosystem_type_2015"] <- "Ecosystem type 2015"
H_imp_var_long$w[H_imp_var_long$w =="Ecosystem_type_2009"] <- "Ecosystem type 2009"
H_imp_var_long$w[H_imp_var_long$w =="LC1"] <- "Intensity gradient 2009-2018"
#set order of variables
H_imp_var_long$group <- factor(H_imp_var_long$group, levels=c("Fungi","Protists", "Tardigrades", "Nematodes","Arthropods","Annelids"))
H_imp_var_long$w <- factor(H_imp_var_long$w, levels=c("Ecosystem type 2018","Ecosystem type 2015","Ecosystem type 2009","Intensity gradient 2009-2018","Plant cover","Aspect","Elevation","Slope","Basal respiration","rotibial biomass carbon","Respiration quotient","Sample month","Bulk density","Carbonates","Clay", "C:N ratio", "Coarse fragments","Electrical conductivity","Erosion risk", "Potassium","Phosphorus", "pH","Sand","Soil depth","Surface moisture","Surface moisture abnormality","Soil water content","Aridity in sample month","Precipitation in sample month","Temperature in sample month","Mean temperature 2000-2018","Mean precipitation 2000-2018","Aridity 1970-2000","Isothermality 1970-2000","Max temperature of warmest month 1970-2000","Mean annual temperature 1970-2000","Mean diurnal range 1970-2000","Mean temperature of coldest quarter 1970-2000","Mean temperature of driest quarter 1970-2000","Mean temperature of wettest quarter 1970-2000","Mean temperature of warmest quarter 1970-2000","Min temperature of coldest month 1970-2000","Precipitation of coldest quarter 1970-2000","Precipitation of driest month 1970-2000","Precipitation of driest quarter 1970-2000","Precipitation of warmest quarter 1970-2000","Precipitation of wettest month 1970-2000","Precipitation of wettest quarter 1970-2000","Precipitation seasonality 1970-2000","Temperature annual range 1970-2000","Temperature seasonality 1970-2000","Total annual precipitation 1970-2000"))
H_imp_var_long <- H_imp_var_long %>% tidyr::gather(group, value, -v, -w,-group)

#setting colors
colH <- c("#387338","#387338","#387338","#387338","#387338",
          "#a0a7b0","#a0a7b0","#a0a7b0",
          "#b352a9","#b352a9","#b352a9",
          "#e5e83a",
          "#b38952","#b38952","#b38952","#b38952","#b38952","#b38952","#b38952","#b38952","#b38952","#b38952","#b38952","#b38952","#b38952","#b38952",
          "#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe","#77b5fe") 

#PLOT FIGURE S9
plot_H_imp_var_long <- ggplot(H_imp_var_long,aes(x=v,y=w))
svg("1605_FigS9.svg",width=10,height=8)
plot_H_imp_var_long + geom_bar(aes(fill=w),stat="identity") + 
  facet_wrap(~group,nrow=1, scales = "free_x")+
  ylab("\nVariation explained") +
  ggtitle("Shannon index")+
  xlab("\nVariable importance") +
  scale_fill_manual(values=c(colH)) +  
  theme_jk5()+theme(legend.position="none",
                    axis.text.x = element_text(size=10,angle=90, vjust = 0.5, hjust=1),
                    axis.text.y = element_text(size = 9),
                    strip.text = element_text(face = "bold",size=9))

dev.off() # Close the graphics device



