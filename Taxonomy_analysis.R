#############################################
#Taxonomy analysis, Figure S1-S6
#Set graphical themes
theme_jk <- function(){
  theme_bw() +
    theme(text = element_text(size=14, family = "sans"),
          axis.text = element_text(size = 14), 
          axis.title = element_text(size = 16, vjust = -3),
          plot.margin = unit(c(1, 1, 1, 1), units =, "cm"),
          plot.title = element_text(size = 16, vjust = 1, hjust = 0),
          legend.text = element_text(size = 15), 
          legend.key.size = unit(0,"cm"),
          legend.key.width= unit(0.2, "cm"),
          legend.title = element_text(size=20),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 0.5, linetype = "blank"))
}

#for species
theme_jk3 <- function(){
  theme_bw() +
    theme(text = element_text(size=14, family = "sans"),
          axis.text = element_text(size = 14), 
          axis.title = element_text(size = 16, vjust = -3),
          plot.margin = unit(c(1, 1, 1, 1), units =, "cm"),
          plot.title = element_text(size = 12, vjust = 1, hjust = 0),
          legend.text = element_text(size = 10), 
          legend.key.size = unit(0,"cm"),
          legend.key.width= unit(0.2, "cm"),
          legend.title = element_blank(),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 0.5, linetype = "blank"))
}

mydf <- data.frame(group=paste0('gr',1:10), var=paste('some long text -', 1:50), value=runif(500, 0, 100)) 
width_scale <- 16 * 27 / length(unique(mydf$var))

#rotifers
#CLASS LEVEL TOP 10
lc_a_roti_class <- phyloseq::tax_glom(lc_a_roti_rel_abund, "class")
melt_lc_a_roti_class <- phyloseq::psmelt(lc_a_roti_class) 
melt_lc_a_roti_class$LC1_2018 <- factor(melt_lc_a_roti_class$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved"), "WL_coniferous")
#PLOT
melt_lc_a_roti_class_plot <-  ggplot(data = melt_lc_a_roti_class, aes(x = LC1_2018, y=Abundance,fill= class, color=class)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Class", color="Class") +
  ggtitle("Rotifers") + 
  scale_fill_manual(values = c("#E47250",  "#EBB261","#2D4766")) +
  scale_colour_manual(values = c("#E47250",  "#EBB261","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.title=element_text(size=16, hjust=0.01))

#ORDER level TOP 10
top.order_roti <- sort(tapply(taxa_sums(lc_a_roti_rel_abund), tax_table(lc_a_roti_rel_abund)[, "order"], sum), TRUE)
top.order_roti <- top.order_roti[1:10]
# Prune to just the most-abundant 10 order
lc_roti_subset_order <- subset_taxa(lc_a_roti_rel_abund, order %in% names(top.order_roti))
lc_roti_order <- phyloseq::tax_glom(lc_roti_subset_order, "order")
roti_order <- psmelt(lc_roti_order)
unique(roti_order$order)
roti_order$order <- factor(roti_order$order, levels=c("Adinetida","Flosculariacea","Ploimida","Philodinida","Unknown"))
roti_order$LC1_2018 <- factor(roti_order$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
#PLOT
melt_lc_a_roti_order_plot <- roti_order %>%
  ggplot(data = ., aes(x = LC1_2018, y=Abundance,fill= order, color=order)) +
  geom_bar(stat="identity",position="fill") +
  labs( y = "Relative abundance\n",fill="Order", color="Order") +
  ggtitle("Rotifers") + 
  scale_fill_manual(values = c("#900C3F", "#508578", "#E47250",  "#EBB261","#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578", "#E47250",  "#EBB261","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.title=element_text(size=16, hjust=0.01))


#FAMILY LEVEL TOP 10
top.family_roti <- sort(tapply(taxa_sums(lc_a_roti_rel_abund), tax_table(lc_a_roti_rel_abund)[, "family"], sum), TRUE)
top.family_roti <- top.family_roti[1:10]
# Prune to just the most-abundant 10 family
lc_roti_subset_family <- subset_taxa(lc_a_roti_rel_abund, family %in% names(top.family_roti))
lc_roti_family <- phyloseq::tax_glom(lc_roti_subset_family, "family")
roti_family <- psmelt(lc_roti_family)
unique(roti_family$family)
roti_family$family[roti_family$family =="No rank"] <- "no rank"
roti_family$family <- factor(roti_family$family, levels=c("Flosculariidae","Trichodoridae","no rank","Unknown"))
roti_family$LC1_2018 <- factor(roti_family$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved","WL_coniferous"))
#PLOT
melt_lc_roti_family_plot <- roti_family %>%
  ggplot(data = ., aes(x = LC1_2018, y=Abundance,fill= family, color=family)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Family", color="Family") +
  ggtitle("Rotifers") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#91765e","#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#91765e","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.title=element_text(size=16, hjust=0.01))

#GENUS LEVEL TOP 10
top.genus_roti <- sort(tapply(taxa_sums(lc_a_roti_rel_abund), tax_table(lc_a_roti_rel_abund)[, "genus"], sum), TRUE)
top.genus_roti <- top.genus_roti[1:10]
# Prune to just the most-abundant 10 genus
lc_roti_subset_genus <- subset_taxa(lc_a_roti_rel_abund, genus %in% names(top.genus_roti))
lc_roti_genus <- phyloseq::tax_glom(lc_roti_subset_genus, "genus")
roti_genus <- phyloseq::psmelt(lc_roti_genus) 
roti_genus$LC1_2018 <- factor(roti_genus$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
#PLOT
melt_lc_roti_genus_plot <-  ggplot(data = roti_genus, aes(x = LC1_2018, y=Abundance,fill= genus, color=genus)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Genus", color="Genus") +
  ggtitle("Rotifers") + 
  scale_fill_manual(values = c("#900C3F", "#508578", "#E47250",  "#EBB261", "#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#E47250",  "#EBB261","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.text=element_text(size=14,face="italic"),
        legend.title=element_text(size=16, hjust=0.01))


#SPECIES LEVEL TOP 10
top.species_roti <- sort(tapply(taxa_sums(lc_a_roti_rel_abund), tax_table(lc_a_roti_rel_abund)[, "species"], sum), TRUE)
top.species_roti <- top.species_roti[1:10]
# Prune to just the most-abundant 10 species
lc_roti_subset_species <- subset_taxa(lc_a_roti_rel_abund, species %in% names(top.species_roti))
lc_roti_species <- phyloseq::tax_glom(lc_roti_subset_species, "species")
mydf <- data.frame(group=paste0('gr',1:10), var=paste('some long text -', 1:50), value=runif(500, 0, 100)) 
width_scale <- 16 * 27 / length(unique(mydf$var))
legend_title <-"Species"
roti_species <- psmelt(lc_roti_species)
unique(roti_species$species)
roti_species$species[roti_species$species =="Rotaria_rotatoria"] <- "Rotaria rotatoria"
roti_species$species[roti_species$species =="Collotheca_campanulata"] <- "Collotheca campanulata"
roti_species$species[roti_species$species =="Trichocerca_rattus"] <- "Trichocerca rattus"
roti_species$species[roti_species$species =="Lacinularia_flosculosa"] <- "Lacinularia flosculosa"

roti_species$LC1_2018 <- factor(roti_species$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
roti_species$species<- factor(roti_species$species, levels=c("Collotheca campanulata","Lacinularia flosculosa","Rotaria rotatoria","Trichocerca rattus","No rank","Unknown"))
#PLOT
melt_lc_roti_species_plot <- roti_species %>%
  ggplot(data = ., aes(x = LC1_2018, y=Abundance,fill= species,color=species)) +
  geom_bar(stat="identity",position="fill") +
  labs( y = "Relative abundance\n",fill="Species", color="Species") +
  ggtitle("Rotifers") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#E47250",  "#EBB261", "#91765e","#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#E47250",  "#EBB261","#91765e","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.title=element_text(size=16, hjust=0.01),
        legend.text=element_text(size=14,face="italic"))

#nematodes
#CLASS LEVEL TOP 10
lc_a_nema_class <- phyloseq::tax_glom(lc_a_nema_rel_abund, "class")
nema_class <- phyloseq::psmelt(lc_a_nema_class)
nema_class$LC1_2018 <- factor(nema_class$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
#PLOT
melt_lc_a_nema_class_plot <-  ggplot(data = nema_class, aes(x = LC1_2018, y=Abundance,fill= class, color=class)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Class", color="Class") +
  ggtitle("Nematodes") + 
  scale_fill_manual(values = c("#E47250",  "#EBB261","#2D4766")) +
  scale_colour_manual(values = c("#E47250",  "#EBB261","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.title=element_text(size=16, hjust=0.01))

#ORDER LEVEL TOP 10
top.order_nema <- sort(tapply(taxa_sums(lc_a_nema_rel_abund), tax_table(lc_a_nema_rel_abund)[, "order"], sum), TRUE)
top.order_nema <- top.order_nema[1:10]
lc_nema_subset_order <- subset_taxa(lc_a_nema_rel_abund, order %in% names(top.order_nema))# Prune to just the most-abundant 10 order
lc_nema_order <- phyloseq::tax_glom(lc_nema_subset_order, "order")
nema_order <- phyloseq::psmelt(lc_nema_order) 
unique(nema_order$order)
nema_order$order <- factor(nema_order$order, levels=c("Desmodorida","Diplogasterida","Dorylaimida","Mononchida" ,"Oxyurida","Rhabditida","Triplonchida","Tylenchida","No rank","Unknown"))
nema_order$LC1_2018 <- factor(nema_order$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
#PLOT
melt_lc_nema_order_plot <-  ggplot(data = nema_order, aes(x = LC1_2018, y=Abundance,fill= order, color=order)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Order", color="Order") +
  ggtitle("Nematodes") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.title=element_text(size=16, hjust=0.01),
        legend.text=element_text(size=14))

#FAMILY LEVEL TOP 10
top.family_nema <- sort(tapply(taxa_sums(lc_a_nema_rel_abund), tax_table(lc_a_nema_rel_abund)[, "family"], sum), TRUE)
top.family_nema <- top.family_nema[1:10]
# Prune to just the most-abundant 10 family
lc_nema_subset_family <- subset_taxa(lc_a_nema_rel_abund, family %in% names(top.family_nema))
lc_nema_family <- phyloseq::tax_glom(lc_nema_subset_family, "family")
nema_family <- psmelt(lc_nema_family)
unique(nema_family$family)
nema_family$family <- factor(nema_family$family, levels=c("Heteroderidae","Hoplolaimidae","Mononchidae","Mylonchulidae","Pratylenchidae","Prismatolaimidae","Rhabditidae","Tylenchulidae","No rank","Unknown"))
nema_family$LC1_2018 <- factor(nema_family$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
#PLOT
melt_lc_nema_family_plot <- nema_family %>%
  ggplot(data = ., aes(x = LC1_2018, y=Abundance,fill= family, color=family)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Family", color="Family") +
  ggtitle("Nematodes") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.title=element_text(size=16, hjust=0.01))

#GENUS LEVEL TOP 10
top.genus_nema <- sort(tapply(taxa_sums(lc_a_nema_rel_abund), tax_table(lc_a_nema_rel_abund)[, "genus"], sum), TRUE)
top.genus_nema <- top.genus_nema[1:10]
# Prune to just the most-abundant 10 genus
lc_nema_subset_genus <- subset_taxa(lc_a_nema_rel_abund, genus %in% names(top.genus_nema))
lc_nema_genus <- phyloseq::tax_glom(lc_nema_subset_genus, "genus")
nema_genus <- phyloseq::psmelt(lc_nema_genus) 
nema_genus$LC1_2018 <- factor(nema_genus$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
#PLOT
melt_lc_nema_genus_plot <-  ggplot(data = nema_genus, aes(x = LC1_2018, y=Abundance,fill= genus, color=genus)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Genus", color="Genus") +
  ggtitle("Nematodes") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.text=element_text(size=14,face="italic"),
        legend.title=element_text(size=16, hjust=0.01))

#SPECIES LEVEL TOP 10
top.species_nema <- sort(tapply(taxa_sums(lc_a_nema_rel_abund), tax_table(lc_a_nema_rel_abund)[, "species"], sum), TRUE)
top.species_nema <- top.species_nema[1:10]
# Prune to just the most-abundant 10 species
lc_nema_subset_species <- subset_taxa(lc_a_nema_rel_abund, species %in% names(top.species_nema))
lc_nema_species <- phyloseq::tax_glom(lc_nema_subset_species, "species")
nema_species <- psmelt(lc_nema_species)
unique(nema_species$species)
nema_species$species[nema_species$species =="Clarkus_sp._PDL-2005"] <- "Clarkus sp. PDL-2005"
nema_species$species[nema_species$species =="Oscheius_tipulae"] <- "Oscheius tipulae"
nema_species$species[nema_species$species =="Paratylenchus_dianthus"] <- "Paratylenchus dianthus"
nema_species$species[nema_species$species =="Phasmarhabditis_sp._EM434"] <- "Phasmarhabditis sp. EM434"
nema_species$species[nema_species$species =="Prismatolaimus_cf._dolichurus_JH-2004"] <- "Prismatolaimus cf. dolichurus JH-2004"
nema_species$species[nema_species$species =="Prismatolaimus_intermedius"] <- "Prismatolaimus intermedius"
nema_species$species[nema_species$species =="Rhabditis_cf._terricola_JH-2004"] <- "Rhabditis cf. terricola JH-2004"
nema_species$species[nema_species$species =="Rhabditis_sp._DF5059"] <- "Rhabditis sp. DF5059"
nema_species$LC1_2018 <- factor(nema_species$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
nema_species$species<- factor(nema_species$species, levels=c("Clarkus sp. PDL-2005","Oscheius tipulae","Paratylenchus dianthus","Phasmarhabditis sp. EM434","Prismatolaimus cf. dolichurus JH-2004","Prismatolaimus intermedius","Rhabditis cf. terricola JH-2004","Rhabditis sp. DF5059","No rank","Unknown"))


#PLOT
melt_lc_nema_species_plot <- nema_species %>%
  ggplot(data = ., aes(x = LC1_2018, y=Abundance,fill= species,color=species)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Species", color="Species") +
  ggtitle("Nematodes") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.title=element_text(size=16, hjust=0.01),
        legend.text=element_text(size=14,face="italic"))


#tardigrades
#CLASS LEVEL TOP 
lc_a_tardi_class <- phyloseq::tax_glom(lc_a_tardi_rel_abund, "class")
tardi_class <- phyloseq::psmelt(lc_a_tardi_class) 
tardi_class$LC1_2018 <- factor(tardi_class$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
#PLOT
melt_lc_a_tardi_class_plot <-   ggplot(data = tardi_class, aes(x = LC1_2018, y=Abundance,fill= class, color=class)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Class", color="Class") +
  ggtitle("Tardigrades") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.title=element_text(size=16, hjust=0.01))


#ORDER LEVEL TOP 10
top.order_tardi <- sort(tapply(taxa_sums(lc_a_tardi_rel_abund), tax_table(lc_a_tardi_rel_abund)[, "order"], sum), TRUE)
top.order_tardi <- top.order_tardi[1:10]
lc_tardi_subset_order <- subset_taxa(lc_a_tardi_rel_abund, order %in% names(top.order_tardi))# Prune to just the most-abundant 10 order
lc_tardi_order <- phyloseq::tax_glom(lc_tardi_subset_order, "order")
tardi_order <- psmelt(lc_tardi_order)
tardi_order$LC1_2018 <- factor(tardi_order$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
#PLOT
melt_lc_tardi_order_plot <- tardi_order %>%
  ggplot(data = ., aes(x = LC1_2018, y=Abundance,fill= order, color=order)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Order", color="Order") +
  ggtitle("Tardigrades") + 
  scale_fill_manual(values = c("#900C3F", "#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.title=element_text(size=16, hjust=0.01))


#FAMILY LEVEL TOP 10
top.family_tardi <- sort(tapply(taxa_sums(lc_a_tardi_rel_abund), tax_table(lc_a_tardi_rel_abund)[, "family"], sum), TRUE)
top.family_tardi <- top.family_tardi[1:10]
lc_tardi_subset_family <- subset_taxa(lc_a_tardi_rel_abund, family %in% names(top.family_tardi))# Prune to just the most-abundant 10 family
lc_tardi_family <- phyloseq::tax_glom(lc_tardi_subset_family, "family")
tardi_family <- psmelt(lc_tardi_family)
tardi_family$family <- factor(tardi_family$family, levels=c("Eohypsibiidae","Hypsibiidae","Isohypsibiidae","Macrobiotidae","Murrayidae" ,"No rank","Unknown"))
tardi_family$LC1_2018 <- factor(tardi_family$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
#PLOT
melt_lc_tardi_family_plot <- tardi_family %>%
  ggplot(data = ., aes(x = LC1_2018, y=Abundance,fill= family, color=family)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Family", color="Family") +
  ggtitle("Tardigrades") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#91765e","#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#91765e","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.title=element_text(size=16, hjust=0.01))


#GENUS LEVEL TOP 10
top.genus_tardi <- sort(tapply(taxa_sums(lc_a_tardi_rel_abund), tax_table(lc_a_tardi_rel_abund)[, "genus"], sum), TRUE)
top.genus_tardi <- top.genus_tardi[1:10]
# Prune to just the most-abundant 10 genus
lc_tardi_subset_genus <- subset_taxa(lc_a_tardi_rel_abund, genus %in% names(top.genus_tardi))
lc_tardi_genus <- phyloseq::tax_glom(lc_tardi_subset_genus, "genus")
tardi_genus <- psmelt(lc_tardi_genus)
tardi_genus$genus<- factor(tardi_genus$genus, levels=c("Astatumen","Diphascon","Eremobiotus","Hypsibius","Macrobiotus","Mesobiotus","Mesocrista","Murrayon","Paramacrobiotus","Unknown"))
tardi_genus$LC1_2018 <- factor(tardi_genus$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
#PLOT
melt_lc_tardi_genus_plot <- tardi_genus %>%
  ggplot(data = ., aes(x = LC1_2018, y=Abundance,fill= genus, color=genus)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Genus", color="Genus") +
  ggtitle("Tardigrades") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.text=element_text(size=14,face="italic"),
        legend.title=element_text(size=16, hjust=0.01),
        legend.text.align = 0)


#SPECIES LEVEL TOP 10
top.species_tardi <- sort(tapply(taxa_sums(lc_a_tardi_rel_abund), tax_table(lc_a_tardi_rel_abund)[, "species"], sum), TRUE)
top.species_tardi <- top.species_tardi[1:10]
# Prune to just the most-abundant 10 species
lc_tardi_subset_species <- subset_taxa(lc_a_tardi_rel_abund, species %in% names(top.species_tardi))
lc_tardi_species <- phyloseq::tax_glom(lc_tardi_subset_species, "species")
tardi_species <- psmelt(lc_tardi_species)
unique(tardi_species$species)
tardi_species$species[tardi_species$species =="Paramacrobiotus_richtersi"] <- "Paramacrobiotus richtersi"
tardi_species$species[tardi_species$species =="Eremobiotus_alicatai"] <- "Eremobiotus alicatai"
tardi_species$species[tardi_species$species =="Astatumen_trinacriae"] <- "Astatumen trinacriae"
tardi_species$species[tardi_species$species =="Mesocrista_revelata"] <- "Mesocrista revelata"
tardi_species$species[tardi_species$species =="Murrayon_dianeae"] <- "Murrayon dianeae"
tardi_species$species[tardi_species$species =="Hypsibius_dujardini"] <- "Hypsibius dujardini"
tardi_species$species[tardi_species$species =="Macrobiotus_hufelandi_group_sp._NG-2008"] <- "Macrobiotus hufelandi group sp. NG-2008"
tardi_species$species[tardi_species$species =="Mesobiotus_philippinicus"] <- "Mesobiotus philippinicus"
tardi_species$species<- factor(tardi_species$species, levels=c("Astatumen trinacriae","Eremobiotus alicatai","Hypsibius dujardini","Macrobiotus hufelandi group sp. NG-2008","Mesobiotus philippinicus","Mesocrista revelata","Murrayon dianeae","Paramacrobiotus richtersi","No rank","Unknown"))
tardi_species$LC1_2018 <- factor(tardi_species$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged", "WL_coniferous","WL_broadleaved"))
#PLOT
melt_lc_tardi_species_plot <- tardi_species %>%
  ggplot(data = ., aes(x = LC1_2018, y=Abundance,fill= species,color=species)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Species", color="Species") +
  ggtitle("Tardigrades") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.title=element_text(size=16, hjust=0.01),
        legend.text=element_text(size=14,face="italic"))
mydf <- data.frame(group=paste0('gr',1:10), var=paste('some long text -', 1:50), value=runif(500, 0, 100)) 
width_scale <- 18 * 22 / length(unique(mydf$var))






#arthropods
#CLASS LEVEL TOP 
lc_a_arthro_class <- phyloseq::tax_glom(lc_a_arthro_rel_abund, "class")
arthro_class <- phyloseq::psmelt(lc_a_arthro_class) 
arthro_class$LC1_2018 <- factor(arthro_class$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
#PLOT
melt_lc_a_arthro_class_plot <-   ggplot(data = arthro_class, aes(x = LC1_2018, y=Abundance,fill= class, color=class)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Class", color="Class") +
  ggtitle("Arthropods") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.title=element_text(size=16, hjust=0.01))


#ORDER LEVEL TOP 10
top.order_arthro <- sort(tapply(taxa_sums(lc_a_arthro_rel_abund), tax_table(lc_a_arthro_rel_abund)[, "order"], sum), TRUE)
top.order_arthro <- top.order_arthro[1:10]
lc_arthro_subset_order <- subset_taxa(lc_a_arthro_rel_abund, order %in% names(top.order_arthro))# Prune to just the most-abundant 10 order
lc_arthro_order <- phyloseq::tax_glom(lc_arthro_subset_order, "order")
arthro_order <- psmelt(lc_arthro_order)
unique(arthro_order$order)
arthro_order$order<- factor(arthro_order$order, levels=c("Coleoptera","Collembola","Diptera","Embioptera","Hymenoptera","Protura","Sarcoptiformes","Tectocepheidae","Trombidiformes","Unknown"))
arthro_order$LC1_2018 <- factor(arthro_order$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
#PLOT
melt_lc_arthro_order_plot <- arthro_order %>%
  ggplot(data = ., aes(x = LC1_2018, y=Abundance,fill= order, color=order)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Order", color="Order") +
  ggtitle("Arthropods") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.title=element_text(size=16, hjust=0.01))


#FAMILY LEVEL TOP 10
top.family_arthro <- sort(tapply(taxa_sums(lc_a_arthro_rel_abund), tax_table(lc_a_arthro_rel_abund)[, "family"], sum), TRUE)
top.family_arthro <- top.family_arthro[1:10]
lc_arthro_subset_family <- subset_taxa(lc_a_arthro_rel_abund, family %in% names(top.family_arthro))# Prune to just the most-abundant 10 family
lc_arthro_family <- phyloseq::tax_glom(lc_arthro_subset_family, "family")
arthro_family <- psmelt(lc_arthro_family)
arthro_family$family <- factor(arthro_family$family, levels=c("Bourletiellidae","Brachychthoniidae","Camisiidae","Embiidae","Epilohmanniidae","Steganacaridae","Tectocepheidae","Tullbergiidae","No rank","Unknown"))
arthro_family$LC1_2018 <- factor(arthro_family$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
#PLOT
melt_lc_arthro_family_plot <- arthro_family %>%
  ggplot(data = ., aes(x = LC1_2018, y=Abundance,fill= family, color=family)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Family", color="Family") +
  ggtitle("Arthropods") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.title=element_text(size=16, hjust=0.01))


#GENUS LEVEL TOP 10
top.genus_arthro <- sort(tapply(taxa_sums(lc_a_arthro_rel_abund), tax_table(lc_a_arthro_rel_abund)[, "genus"], sum), TRUE)
top.genus_arthro <- top.genus_arthro[1:10]
# Prune to just the most-abundant 10 genus
lc_arthro_subset_genus <- subset_taxa(lc_a_arthro_rel_abund, genus %in% names(top.genus_arthro))
lc_arthro_genus <- phyloseq::tax_glom(lc_arthro_subset_genus, "genus")
arthro_genus <- psmelt(lc_arthro_genus)
arthro_genus$genus<- factor(arthro_genus$genus, levels=c("Epilohmannia","Embia","Heterosminthurus","Oppiella","Platynothrus","Steganacarus","Tectocepheus","Tullbergia","No rank","Unknown"))
arthro_genus$LC1_2018 <- factor(arthro_genus$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
#PLOT
melt_lc_arthro_genus_plot <- arthro_genus %>%
  ggplot(data = ., aes(x = LC1_2018, y=Abundance,fill= genus, color=genus)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Genus", color="Genus") +
  ggtitle("Arthropods") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.text=element_text(size=14,face="italic"),
        legend.title=element_text(size=16, hjust=0.01))


#SPECIES LEVEL TOP 10
top.species_arthro <- sort(tapply(taxa_sums(lc_a_arthro_rel_abund), tax_table(lc_a_arthro_rel_abund)[, "species"], sum), TRUE)
top.species_arthro <- top.species_arthro[1:10]
# Prune to just the most-abundant 10 species
lc_arthro_subset_species <- subset_taxa(lc_a_arthro_rel_abund, species %in% names(top.species_arthro))
lc_arthro_species <- phyloseq::tax_glom(lc_arthro_subset_species, "species")
arthro_species <- psmelt(lc_arthro_species)
unique(arthro_species$species)
arthro_species$species[arthro_species$species =="Brachychthoniidae_gen._sp._3_AD1301"] <- "Brachychthoniidae gen. sp. 3 AD1301"
arthro_species$species[arthro_species$species =="Embia_tyrrhenica"] <- "Embia_tyrrhenica"
arthro_species$species[arthro_species$species =="Epilohmannia_pallida"] <- "Epilohmannia pallida"
arthro_species$species[arthro_species$species =="Tullbergia_yosii"] <- "Tullbergia yosii"
arthro_species$species[arthro_species$species =="Heterosminthurus_sp._SF-2017"] <- "Heterosminthurus sp. SF-2017"
arthro_species$species[arthro_species$species =="Platynothrus_peltifer"] <- "Platynothrus peltifer"
arthro_species$species[arthro_species$species =="Steganacarus_magnus"] <- "Steganacarus magnus"
arthro_species$species[arthro_species$species =="Tectocepheus_velatus_sarekensis"] <- "Tectocepheus velatus sarekensis"
arthro_species$species<- factor(arthro_species$species, levels=c("Brachychthoniidae gen. sp. 3 AD1301","Embia_tyrrhenica","Epilohmannia pallida","Heterosminthurus sp. SF-2017","Platynothrus peltifer","Steganacarus magnus","Tectocepheus velatus sarekensis","Tullbergia yosii" ,"No rank","Unknown"))
arthro_species$LC1_2018 <- factor(arthro_species$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged", "WL_coniferous","WL_broadleaved"))
#PLOT
melt_lc_arthro_species_plot <- arthro_species %>%
  ggplot(data = ., aes(x = LC1_2018, y=Abundance,fill= species,color=species)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Species", color="Species") +
  ggtitle("Arthropods") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.title=element_text(size=16, hjust=0.01),
        legend.text=element_text(size=14,face="italic"))
mydf <- data.frame(group=paste0('gr',1:10), var=paste('some long text -', 1:50), value=runif(500, 0, 100)) 
width_scale <- 18 * 22 / length(unique(mydf$var))





#annelids
#CLASS LEVEL TOP 
lc_a_anneli_class <- phyloseq::tax_glom(lc_a_anneli_rel_abund, "class")
anneli_class <- phyloseq::psmelt(lc_a_anneli_class) 
anneli_class$LC1_2018 <- factor(anneli_class$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
#PLOT
melt_lc_a_anneli_class_plot <-   ggplot(data = anneli_class, aes(x = LC1_2018, y=Abundance,fill= class, color=class)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Class", color="Class") +
  ggtitle("Annelids") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.title=element_text(size=16, hjust=0.01))


#ORDER LEVEL TOP 10
top.order_anneli <- sort(tapply(taxa_sums(lc_a_anneli_rel_abund), tax_table(lc_a_anneli_rel_abund)[, "order"], sum), TRUE)
top.order_anneli <- top.order_anneli[1:10]
lc_anneli_subset_order <- subset_taxa(lc_a_anneli_rel_abund, order %in% names(top.order_anneli))# Prune to just the most-abundant 10 order
lc_anneli_order <- phyloseq::tax_glom(lc_anneli_subset_order, "order")
anneli_order <- psmelt(lc_anneli_order)
anneli_order$order<- factor(anneli_order$order, levels=c("Haplotaxida","No rank"))
anneli_order$LC1_2018 <- factor(anneli_order$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
#PLOT
melt_lc_anneli_order_plot <- anneli_order %>%
  ggplot(data = ., aes(x = LC1_2018, y=Abundance,fill= order, color=order)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Order", color="Order") +
  ggtitle("Annelids") + 
  scale_fill_manual(values = c("#900C3F","#2D4766")) +
  scale_colour_manual(values = c("#900C3F",  "#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.title=element_text(size=16, hjust=0.01))


#FAMILY LEVEL TOP 10
top.family_anneli <- sort(tapply(taxa_sums(lc_a_anneli_rel_abund), tax_table(lc_a_anneli_rel_abund)[, "family"], sum), TRUE)
top.family_anneli <- top.family_anneli[1:10]
lc_anneli_subset_family <- subset_taxa(lc_a_anneli_rel_abund, family %in% names(top.family_anneli))# Prune to just the most-abundant 10 family
lc_anneli_family <- phyloseq::tax_glom(lc_anneli_subset_family, "family")
anneli_family <- psmelt(lc_anneli_family)
anneli_family$family <- factor(anneli_family$family, levels=c("Aeolosomatidae","Enchytraeidae","Megascolecidae","Lumbricidae","Parergodrilidae","No rank"))
anneli_family$LC1_2018 <- factor(anneli_family$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
#PLOT
melt_lc_anneli_family_plot <- anneli_family %>%
  ggplot(data = ., aes(x = LC1_2018, y=Abundance,fill= family, color=family)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Family", color="Family") +
  ggtitle("Annelids") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#E47250",  "#EBB261", "#9D5A6C","#91765e")) +
  scale_colour_manual(values = c("#900C3F", "#508578", "#E47250",  "#EBB261", "#9D5A6C","#91765e")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.title=element_text(size=16, hjust=0.01))


#GENUS LEVEL TOP 10
top.genus_anneli <- sort(tapply(taxa_sums(lc_a_anneli_rel_abund), tax_table(lc_a_anneli_rel_abund)[, "genus"], sum), TRUE)
top.genus_anneli <- top.genus_anneli[1:10]
# Prune to just the most-abundant 10 genus
lc_anneli_subset_genus <- subset_taxa(lc_a_anneli_rel_abund, genus %in% names(top.genus_anneli))
lc_anneli_genus <- phyloseq::tax_glom(lc_anneli_subset_genus, "genus")
anneli_genus <- psmelt(lc_anneli_genus)
unique(anneli_genus$genus)
anneli_genus$genus<- factor(anneli_genus$genus, levels=c("Achaeta","Buchholzia","Cernosvitoviella","Enchytronia","Hrabeiella","Marionina","Mesenchytraeus","Lumbricus","Parergodrilus","Rheomorpha","No rank","Unknown"))
anneli_genus$LC1_2018 <- factor(anneli_genus$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
#PLOT
melt_lc_anneli_genus_plot <- anneli_genus %>%
  ggplot(data = ., aes(x = LC1_2018, y=Abundance,fill= genus, color=genus)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Genus", color="Genus") +
  ggtitle("Annelids") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text=element_text(size=14,face="italic"),
        legend.text.align = 0)


#SPECIES LEVEL TOP 10
top.species_anneli <- sort(tapply(taxa_sums(lc_a_anneli_rel_abund), tax_table(lc_a_anneli_rel_abund)[, "species"], sum), TRUE)
top.species_anneli <- top.species_anneli[1:10]
# Prune to just the most-abundant 10 species
lc_anneli_subset_species <- subset_taxa(lc_a_anneli_rel_abund, species %in% names(top.species_anneli))
lc_anneli_species <- phyloseq::tax_glom(lc_anneli_subset_species, "species")
anneli_species <- psmelt(lc_anneli_species)
unique(anneli_species$species)
anneli_species$species[anneli_species$species =="Achaeta_bibulba"] <- "Achaeta bibulba"
anneli_species$species[anneli_species$species =="Buchholzia_appendiculata"] <- "Buchholzia appendiculata"
anneli_species$species[anneli_species$species =="Cernosvitoviella_cf._atrata_PDW-2010"] <- "Cernosvitoviella cf. atrata PDW-2010"
anneli_species$species[anneli_species$species =="Enchytronia_parva"] <- "Enchytronia parva"
anneli_species$species[anneli_species$species =="Cernosvitoviella_atrata"] <- "Cernosvitoviella atrata"
anneli_species$species[anneli_species$species =="Hrabeiella_periglandulata"] <- "Hrabeiella periglandulata"
anneli_species$species[anneli_species$species =="Marionina_communis"] <- "Marionina communis"
anneli_species$species[anneli_species$species =="Marionina_southerni"] <- "Marionina southerni"
anneli_species$species[anneli_species$species =="Rheomorpha_neiswestonovae"] <- "Rheomorpha neiswestonovae"
anneli_species$species[anneli_species$species =="Cernosvitoviella_atrata"] <- "Cernosvitoviella atrata"
anneli_species$species[anneli_species$species =="Lumbricus polyphemus"] <- "Lumbricus polyphemus"
anneli_species$species<- factor(anneli_species$species, levels=c("Achaeta bibulba","Buchholzia appendiculata","Cernosvitoviella atrata" ,"Cernosvitoviella cf. atrata PDW-2010","Enchytronia parva","Hrabeiella periglandulata","Marionina communis","Marionina southerni","Lumbricus_polyphemus","Rheomorpha neiswestonovae"))
anneli_species$LC1_2018 <- factor(anneli_species$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged", "WL_coniferous","WL_broadleaved"))
#PLOT
melt_lc_anneli_species_plot <- anneli_species %>%
  ggplot(data = ., aes(x = LC1_2018, y=Abundance,fill= species,color=species)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Species", color="Species") +
  ggtitle("Annelids") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.title=element_text(size=16, hjust=0.01),
        legend.text=element_text(size=14,face="italic"))
mydf <- data.frame(group=paste0('gr',1:10), var=paste('some long text -', 1:50), value=runif(500, 0, 100)) 
width_scale <- 18 * 22 / length(unique(mydf$var))











#Fungi
#PHYLUM LEVEL TOP 10
top.phyla_f <- sort(tapply(taxa_sums(lc_f_rel_abund), tax_table(lc_f_rel_abund)[, "phylum"], sum), TRUE)
top.phyla_f <- top.phyla_f[1:10]
lc_f_subset_class <- subset_taxa(lc_f_rel_abund, phylum %in% names(top.phyla_f))
lc_f_phylum <- phyloseq::tax_glom(lc_f_subset_class, "phylum")
f_phylum <- phyloseq::psmelt(lc_f_phylum)
f_phylum$LC1_2018 <- factor(f_phylum$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))

#PLOT
melt_lc_f_phylum_plot <-ggplot(data = f_phylum, aes(x = LC1_2018, y=Abundance,fill= phylum, color=phylum)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Phylum", color="Phylum") +
  ggtitle("Fungi") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.title=element_text(size=16, hjust=0.01),
        legend.text = element_text(size=15))

#CLASS LEVEL TOP 10
top.class_f <- sort(tapply(taxa_sums(lc_f_rel_abund), tax_table(lc_f_rel_abund)[, "class"], sum), TRUE)
top.class_f <- top.class_f[1:10]
# Prune to just the most-abundant 10 phyla
lc_f_subset_class <- subset_taxa(lc_f_rel_abund, class %in% names(top.class_f))
lc_f_class <- phyloseq::tax_glom(lc_f_subset_class, "class")
f_class <- phyloseq::psmelt(lc_f_class)
f_class$LC1_2018 <- factor(f_class$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
#PLOT
melt_lc_f_class_plot <-  ggplot(data = f_class, aes(x = LC1_2018, y=Abundance,fill= class, color=class)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Class", color="Class") +
  ggtitle("Fungi") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.title=element_text(size=16, hjust=0.01))


#ORDER LEVEL TOP 10
top.order_f <- sort(tapply(taxa_sums(lc_f_rel_abund), tax_table(lc_f_rel_abund)[, "order"], sum), TRUE)
top.order_f <- top.order_f[1:10]
lc_f_subset_order <- subset_taxa(lc_f_rel_abund, order %in% names(top.order_f))# Prune to just the most-abundant 10 order
lc_f_order <- phyloseq::tax_glom(lc_f_subset_order, "order")
f_order <- phyloseq::psmelt(lc_f_order) 
f_order$LC1_2018 <- factor(f_order$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
#PLOT
melt_lc_f_order_plot <- ggplot(data = f_order, aes(x = LC1_2018, y=Abundance,fill= order, color=order)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Order", color="Order") +
  ggtitle("Fungi") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.title=element_text(size=16, hjust=0.01))

#FAMILY LEVLE TOP 10
top.family_f <- sort(tapply(taxa_sums(lc_f_rel_abund), tax_table(lc_f_rel_abund)[, "family"], sum), TRUE)
top.family_f <- top.family_f[1:10]
lc_f_subset_family <- subset_taxa(lc_f_rel_abund, family %in% names(top.family_f))# Prune to just the most-abundant 10 family
lc_f_family <- phyloseq::tax_glom(lc_f_subset_family, "family")
f_family <- phyloseq::psmelt(lc_f_family)
f_family$LC1_2018 <- factor(f_family$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
#PLOT
melt_lc_f_family_plot <- ggplot(data = f_family, aes(x = LC1_2018, y=Abundance,fill= family, color=family)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Family", color="Family") +
  ggtitle("Fungi") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.title=element_text(size=16, hjust=0.01),
        legend.text=element_text(size=12))


#GENUS LEVEL TOP 10
top.genus_f <- sort(tapply(taxa_sums(lc_f_rel_abund), tax_table(lc_f_rel_abund)[, "genus"], sum), TRUE)
top.genus_f <- top.genus_f[1:10]
lc_f_subset_genus <- subset_taxa(lc_f_rel_abund, genus %in% names(top.genus_f))# Prune to just the most-abundant 10 genus
lc_f_genus <- phyloseq::tax_glom(lc_f_subset_genus, "genus")
f_genus <- psmelt(lc_f_genus)
f_genus$genus<- factor(f_genus$genus, levels=c("Archaeorhizomyces","Chaetomium","Craterellus","Fusarium","Gymnosporangium","Lactarius","Mortierella","Russula","No rank","Unknown"))
f_genus$LC1_2018 <- factor(f_genus$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
#PLOT
melt_lc_f_genus_plot <- f_genus %>%
  ggplot(data = ., aes(x = LC1_2018, y=Abundance,fill= genus, color=genus)) +
  geom_bar(stat="identity",position="fill") +
  labs( y = "Relative abundance\n",fill="Genus", color="Genus") +
  ggtitle("Fungi") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.text=element_text(size=14,face="italic"),
        legend.title=element_text(size=16, hjust=0.01))

#SPECIES LEVEL TOP 10
top.species_f <- sort(tapply(taxa_sums(lc_f_rel_abund), tax_table(lc_f_rel_abund)[, "species"], sum), TRUE)
top.species_f <- top.species_f[1:10]
lc_f_subset_species <- subset_taxa(lc_f_rel_abund, species %in% names(top.species_f))# Prune to just the most-abundant 10 species
lc_f_species <- phyloseq::tax_glom(lc_f_subset_species, "species")
f_species <- psmelt(lc_f_species)
f_species$species[f_species$species =="Archaeorhizomyces_borealis"] <- "Archaeorhizomyces borealis"
f_species$species[f_species$species =="Boeremia_exigua_var._exigua"] <- "Boeremia exigua var.exigua"
f_species$species[f_species$species =="Gymnosporangium_ellisii"] <- "Gymnosporangium ellisii"
f_species$species[f_species$species =="Mycena_plumbea"] <- "Mycena plumbea"
f_species$species[f_species$species =="Piloderma_fallax"] <- "Piloderma fallax"
f_species$species[f_species$species =="Pseudogymnoascus_destructans"] <- "Pseudogymnoascus destructans"
f_species$species[f_species$species =="Trechispora_alnicola"] <- "Trechispora alnicola"
f_species$species[f_species$species =="uncultured_fungus"] <- "uncultured fungus"
f_species$species[f_species$species =="uncultured_Madurella"] <- "uncultured Madurella"
f_species$species<- factor(f_species$species, levels=c("Archaeorhizomyces borealis","Boeremia exigua var.exigua","Gymnosporangium ellisii","Mycena plumbea","Piloderma fallax","Pseudogymnoascus destructans","Trechispora alnicola","uncultured fungus","uncultured Madurella","Unknown"))
f_species$LC1_2018 <- factor(f_species$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
#PLOT
melt_lc_f_species_plot <- f_species %>%
  ggplot(data = ., aes(x = LC1_2018, y=Abundance,fill= species,color=species)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Species", color="Species") +
  ggtitle("Fungi") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.title=element_text(size=16, hjust=0.01),
        legend.text=element_text(size=14,face="italic"))

#PROTISTS
#PHYLUM
top.phyla_p <- sort(tapply(taxa_sums(lc_p_rel_abund), tax_table(lc_p_rel_abund)[, "phylum"], sum), TRUE)
top.phyla_p <- top.phyla_p[1:10]
lc_p_subset_phylum <- subset_taxa(lc_p_rel_abund, phylum %in% names(top.phyla_p))# Prune to just the most-abundant 10 phyla
lc_p_phylum <- phyloseq::tax_glom(lc_p_subset_phylum, "phylum")
p_phylum <- phyloseq::psmelt(lc_p_phylum)
p_phylum$LC1_2018 <- factor(p_phylum$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
#PLOT
melt_lc_p_phylum_plot <- ggplot(data = p_phylum, aes(x = LC1_2018, y=Abundance,fill= phylum, color=phylum)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Phylum", color="Phylum") +
  ggtitle("Protists") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.title=element_text(size=16, hjust=0.01),
        legend.text = element_text(size=15))


#CLASS LEVEL TOP 10
top.class_p <- sort(tapply(taxa_sums(lc_p_rel_abund), tax_table(lc_p_rel_abund)[, "class"], sum), TRUE)
top.class_p <- top.class_p[1:10]
lc_p_subset_class <- subset_taxa(lc_p_rel_abund, class %in% names(top.class_p))# Prune to just the most-abundant 10 phyla
lc_p_class <- phyloseq::tax_glom(lc_p_subset_class, "class")
p_class <- psmelt(lc_p_class)
unique(p_class$class)
p_class$class <- factor(p_class$class, levels=c("Cercomonadidae","Glissomonadida","Imbricatea","Intramacronucleata","Phytomyxea","Thecofilosea","Tubulinea","Vampyrellidae","No rank", "Unknown"))
p_class$LC1_2018 <- factor(p_class$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
#PLOT
melt_lc_p_class_plot <- p_class %>%
  ggplot(data = ., aes(x = as.factor(LC1_2018), y=Abundance,fill= class, color=class)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Class", color="Class") +
  ggtitle("Protists") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
                axis.text.y= element_text(size=13),
                legend.title=element_text(size=16, hjust=0.01))

#ORDER LEVEL TOP 10
top.order_p <- sort(tapply(taxa_sums(lc_p_rel_abund), tax_table(lc_p_rel_abund)[, "order"], sum), TRUE)
top.order_p <- top.order_p[1:10]
lc_p_subset_order <- subset_taxa(lc_p_rel_abund, order %in% names(top.order_p))# Prune to just the most-abundant 10 order
lc_p_order <- phyloseq::tax_glom(lc_p_subset_order, "order")
p_order <- psmelt(lc_p_order)
p_order$order <- factor(p_order$order, levels=c("Cercomonadida","Conthreep","Eucoccidiorida","Litostomatea","Plasmodiophorida","Pythiales","Spirotrichea","Vampyrellida","No rank","Unknown"))
p_order$LC1_2018 <- factor(p_order$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
#PLOT
melt_lc_p_order_plot <- p_order %>%
  ggplot(data = ., aes(x = LC1_2018, y=Abundance,fill= order, color=order)) +
  geom_bar(stat="identity",position="fill") +
  labs( y = "Relative abundance\n",fill="Order", color="Order") +
  ggtitle("Protists") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.title=element_text(size=16, hjust=0.01))

#FAMILY level TOP 10
top.family_p <- sort(tapply(taxa_sums(lc_p_rel_abund), tax_table(lc_p_rel_abund)[, "family"], sum), TRUE)
top.family_p <- top.family_p[1:10]
lc_p_subset_family <- subset_taxa(lc_p_rel_abund, family %in% names(top.family_p))# Prune to just the most-abundant 10 family
lc_p_family <- phyloseq::tax_glom(lc_p_subset_family, "family")
p_family <- psmelt(lc_p_family)
p_family$family <- factor(p_family$family, levels=c("Cercomonadidae","Colpodea","Eimeriidae","Haptoria","Hypotrichia","Leptophryidae","Plasmodiophoridae","Pythiaceae","No rank","Unknown"))
p_family$LC1_2018 <- factor(p_family$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
#PLOT
melt_lc_p_family_plot <- p_family %>%
  ggplot(data = ., aes(x = LC1_2018, y=Abundance,fill= family, color=family)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Family", color="Family") +
  ggtitle("Protists") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.title=element_text(size=16, hjust=0.01))

#GENUS LEVEL TOP 10
top.genus_p <- sort(tapply(taxa_sums(lc_p_rel_abund), tax_table(lc_p_rel_abund)[, "genus"], sum), TRUE)
top.genus_p <- top.genus_p[1:10]
lc_p_subset_genus <- subset_taxa(lc_p_rel_abund, genus %in% names(top.genus_p))# Prune to just the most-abundant 10 genus
lc_p_genus <- phyloseq::tax_glom(lc_p_subset_genus, "genus")
p_genus <- psmelt(lc_p_genus)
p_genus$genus<- factor(p_genus$genus, levels=c("Arachnula","Cercomonas","Colpoda","Globisporangium","Gonostomum","Polymyxa","Rhogostoma","Uroleptoides","No rank","Unknown"))
p_genus$LC1_2018 <- factor(p_genus$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
#PLOT
melt_lc_p_genus_plot <- p_genus %>%
  ggplot(data = ., aes(x = LC1_2018, y=Abundance,fill= genus, color=genus)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Genus", color="Genus") +
  ggtitle("Protists") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.title=element_text(size=16, hjust=0.01),
        legend.text=element_text(size=14,face="italic"))

#SPECIES LEVEL TOP 10
top.species_p <- sort(tapply(taxa_sums(lc_p_rel_abund), tax_table(lc_p_rel_abund)[, "species"], sum), TRUE)
top.species_p <- top.species_p[1:10]
lc_p_subset_species <- subset_taxa(lc_p_rel_abund, species %in% names(top.species_p))# Prune to just the most-abundant 10 species
lc_p_species <- phyloseq::tax_glom(lc_p_subset_species, "species")
p_species <- psmelt(lc_p_species)
p_species$species
p_species$species[p_species$species =="Arachnula_impatiens"] <- "Arachnula impatiens"
p_species$species[p_species$species =="Heteromita_globosa"] <- "Heteromita globosa"
p_species$species[p_species$species =="Paragonostomoides_xianicum"] <- "Paragonostomoides xianicum"
p_species$species[p_species$species =="Polymyxa_graminis"] <- "Polymyxa graminis"
p_species$species[p_species$species =="uncultured_Eimeriidae"] <- "uncultured Eimeriidae"
p_species$species[p_species$species =="uncultured_Cercozoa"] <- "uncultured Cercozoa"
p_species$species[p_species$species =="uncultured_eukaryote"] <- "uncultured eukaryote"
p_species$species[p_species$species =="uncultured_Oxytrichidae"] <- "uncultured Oxytrichidae"
p_species$species[p_species$species =="Uroleptoides_longiseries"] <- "Uroleptoides longiseries"
p_species$species<- factor(p_species$species, levels=c("Arachnula impatiens","Heteromita globosa","Paragonostomoides xianicum","Polymyxa graminis","uncultured Eimeriidae","uncultured Cercozoa","uncultured eukaryote","uncultured Oxytrichidae","Uroleptoides longiseries","Unknown"))
p_species$LC1_2018 <- factor(p_species$LC1_2018, levels=c("CL_annual", "CL_permanent", "GL_managed", "GL_unmanaged","WL_broadleaved", "WL_coniferous"))
#PLOT
melt_lc_p_species_plot <- p_species %>%
  ggplot(data = ., aes(x = LC1_2018, y=Abundance,fill= species,color=species)) +
  geom_bar(stat="identity",position="fill") +
  labs(y = "Relative abundance\n",fill="Species", color="Species") +
  ggtitle("Protists") + 
  scale_fill_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_colour_manual(values = c("#900C3F", "#508578","#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C","#708238","#BAD6E5","#91765e","#2D4766")) +
  scale_x_discrete(name="Ecosystem type",labels=c("CL1","CL2","GL1","GL2","WL1","WL2"))+
  theme_jk()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=13),
        axis.text.y= element_text(size=13),
        legend.title=element_text(size=16, hjust=0.01),
        legend.text=element_text(size=14,face="italic"))


#PLOT ALL EUKARYOTIC GROUPS
#PHYLUM Figure S1
svg("1605_FigS1.svg",width=12,height=5)
#plot
ggarrange(melt_lc_f_phylum_plot,melt_lc_p_phylum_plot,
           labels = c("(a)", "(b)"),align="v",common.legend = F,
           ncol = 2, nrow = 1)
dev.off() # Close the graphics device

png("1905_FigS1.png",width=1500,height=800)
ggarrange(melt_lc_f_phylum_plot,melt_lc_p_phylum_plot,
          labels = c("(a)", "(b)"),align="v",common.legend = F,
          ncol = 2, nrow = 1)
dev.off() # Close the graphics device



#CLASS Figure S2
svg("1805_FigS2.svg",width=14,height=13)
ggarrange(melt_lc_f_class_plot,melt_lc_p_class_plot,melt_lc_a_roti_class_plot,melt_lc_a_tardi_class_plot,melt_lc_a_nema_class_plot,melt_lc_a_arthro_class_plot,melt_lc_a_anneli_class_plot,
          labels = c("(a)", "(b)","(c)","(d)","(e)","(f)","(g)"),align="v",common.legend = F,
           ncol = 2, nrow = 4)
dev.off() # Close the graphics device


png("1905_FigS2.png",width=1500,height=1300)
ggarrange(melt_lc_f_class_plot,melt_lc_p_class_plot,melt_lc_a_roti_class_plot,melt_lc_a_tardi_class_plot,melt_lc_a_nema_class_plot,melt_lc_a_arthro_class_plot,melt_lc_a_anneli_class_plot,
          labels = c("(a)", "(b)","(c)","(d)","(e)","(f)","(g)"),align="v",common.legend = F,
          ncol = 2, nrow = 4)
dev.off() # Close the graphics device




#ORDER Figure S3
svg("1805_FigS3.svg",width=14,height=13)
ggarrange(melt_lc_f_order_plot,melt_lc_p_order_plot,melt_lc_a_roti_order_plot,melt_lc_tardi_order_plot,melt_lc_nema_order_plot,melt_lc_arthro_order_plot,melt_lc_anneli_order_plot,
           labels = c("(a)", "(b)","(c)","(d)","(e)","(f)","(g)"),align="v",common.legend = F,
          ncol = 2, nrow = 4)
dev.off() # Close the graphics device


png("1905_FigS3.png",width=1500,height=1300)
ggarrange(melt_lc_f_order_plot,melt_lc_p_order_plot,melt_lc_a_roti_order_plot,melt_lc_tardi_order_plot,melt_lc_nema_order_plot,melt_lc_arthro_order_plot,melt_lc_anneli_order_plot,
          labels = c("(a)", "(b)","(c)","(d)","(e)","(f)","(g)"),align="v",common.legend = F,
          ncol = 2, nrow = 4)
dev.off() # Close the graphics device


#FAMILY Figure S4
svg("1805_FigS4.svg",width=14,height=13)
ggarrange(melt_lc_f_family_plot,melt_lc_p_family_plot,melt_lc_roti_family_plot,melt_lc_tardi_family_plot,melt_lc_nema_family_plot,melt_lc_arthro_family_plot,melt_lc_anneli_family_plot,
          labels = c("A", "B", "C","D","E","F","G"),align="v",common.legend = F,
          ncol = 2, nrow = 4)
dev.off() # Close the graphics device

png("1905_FigS4.png",width=1500,height=1300)
ggarrange(melt_lc_f_family_plot,melt_lc_p_family_plot,melt_lc_roti_family_plot,melt_lc_tardi_family_plot,melt_lc_nema_family_plot,melt_lc_arthro_family_plot,melt_lc_anneli_family_plot,
          labels = c("A", "B", "C","D","E","F","G"),align="v",common.legend = F,
          ncol = 2, nrow = 4)
dev.off() # Close the graphics device




#GENUS Figure S5
svg("1805_FigS5.svg",width=14,height=13)
ggarrange(melt_lc_f_genus_plot,melt_lc_p_genus_plot,melt_lc_roti_genus_plot,melt_lc_tardi_genus_plot,melt_lc_nema_genus_plot,melt_lc_arthro_genus_plot,melt_lc_anneli_genus_plot,
          labels = c("A", "B", "C","D","E","F","G"),align="v",common.legend = F,
          ncol = 2, nrow = 4)
dev.off() # Close the graphics device

png("1905_FigS5.png",width=1500,height=1300)
ggarrange(melt_lc_f_genus_plot,melt_lc_p_genus_plot,melt_lc_roti_genus_plot,melt_lc_tardi_genus_plot,melt_lc_nema_genus_plot,melt_lc_arthro_genus_plot,melt_lc_anneli_genus_plot,
          labels = c("(a)", "(b)","(c)","(d)","(e)","(f)","(g)"),align="v",common.legend = F,
          ncol = 2, nrow = 4)
dev.off() # Close the graphics device



#SPECIES Figure S6
svg("1805_FigS6.svg",width=14,height=13)
ggarrange(melt_lc_f_species_plot,melt_lc_p_species_plot,melt_lc_roti_species_plot,melt_lc_tardi_species_plot,melt_lc_nema_species_plot,melt_lc_arthro_species_plot,melt_lc_anneli_species_plot,
          labels = c("(a)", "(b)","(c)","(d)","(e)","(f)","(g)"),align="v",common.legend = F,ncol = 2, nrow =4)
dev.off() # Close the graphics device


png("1905_FigS6.png",width=1500,height=1300)
ggarrange(melt_lc_f_species_plot,melt_lc_p_species_plot,melt_lc_roti_species_plot,melt_lc_tardi_species_plot,melt_lc_nema_species_plot,melt_lc_arthro_species_plot,melt_lc_anneli_species_plot,
          labels = c("(a)", "(b)","(c)","(d)","(e)","(f)","(g)"),align="v",common.legend = F,ncol = 2, nrow =4)
dev.off() # Close the graphics device








