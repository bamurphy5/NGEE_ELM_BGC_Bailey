#5/29/24

#looking at the input data for our BEO ELM-PFLOTRAN simulations so I know what observational data to try and find from the list of datasets
#that Neslihan sent, to make sure our inputs all look good/as close to reality as possible before moving forward w/trying to fix chemistry issues


#----------------------------------------------------------------------------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ncdf4)
library(ggplot2)
library(lemon)



#----------------
#Hard-coded variables
#----------------

site_name <- "BEO" 

run_type <- "Default"

dat.base <- paste0("~/NGEE_modeling/ELM_PFLOTRAN_BEO/1_Runs/", run_type, "_PFT")
veg_data <- "~/NGEE_modeling/ELM_PFLOTRAN_BEO/2_Analysis/observational_data/vegetation/community_comp/data"
veg_save <- "~/NGEE_modeling/ELM_PFLOTRAN_BEO/2_Analysis/observational_data/vegetation"
soil_phys_chem_path <- "~/NGEE_modeling/ELM_PFLOTRAN_BEO/2_Analysis/observational_data/soil/physical_chem_hydro/data"
setwd(dat.base)

#----------------
#Load & process data 
#----------------

surf_data <- ncdf4::nc_open("BEO_surfdata_multicell.nc")

attributes(surf_data$dim)
names(surf_data$var)


print(surf_data) #this shows all the variables, their long name, units, data type, etc.

# Save the print(nc) dump to a text file, then can quickly access info for all vars
{
  sink('BEO_surfdata_multicell.txt')
  print(surf_data)
  sink()
}

# {
#   sink('clm_params.txt')
#   print(params)
#   sink()
# }

lon <- ncvar_get(surf_data, "LONGXY")
lat <- ncvar_get(surf_data, "LATIXY")

PFT <- ncvar_get(surf_data, "PCT_NAT_PFT")
LAI <- ncvar_get(surf_data, "MONTHLY_LAI")
SAI <- ncvar_get(surf_data, "MONTHLY_SAI")
ORG <- ncvar_get(surf_data, "ORGANIC")



#pull in obs data of fractional coverage to compare to what we have in the current model setup------------------------------------------------------
#data are from 2012 intensive field survey along the BEO transects, polygon types specified
#PFT's associated w/species are described in "V_1_2_species_list_plant_community_survey_Barrow_2012_v1.pdf" (accompanies dataset)
obs_veg_cover <- read.csv(file.path(veg_data, "V_1_2_plant_community_composition_Barrow_2012_v1.csv"), header = TRUE)

#veg is listed by species, bare ground also included, everything reported as fractional cover
#polygon types are split into two cols, "polygon_sub_unit" says whether its center, edge, or trough and "polygon_type"
#says whether high or low centered 

#first create a new polygon_type col that is same format I've been using
obs_veg_cover <- obs_veg_cover %>%
  mutate(polygon_type = case_when(polygon_type == "low-center" ~ "LCP",
                                  polygon_type == "high-center" ~ "HCP",
                                  polygon_type == "transitional" ~ "FCP")) #"transition"
#edge is same as rim
obs_veg_cover$polygon_sub_unit <- ifelse(obs_veg_cover$polygon_sub_unit == "edge", "rim", obs_veg_cover$polygon_sub_unit)


obs_veg_cover <- obs_veg_cover %>%
  mutate(polygon_type = paste0(polygon_type, polygon_sub_unit))

#make "transitioncenter" "HCPtransition" and "transitiontrough" "LCPtransition"
# obs_veg_cover <- obs_veg_cover %>%
#   mutate(polygon_type = case_when(polygon_type == "transitioncenter" ~ "HCPtransition",
#                                   polygon_type == "transitiontrough" ~ "LCPtransition",
#                                   .default = as.character(polygon_type)))


# #now fix it so anything that has 'trough' in the string just says 'Trough' (instead of "LCPtrough" for example)
# obs_veg_cover <- obs_veg_cover %>%
#   mutate(polygon_type = str_replace_all(polygon_type, ".*trough.*", "Trough"))
# 
# #do same thing so anything with 'edge' becomes 'Rim'
# obs_veg_cover <- obs_veg_cover %>%
#   mutate(polygon_type = str_replace_all(polygon_type, ".*edge.*", "Rim"))


#now handle converting the species into PFTs
#species are currently organized as individual cols, so convert to long format first
names(obs_veg_cover)

# [13] "Carex_aquatilis"             "Dupontia_fisheri"            "Eriophorum_angustifolium"    "Eriophorum_sp."             
# [17] "Juncus_biglumis"             "Arctoagrostis_latifolia"     "Alopercurus_alpinus"         "Poa_arctica"                
# [21] "Poa_vivaparum"               "Luzula_arctica"              "Luzula_confusa"              "Salix_rotundifolia"         
# [25] "Vaccinium_vitis_idaea"       "Cardamine_pratensis"         "Cochlearia_officinalis"      "Petasites_frigidus"         
# [29] "Ranunculus_nivalis"          "Ranunculus_pallasii"         "Saxifraga_cernua"            "Saxifraga_foliolosa"        
# [33] "Saxifraga_hirculus"          "Stellaria_sp."               "Alectoria_sp."               "Cetraria_cucullata"         
# [37] "Cetraria_islandica"          "Cetraria_laevigata"          "Cetraria_nivalis"            "Cetraria_richardsonii"      
# [41] "Cladonia_cornuta"            "Cladina_mitis"               "Cladonia_pyxidata"           "Cladonia_squamosa"          
# [45] "Cladonia_uncialis"           "Cladonia_amaurocraea"        "Dactylina_arctica"           "Peltigera_apthosa"          
# [49] "Peltigera_canina"            "Sphaerophorus_globosus"      "Thamnolia_subuliformis"      "White_crust"                
# [53] "Aulacomnion_turgidum"        "Bryum_sp."                   "Calligerion_sarmentosum"     "Cinclidium_sp."             
# [57] "Drepanocladus_sp."           "Pogonatum_sp."               "Polytrichum_sp."             "Sphagnum_sp."               
# [61] "Mixed_mosses_and_liverworts" "Leafy_liverworts"            "Thallose_liverworts"         "Bare_ground" 

obs_veg_cover_long <- obs_veg_cover %>%
  pivot_longer(cols = Carex_aquatilis:Bare_ground, names_to = "species", values_to = "fractional_cover")

#since current simulations only have default PFTs, grouping into these instead of the more detailed PFTs provided in the 
#pdf accompanying the data, this means some things won't really be counted...like lichen and mosses

#with the default PFT breakdown......---------------------------------------------

#1 (not vegetated)
#10 (broadleaf_evergreen_shrub)
#12 (broadleaf_deciduous_boreal_shrub)
#13 (c3_arctic_grass) 

PFT1 <- c("Bare_ground", "Alectoria_sp.", "Cetraria_cucullata", "Cetraria_islandica", "Cetraria_laevigata",
          "Cetraria_nivalis", "Cetraria_richardsonii", "Cladonia_cornuta", "Cladina_mitis", "Cladonia_pyxidata", "Cladonia_squamosa",
          "Cladonia_uncialis", "Cladonia_amaurocraea", "Dactylina_arctica", "Peltigera_apthosa",          
          "Peltigera_canina", "Sphaerophorus_globosus", "Thamnolia_subuliformis", "White_crust",                
          "Aulacomnion_turgidum", "Bryum_sp.", "Calligerion_sarmentosum", "Cinclidium_sp.",             
          "Drepanocladus_sp.", "Pogonatum_sp.", "Polytrichum_sp.", "Sphagnum_sp.",               
          "Mixed_mosses_and_liverworts", "Leafy_liverworts", "Thallose_liverworts") #grouping moss and lichen in w/bare ground

#9/10/24: trying a version where lichen is grouped into bare ground, but moss is grouped into grasses (PFT13)
PFT1 <- c("Bare_ground", "Alectoria_sp.", "Cetraria_cucullata", "Cetraria_islandica", "Cetraria_laevigata",
          "Cetraria_nivalis", "Cetraria_richardsonii", "Cladonia_cornuta", "Cladina_mitis", "Cladonia_pyxidata", "Cladonia_squamosa",
          "Cladonia_uncialis", "Cladonia_amaurocraea", "Dactylina_arctica", "Peltigera_apthosa",          
          "Peltigera_canina", "Sphaerophorus_globosus", "Thamnolia_subuliformis", "White_crust") #bare ground and lichen

PFT10 <- c("Vaccinium_vitis_idaea") #evergreen shrub

PFT12 <- c("Salix_rotundifolia") #deciduous shrub

PFT13 <- c("Carex_aquatilis", "Dupontia_fisheri", "Eriophorum_angustifolium", "Eriophorum_sp.", "Juncus_biglumis",
           "Arctoagrostis_latifolia", "Alopercurus_alpinus", "Poa_arctica", "Poa_vivaparum", "Luzula_arctica",
           "Luzula_confusa", "Cardamine_pratensis", "Cochlearia_officinalis", "Petasites_frigidus", "Ranunculus_nivalis", 
           "Ranunculus_pallasii",
           "Saxifraga_cernua", "Saxifraga_foliolosa", "Saxifraga_hirculus", "Stellaria_sp.") #c3_arctic_grass 

#version w/mosses grouped into grasses
PFT13 <- c("Carex_aquatilis", "Dupontia_fisheri", "Eriophorum_angustifolium", "Eriophorum_sp.", "Juncus_biglumis",
           "Arctoagrostis_latifolia", "Alopercurus_alpinus", "Poa_arctica", "Poa_vivaparum", "Luzula_arctica",
           "Luzula_confusa", "Cardamine_pratensis", "Cochlearia_officinalis", "Petasites_frigidus", "Ranunculus_nivalis", 
           "Ranunculus_pallasii",
           "Saxifraga_cernua", "Saxifraga_foliolosa", "Saxifraga_hirculus", "Stellaria_sp.",
           "Aulacomnion_turgidum", "Bryum_sp.", "Calligerion_sarmentosum", "Cinclidium_sp.",             
           "Drepanocladus_sp.", "Pogonatum_sp.", "Polytrichum_sp.", "Sphagnum_sp.",               
           "Mixed_mosses_and_liverworts", "Leafy_liverworts", "Thallose_liverworts") #c3_arctic_grass and mosses


obs_veg_cover_long <- obs_veg_cover_long %>%
  mutate(PFT = case_when(species %in% PFT1 ~ "PFT1",
                         species %in% PFT10 ~ "PFT10",
                         species %in% PFT12 ~ "PFT12",
                         species %in% PFT13 ~ "PFT13"))



#with the default PFT breakdown still, but separating moss and lichen out from bare ground---------------------------

#1 (not vegetated)
#10 (broadleaf_evergreen_shrub)
#12 (broadleaf_deciduous_boreal_shrub)
#13 (c3_arctic_grass) 

Moss <- c("Aulacomnion_turgidum", "Bryum_sp.", "Calligerion_sarmentosum", "Cinclidium_sp.",             
          "Drepanocladus_sp.", "Pogonatum_sp.", "Polytrichum_sp.", "Sphagnum_sp.",               
          "Mixed_mosses_and_liverworts", "Leafy_liverworts", "Thallose_liverworts")

Lichen <- c("Alectoria_sp.", "Cetraria_cucullata", "Cetraria_islandica", "Cetraria_laevigata",
            "Cetraria_nivalis", "Cetraria_richardsonii", "Cladonia_cornuta", "Cladina_mitis", "Cladonia_pyxidata", "Cladonia_squamosa",
            "Cladonia_uncialis", "Cladonia_amaurocraea", "Dactylina_arctica", "Peltigera_apthosa",          
            "Peltigera_canina", "Sphaerophorus_globosus", "Thamnolia_subuliformis", "White_crust")

PFT1 <- c("Bare_ground") 

PFT10 <- c("Vaccinium_vitis_idaea") #evergreen shrub

PFT12 <- c("Salix_rotundifolia") #deciduous shrub

PFT13 <- c("Carex_aquatilis", "Dupontia_fisheri", "Eriophorum_angustifolium", "Eriophorum_sp.", "Juncus_biglumis",
           "Arctoagrostis_latifolia", "Alopercurus_alpinus", "Poa_arctica", "Poa_vivaparum", "Luzula_arctica",
           "Luzula_confusa", "Cardamine_pratensis", "Cochlearia_officinalis", "Petasites_frigidus", "Ranunculus_nivalis", 
           "Ranunculus_pallasii",
           "Saxifraga_cernua", "Saxifraga_foliolosa", "Saxifraga_hirculus", "Stellaria_sp.") #c3_arctic_grass




obs_veg_cover_long <- obs_veg_cover_long %>%
  mutate(PFT = case_when(species %in% PFT1 ~ "PFT1",
                         species %in% PFT10 ~ "PFT10",
                         species %in% PFT12 ~ "PFT12",
                         species %in% PFT13 ~ "PFT13",
                         species %in% Moss ~ "Moss",
                         species %in% Lichen ~ "Lichen"))


#with the arctic PFT breakdown ---------------------------------------------------------------

# ELM_PFT_arctic_number	ELM_PFT_arctic_name
# 1	not_vegetated
# 2	arctic_lichen
# 3	arctic_bryophyte
# 4	arctic_needleleaf_tree 
# 5	arctic_broadleaf_tree  
# 6	arctic_evergreen_shrub_dwarf 
# 7	arctic_evergreen_shrub_tall  
# 8	arctic_deciduous_shrub_dwarf 
# 9	arctic_deciduous_shrub_low 
# 10	arctic_deciduous_shrub_tall  
# 11	arctic_deciduous_shrub_alder 
# 12	arctic_forb  
# 13	arctic_dry_graminoid 
# 14	arctic_wet_graminoid


PFT1 <- c("Bare_ground") 

PFT2 <- c("Alectoria_sp.", "Cetraria_cucullata", "Cetraria_islandica", "Cetraria_laevigata",
            "Cetraria_nivalis", "Cetraria_richardsonii", "Cladonia_cornuta", "Cladina_mitis", "Cladonia_pyxidata", "Cladonia_squamosa",
            "Cladonia_uncialis", "Cladonia_amaurocraea", "Dactylina_arctica", "Peltigera_apthosa",          
            "Peltigera_canina", "Sphaerophorus_globosus", "Thamnolia_subuliformis", "White_crust") #arctic_lichen

PFT3 <- c("Aulacomnion_turgidum", "Bryum_sp.", "Calligerion_sarmentosum", "Cinclidium_sp.",             
          "Drepanocladus_sp.", "Pogonatum_sp.", "Polytrichum_sp.", "Sphagnum_sp.",               
          "Mixed_mosses_and_liverworts", "Leafy_liverworts", "Thallose_liverworts") #arctic_bryophyte


PFT6 <- c("Vaccinium_vitis_idaea") #arctic_evergreen_shrub_dwarf

PFT8 <- c("Salix_rotundifolia") #arctic_deciduous_shrub_dwarf

PFT12 <- c("Cardamine_pratensis", "Cochlearia_officinalis", "Petasites_frigidus", "Ranunculus_nivalis", "Ranunculus_pallasii",
           "Saxifraga_cernua", "Saxifraga_foliolosa", "Saxifraga_hirculus", "Stellaria_sp.") #arctic forb


PFT13 <- c("Arctoagrostis_latifolia", "Alopercurus_alpinus", "Poa_arctica", "Poa_vivaparum", "Luzula_arctica",
           "Luzula_confusa") #arctic_dry_graminoid

PFT14 <- c("Carex_aquatilis", "Dupontia_fisheri", "Eriophorum_angustifolium", "Eriophorum_sp.", "Juncus_biglumis") #arctic_wet_graminoid



obs_veg_cover_long <- obs_veg_cover_long %>%
  mutate(PFT = case_when(species %in% PFT1 ~ "PFT1",
                         species %in% PFT2 ~ "PFT2",
                         species %in% PFT3 ~ "PFT3",
                         species %in% PFT6 ~ "PFT6",
                         species %in% PFT8 ~ "PFT8",
                         species %in% PFT12 ~ "PFT12",
                         species %in% PFT13 ~ "PFT13",
                         species %in% PFT14 ~ "PFT14"))


#now group by plot_ID and PFT and sum to get the total fractional coverage of each PFT for each observation 
#have to sum before avg'ing b/c going from a ton of veg cover options down to only 4
obs_veg_cover_long <- obs_veg_cover_long %>%
  group_by(plot_ID, PFT) %>%
  mutate(fractional_cover_BAM = sum(fractional_cover, na.rm = TRUE))


obs_veg_cover_long <- obs_veg_cover_long %>%
  group_by(polygon_type, PFT) %>%
  summarize(fractional_cover_avg = mean(fractional_cover_BAM, na.rm = TRUE))

#coverages for polygon types from the obs data are all slightly over 100%, since they were summed and avg'd across multiple polygons 
#that were evaluated for each polygon type
#scale these to sum to 100 by calculating an adjustment factor
obs_veg_cover_long <- obs_veg_cover_long %>%
  group_by(polygon_type) %>%
  mutate(total_coverage = sum(fractional_cover_avg),
         adjustment_factor = 100 / total_coverage,
         adjusted_fractional_coverage = fractional_cover_avg * adjustment_factor) %>%
  select(-total_coverage, -adjustment_factor) # Remove intermediate columns


#save
write.csv(obs_veg_cover_long, file.path(veg_save, "Default_PFT_obs_veg_fractional_coverage_datasetpolygondef_moss_lichen_seperate.csv"), row.names = FALSE)
write.csv(obs_veg_cover_long, file.path(veg_save, "Arctic_PFT_obs_veg_fractional_coverage_datasetpolygondef.csv"), row.names = FALSE)
write.csv(obs_veg_cover_long, file.path(veg_save, "Default_PFT_obs_veg_fractional_coverage_datasetpolygondef_lichenbareground.csv"), row.names = FALSE) #this version has lichen grouped w/bare ground and moss grouped w/grass

#pasted the model and obs PFT fractional coverage datasets together in excel since they were small, pull in now and plot
veg_cover_compare <- read.csv(file.path(veg_save, "model_obs_PFT_fractional_cover_compare_datasetpolygondef.csv"), header = TRUE)
veg_cover_compare <- read.csv(file.path(veg_save, "model_obs_PFT_fractional_cover_compare_datasetpolygondef_moss_lichen_seperate.csv"), header = TRUE)
veg_cover_compare <- read.csv(file.path(veg_save, "model_obs_PFT_fractional_cover_compare_datasetpolygondef_lichenbareground.csv"), header = TRUE)



ggplot(veg_cover_compare, aes(x = factor(polygon_type), y = fractional_cover)) +
  facet_grid(~data_type) +
  geom_bar(stat = "identity", aes(fill = PFT)) +
  #scale_fill_manual(values = pft_colors) +
  #scale_x_discrete(labels=c("Trough", "LCPcenter", "LCPtransition", "HCPcenter", "HCPtransition")) +
  labs(x = NULL, y = "Fractional Coverage", title = "PFT Fractional Coverage by Polygon Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        legend.title = element_blank(), axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 15), axis.title = element_text(size = 15)) 

# pft_colors <- c("#332288", "#6699cc", "#88ccee", "#44aa99", "#117733", "#78b33e", "#999933", "#888888",
#                  "#dcd300", "#ddcc77", "#c4451c", "#661100", "#cc6677", "#aa4499", "#967acc")

pft_colors <- c("Lichen" = "#ddcc77", "Moss" = "#999933", "PFT1" = "#888888", "PFT10" = "#117733",
                "PFT12" = "#967acc", "PFT13" = "#6699cc")

#when plotting w/different # of polygon types, need to make two plots otherwise looks shit
veg_cover_compare_model <- subset(veg_cover_compare, data_type == "model")
p1 <- ggplot(veg_cover_compare_model, aes(x = factor(polygon_type), y = fractional_cover)) +
  #facet_grid(~data_type) +
  geom_bar(stat = "identity", aes(fill = PFT)) +
  scale_fill_manual(values = pft_colors) +
  #scale_x_discrete(labels=c("Trough", "LCPcenter", "LCPtransition", "HCPcenter", "HCPtransition")) +
  labs(x = NULL, y = "Fractional Coverage", title = "Model") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        legend.title = element_blank(), axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 15), axis.title = element_text(size = 15)) 


veg_cover_compare_obs <- subset(veg_cover_compare, data_type == "obs")
p2 <- ggplot(veg_cover_compare_obs, aes(x = factor(polygon_type), y = fractional_cover)) +
  #facet_grid(~data_type) +
  geom_bar(stat = "identity", aes(fill = PFT)) +
  scale_fill_manual(values = pft_colors) +
  #scale_x_discrete(labels=c("Trough", "LCPcenter", "LCPtransition", "HCPcenter", "HCPtransition")) +
  labs(x = NULL, y = NULL, title = "Observations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        legend.title = element_blank(), axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 15), axis.title = element_text(size = 15)) 


plot_list <- list(p1, p2)

plot_final <- grid_arrange_shared_legend(plot_list[[1]], plot_list[[2]], 
                                         ncol = 2, nrow = 1, position='right')



#plot just the obs for arctic PFTs--------------

pft_colors <- c("PFT2" = "#ddcc77", "PFT3" = "#999933", "PFT1" = "#888888", "PFT6" = "#117733",
                "PFT8" = "#967acc", "PFT12" = "#6699cc", "PFT13" = "#78b33e", "PFT14" = "#cc6677")

ggplot(obs_veg_cover_long, aes(x = factor(polygon_type), y = adjusted_fractional_coverage)) +
  #facet_grid(~data_type) +
  geom_bar(stat = "identity", aes(fill = PFT)) +
  scale_fill_manual(values = pft_colors) +
  #scale_x_discrete(labels=c("Trough", "LCPcenter", "LCPtransition", "HCPcenter", "HCPtransition")) +
  labs(x = NULL, y = NULL, title = "Observations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        legend.title = element_blank(), axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 15), axis.title = element_text(size = 15)) 




#Check out organic matter by soil layer-------------------------------------------------------------------------------------------

#start by pulling in what the current model setup is (ORG)
polygon_type <- c("Trough", "LCPcenter", "LCPtransition", "HCPcenter", "HCPtransition", "Rim", "Mean")
org_matter <- data.frame(polygon_type, ORG)
names(org_matter) <- c("polygon_type", "SL1", "SL2", "SL3", "SL4", "SL5", "SL6", "SL7", "SL8", "SL9", "SL10")

#convert to long format
org_matter_long <- org_matter %>%
  pivot_longer(cols = SL1:SL10, names_to = "soil_layer", values_to = "value")

#drop the "SL" before soil layer values
org_matter_long <- org_matter_long %>%
  mutate(soil_layer = str_replace(soil_layer, "^[A-Za-z]{2}", ""))

org_matter_long$soil_layer <- as.numeric(org_matter_long$soil_layer)
#org_matter_long$soil_layer <- factor(org_matter_long$soil_layer, levels = c("1.04", "0.62", "0.37", "0.21", "0.12", "0.06", "0.03", "0.01"))

fill_name <- "OM Density \n(kg m-3)"
#need to convert soil depth to factor so it extends the width of the line and fills out the plot
g1 <- ggplot(data = org_matter_long, aes(x = polygon_type, y = soil_layer, z = value))  +
  #scale_y_continuous(breaks = seq(0, 10, 1)) + #breaks = seq(0, 3, 0.5)
  #scale_y_reverse() +
  #facet_grid(polygon_type~.) +
  scale_y_continuous(trans = "reverse", breaks = seq(1, 10, 1)) +
  #scale_y_discrete(breaks = levels(subset_data$soil_depth)[seq(1, length(levels(subset_data$soil_depth)), by = 2)]) +
  #scale_x_date(date_labels = "%b") +  # Format to display only month
  labs(x = NULL, y = "Soil Layer",
       fill = fill_name, title = "Model")

p1 <- g1 + geom_tile(aes(fill = value)) + #color= "white",size=0.1, alpha = .7, 
  #scale_fill_viridis_c(option = "mako") +
  scale_fill_viridis(option ="C", limits=c(rng[1], rng[2])) +
  theme_minimal() +
  theme(axis.text = element_text(size = 12), #axis.title.x = element_blank(), axis.title.y = element_blank(), 
        legend.text = element_text(size = 13), legend.title = element_text(size = 13),
        title = element_text(size = 13), strip.text.y = element_text(size = 10), legend.position = "right")


p1



#pull in the obs data that has organic matter and soil bulk density by layer and polygon type---------------------
soil_phys_chem <- read.csv(file.path(soil_phys_chem_path, "core_physical_chem_data_20180321.csv"), skip = 5, header = TRUE) #has a few metadata lines at top to be skipped
soil_phys_chem <- soil_phys_chem[-1,] #drop the first row, units

#has column called "Geomorph_Feature" that is polygon type (options HCP, LCP, FCP) and col "Polygon_feature" that has options "Center", "Rim", "Trough" for each polygon type
#change these to lower case, then combine the two cols
soil_phys_chem <- soil_phys_chem %>%
  mutate(Polygon_feature = case_when(Polygon_feature == "Center" ~ "center",
                                     Polygon_feature == "Rim" ~ "rim",
                                     Polygon_feature == "Trough" ~ "trough")) 

soil_phys_chem <- soil_phys_chem %>%
  mutate(polygon_type = paste0(Geomorph_Feature, Polygon_feature))

unique(soil_phys_chem$polygon_type)
#hmmmm ok so for the chemical/physical data there's only measurements for "LCPcenter", "FCPcenter", "FCPrim", "HCPcenter", "HCPtrough"





org_matter_obs <- data.frame(soil_phys_chem$polygon_type, soil_phys_chem$Mid_Depth_Core_Section, soil_phys_chem$Dry_Bulk_Density, 
                             soil_phys_chem$OM, soil_phys_chem$Organic_Matter_Content)
names(org_matter_obs) <- c("polygon_type", "soil_layer", "bulk_dens", "OM_vol", "OM_wt")
org_matter_obs[org_matter_obs == -9999.00] <- NA
org_matter_obs$soil_layer <- as.numeric(org_matter_obs$soil_layer)
org_matter_obs$bulk_dens <- as.numeric(org_matter_obs$bulk_dens)
org_matter_obs$OM_vol <- as.numeric(org_matter_obs$OM_vol)
org_matter_obs$OM_wt <- as.numeric(org_matter_obs$OM_wt)
#soil layer is really still soil depth, units are cm, bulk_dens is dry bulk density in units of g/cm3, OM_vol is organic matter content as volume %
#convert from cm to m
org_matter_obs$soil_layer <- org_matter_obs$soil_layer/100
length(unique(org_matter_obs$soil_layer)) #102 unique soil depth values, ranging from 0.025m - 3.835m

#model soil depths corresponding to the first 10 layers are "2.86", "1.73", "1.04", "0.62", "0.37", "0.21", "0.12", "0.06", "0.03", "0.01" (listed from deepest to shallowest)
#start by cutting off anything deeper than 3m
org_matter_obs <- subset(org_matter_obs, soil_layer <= 3)

#convert bulk density from g/cm3 to kg/m3
org_matter_obs$bulk_dens <- org_matter_obs$bulk_dens * 1000 #1 g/cm3 is equal to 1000 kg/m3

# Convert volume percent of organic matter to kg/m3
#org_matter_obs$OM_density <- (org_matter_obs$OM_vol/100) * org_matter_obs$bulk_dens
#see slides, multiply obs OM as %vol by 130 kg/m3 (model default density value it uses internally to get to a volume fraction)
org_matter_obs$OM_density <- (org_matter_obs$OM_wt/100) * 130

#looks like a bunch of the FCPcenter values are NA, remove any rows where OM_density is NA
org_matter_obs <- org_matter_obs %>%
  filter(!is.na(OM_density))


#next group depths into layer bins, following how the model layer to depth relationships are set up
# org_matter_obs <- org_matter_obs %>%
#   mutate(soil_layer_test = case_when(soil_layer <= 0.03 ~ 1,
#                                      soil_layer > 0.03 & soil_layer <= 0.06 ~ 2,
#                                      soil_layer > 0.06 & soil_layer <= 0.12 ~ 3,
#                                      soil_layer > 0.12 & soil_layer <= 0.21 ~ 4,
#                                      soil_layer > 0.21 & soil_layer <= 0.37 ~ 5,
#                                      soil_layer > 0.37 & soil_layer <= 0.62 ~ 6,
#                                      soil_layer > 0.62 & soil_layer <= 1.04 ~ 7,
#                                      soil_layer > 1.04 & soil_layer <= 1.73 ~ 8,
#                                      soil_layer > 1.73 & soil_layer <= 2.86 ~ 9,
#                                      soil_layer > 2.86 ~ 10))

#using the above definitions there is no layer 1...
# org_matter_obs <- org_matter_obs %>%
#   mutate(soil_layer_test = case_when(soil_layer <= 0.03 ~ 1,
#                                      soil_layer > 0.03 & soil_layer <= 0.06 ~ 2,
#                                      soil_layer > 0.06 & soil_layer <= 0.12 ~ 3,
#                                      soil_layer > 0.12 & soil_layer <= 0.21 ~ 4,
#                                      soil_layer > 0.21 & soil_layer <= 0.37 ~ 5,
#                                      soil_layer > 0.37 & soil_layer <= 0.62 ~ 6,
#                                      soil_layer > 0.62 & soil_layer <= 1.04 ~ 7,
#                                      soil_layer > 1.04 & soil_layer <= 1.73 ~ 8,
#                                      soil_layer > 1.73 & soil_layer <= 2.86 ~ 9,
#                                      soil_layer > 2.86 ~ 10))

plot(org_matter_obs$soil_layer_test, org_matter_obs$soil_layer)
#using the above definitions there's no layer 2...very few obs at surface layers, lots at deep layers
#in the model the same values are used for layer 1 and 2, just do this, set up below then copy the layer 1 values and call it layer 2
# org_matter_obs <- org_matter_obs %>%
#   mutate(soil_layer_test = case_when(soil_layer <= 0.06 ~ 1,
#                                      soil_layer > 0.06 & soil_layer <= 0.12 ~ 3,
#                                      soil_layer > 0.12 & soil_layer <= 0.21 ~ 4,
#                                      soil_layer > 0.21 & soil_layer <= 0.37 ~ 5,
#                                      soil_layer > 0.37 & soil_layer <= 0.62 ~ 6,
#                                      soil_layer > 0.62 & soil_layer <= 1.04 ~ 7,
#                                      soil_layer > 1.04 & soil_layer <= 1.73 ~ 8,
#                                      soil_layer > 1.73 & soil_layer <= 2.86 ~ 9,
#                                      soil_layer > 2.86 ~ 10))
# unique(org_matter_obs$soil_layer_test)
# 
# layer2 <- subset(org_matter_obs, soil_layer_test == 1) #only LCP center has a value for layer 1...
# 
# org_matter_obs <- org_matter_obs %>%
#   mutate(soil_layer_test = case_when(soil_layer <= 0.12 ~ 1, #basically going to have to repeat values for layers 1-3
#                                      soil_layer > 0.12 & soil_layer <= 0.21 ~ 4,
#                                      soil_layer > 0.21 & soil_layer <= 0.37 ~ 5,
#                                      soil_layer > 0.37 & soil_layer <= 0.62 ~ 6,
#                                      soil_layer > 0.62 & soil_layer <= 1.04 ~ 7,
#                                      soil_layer > 1.04 & soil_layer <= 1.73 ~ 8,
#                                      soil_layer > 1.73 & soil_layer <= 2 ~ 9, #making layer 9 a narrower interval
#                                      soil_layer > 2 ~ 10))
# 
# layer1 <- subset(org_matter_obs, soil_layer_test == 1)
# #layer 1 still only has data for LCPcenter, FCPcenter, FCPrim, and HCPcenter...pick which the other polygon type/feature combos are most similar to and use that?
# #LCPrim, LCPtrough, FCPtrough, HCPrim, HCPtrough
# #k no don't do that, only "LCPcenter, "FCPcenter“, "FCPrim“,  "HCPcenter“, and "HCPtrough“ have any data in the whole dataset
# 
# #try just breaking up into 10 even layers
# org_matter_obs <- org_matter_obs %>%
#   mutate(soil_layer_test = case_when(soil_layer <= 0.3 ~ 1, 
#                                      soil_layer > 0.3 & soil_layer <= 0.6 ~ 2,
#                                      soil_layer > 0.6 & soil_layer <= 0.9 ~ 3,
#                                      soil_layer > 0.9 & soil_layer <= 1.20 ~ 4,
#                                      soil_layer > 1.20 & soil_layer <= 1.50 ~ 5,
#                                      soil_layer > 1.50 & soil_layer <= 1.80 ~ 6,
#                                      soil_layer > 1.80 & soil_layer <= 2.10 ~ 7,
#                                      soil_layer > 2.10 & soil_layer <= 2.40 ~ 8,
#                                      soil_layer > 2.40 & soil_layer <= 2.70 ~ 9, #making layer 9 a narrower interval
#                                      soil_layer > 2.70 ~ 10))

#shitty, trying breaking up based on difference in size of depth intervals (see notebook)
# org_matter_obs <- org_matter_obs %>%
#   mutate(soil_layer_test = case_when(soil_layer <= 0.07 ~ 1, 
#                                      soil_layer > 0.07 & soil_layer <= 0.13 ~ 2,
#                                      soil_layer > 0.13 & soil_layer <= 0.21 ~ 3,
#                                      soil_layer > 0.21 & soil_layer <= 0.33 ~ 4,
#                                      soil_layer > 0.33 & soil_layer <= 0.60 ~ 5,
#                                      soil_layer > 0.60 & soil_layer <= 1.09 ~ 6,
#                                      soil_layer > 1.09 & soil_layer <= 1.60 ~ 7,
#                                      soil_layer > 1.60 & soil_layer <= 2.08 ~ 8,
#                                      soil_layer > 2.08 & soil_layer <= 2.70 ~ 9, #making layer 9 a narrower interval
#                                      soil_layer > 2.70 ~ 10))

#same approach as above, redoing after switching to determining OM density from weight% using model default density value
org_matter_obs <- org_matter_obs %>%
  mutate(soil_layer_test = case_when(soil_layer <= 0.03 ~ 1, 
                                     soil_layer > 0.03 & soil_layer <= 0.08 ~ 2,
                                     soil_layer > 0.08 & soil_layer <= 0.18 ~ 3,
                                     soil_layer > 0.18 & soil_layer <= 0.30 ~ 4,
                                     soil_layer > 0.30 & soil_layer <= 0.48 ~ 5,
                                     soil_layer > 0.48 & soil_layer <= 0.90 ~ 6,
                                     soil_layer > 0.90 & soil_layer <= 1.60 ~ 7,
                                     soil_layer > 1.60 & soil_layer <= 2.08 ~ 8,
                                     soil_layer > 2.08 & soil_layer <= 2.70 ~ 9, #making layer 9 a narrower interval
                                     soil_layer > 2.70 ~ 10))

#hmmmm try just seeing what the OM density by depth looks like for obs

fill_name <- "OM Density \n(kg m-3)"
#need to convert soil depth to factor so it extends the width of the line and fills out the plot
g2 <- ggplot(data = org_matter_obs, aes(x = polygon_type, y = soil_layer, color = OM_density))+ #, z = OM_density))  +
  #scale_y_continuous(breaks = seq(0, 10, 1)) + #breaks = seq(0, 3, 0.5)
  scale_y_reverse() +
  #facet_grid(polygon_type~.) +
  #scale_y_continuous(trans = "reverse", breaks = seq(1, 10, 1)) +
  #scale_y_discrete(breaks = levels(subset_data$soil_depth)[seq(1, length(levels(subset_data$soil_depth)), by = 2)]) +
  #scale_x_date(date_labels = "%b") +  # Format to display only month
  labs(x = NULL, y = "Soil Depth",
       fill = fill_name, title = "Observations")

p2 <- g2 + geom_point(size = 3) + #color= "white",size=0.1, alpha = .7, 
  #scale_fill_viridis_c(option = "mako") +
  scale_color_viridis(option ="C") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12), #axis.title.x = element_blank(), axis.title.y = element_blank(), 
        legend.text = element_text(size = 13), legend.title = element_text(size = 13),
        title = element_text(size = 13), strip.text.y = element_text(size = 10), legend.position = "right")


p2


#geom tile version
fill_name <- "OM Density \n(kg m-3)"
#need to convert soil depth to factor so it extends the width of the line and fills out the plot
g2 <- ggplot(data = org_matter_obs, aes(x = polygon_type, y = soil_layer_test, fill = OM_density))+ #, z = OM_density))  +
  #scale_y_continuous(breaks = seq(0, 10, 1)) + #breaks = seq(0, 3, 0.5)
  #scale_y_reverse() +
  #facet_grid(polygon_type~.) +
  scale_y_continuous(trans = "reverse", breaks = seq(1, 10, 1)) +
  #scale_y_discrete(breaks = levels(subset_data$soil_depth)[seq(1, length(levels(subset_data$soil_depth)), by = 2)]) +
  #scale_x_date(date_labels = "%b") +  # Format to display only month
  labs(x = NULL, y = "Soil Layer",
       fill = fill_name, title = "Observations")

p2 <- g2 + geom_tile(aes(fill = OM_density)) + #color= "white",size=0.1, alpha = .7, 
  #scale_fill_viridis_c(option = "mako") +
  scale_fill_viridis(option ="C", limits=c(rng[1], rng[2])) +
  theme_minimal() +
  theme(axis.text = element_text(size = 12), #axis.title.x = element_blank(), axis.title.y = element_blank(), 
        legend.text = element_text(size = 13), legend.title = element_text(size = 13),
        title = element_text(size = 13), strip.text.y = element_text(size = 10), legend.position = "right")


p2

#set shared colorbar range and plot together so comparable (just go back up and edit geom_tile plots to add this)
rng <- range(c((org_matter_obs$OM_density), (org_matter_long$value)))
#rng <- range(c((org_matter_obs$OM_vol), (org_matter_long$value)))


plot_list <- list(p1, p2)

plot_final <- grid_arrange_shared_legend(plot_list[[1]], plot_list[[2]], 
                                         ncol = 2, nrow = 1, position='right')




#look at canopy height data...hand measured obs are only from July 2012, while model needs monthly avg values, but can check if in reasonable range------------------------------

#look at how model is set up 
plant_height <- ncvar_get(surf_data, "MONTHLY_HEIGHT_TOP")
#formatted as a slice for each month, rows are gridcells cols are PFTs
#says monthly but looks like the same values are repeated for each month slice
#same values also repeated for each grid cell
#guessing these are just default values b/c most of these PFTs aren't even present in these sims, only 10, 12, 13 exist in these sims
#[,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15] [,16] [,17]
#0   17   17   14   35   35   18   20   20   0.5   0.5   0.5   0.5   0.5   0.5   0.5   0.5

#so PFT 10 = 0.5m, PFT12 = 0.5m, PFT13 = 0.5m

#if looking just at "MONTHLY_HEIGHT_TOP" compared to obs measurements of tallest leaf height, "MONTHLY_HEIGHT_TOP" are way higher, think need to subtract
#bottom height and get actual plant height 
plant_height_b <- ncvar_get(surf_data, "MONTHLY_HEIGHT_BOT")
#k seems to be 0.1m for PFT10 and 12, 0.01m for PFT13, so won't make huge diff but doing anyways, so then PFT10 and PFT 12 have height 0.4m, PFT13 has height 0.49m


#obs canopy height data from the 2012 intensive surveys (V_1_2_vegetation_canopy_height_Barrow_2012_v1.csv) aren't broken up by PFT, but the leaf height
#file (V_1_2_species_leaf_height_Barrow_2012_v1.csv) are broken up by species, so use that
#"Height of tallest leaf of species within canopy in ten random locations within the 1 x 1 m vegetation survey plots.
#Where a species was insufficiently abundant for ten measurements, five were made."

plant_height_obs <- read.csv(file.path(veg_data, "V_1_2_species_leaf_height_Barrow_2012_v1.csv"), header = TRUE)
#heights are in cm
plant_height_obs <- plant_height_obs[-1,] #drop first row, units
#convert from cm to m
plant_height_obs$leaf_height <- as.numeric(plant_height_obs$leaf_height)
plant_height_obs$leaf_height <- plant_height_obs$leaf_height/100

#first create a new polygon_type col that is same format I've been using
#check polygon type names first
unique(plant_height_obs$polygon_type)
#[1] "low-centered"  "high-centered" "transitional"
plant_height_obs <- plant_height_obs %>%
  mutate(polygon_type = case_when(polygon_type == "low-centered" ~ "LCP",
                                  polygon_type == "high-centered" ~ "HCP",
                                  polygon_type == "transitional" ~ "FCP")) #"transition"

unique(plant_height_obs$polygon_sub_unit)
#polygon subunits are labelled differently here
plant_height_obs <- plant_height_obs %>%
  mutate(polygon_sub_unit = case_when(polygon_sub_unit == "C" ~ "center",
                                      polygon_sub_unit == "E" ~ "rim", #edge is same as rim
                                      polygon_sub_unit == "T" ~ "trough"))


plant_height_obs <- plant_height_obs %>%
  mutate(polygon_type = paste0(polygon_type, polygon_sub_unit))


#use PFT-species groupings from above to add a col for PFT
#jk there are way fewer and something is throwing an error...one must not exist in the previous groupings...mmm nvm notation is diff, spaces instead of _
PFT1 <- c() #grouping moss and lichen in w/bare ground

PFT10 <- c("Vaccinium vitis-idaea", "Salix rotundifolia") #shrub

PFT12 <- c("Petasites frigidus", "Stellaria sp. ", "Saxifraga cernua", "Saxifraga foliolosa", "Ranunculus pallasii") #arctic forb

PFT13 <- c("Carex aquatilis", "Eriophorum sp.", "Dupontia fisheri", "Eriophorum angustifolium", "Eriophorum russeolum", "Luzula arctica",          
           "Arctagrostis latifolia", "Luzula confusa", "Poa arctica", "Eriophorum vaginatum", "Juncus biglumis") #arctic_dry_graminoid


plant_height_obs <- plant_height_obs %>%
  mutate(PFT = case_when(species %in% PFT1 ~ "PFT1",
                         species %in% PFT10 ~ "PFT10",
                         species %in% PFT12 ~ "PFT12",
                         species %in% PFT13 ~ "PFT13"))


#group by PFT and avg
plant_height_obs_summary <- plant_height_obs %>%
  group_by(polygon_type, PFT) %>%
  summarise(height_mean = mean(leaf_height, na.rm = TRUE))

#save, then open as csv and paste in the model input values
write.csv(plant_height_obs_summary, file.path(veg_save, "Default_PFT_obs_veg_height_datasetpolygondef.csv"), row.names = FALSE)

plant_height_compare <- read.csv(file.path(veg_save, "Default_PFT_obs_veg_height_compare.csv"), header = TRUE)



ggplot(plant_height_compare, aes(x = factor(polygon_type), y = height_mean, color = PFT, shape = data_type)) +
  #facet_grid(~data_type) +
  geom_point(size = 3) +
  #scale_fill_manual(values = pft_colors) +
  #scale_x_discrete(labels=c("Trough", "LCPcenter", "LCPtransition", "HCPcenter", "HCPtransition")) +
  labs(x = NULL, y = "Canopy Height (m)", title = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        legend.title = element_blank(), axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 15), axis.title = element_text(size = 15)) 




