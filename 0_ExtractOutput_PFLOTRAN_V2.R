
#Created: 4/26/24
#Author: Bailey Murphy

#pulling in output from ELM-PFLOTRAN simulations Ben ran at Barrow, did a 7 cell transect to capture the different polygonal tundra landforms
#uses the default ELM arctic PFTs, not set up w/the new PFTs yet

#5/1/24: updating to only pull first 10 soil layers, Ben said beyond that soil biogeochemistry isn't really defined, only hydrology/temp stuff is defined


#-----------------------------------------------------------------------------------------------------------------------------

#Code inputs:


#Code outputs:

#---------------------------------------------------------------------------------------------------------------------------

#----------------
#install/load required packages
#----------------


packages <- c("ncdf4", "dplyr", "stringr", "readbulk", "lubridate", "tidyr", "magrittr", "ggplot2", "ggridges",
              "reshape2", "lemon") 
#"magrittr" is for pipes

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages],
                   repos="http://cran.us.r-project.org")
}

invisible(lapply(packages, library, character.only = TRUE))

#----------------
#Hard-coded variables
#----------------

site_name <- "BEO" 

run_type <- "Arctic"#specify whether using model output from simulations with the new Arctic PFTs or w/default ELM PFTs
output_yr <- "1997"
start_yr <- "1995" 
end_yr <- "1999" 

dat.base_main <- "~/NGEE_modeling/ELM_PFLOTRAN_BEO"
dat.base <- paste0("~/NGEE_modeling/ELM_PFLOTRAN_BEO/1_Runs/", run_type, "_PFT/BEO_7cell_RAW_h0") #this is where the raw model output is stored
path.out <- paste0("~/NGEE_modeling/ELM_PFLOTRAN_BEO/1_Runs/", run_type, "_PFT/BEO_7cell_EXTRACTED") 
run_list_path <- paste0("~/NGEE_modeling/ELM_PFLOTRAN_BEO/1_Runs/", run_type, "_PFT") 

setwd(run_list_path)
run.list <- read.csv(paste0("RUNID_list_", site_name, "_7cell_h0_", run_type, ".csv"), header = FALSE) #UPDATE THIS WHEN USING FULL RUNID LIST


#open connection to the file, want to poke around and see what variables exist/are of interest to extract
run.list <- c(run.list$V1)
RUNID <- run.list[3]
#model_data <- ncdf4::nc_open(file.path(dat.base, RUNID))

attributes(model_data$dim)
names(model_data$var)  #this will list the names of all the variables
#summary(model_data$var) # This will list what variables we can extract, redundant to 'names', but tells you length
#nc.close(model_data)

#only run this the first time through---------
# print(model_data) #this shows all the variables, their long name, units, data type, etc.
# 
# # Save the print(nc) dump to a text file, then can quickly access info for all vars
# {
#   sink('Alaska_alquimia_7cell_AK-BEO.h0_metadata.txt')
#   print(model_data)
#   sink()
# }

#---------------------------------------------

#some variables are vertically resolved by soil layer so they have a third dimension to the data, extract these separately
#read list of 2D variables to extract:
variable_list_2D <- read.csv(file.path(run_list_path, "ELM-PFLOTRAN_2D_variables_h0.csv"), header = FALSE)
variable_list_2D <- c(variable_list_2D$V1)

variable_list_3D <- read.csv(file.path(run_list_path, "ELM-PFLOTRAN_3D_variables_h0.csv"), header = FALSE) 
variable_list_3D <- c(variable_list_3D$V1)


#----------------
#Load & process data 
#----------------


print(paste0("now processing: ", RUNID))

run.splt <- stringr::str_split(RUNID, ".h0")[[1]] # Splits apart the run name to extract different bits of info
run.splt <- stringr::str_split(run.splt[2], "-")[[1]]

model_data <- ncdf4::nc_open(file.path(dat.base, RUNID)) 

#first pull out lat, lon, and time
lon <- ncvar_get(model_data, "lon") #length is 7
lat <- ncvar_get(model_data, "lat", verbose = F)
date <- ncvar_get(model_data, "mcdate") #format is YYYYMMDD, length 8760, so hourly data for a year 
sec <- ncvar_get(model_data, "mcsec") #current seconds of current day, starts at 0 goes to 82800, 24 entries per day

#Create 2D matrix of long and date (doing lon, lat, AND date will take 10 million years)
lon_date <- as.matrix(expand.grid(lon, date)) 

#then do the same for lat and sec
lat_date <- as.matrix(expand.grid(lat, sec))

#combine the two and remove the duplicate date col
lonlatdate <- data.frame(cbind(lon_date[,1], lat_date[,1], 
                               lon_date[,2], lat_date[,2]))
#clean up
rm(lon_date)
rm(lat_date)
#rm(lon)
#rm(lat)
#rm(date)


#set up empty df to store variable vectors, needs to match dimensions of what the variable vectors will be
var_df_2D <- data.frame(matrix(ncol = 1, nrow = 61320)) #61320 comes from 7 grid cells, 8760 entries (hourly)

print("Now looping through variables")

#LOOP THROUGH EACH VARIABLE WITHIN EACH FILE STARTS HERE--------------------------------
for(variable_name in variable_list_2D){
  #Read in the data from the ecosystem variable  
  var.array <- ncvar_get(model_data, variable_name) # store the data in a 2-dimensional array...lat/lon aren't own dimensions here, its time and land unit grid cell as dims
  
  fillvalue <- ncatt_get(model_data, variable_name, "_FillValue")
  
  #replace fill values w/NA
  var.array[var.array == fillvalue$value] <- NA
  
  #collapse the variable into a single vector
  var_vec <- as.vector(var.array)
  
  #add to var_df
  var_df_2D <- cbind(var_df_2D, var_vec)
  
}

#remove 1st dummy col of df
var_df_2D[1] <- NULL

#now combine everything into 1 df
full_df_2D <- data.frame(cbind(lonlatdate, var_df_2D))

#clean up
rm(lonlatdate)
rm(var_df_2D)
rm(var.array)
rm(var_vec)

#update column names
colnames_list <- c("lon", "lat", "date", "sec")
colnames_list <- append(colnames_list, variable_list_2D)
#colnames_list <- append(colnames_list, "WATSAT") #add variable name for the porosity vector that got added on
colnames(full_df_2D) <- colnames_list

#THIS IS WRONG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# #this time forgot to add H2OSFC_TIDE to the list of variables to output...should be in order though, 7 grid cells 
# polygon_type <- c("LCPtrough", "LCPcenter", "FCPtrough", "FCPcenter", "LCPrim", "HCPtrough", "HCPcenter")
# 
# #full_df_2D has 61320 rows, so polygon types need to be repeated 8760 times
# polygon_type <- rep(polygon_type, length.out = 8760)
# 
# full_df_2D$polygon_type <- polygon_type


#12/5/24: edit this once I actually output H2OSFC_TIDE
# #add polygon type column, based on tide height above soil surface (7 unique)
# full_df_2D$H2OSFC_TIDE <- as.character(full_df_2D$H2OSFC_TIDE) #convert to char first so rounding doesn't mess up classification
# 
# full_df_2D <- full_df_2D %>%
#   mutate(polygon_type = case_when(H2OSFC_TIDE == "100" ~ "Trough",
#                                   H2OSFC_TIDE == "59.1641693115234" ~ "LCPcenter",
#                                   H2OSFC_TIDE == "62.7643089294434" ~ "LCPtransition",
#                                   H2OSFC_TIDE == "-115.806213378906" ~ "HCPcenter",
#                                   H2OSFC_TIDE == "-121.924034118652" ~ "HCPtransition",
#                                   H2OSFC_TIDE == "-18.96484375" ~ "Rim",
#                                   H2OSFC_TIDE == "20.1018199920654" ~ "Mean"))
# 
# full_df_2D$H2OSFC_TIDE <- as.numeric(full_df_2D$H2OSFC_TIDE) #switch back to numeric






#now repeat the same thing for the 3D variables (the ones that are vertically resolved)------------------------------------------------------
#set up empty df to store variable vectors, needs to match dimensions of what the variable vectors will be
#var_df_3D <- data.frame(matrix(ncol = 1, nrow = 919800)) #919800 comes from 7 grid cells, 8760 entries (hourly), and 15 soil layers

#only doing first 10 soil layers now, but being lazy and just chopping off the additional layers after

# Generate repetitions of location, soil layer, and polygon type
#polygon_type <- c("Trough", "LCPcenter", "LCPtransition", "HCPcenter", "HCPtransition", "Rim", "Mean")
polygon_type <- c("LCPtrough", "LCPcenter", "FCPtrough", "FCPcenter", "LCPrim", "HCPtrough", "HCPcenter")
#polygon_type <- rep(polygon_type, each = 15)
#lat <- rep(lat, each = 15)
#lon <- rep(lon, each = 15)
# soil_layers <- rep(1:10)
# #should all be length 105 if looking at all 15 layers
# length(polygon_type)
# length(lat)
# length(lon)
# length(soil_layers)

var_df_3D <- data.frame(matrix(ncol = 8, nrow = 0)) #set up empty df to store the vertically resolved variables once extracted

#set column names, so can rbind df generated for each variable that gets looped through
#colnames(var_df_3D) <- c("dummy", "value", "variable", "lat", "lon", "polygon_type", "soil_layer", "date", "sec")
colnames(var_df_3D) <- c("polygon_type_drop", "soil_layer", "time_drop", "value", "variable", "polygon_type", "date", "sec")


for(variable_name in variable_list_3D){ 
  array_data <- ncvar_get(model_data, variable_name) # store the data in an array
  
  # Subset the array to the first 10 soil layers
  array_data <- array_data[, 1:10, ]
  
  fillvalue <- ncatt_get(model_data, variable_name, "_FillValue")
  
  #replace fill values w/NA
  array_data[array_data == fillvalue$value] <- NA
  
  # Reshape the array to long format
  long_data <- melt(array_data, varnames = c("polygon_type_drop", "soil_layer", "time_drop"))
  long_data$variable <- variable_name
  
  long_data <- long_data %>%
    mutate(
      polygon_type = rep(polygon_type, times = 87600), #repeats 1:7 87600 times, format "1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4...", total length 613200
      date = rep(date, each = 70), #repeats first date 70 times before moving on to second date, total length 613200
      sec = rep(sec, each = 70)) #repeats first time 70 times before moving on to second time, total length 613200

#add to var_df
var_df_3D <- rbind(var_df_3D, long_data)

}

var_df_3D <- var_df_3D[,-1]
var_df_3D <- var_df_3D[,-2]


rm(long_data)
rm(array_data)
rm(date)
rm(sec)
rm(lat)
rm(lon)

#fix later to pull this in earlier, but pulling in the actual soil depth info instead of just using layer number
soil_layer_depth <- model_data[["dim"]][["levdcmp"]]$vals #depths for each layer are in m, same values for "levdcmp" and "levgrnd"
#round them so they're cleaner
soil_layer_depth <- round(soil_layer_depth, 2)

print("Closing netcdf file")
ncdf4::nc_close(model_data) # Closing the connection to the netcdf file; this is important otherwise your computer will get confused
rm(model_data)



#----------------------------------------------------------------------------------

#do unit conversions for 2D variables, pretty much just converting mm water to m
#jk only convert the water height variables, runoff and drainage numbers are extremely tiny if in m
full_df_2D$H2OSFC <- full_df_2D$H2OSFC / 1000
full_df_2D$H2OSFC_TIDE <- full_df_2D$H2OSFC_TIDE / 1000 #commenting out for now...

#full_df_2D$QDRAI <- full_df_2D$QDRAI / 1000
#full_df_2D$QFLX_EVAP_TOT <- full_df_2D$QFLX_EVAP_TOT / 1000
#full_df_2D$QFLX_LAT_AQU <- full_df_2D$QFLX_LAT_AQU / 1000
#full_df_2D$QVEGT <- full_df_2D$QVEGT / 1000

#unit conversions for 3D variables
var_df_3D$value <- ifelse(var_df_3D$variable == "TSOI", var_df_3D$value - 273.15, var_df_3D$value) #convert temp from K to C
var_df_2D$value <- ifelse(var_df_2D$variable == "TSA", var_df_2D$value - 273.15, var_df_2D$value) #convert temp from K to C, came back to add later, forgot

#DON'T do this!! think this is what made VWC so small, just leave it in mm
#water_vars <- c("H2OSOI", "QDRAI_VR", "QFLX_ADV")
#var_df_3D$value <- ifelse(var_df_3D$variable %in% water_vars, var_df_3D$value / 1000, var_df_3D$value) #convert water related variables from mm to m

#want to have versions of the chem variables that are gCm-3 AND mmolL-1 water so can have values that are useful for comparison to obs, and vals that are useful for thinking ~C balances
chem_vars <- c("soil_Fe2", "soil_FeOxide", "soil_FeS", "soil_O2", "soil_sulfate", "soil_sulfide", "soil_acetate")
#var_df_3D$value <- ifelse(var_df_3D$variable %in% chem_vars, var_df_3D$value * 1000, var_df_3D$value) #convert chemistry related variables from mol/m3 to mmol/m3
#convert from mol/m-3 of soil to mmol L-1 of water, #hold off on this:so can express as micromolar concentration (uM of substance per L of water)
#to do this, divide the value expressed as mol/m-3 of soil by the water volume (the porosity, so 'watsat')
#(mol/m3 soil) / (m3 water/m3 soil) = mol/m3 water = mmol/L water
WATSAT <- subset(var_df_3D, variable == "watsat")
porosity <- WATSAT$value
rm(WATSAT)
#var_df_3D$value <- ifelse(var_df_3D$variable %in% chem_vars, var_df_3D$value / porosity, var_df_3D$value) #now chemistry variables expressed as mmol/L water
soil_Fe2_aq <- subset(var_df_3D, variable == "soil_Fe2")
soil_Fe2_aq$value <- soil_Fe2_aq$value / porosity
#change variable name and append to the main df
soil_Fe2_aq$variable <- "soil_Fe2_aq"
var_df_3D <- rbind(var_df_3D, soil_Fe2_aq)

soil_FeOxide_aq <- subset(var_df_3D, variable == "soil_FeOxide")
soil_FeOxide_aq$value <- soil_FeOxide_aq$value / porosity
#change variable name and append to the main df
soil_FeOxide_aq$variable <- "soil_FeOxide_aq"
var_df_3D <- rbind(var_df_3D, soil_FeOxide_aq)

soil_FeS_aq <- subset(var_df_3D, variable == "soil_FeS")
soil_FeS_aq$value <- soil_FeS_aq$value / porosity
#change variable name and append to the main df
soil_FeS_aq$variable <- "soil_FeS_aq"
var_df_3D <- rbind(var_df_3D, soil_FeS_aq)

soil_O2_aq <- subset(var_df_3D, variable == "soil_O2")
soil_O2_aq$value <- soil_O2_aq$value / porosity
#change variable name and append to the main df
soil_O2_aq$variable <- "soil_O2_aq"
var_df_3D <- rbind(var_df_3D, soil_O2_aq)

soil_sulfate_aq <- subset(var_df_3D, variable == "soil_sulfate")
soil_sulfate_aq$value <- soil_sulfate_aq$value / porosity
#change variable name and append to the main df
soil_sulfate_aq$variable <- "soil_sulfate_aq"
var_df_3D <- rbind(var_df_3D, soil_sulfate_aq)

soil_sulfide_aq <- subset(var_df_3D, variable == "soil_sulfide")
soil_sulfide_aq$value <- soil_sulfide_aq$value / porosity
#change variable name and append to the main df
soil_sulfide_aq$variable <- "soil_sulfide_aq"
var_df_3D <- rbind(var_df_3D, soil_sulfide_aq)

soil_acetate_aq <- subset(var_df_3D, variable == "soil_acetate")
soil_acetate_aq$value <- soil_acetate_aq$value / porosity
#change variable name and append to the main df
soil_acetate_aq$variable <- "soil_acetate_aq"
var_df_3D <- rbind(var_df_3D, soil_acetate_aq)


rm(soil_Fe2_aq)
rm(soil_FeOxide_aq)
rm(soil_FeS_aq)
rm(soil_O2_aq)
rm(soil_sulfate_aq)
rm(soil_sulfide_aq)
rm(soil_acetate_aq)

#normalize volumetric water content by porosity
H2OSOI <- subset(var_df_3D, variable == "H2OSOI")
H2OSOI$value <- H2OSOI$value / porosity #now this is in % form, the water content/total potential water content, 1 = fully saturated, 0.7 = 70% sat
#change variable name and append to the main df
H2OSOI$variable <- "VWC"
var_df_3D <- rbind(var_df_3D, H2OSOI)

#calc frozen fraction of soil
SOILICE <- subset(var_df_3D, variable == "SOILICE")
SOILLIQ <- subset(var_df_3D, variable == "SOILLIQ")

frozen_frac <- SOILICE$value / (SOILICE$value + SOILLIQ$value)
#change variable name and append to the main df
SOILICE$value <- frozen_frac
SOILICE$variable <- "frozen_frac"
var_df_3D <- rbind(var_df_3D, SOILICE)

rm(H2OSOI)
rm(frozen_frac)
rm(SOILICE)
rm(SOILLIQ)

#convert vertically resolved DOC, DIC, CH4 from gC/m^3 to mmol C/L (also keep in gC form though)
#carbon_vr_vars <- c("DIC_vr", "CH4_vr", "DOC_vr")
#var_df_3D$value <- ifelse(var_df_3D$variable %in% carbon_vr_vars, var_df_3D$value / porosity / 12.011, var_df_3D$value) #now carbon variables expressed as mmol C/L water
DIC_vr_aq <- subset(var_df_3D, variable == "DIC_vr")
DIC_vr_aq$value <- (DIC_vr_aq$value/12.01) / porosity
#change variable name and append to the main df
DIC_vr_aq$variable <- "DIC_vr_aq"
var_df_3D <- rbind(var_df_3D, DIC_vr_aq)

CH4_vr_aq <- subset(var_df_3D, variable == "CH4_vr")
CH4_vr_aq$value <- (CH4_vr_aq$value/12.01) / porosity
#change variable name and append to the main df
CH4_vr_aq$variable <- "CH4_vr_aq"
var_df_3D <- rbind(var_df_3D, CH4_vr_aq)

DOC_vr_aq <- subset(var_df_3D, variable == "DOC_vr")
DOC_vr_aq$value <- (DOC_vr_aq$value/12.01) / porosity
#change variable name and append to the main df
DOC_vr_aq$variable <- "DOC_vr_aq"
var_df_3D <- rbind(var_df_3D, DOC_vr_aq)

rm(DIC_vr_aq)
rm(CH4_vr_aq)
rm(DOC_vr_aq)


#calculate total soil carbon concentration, express as kg C m-3
SOIL1C_vr <- subset(var_df_3D, variable == "SOIL1C_vr")
SOIL2C_vr <- subset(var_df_3D, variable == "SOIL2C_vr")
SOIL3C_vr <- subset(var_df_3D, variable == "SOIL3C_vr")
SOIL4C_vr <- subset(var_df_3D, variable == "SOIL4C_vr")
LITR1C_vr <- subset(var_df_3D, variable == "LITR1C_vr")
LITR2C_vr <- subset(var_df_3D, variable == "LITR2C_vr")
LITR3C_vr <- subset(var_df_3D, variable == "LITR3C_vr")

soilC_concen <- (SOIL1C_vr$value + SOIL2C_vr$value + SOIL3C_vr$value + SOIL4C_vr$value + LITR1C_vr$value + LITR2C_vr$value + LITR3C_vr$value)/1000
#change variable name and append to the main df
SOIL1C_vr$value <- soilC_concen
SOIL1C_vr$variable <- "soilC_concen"
var_df_3D <- rbind(var_df_3D, SOIL1C_vr)

rm(soilC_concen)
rm(SOIL1C_vr)
rm(SOIL2C_vr)
rm(SOIL3C_vr)
rm(SOIL4C_vr)
rm(LITR1C_vr)
rm(LITR2C_vr)
rm(LITR3C_vr)

rm(porosity)

#a few more units could be changed...leaving for now (look at https://github.com/bsulman/REDOX-PFLOTRAN/blob/master/plot_ELM_alquimia_result.py)

var_df_3D <- var_df_3D %>%
  mutate(soil_depth = case_when(soil_layer == 1 ~ soil_layer_depth[1],
                                soil_layer == 2 ~ soil_layer_depth[2],
                                soil_layer == 3 ~ soil_layer_depth[3],
                                soil_layer == 4 ~ soil_layer_depth[4],
                                soil_layer == 5 ~ soil_layer_depth[5],
                                soil_layer == 6 ~ soil_layer_depth[6],
                                soil_layer == 7 ~ soil_layer_depth[7],
                                soil_layer == 8 ~ soil_layer_depth[8],
                                soil_layer == 9 ~ soil_layer_depth[9],
                                soil_layer == 10 ~ soil_layer_depth[10]))

rm(soil_layer_depth)



#convert time into date format
full_df_2D$date <- as.Date(as.character(full_df_2D$date), format='%Y%m%d') 

#add day, month, year, hour cols
full_df_2D$year <- format(full_df_2D$date, '%Y')
full_df_2D$month <- format(full_df_2D$date, '%m')
full_df_2D$day <- format(full_df_2D$date, '%d')
full_df_2D$hour <- full_df_2D$sec / 3600

#save the 2D df
write.csv(full_df_2D, file.path(path.out, paste0(site_name, "_hourly_polygon_2D_", run_type, ".", output_yr, ".csv")), row.names = FALSE)

#now do the vertically resolved df-------
#convert time into date format
var_df_3D$date <- as.Date(as.character(var_df_3D$date), format='%Y%m%d') #this takes a min

#do these later, super slow
# #add day, month, year, hour cols
# var_df_3D$year <- format(var_df_3D$date, '%Y')
# var_df_3D$month <- format(var_df_3D$date, '%m')
# var_df_3D$day <- format(var_df_3D$date, '%d')
# var_df_3D$hour <- var_df_3D$sec / 3600

#save the 3D df
write.csv(var_df_3D, file.path(path.out, paste0(site_name, "_hourly_polygon_3D_", run_type, ".", output_yr, ".csv")), row.names = FALSE)

# #STOPPED HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# #pulling back in b/c R restarted
# full_df_2D <- read.csv("~/NGEE_modeling/ELM_PFLOTRAN_BEO/1_Runs/Arctic_PFT/BEO_7cell_EXTRACTED/BEO_hourly_polygon_2D_Arctic.1996.csv", header = T)
# var_df_3D <- read.csv("~/NGEE_modeling/ELM_PFLOTRAN_BEO/1_Runs/Arctic_PFT/BEO_7cell_EXTRACTED/BEO_hourly_polygon_3D_Arctic.1996.csv", header = T)

#now calculate daily scale for both 2D and 3D------

#2D df variables that should be grouped by day and summed (all per second fluxes): 
#rain is in mm/sec, so also want to sum it to get mm/day
sum_vars <- c("CH4FLUX_ALQUIMIA", "DIC_RUNOFF", "DOC_RUNOFF", "GPP", "HR", "NEE", "NPP", "QDRAI", "QFLX_EVAP_TOT", "QVEGT", #"QFLX_LAT_AQU",
                 "SMINN_TO_PLANT", "SMIN_NO3_RUNOFF", "FCH4", "RAIN", "SMINN") # just moved "SMINN" here, should be summed to get from gN/m^2/s to gN/m^2/day, forgot!

#2D df variables that should be grouped by day and averaged: 
avg_vars <- c("H2OSFC", "H2OSFC_TIDE", "LEAFC", "TOTLITC", "TOTSOMC", "TOTVEGC", "ZWT_PERCH", "chem_dt", "ALT", "FSAT", "TSA") #"SALINITY", 

#first convert 2D df into long format
full_df_2D_long <- full_df_2D %>%
  pivot_longer(
    cols = CH4FLUX_ALQUIMIA:FSAT, #MANUALLY CHANGE HERE IF VARIABLE LIST CHANGES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    names_to = "variable",
    values_to = "value")

sum_vars_sub <- subset(full_df_2D_long, variable %in% sum_vars)
#fluxes/rates are per second, reported every hour. Multiply by 3600 to get total flux per hour, then can group by day and sum to get 
#total flux per day
sum_vars_sub$value <- sum_vars_sub$value * 3600

sum_vars_sub <- sum_vars_sub %>%
  group_by(polygon_type, variable, year, month, day) %>%
  summarise(value_daily = sum(value, na.rm = TRUE))

avg_vars_sub <- subset(full_df_2D_long, variable %in% avg_vars)

avg_vars_sub <- avg_vars_sub %>%
  group_by(polygon_type, variable, year, month, day) %>%
  summarise(value_daily = mean(value, na.rm = TRUE))

#combine them back into single df
full_df_2D_long_daily <- rbind(sum_vars_sub, avg_vars_sub)

rm(sum_vars_sub)
rm(avg_vars_sub)

#using Sigrid's season approach
snowmelt <- c(4, 5)
growing <- c(6, 7, 8, 9)
freeze <- c(10, 11)
winter <- c(1, 2, 3, 12)

full_df_2D_long_daily <- full_df_2D_long_daily %>%
  mutate(season = case_when(month %in% snowmelt ~ "Snowmelt",
                            month %in% growing ~ "Growing Season",
                            month %in% freeze ~ "Freeze Up",
                            month %in% winter ~ "Winter Dormant"))

#save 2D variable, daily timescale
write.csv(full_df_2D_long_daily, file.path(path.out, paste0(site_name, "_daily_polygon_2D_", run_type, ".", output_yr, ".csv")), row.names = FALSE)


#now bring the 3D df to daily resolution---------

#3D df variables that should be grouped by day and summed (all per second fluxes): 
#"QDRAI_VR", "QFLX_ADV"
sum_vars <- c("QDRAI_VR", "QFLX_ADV")

#3D df variables that should be grouped by day and averaged:
#"CH4_vr", "DIC_vr", "DOC_vr", "H2OSOI", "LITR1C_vr", "LITR2C_vr", "LITR3C_vr", "SIC_vr", "SOIL1C_vr", "SOIL2C_vr", "SOIL3C_vr", "SOIL4C_vr", 
#"SOILICE", "SOILLIQ", "TSOI", "soil_Fe2", "soil_FeOxide", "soil_FeS", "soil_O2", "soil_pH", "soil_salinity", "soil_sulfate", "soil_sulfide", "watsat"
#"soilC_concen", "DIC_vr_aq", "CH4_vr_aq", "DOC_vr_aq", "frozen_frac", "VWC", "soil_Fe2_aq", "soil_FeOxide_aq", "soil_FeS_aq", 
#"soil_O2_aq", "soil_sulfate_aq", "soil_sulfide_aq"

#....basically everything in the 3D df except "QDRAI_VR" and "QFLX_ADV" should just be grouped by day and avg'd

var_df_3D_sum <- subset(var_df_3D, variable %in% sum_vars)
var_df_3D_sum$soil_depth <- as.character(var_df_3D_sum$soil_depth) #for grouping
#fluxes/rates are per second, reported every hour. Multiply by 3600 to get total flux per hour, then can group by day and sum to get 
#total flux per day
var_df_3D_sum$value <- var_df_3D_sum$value * 3600

var_df_3D_sum <- var_df_3D_sum %>%
  group_by(variable, polygon_type, soil_depth, date) %>%
  summarise(value_daily = sum(value, na.rm = TRUE))

var_df_3D_avg <- subset(var_df_3D, variable != "QDRAI_VR" & variable != "QFLX_ADV")
var_df_3D_avg$soil_depth <- as.character(var_df_3D_avg$soil_depth) #for grouping

var_df_3D_avg <- var_df_3D_avg %>%
  group_by(variable, polygon_type, soil_depth, date) %>%
  summarise(value_daily = mean(value, na.rm = TRUE))


#combine them back into single df
full_df_3D_long_daily <- rbind(var_df_3D_sum, var_df_3D_avg)

rm(var_df_3D_sum)
rm(var_df_3D_avg)
rm(var_df_3D)

#convert soil depth back to a number
full_df_3D_long_daily$soil_depth <- as.numeric(full_df_3D_long_daily$soil_depth)

full_df_3D_long_daily$year <- format(full_df_3D_long_daily$date, '%Y')
full_df_3D_long_daily$month <- format(full_df_3D_long_daily$date, '%m')
full_df_3D_long_daily$day <- format(full_df_3D_long_daily$date, '%d')


full_df_3D_long_daily <- full_df_3D_long_daily %>%
  mutate(season = case_when(month %in% snowmelt ~ "Snowmelt",
                            month %in% growing ~ "Growing Season",
                            month %in% freeze ~ "Freeze Up",
                            month %in% winter ~ "Winter Dormant"))

#save 3D variable, daily timescale
write.csv(full_df_3D_long_daily, file.path(path.out, paste0(site_name, "_daily_polygon_3D_", run_type, ".", output_yr, ".csv")), row.names = FALSE)


#clean up, then go back to top and extract second year of output
rm(full_df_2D)
rm(full_df_2D_long)
rm(full_df_2D_long_daily)
rm(full_df_3D_long_daily)

#THIS IS WHERE LOOPING THROUGH INDIVIDUAL ANNUAL FILES ENDS!!!!!!!!!!!!!!!!

#--------------------------------------------------------------------------------------------------------------


#now pull in all years of model output, combine 
full_df_2D_long_daily <- readbulk::read_bulk(directory = path.out, extension = paste0("BEO_daily_polygon_2D_", run_type, ".*.csv"), header = TRUE)
write.csv(full_df_2D_long_daily, file.path(path.out, paste0(site_name, "_daily_polygon_2D_", run_type, "_", start_yr, "_", end_yr,".csv")), row.names=F)

full_df_3D_long_daily <- readbulk::read_bulk(directory = path.out, extension = paste0("BEO_daily_polygon_3D_", run_type, ".*.csv"), header = TRUE)
write.csv(full_df_3D_long_daily, file.path(path.out, paste0(site_name, "_daily_polygon_3D_", run_type, "_", start_yr, "_", end_yr,".csv")), row.names=F)


#trying w/single yr, seasonal pattern for CH4 weird...
#full_df_3D_long_daily <- subset(full_df_3D_long_daily, year == 1999)
#full_df_3D_long_daily$File <- NULL

#don't know why this is added to avg again, should already be daily avg, not running again
# full_df_2D_long_daily <- full_df_2D_long_daily %>%
#   group_by(polygon_type, variable, month, day, season) %>%
#   summarise(value_daily = mean(value_daily, na.rm = T))
# 
# full_df_3D_long_daily <- full_df_3D_long_daily %>%
#   group_by(polygon_type, variable, soil_depth, month, day, season) %>%
#   summarise(value_daily = mean(value_daily, na.rm = T))

#but DO want to drop the file col, don't need
full_df_2D_long_daily$File <- NULL
full_df_3D_long_daily$File <- NULL
#hmm looks like season got messed up on the 3D file...
#try running again
full_df_3D_long_daily <- full_df_3D_long_daily %>%
  mutate(season = case_when(month %in% snowmelt ~ "Snowmelt",
                            month %in% growing ~ "Growing Season",
                            month %in% freeze ~ "Freeze Up",
                            month %in% winter ~ "Winter Dormant"))


#just realized CH4FLUX_ALQUIMIA is natively in gC/m^2/s and I converted it to gC/m^2/day, but FCH4 (default methane model) is in kgC/m2/s
#I already converted it to daily, but need to go from kgC to gC
full_df_2D_long_daily$value_daily <- ifelse(full_df_2D_long_daily$variable == "FCH4", full_df_2D_long_daily$value_daily*1000, full_df_2D_long_daily$value_daily)

#also forgot to convert TSA (2m air temp) from K to C...
full_df_2D_long_daily$value_daily <- ifelse(full_df_2D_long_daily$variable == "TSA", full_df_2D_long_daily$value_daily-273.15, full_df_2D_long_daily$value_daily)

#before plotting want to remove outliers by screening out anything thats > +- 5 sd
#I think should group by variable, polygon type, year, and season for outlier screening
#test <- full_df_2D_long_daily[1:200,]
full_df_2D_long_daily <- full_df_2D_long_daily %>%
  group_by(variable, polygon_type, year, season) %>%
  mutate(mean_value = mean(value_daily, na.rm = TRUE),
         sd_value = sd(value_daily, na.rm = TRUE),
         is_outlier = abs(value_daily - mean_value) > 5 * sd_value,
         value_daily_OR = ifelse(is_outlier, NA, value_daily)) %>%
  select(-mean_value, -sd_value, -is_outlier) # Remove temporary columns used in calculation

#hmm don't appear to be any outliers in the daily 2D data...

full_df_3D_long_daily <- full_df_3D_long_daily %>%
  group_by(variable, polygon_type, year, season) %>%
  mutate(mean_value = mean(value_daily, na.rm = TRUE),
         sd_value = sd(value_daily, na.rm = TRUE),
         is_outlier = abs(value_daily - mean_value) > 5 * sd_value,
         value_daily_OR = ifelse(is_outlier, NA, value_daily)) %>%
  select(-mean_value, -sd_value, -is_outlier) # Remove temporary columns used in calculation

#nor in the 3D data...

#alquimia ch4 looks like it still has outliers..
test <- subset(full_df_2D_long_daily, variable == "CH4FLUX_ALQUIMIA")
#ok think first need to remove implausible values THEN do the outlier screening otherwise the really crazy values throw things off
test <- test %>%
  group_by(variable, polygon_type, year, season) %>% #trying w/o season
  mutate(mean_value = mean(value_daily, na.rm = TRUE),
         sd_value = sd(value_daily, na.rm = TRUE),
         is_outlier = abs(value_daily - mean_value) > (5 * sd_value),
         value_daily_OR = ifelse(is_outlier, NA, value_daily)) #%>%
  #select(-mean_value, -sd_value, -is_outlier) 


#h20sfc_tide should NOT be variable in time!!! I think I might've messed up assigning the polygon types...
#which means a whole bunch of stuff is messed up because the wrong polygon types are linked with data
test <- subset(full_df_2D_long_daily, variable == "H2OSFC_TIDE")

#----------------
#PLOTS 
#----------------

#violin plots aren't super useful, dropping for now----------------------------------------

# #violin plots of 2D variables, daily data across whole year (1999 only right now)
# c_vars_1 <- c("CH4FLUX_ALQUIMIA", "GPP", "HR", "NEE", "NPP")
# c_vars_2 <- c("LEAFC", "TOTLITC", "TOTSOMC", "TOTVEGC")
# water_vars_1 <- c("DIC_RUNOFF", "DOC_RUNOFF", "QDRAI", "QFLX_EVAP_TOT", "QFLX_LAT_AQU")
# water_vars_2 <- c("QVEGT", "H2OSFC", "H2OSFC_TIDE", "SALINITY", "ZWT")
# other_vars <- c("SMINN_TO_PLANT", "SMIN_NO3_RUNOFF", "SMINN", "chem_dt")
# 
# var_list <- c_vars_1 #pick which var list to plot, CHANGE HERE!!!
# 
# full_df_2D_long_daily_sub <- subset(full_df_2D_long_daily, variable %in% var_list)
# 
# #idk that these violin plots are really that helpful...
# ggplot(full_df_2D_long_daily_sub, aes(x = factor(polygon_type), y = value_daily, fill = polygon_type, color = polygon_type)) +
#   geom_violin(trim = 0.05) + # Trim 5% from each end of distribution (so only plotting up to the 95th percentile) to reduce visual impact of outliers
#   #stat_summary(fun = mean, geom = "crossbar", color = "black", linetype = "dashed") + #width = 0.5, linewidth = 1,
#   stat_summary(fun = median, geom = "crossbar", width = 0.5, linewidth = 0.5, color = "black", show.legend = FALSE) +
#   facet_wrap(~variable, scales = "free") +
#   scale_fill_manual(values = c("FCPcenter" = "#abba8f", "LCPcenter" = "#2a7c6e", "LCPtrough" = "#2aaeae", "HCPcenter" = "#c85c35", 
#                                "HCPtrough" = "#c69a51", "LCPrim" = "#b17843", "FCPtrough" = "#cdd0bb")) +
#   scale_color_manual(values = c("FCPcenter" = "#abba8f", "LCPcenter" = "#2a7c6e", "LCPtrough" = "#2aaeae", "HCPcenter" = "#c85c35", 
#                                 "HCPtrough" = "#c69a51", "LCPrim" = "#b17843", "FCPtrough" = "#cdd0bb")) +
#   theme_minimal() +
#   labs(x = NULL, y = "ADD UNITS", fill = "Polygon Type") +
#   guides(color = "none") +
#   theme(axis.text.x = element_blank(), strip.text.x = element_text(size = 12)) 
#   #ggtitle("Violin Plots of Variable1 by Soil Type, Faceted by Variable2")
# 




#ridgeplots of 2D variables by season---------------------------------------

#season_list <- c("Spring", "Summer", "Winter", "Fall") #old ones
season_list <- c("Snowmelt", "Growing Season", "Freeze Up", "Winter Dormant")
length(season_list)
var_list <- unique(full_df_2D_long_daily$variable)
#var_list <- var_list[1:2]
#units_list <- c("gC m-2 day-1", )
#units_plot <- grid::textGrob(var_plot, gp=grid::gpar(fontsize=14))
#title_plot <- grid::textGrob(var_plot, gp=grid::gpar(fontsize=14, fontface='bold'))
#full_df_2D_long_daily_sub <- subset(full_df_2D_long_daily, variable == var_plot)
#range_plot <- range(full_df_2D_long_daily_sub$value_daily) #looks shitty when they all have same x scale...can't see variability in summer
for(var_plot in unique(var_list)) {
  plot_list <- list()
  title_plot <- grid::textGrob(var_plot, gp=grid::gpar(fontsize=14, fontface='bold'))
  
for(season_plot in unique(season_list)) {
  
full_df_2D_long_daily_sub <- subset(full_df_2D_long_daily, variable == var_plot & season == season_plot)

# Calculate median for each group
medians <- full_df_2D_long_daily_sub %>%
  group_by(variable, polygon_type, season) %>%
  summarise(median = median(value_daily, na.rm = TRUE))

p <- ggplot(data = full_df_2D_long_daily_sub, aes(x = value_daily, y = season, fill = polygon_type)) +
  geom_density_ridges(alpha = .7, rel_min_height = .01, color = "white") +
  #facet_wrap(~ season, ncol = 2, nrow = 2) + #scales = "free",
  scale_fill_manual(values = c("FCPcenter" = "#abba8f", "LCPcenter" = "#2a7c6e", "LCPtrough" = "#2aaeae", "HCPcenter" = "#781e1f", 
                               "HCPtrough" = "#c85c35", "LCPrim" = "#b4f0f8", "FCPtrough" = "#cdd0bb")) +
  
  #using scale_fill_cyclical was messing up the group to color mapping!!!!!!!!!!!!
  #scale_fill_cyclical(labels = c(Trough = "Trough", LCPcenter = "LCPcenter", LCPtransition = "LCPtransition",
                                 #HCPcenter = "HCPcenter", HCPtransition = "HCPtransition", Rim = "Rim", Mean = "Mean"),
                      #values = c("#abba8f", "#2a7c6e", "#2aaeae", "#c85c35", "#c69a51", "#b17843", "#cdd0bb"), name = NULL, guide = "legend") +
  
  
  geom_vline(data = medians, aes(xintercept = median, color = polygon_type), linetype = "dashed", linewidth = 1.5) +
  scale_color_manual(values = c("FCPcenter" = "#abba8f", "LCPcenter" = "#2a7c6e", "LCPtrough" = "#2aaeae", "HCPcenter" = "#781e1f", 
                                "HCPtrough" = "#c85c35", "LCPrim" = "#b4f0f8", "FCPtrough" = "#cdd0bb"), name = NULL, guide = "none") +
  theme_ridges(font_size = 13, grid = T) +
  #scale_x_continuous(limits = c(range_plot[1], 0.02)) +  
  #stat_summary(fun = median, geom = "crossbar", width = 0.5, linewidth = 0.5, color = "black", show.legend = FALSE) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "top", strip.background = element_blank(), strip.text = element_blank(), legend.text = element_text(size = 15),
        axis.text.y = element_text(size = 15), axis.title = element_text(size = 15), axis.text.x = element_text(size = 13),
        legend.title = element_blank()) 

plot_list[[season_plot]] <- p

}

grid_arrange_shared_legend(plot_list[[1]], plot_list[[2]], plot_list[[3]],
                           plot_list[[4]], #plot_list[[5]], plot_list[[6]], #plot_list[[7]], #plot_list[[8]], 
                           #plot_list[[9]],
                           ncol = 2, nrow = 2, position='top', top = title_plot) #bottom = units_plot)

}



#timeseries of a few of the 2D variables----------------------------
#"CH4FLUX_ALQUIMIA", "GPP", "HR", "NEE", "ZWT"

#set up function for facet tags
tag_facet <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf, 
                      hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
  
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
  p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE) 
}


#want to focus in on just net carbon balance related variables for now, so "FCH4", "GPP", "ER", where for now NECB = GPP - ER - FCH4 (eventually will include leaching)
variable_list <- c("CH4FLUX_ALQUIMIA", "FCH4", "ZWT_PERCH")
variable_list <- c("GPP", "HR", "NEE")
variable_list <- c("ALT", "H2OSFC_TIDE", "FSAT")

subset_data <- subset(full_df_2D_long_daily, variable %in% variable_list)
#subset_data$month <- as.numeric(subset_data$month)
#subset_data$day <- as.numeric(subset_data$day)
# units <- expression(kgC ~m^{-2}~ month ^{-1})
# #units <- expression(atop(GPP, paste(~kgC ~m^{-2}~ month ^{-1}))) #atop() basically fakes a line break

p <- ggplot(data = subset_data) + 
  facet_grid(variable ~ ., scales = "free_y") + #, scales = "free_y"allows y axis to vary, but same xaxis for all
  geom_smooth(aes(x = month, y= value_daily, color = polygon_type), linewidth = 0.8) + #linetype = "dotdash") + #add linear fit line
  # annotate("text", label=(paste0("slope==", coef(lm(function_sub_single$MEAN~function_sub_single$year))[2])),
  #          parse=TRUE) +
  # geom_ribbon(aes(x = month, y= value_daily, ymin = value - SD, ymax = value + SD, fill = run_type),
  #             color = "transparent", alpha = .3) +
  # scale_fill_manual(values = c("Arctic" = "#8c86a0", "Default" = "#7c544c", "obs" = "#cf932c")) +
  #scale_color_manual(values = c("Arctic" = "#8c86a0", "Default" = "#7c544c", "obs" = "#cf932c"), labels = c("Model (Arctic PFT)", "Model (Default PFT)", "Flux Tower")) + 
  #labs(x = "Year", y = units, title = title_label) + 
  scale_color_manual(values = c("FCPcenter" = "#abba8f", "LCPcenter" = "#2a7c6e", "LCPtrough" = "#2aaeae", "HCPcenter" = "#781e1f", 
                                "HCPtrough" = "#c85c35", "LCPrim" = "#b4f0f8", "FCPtrough" = "#cdd0bb")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_continuous(breaks=seq(1,12,1)) +
  labs(x = "Month", y = NULL) + 
  #guides(fill = "none") + #don't display fill legend
  theme_bw() +
  theme(legend.title = element_blank(), axis.text = element_text(size = 10), axis.title.x = element_text(size = 12), 
        legend.position = "top", legend.text = element_text(size = 12), strip.text.x = element_text(size = 12), 
        strip.text.y = element_text(size = 12), strip.background = element_blank()) 

tag_facet(p) 


#time series but broken up by polygon type and not smoothed-------------
#need a dummy year to set up date object, but avg'd across 1999-2000, just say 2000
subset_data$year <- 1998

# Function to convert integer to character with leading zero if single digit
add_leading_zero <- function(num) {
  formatted_num <- sprintf("%02d", num)
  return(formatted_num)
}

# Convert month/day to have leading zeros
subset_data$month <- sapply(subset_data$month, add_leading_zero)
subset_data$day <- sapply(subset_data$day, add_leading_zero)

subset_data$date <- paste0(subset_data$year, subset_data$month, subset_data$day)
subset_data$date <- as.Date(subset_data$date, format='%Y%m%d')

#subset_data$date <- as.Date(as.character(subset_data$date), format='%Y%m%d')
#summer_months <- c("06", "07", "08", "09")
#subset_data <- subset(subset_data, month %in% summer_months)

#want to see whats going on in 1998
test <- subset(subset_data, year == 1998)

p <- ggplot(data = test) + 
  facet_grid(variable ~ ., scales = "free_y") + #, scales = "free_y"allows y axis to vary, but same xaxis for all
  geom_line(aes(x = date, y= value_daily_OR, color = polygon_type), linewidth = 0.8) + #linetype = "dotdash") + #add linear fit line
  # annotate("text", label=(paste0("slope==", coef(lm(function_sub_single$MEAN~function_sub_single$year))[2])),
  #          parse=TRUE) +
  # geom_ribbon(aes(x = month, y= value_daily, ymin = value - SD, ymax = value + SD, fill = run_type),
  #             color = "transparent", alpha = .3) +
  # scale_fill_manual(values = c("Arctic" = "#8c86a0", "Default" = "#7c544c", "obs" = "#cf932c")) +
  #scale_color_manual(values = c("Arctic" = "#8c86a0", "Default" = "#7c544c", "obs" = "#cf932c"), labels = c("Model (Arctic PFT)", "Model (Default PFT)", "Flux Tower")) + 
  #labs(x = "Year", y = units, title = title_label) + 
  scale_color_manual(values = c("FCPcenter" = "#abba8f", "LCPcenter" = "#2a7c6e", "LCPtrough" = "#2aaeae", "HCPcenter" = "#781e1f", 
                                "HCPtrough" = "#c85c35", "LCPrim" = "#b4f0f8", "FCPtrough" = "#cdd0bb")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #scale_x_continuous(breaks=seq(1,12,1)) +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "1 month") + # "%b" Display month names, "%m" is month numbers, expand = c(0, 0) removes extra space on the x-axis
  labs(x = "Year-Month", y = NULL) + 
  #guides(fill = "none") + #don't display fill legend
  theme_bw() +
  theme(legend.title = element_blank(), axis.text = element_text(size = 10), axis.title.x = element_text(size = 12), 
        legend.position = "top", legend.text = element_text(size = 12), strip.text.x = element_text(size = 12), 
        strip.text.y = element_text(size = 12), strip.background = element_blank()) 

#tag_facet(p) #not working for some reason, just print out plot w/o facet tags


#3D variable plots----------------------------------------------------------------------------------------------

# #quick plot of hourly data to see how things look
# summer_months <- c("06", "07", "08")
# test <- subset(full_df_3D_long_daily, variable == "CH4_vr" & season == "Summer") #& month == "08"
# 
# ggplot(test, aes(x = value_daily, y = soil_depth, color = polygon_type)) +
#   geom_point() +
#   geom_line() +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   labs(x = "Dissolved Methane (g C m-3)", y = "Soil Depth (m)", color = "Polygon Type") +
#   scale_y_reverse() + #reverse y axis order so soil surface is at top
#   theme_minimal()



#tile plots of variables by depth and time, faceted by month and polygon type------------------------------------------

subset_data <- subset(full_df_3D_long_daily, variable == "TSOI") #& polygon_type == "HCPcenter")

#fill_name <- expression(atop("Soil Temp.", paste(~gC ~m^{-3})))
fill_name <- "Temp. (C)"

subset_data$day <- as.numeric(subset_data$day)

subset_data$soil_depth <- factor(subset_data$soil_depth, levels = c("2.86", "1.73", "1.04", "0.62", "0.37", "0.21", "0.12", "0.06", 
                                                                    "0.03", "0.01"))

#need to convert soil depth to factor so it extends the width of the line and fills out the plot
g1 <- ggplot(data = subset_data, aes(x = day, y = soil_depth, z = value_daily))  +
  #scale_y_continuous(breaks = seq(0, 3, 0.5)) + #breaks = seq(0, 3, 0.5)
  #scale_y_reverse() +
  facet_grid(polygon_type~month) +
  #scale_y_continuous(trans = "reverse", breaks = seq(0, 3, 1)) +
  scale_y_discrete(breaks = levels(subset_data$soil_depth)[seq(1, length(levels(subset_data$soil_depth)), by = 2)]) +
  scale_x_continuous(breaks =c(1,10,20,31)) +
  labs(x = "Day", y = "Soil Depth (m)",
       fill = fill_name, title = "Soil Temperature")

p1 <- g1 + geom_tile(aes(fill = value_daily)) + #color= "white",size=0.1, alpha = .7, 
  #scale_fill_viridis_c(option = "mako") +
  scale_fill_viridis(option ="C") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12), #axis.title.x = element_blank(), axis.title.y = element_blank(), 
        legend.text = element_text(size = 13), legend.title = element_text(size = 13),
        title = element_text(size = 13))

p1


#simpler version of the above plot, not broken up by month
#tile plots of variables by depth and time, faceted by polygon type------------------------------------------

#need a dummy year to set up date object, but avg'd across 1999-2000, just say 2000
full_df_3D_long_daily$year <- 2000
full_df_3D_long_daily$month <- sapply(full_df_3D_long_daily$month, add_leading_zero)
full_df_3D_long_daily$day <- sapply(full_df_3D_long_daily$day, add_leading_zero)
full_df_3D_long_daily$date <- paste0(full_df_3D_long_daily$year, full_df_3D_long_daily$month, full_df_3D_long_daily$day)
full_df_3D_long_daily$date <- as.Date(full_df_3D_long_daily$date, format='%Y%m%d')

variable_list <- c("CH4_vr_aq", "DIC_vr_aq", "DOC_vr_aq", "VWC", "TSOI", "frozen_frac", "soilC_concen",
                   "soil_Fe2_aq", "soil_FeOxide_aq", "soil_FeS_aq", "soil_O2_aq", 
                   "soil_pH", "soil_sulfate_aq", "soil_sulfide_aq", "soil_acetate")

variable_plot <- variable_list[1]
print(variable_plot)

#subset_data <- subset(full_df_3D_long_daily, variable == "TSOI") #& polygon_type == "HCPcenter")
subset_data <- subset(full_df_3D_long_daily, variable == variable_plot)

#fill_name <- expression(atop("Soil Temp.", paste(~gC ~m^{-3})))
fill_name <- variable_plot
variable_names <- c("Soil Dissolved CH4 (mmol L-1)", "Soil DIC (mmol L-1)", "Soil DOC (mmol L-1)", "VWC",
                    "T Soil (C)", "Frozen Fraction", "Soil C Concentration (gC m-3)", "Soil Fe2 (mmol L-1)",
                    "Soil Fe Oxide (mmol L-1)", "Soil FeS (mmol L-1)", "Soil O2 (mmol L-1)",
                    "Soil pH", "Soil Sulfate (mmol L-1)", "Soil Sulfide (mmol L-1)", "Soil Acetate (mmol L-1)", 
                    "Soil Inorganic C (gC m-3)")
plot_title <- variable_names[1]

#something weird going on w/CH4 in deepest level...cutting off to see if 0-1m look normal or if I made units error
#subset_data <- subset(subset_data, soil_depth <= 1.04) #1.04

#full depth profile
subset_data$soil_depth <- factor(subset_data$soil_depth, levels = c("2.86", "1.73", "1.04", "0.62", "0.37", "0.21", "0.12", "0.06", 
                                                                    "0.03", "0.01"))

#this one is w/o deepest 2 soil layers
subset_data$soil_depth <- factor(subset_data$soil_depth, levels = c("0.62", "0.37", "0.21", "0.12", "0.06", 
                                                                    "0.03", "0.01"))

#top meter-ish
subset_data$soil_depth <- factor(subset_data$soil_depth, levels = c("1.04", "0.62", "0.37", "0.21", "0.12", "0.06", 
"0.03", "0.01"))

#need to convert soil depth to factor so it extends the width of the line and fills out the plot
g1 <- ggplot(data = subset_data, aes(x = date, y = soil_depth, z = value_daily))  +
  #scale_y_continuous(breaks = seq(0, 3, 0.5)) + #breaks = seq(0, 3, 0.5)
  #scale_y_reverse() +
  facet_grid(polygon_type~.) +
  #scale_y_continuous(trans = "reverse", breaks = seq(0, 3, 1)) +
  scale_y_discrete(breaks = levels(subset_data$soil_depth)[seq(1, length(levels(subset_data$soil_depth)), by = 2)]) +
  scale_x_date(date_labels = "%b") +  # Format to display only month
  labs(x = "Time", y = "Soil Depth (m)",
       fill = fill_name, title = plot_title)

p1 <- g1 + geom_tile(aes(fill = value_daily)) + #color= "white",size=0.1, alpha = .7, 
  #scale_fill_viridis_c(option = "mako") +
  scale_fill_viridis(option ="C") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12), #axis.title.x = element_blank(), axis.title.y = element_blank(), 
        legend.text = element_text(size = 13), legend.title = element_text(size = 13),
        title = element_text(size = 13), strip.text.y = element_text(size = 10), legend.position = "right")

#tag_facet(p1) #not working
p1





#plot of profiles grouped by season and avg'd---------------------------------------------------------------

variable_list <- c("CH4_vr_aq", "DIC_vr_aq", "DOC_vr_aq", "VWC", "TSOI", "frozen_frac", "soilC_concen",
                   "soil_Fe2_aq", "soil_FeOxide_aq", "soil_FeS_aq", "soil_O2_aq", 
                   "soil_pH", "soil_sulfate_aq", "soil_sulfide_aq", "soil_acetate_aq", "SIC_vr")

variable_names <- c("Soil Dissolved CH4 (mmol L-1)", "Soil DIC (mmol L-1)", "Soil DOC (mmol L-1)", "VWC",
                    "T Soil (C)", "Frozen Fraction", "Soil C Concentration (gC m-3)", "Soil Fe2 (mmol L-1)",
                    "Soil Fe Oxide (mmol L-1)", "Soil FeS (mmol L-1)", "Soil O2 (mmol L-1)",
                    "Soil pH", "Soil Sulfate (mmol L-1)", "Soil Sulfide (mmol L-1)", "Soil Acetate (mmol L-1)", 
                    "Soil Inorganic C (gC m-3)")


variable_plot <- variable_list[16]
print(variable_plot)

plot_title <- variable_names[16] #"Soil C Concentration (gC m-3)" #"Soil Inorganic C (gC m-3)"

subset_data <- subset(full_df_3D_long_daily, variable == variable_plot)

subset_data <- subset_data %>%
  group_by(variable, polygon_type, soil_depth, season) %>%
  summarise(value_daily = mean(value_daily, na.rm = T))


subset_data <- subset(subset_data, soil_depth <= 1.04) #1.04 #0.5


subset_data$season <- factor(subset_data$season, levels = c("Winter", "Spring", "Summer", "Fall"))



ggplot(subset_data, aes(x = soil_depth, y = value_daily, color = polygon_type)) +
  facet_grid(.~season) +
  #geom_point() +
  geom_line(aes(group = polygon_type), linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("FCPcenter" = "#abba8f", "LCPcenter" = "#2a7c6e", "LCPtrough" = "#2aaeae", "HCPcenter" = "#781e1f", 
                                "HCPtrough" = "#c85c35", "LCPrim" = "#b4f0f8", "FCPtrough" = "#cdd0bb")) +
  #scale_y_discrete(breaks = levels(subset_data$soil_depth)) +
  scale_x_reverse() + #reverse y axis order so soil surface is at top
  coord_flip() + #have to do it this way otherwise it doesn't get plotted shallowest to deepest!!!
  labs(y = plot_title, x = "Soil Depth (m)", color = "Polygon Type") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"), strip.text.x = element_text(size = 14),
        legend.title = element_blank(), axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14), axis.title = element_text(size = 14), legend.position = "top")







#correlation plots for chemistry variables, separate plot for each polygon type-------------------------------------------

library(reshape2)

polygon_plot <- "HCPcenter"

#3D variables------------------------------------

#pulled out for now: "SOIL2C_vr","SOIL3C_vr","SOIL4C_vr"
variables_compare <- c("QDRAI_VR","QFLX_ADV","CH4_vr","DIC_vr","DOC_vr",         
                       "LITR1C_vr","LITR2C_vr","LITR3C_vr","SOIL1C_vr",         
                       "TSOI","VWC","frozen_frac","soilC_concen","soil_Fe2","soil_FeOxide",    
                       "soil_FeS","soil_O2","soil_pH","soil_salinity","soil_sulfate","soil_sulfide", 
                       "soil_acetate")


polygon_compare <- subset(full_df_3D_long_daily, polygon_type == polygon_plot & variable %in% variables_compare)

polygon_compare <- polygon_compare %>%
  pivot_wider(names_from = variable, values_from = value_daily)

cor_polygon_compare <- round(cor(polygon_compare[,8:29], method = "spearman"),2)

# Melt the correlation matrix
cor_melted <- melt(cor_polygon_compare)
# Convert factors to character vectors
cor_melted$Var1 <- as.character(cor_melted$Var1)
cor_melted$Var2 <- as.character(cor_melted$Var2)

cor_melted <- cor_melted[cor_melted$Var1 != cor_melted$Var2 & cor_melted$Var1 < cor_melted$Var2, ] # Filter out lower triangle and diagonal


# Plotting
p1 <- ggplot(cor_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#ef8a62", high = "#67a9cf", mid = "#F0EAD6", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab",
                       name="Spearman \nCorrelation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 13), axis.text.y = element_text(size = 13),
        legend.title = element_text(size = 13)) +
  labs(title = polygon_plot, x = "", y = "") +
  coord_fixed()

p1


#2D variables------------------------------------

variables_compare <- c("CH4FLUX_ALQUIMIA", "DIC_RUNOFF","DOC_RUNOFF","GPP","HR","NEE","NPP","QDRAI","QFLX_EVAP_TOT",   
"QFLX_LAT_AQU","QVEGT","SMINN_TO_PLANT","SMIN_NO3_RUNOFF","H2OSFC","LEAFC","SMINN",           
"TOTLITC","TOTSOMC","TOTVEGC","ZWT")  

polygon_plot <- "HCPcenter"

polygon_compare <- subset(full_df_2D_long_daily, polygon_type == polygon_plot & variable %in% variables_compare)

polygon_compare <- polygon_compare %>%
  pivot_wider(names_from = variable, values_from = value_daily)

cor_polygon_compare <- round(cor(polygon_compare[,5:24], method = "spearman"),2)

# Melt the correlation matrix
cor_melted <- melt(cor_polygon_compare)
# Convert factors to character vectors
cor_melted$Var1 <- as.character(cor_melted$Var1)
cor_melted$Var2 <- as.character(cor_melted$Var2)

cor_melted <- cor_melted[cor_melted$Var1 != cor_melted$Var2 & cor_melted$Var1 < cor_melted$Var2, ] # Filter out lower triangle and diagonal


# Plotting
p1 <- ggplot(cor_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#ef8a62", high = "#67a9cf", mid = "#F0EAD6", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab",
                       name="Spearman \nCorrelation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 13), axis.text.y = element_text(size = 13),
        legend.title = element_text(size = 13)) +
  labs(title = polygon_plot, x = "", y = "") +
  coord_fixed()

p1



#barplot of CH4 flux relative to cO2 flux for the different polygon types--------------------------------------------------------------

variables_sub <- c("NEE", "CH4FLUX_ALQUIMIA", "HR", "GPP") #both are in gC/m^2/s

full_df_2D_long_daily_sub <- subset(full_df_2D_long_daily, variable %in% variables_sub)

#try dropping the last year and replotting, since cH4 and HR are weird for 2000
full_df_2D_long_daily_sub <- subset(full_df_2D_long_daily_sub, year != 2000)

#want to average all years together first for these figs
full_df_2D_long_daily_sub <- full_df_2D_long_daily_sub %>%
  group_by(polygon_type, variable, season, month, day) %>%
  summarize(value_daily_sd = sd(value_daily, na.rm = TRUE), value_daily = mean(value_daily, na.rm = TRUE))
  

#now group by polygon type and find the flux associated with each polygon type
match_summary <- full_df_2D_long_daily_sub %>%
  group_by(polygon_type, variable, season) %>%
  summarize(value_daily_avg = mean(value_daily, na.rm = TRUE), value_daily_sd = sd(value_daily, na.rm = TRUE))

match_summary$season <- factor(match_summary$season, levels = c("Winter", "Spring", "Summer", "Fall"))

#polygon_colors <- c("#967acc")
polygon_colors <- c("#abba8f", "#cdd0bb", "#781e1f", "#c85c35", "#2a7c6e", "#b4f0f8", "#2aaeae")

match_summary_CH4 <- subset(match_summary, variable == "CH4FLUX_ALQUIMIA")

ggplot(match_summary_CH4, aes(x= reorder(polygon_type, value_daily_avg), y = value_daily_avg, fill = factor(polygon_type))) +
  #geom_bar(binwidth = 1, fill = PFT, alpha = 0.7) +
  facet_grid(.~season) +
  geom_bar(stat="identity")+
  #geom_bar(aes(fill = factor(PFT_name))+
  #scale_fill_manual(values = pft_colors, labels = names) +
  #scale_fill_manual(values = pft_colors, labels = c("CH4 (gC m-2 s-1)")) +
  scale_fill_manual(values = polygon_colors) +
  geom_errorbar(aes(ymin = value_daily_avg-value_daily_sd, ymax = value_daily_avg+value_daily_sd), width=.2,
                position=position_dodge(.9)) +
  coord_flip() +
  labs(title = "Polygon Type Avg. CH4 Flux", x = NULL,
       y = "CH4 (gC m-2 s-1)") +
  geom_hline(yintercept = 0, linetype = 2, color = "black", linewidth = 1) +
  #theme_minimal() +
  theme_bw() +
  theme(legend.title = element_blank(), axis.text = element_text(size = 12),  
        legend.text = element_text(size = 12), legend.position="bottom", strip.background = element_rect(fill = "white")) 




match_summary_CO2 <- subset(match_summary, variable == "NEE")

ggplot(match_summary_CO2, aes(x= reorder(polygon_type, -value_daily_avg), y = value_daily_avg, fill = factor(polygon_type))) +
  geom_bar(stat="identity")+
  facet_grid(.~season) +
  scale_fill_manual(values = polygon_colors) +
  geom_errorbar(aes(ymin = value_daily_avg-value_daily_sd, ymax = value_daily_avg+value_daily_sd), width=.2,
                position=position_dodge(.9)) +
  coord_flip() +
  labs(title = "Polygon Type Avg. NEE", x = NULL,
       y = "NEE (gC m-2 s-1)") +
  geom_hline(yintercept = 0, linetype = 2, color = "black", linewidth = 1) +
  #theme_minimal() +
  theme_bw() +
  theme(legend.title = element_blank(), axis.text = element_text(size = 12),  
        legend.text = element_text(size = 12), legend.position="bottom", strip.background = element_rect(fill = "white")) 


match_summary_HR <- subset(match_summary, variable == "HR")

ggplot(match_summary_HR, aes(x= reorder(polygon_type, value_daily_avg), y = value_daily_avg, fill = factor(polygon_type))) +
  geom_bar(stat="identity")+
  facet_grid(.~season) +
  scale_fill_manual(values = polygon_colors) +
  geom_errorbar(aes(ymin = value_daily_avg-value_daily_sd, ymax = value_daily_avg+value_daily_sd), width=.2,
                position=position_dodge(.9)) +
  coord_flip() +
  labs(title = "Polygon Type Avg. Respiration", x = NULL,
       y = "HR (gC m-2 s-1)") +
  geom_hline(yintercept = 0, linetype = 2, color = "black", linewidth = 1) +
  #theme_minimal() +
  theme_bw() +
  theme(legend.title = element_blank(), axis.text = element_text(size = 12),  
        legend.text = element_text(size = 12), legend.position="bottom", strip.background = element_rect(fill = "white")) 



match_summary_GPP <- subset(match_summary, variable == "GPP")

ggplot(match_summary_GPP, aes(x= reorder(polygon_type, value_daily_avg), y = value_daily_avg, fill = factor(polygon_type))) +
  geom_bar(stat="identity")+
  facet_grid(.~season) +
  scale_fill_manual(values = polygon_colors) +
  geom_errorbar(aes(ymin = value_daily_avg-value_daily_sd, ymax = value_daily_avg+value_daily_sd), width=.2,
                position=position_dodge(.9)) +
  coord_flip() +
  labs(title = "Polygon Type Avg. GPP", x = NULL,
       y = "GPP (gC m-2 s-1)") +
  geom_hline(yintercept = 0, linetype = 2, color = "black", linewidth = 1) +
  #theme_minimal() +
  theme_bw() +
  theme(legend.title = element_blank(), axis.text = element_text(size = 12),  
        legend.text = element_text(size = 12), legend.position="bottom", strip.background = element_rect(fill = "white")) 








# p2 <- ggplot(cor_melted, aes(Var1, Var2, fill = value)) +
#   geom_tile() +
#   scale_fill_gradient2(low = "#ef8a62", high = "#67a9cf", mid = "#F0EAD6", 
#                        midpoint = 0, limit = c(-1, 1), space = "Lab",
#                        name="Spearman \nCorrelation") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 13), axis.text.y = element_text(size = 13),
#         legend.title = element_text(size = 13)) +
#   labs(title = "b) Great Lakes", x = "", y = "") +
#   coord_fixed()
# 
# 
# 
# plot_list <- list(p1, p2)
# 
# plot_final <- grid_arrange_shared_legend(plot_list[[1]], plot_list[[2]], 
#                                          ncol = 2, nrow = 1, position='right')


#scatterplot matrix of how the different chemistry variables relate----------------
#install.packages("GGally")
library(GGally)

variables_compare <- c("CH4_vr","DIC_vr","DOC_vr",         
                       "soilC_concen",         
                       "TSOI","VWC","soil_Fe2","soil_FeOxide",    
                       "soil_FeS","soil_O2","soil_salinity","soil_sulfate","soil_sulfide", 
                       "soil_acetate")

variables_compare <- c("CH4_vr","soil_Fe2","soil_FeOxide",    
                       "soil_FeS","soil_O2", "soil_sulfate","soil_sulfide", 
                       "soil_acetate")

variables_compare <- c("CH4_vr", "soil_Fe2","soil_FeOxide",    
                       "soil_FeS", "soil_pH")


polygon_compare <- subset(full_df_3D_long_daily, variable %in% variables_compare)
#focus just on summer for now
polygon_compare <- subset(polygon_compare, season == "Summer" & polygon_type != "Mean")

polygon_compare <- polygon_compare %>%
  pivot_wider(names_from = variable, values_from = value_daily)

polygon_compare$polygon_type <- factor(polygon_compare$polygon_type)



ggpairs(polygon_compare, columns = 8:12, aes(color = polygon_type, alpha = 0.7)) +
  scale_color_manual(values = c("HCPcenter" = "#c85c35", "HCPtransition" = "#c69a51", 
                                "LCPcenter" = "#2a7c6e", "LCPtransition"  = "#2aaeae", 
                                "Rim" = "#b17843", "Trough" = "#abba8f")) +
  scale_fill_manual(values = c("HCPcenter" = "#c85c35", "HCPtransition" = "#c69a51", 
                                "LCPcenter" = "#2a7c6e", "LCPtransition"  = "#2aaeae", 
                                "Rim" = "#b17843", "Trough" = "#abba8f")) +
  theme_bw() +
  ggtitle("Soil Chemistry")







