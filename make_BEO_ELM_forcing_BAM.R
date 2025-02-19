
#2/18/25: editing to updated the arctic PFT version of the surface data files

#10/4/24

#converting Ben's make_BEO_ELM_forcing.py to an R script

#-------------------------------------------------------------------------------------------------

# Install necessary libraries (if not already installed)
# install.packages("ncdf4")
#install.packages("arrayhelpers")

library(ncdf4)
library(arrayhelpers)
library(ggplot2)
library(viridis)

#keeping the 7-cell structure for now
# Define polygon types
landcover_types <- c(
  'LCPtrough',
  'LCPcenter',
  'FCPtrough',
  'FCPcenter',
  'LCPrim',
  'HCPtrough',
  'HCPcenter'
)

# Hourly time series, just do one year and the model will repeat it. Currently using water level of zero. 
ntimes <- 365 * 24  
num_grids <- length(landcover_types)  # number of grid cells

site <- 'beo'

# Open NetCDF files, pull in the domain and surface files that Ben already set up for BEO w/this code, don't need to 
#start from Fengming's original stuff again
#surfdata_multicell <- nc_open('~/GitHub/NGEE_ELM_BGC_Bailey/BEO_surfdata_multicell.nc', write = TRUE) #default PFT version
surfdata_multicell <- nc_open('~/GitHub/NGEE_ELM_BGC_Bailey/BEO_surfdata_multicell_arcticpfts.nc', write = TRUE) #arctic PFT version
# #save the surfdata contents to a text file for easy access
# {
#   sink('BEO_surfdata_multicell.txt')
#   print(surfdata_multicell)
#   sink()
# }

domain_multicell <- nc_open('~/GitHub/NGEE_ELM_BGC_Bailey/BEO_domain_multicell.nc', write = TRUE)

#quickly looking at the current values so I can make sure they get updated correctly
names(surfdata_multicell$var)
test <- ncvar_get(surfdata_multicell, "ht_above_stream")
test <- ncvar_get(surfdata_multicell, "PCT_NAT_PFT") #confirmed that there are 12 PFT options, 7 rows for the 7 gridcells being simulated

# Define new inundation values
# ht_above_stream in meters units
new_inundation_values <- c(
  0,          # gridcell 1 (LCPtrough)
  0.02991,    # gridcell 2 (LCPcenter)
  0.008905,   # gridcell 3 (FCPtrough)
  0.1228581,  # gridcell 4 (FCPcenter)
  0.11896484, # gridcell 5 (LCPrim)
  0.01781,    # gridcell 6 (HCPtrough)
  0.21580621  # gridcell 7 (HCPcenter)
)

# Update the ht_above_stream variable with new values
ncvar_put(surfdata_multicell, "ht_above_stream", new_inundation_values)

#check that it worked
test <- ncvar_get(surfdata_multicell, "ht_above_stream") #yep! adjusted them all to be 7 digits past decimal

# Set dist_from_stream to 1.0 (same size as ht_above_stream)
dist_from_stream <- rep(1.0, length(new_inundation_values))
ncvar_put(surfdata_multicell, "dist_from_stream", dist_from_stream)

#check original ORGANIC values
test <- ncvar_get(surfdata_multicell, "ORGANIC")

# Define new organic matter values (7x10 matrix)
new_organic_values <- matrix(c(
  102.973, 56.927, 68.51975, 43.667, 42.3306, 45.864, 22.02525, 6.24, 6.24, 6.24,  # gridcell 1 (LCPtrough)
  102.973, 56.927, 68.51975, 43.667, 42.3306, 45.864, 22.02525, 6.24, 6.24, 6.24,  # gridcell 2 (LCPcenter)
  101.452, 101.452, 101.452, 29.237, 28.886, 50.895, 13.61285714, 13.624, 5.27475, 5.148,  # gridcell 3 (FCPtrough)
  101.452, 101.452, 101.452, 29.237, 28.886, 50.895, 13.61285714, 13.624, 5.27475, 5.148,  # gridcell 4 (FCPcenter)
  102.973, 56.927, 68.51975, 43.667, 42.3306, 45.864, 22.02525, 6.24, 6.24, 6.24,  # gridcell 5 (LCPrim)
  49.621, 49.621, 49.621, 80.314, 62.16166667, 73.53666667, 52.71825, 10.647, 6.097, 6.097,  # gridcell 6 (HCPtrough)
  120.081, 120.081, 11.895, 59.0525, 62.16166667, 73.53666667, 52.71825, 10.647, 6.56175, 4.706  # gridcell 7 (HCPcenter)
), nrow = 7, byrow = TRUE)

# Update the ORGANIC variable with new values
ncvar_put(surfdata_multicell, "ORGANIC", new_organic_values)

# Set fdrain and slope variables
fdrain <- rep(100.0, num_grids)
ncvar_put(surfdata_multicell, "fdrain", fdrain)

slope <- rep(0.05, num_grids)
ncvar_put(surfdata_multicell, "SLOPE", slope)

# change soil order from 2 (andisols) to 3 (gelisols)
#see https://passel2.unl.edu/view/lesson/2eafec8dd762/3#:~:text=This%20lesson%20will%20examine%20each,Histosols%2C%20Aridisols%2C%20and%20Vertisols. 
order_values <- rep(3, num_grids)
ncvar_put(surfdata_multicell, "SOIL_ORDER", order_values)

#check original veg fractional coverage values
#test <- ncvar_get(surfdata_multicell, "PCT_NAT_PFT")

# Define new vegetation fractional cover values (7x17 matrix)
#these get rounded to different positions past the decimal, which makes rows not add to 0...
#try rounding them all to 6 points past decimal
# new_veg_cover_values <- matrix(c(
#   5.241799121, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 94.75820088, 0, 0, 0, 0,  # gridcell 1 (LCPtrough)
#   43.04960196, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 56.95039804, 0, 0, 0, 0,  # gridcell 2 (LCPcenter)
#   2.781795927, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 97.21820407, 0, 0, 0, 0,  # gridcell 3 (FCPtrough)
#   27.09212377, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 72.90787623, 0, 0, 0, 0,  # gridcell 4 (FCPcenter)
#   11.87843147, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.040826416, 87.08074211, 0, 0, 0, 0,  # gridcell 5 (LCPrim)
#   15.46722699, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 84.53277301, 0, 0, 0, 0,  # gridcell 6 (HCPtrough)
#   45.68772594, 0, 0, 0, 0, 0, 0, 0, 0, 11.4764446, 0, 0, 42.83582946, 0, 0, 0, 0  # gridcell 7 (HCPcenter)
# ), nrow = 7, byrow = TRUE)
# 
# new_veg_cover_values <- round(new_veg_cover_values, 6)

#THIS ONE IS FOR DEFAULT PFTS!!!!!!!!!!!!!!!!!!
#below is version rounded to 6 decimal places, still didn't add to 100 perfectly so need to manually adjust some rows
#its looking out to 6 decimal places, so values are super close but for example row 1 sums to 99.999999
new_veg_cover_values <- matrix(c(
  5.241800, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 94.75820, 0, 0, 0, 0,  # gridcell 1 (LCPtrough)
  43.049600, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 56.95040, 0, 0, 0, 0,  # gridcell 2 (LCPcenter)
  2.781800, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 97.21820, 0, 0, 0, 0,  # gridcell 3 (FCPtrough)
  27.092120, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 72.90788, 0, 0, 0, 0,  # gridcell 4 (FCPcenter)
  11.878434, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.040826, 87.08074, 0, 0, 0, 0,  # gridcell 5 (LCPrim)
  15.467230, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 84.53277, 0, 0, 0, 0,  # gridcell 6 (HCPtrough)
  45.687720, 0, 0, 0, 0, 0, 0, 0, 0, 11.47645, 0, 0, 42.83583, 0, 0, 0, 0  # gridcell 7 (HCPcenter)
), nrow = 7, byrow = TRUE)

#this is the version for the arctic PFTs!!!!!!!!!!!!!!!!!
#12 cols (PFT options)...but shouldn't there be 14? which are not used?
#check this file to find out the PFT names associated w/each number
#arctic_params <- nc_open('~/GitHub/NGEE_ELM_BGC_Bailey/clm_params.c141105.arctic.nc', write = TRUE)...NOPE! look at how Ben set it up: https://github.com/bsulman/NGEE_ELM_BGC/blob/master/make_BEO_ELM_forcing.py 
#arctic model setup is as follows: 
# [1] "not_vegetated                           
# [2] "arctic_lichen                           
# [3] "arctic_bryophyte                        
# [4] "arctic_evergreen_shrub_dwarf            
# [5] "arctic_evergreen_shrub_tall             
# [6] "arctic_deciduous_shrub_dwarf            
# [7] "arctic_deciduous_shrub_low              
# [8] "arctic_deciduous_shrub_tall             
# [9] "arctic_deciduous_shrub_alder            
# [10] "arctic_forb                             
# [11] "arctic_dry_graminoid                    
# [12] "arctic_wet_graminoid                    

#these are the adjusted values rounded to 6 decimal places
new_veg_cover_values <- matrix(c(
  5.213617,	0.028182,	61.185323,	0,	0,	0,	0,	0,	0,	1.186450,	0.084545,	32.301883,  # gridcell 1 (LCPtrough)
  43.018983,	0.030618,	25.541947,	0,	0,	0,	0,	0,	0,	1.108389,	1.258420,	29.041643,  # gridcell 2 (LCPcenter)
  2.781796,	0,	61.650161,	0,	0,	0,	0,	0,	0,	15.639257,	1.235117,	18.693669,  # gridcell 3 (FCPtrough)
  24.027191,	3.064932,	24.572199,	0,	0,	0,	0,	0,	0,	3.768167,	9.851149,	34.716362,  # gridcell 4 (FCPcenter)
  3.642892,	8.235540,	46.342796,	0,	0,	0,	1.040826,	0,	0,	8.381255,	3.622076,	28.734615,  # gridcell 5 (LCPrim)
  14.299889,	1.167338,	49.728594,	0,	0,	0,	0,	0,	0,	10.290083,	1.926108,	22.587988,  # gridcell 6 (HCPtrough)
  13.197911,	32.489813,	17.386814,	11.476445,	0,	0,	0,	0,	0,	1.382912,	21.650313,	2.415792  # gridcell 7 (HCPcenter)
), nrow = 7, byrow = TRUE)


# Update PCT_NAT_PFT variable with new values
ncvar_put(surfdata_multicell, "PCT_NAT_PFT", new_veg_cover_values)

# Ensure that PFT sums to 100
veg_sum <- rowSums(new_veg_cover_values)
if (all(veg_sum == 100)) {
  print("all good, PFT rows sum to 100")
} else {
  print("ERROR, PFT rows don't sum to 100")
}


# close the old NetCDF files...nc_close() saves the changes to the file specified in nc_open()
#so file name doesn't get updated, but this is actually probably easier...
nc_close(surfdata_multicell)
nc_close(domain_multicell)






#generate plots of the inundation to soil height level, veg fractional coverage, and soil OM------------------------------------


#first plot OM by soil layer------------------------------------------------------
OM <- ncvar_get(surfdata_multicell, "ORGANIC")

OM <- data.frame(OM)
names(OM) <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
OM$polygon_type <- landcover_types

# Convert from wide to long format
OM <- OM %>%
  pivot_longer(
    cols = 1:10,  # Specify columns to pivot
    names_to = "soil_layer",            
    values_to = "value"   
  )


fill_name <- "OM Density \n(kg m-3)"
#need to convert soil depth to factor so it extends the width of the line and fills out the plot
g1 <- ggplot(data = OM, aes(x = polygon_type, y = soil_layer, fill = value))+ 
  #scale_y_continuous(trans = "reverse") + #breaks = seq(1, 10, 1)) +
  scale_y_discrete(limits = rev(unique(OM$soil_layer))) + #for characters
  labs(x = NULL, y = "Soil Layer",
       fill = fill_name, title = "Updated Model Input")

p1 <- g1 + geom_tile(aes(fill = value)) + #color= "white",size=0.1, alpha = .7, 
  scale_fill_viridis(option ="C") + #limits=c(rng[1], rng[2])) +
  theme_minimal() +
  theme(axis.text = element_text(size = 12), #axis.title.x = element_blank(), axis.title.y = element_blank(), 
        legend.text = element_text(size = 13), legend.title = element_text(size = 13),
        title = element_text(size = 13), strip.text.y = element_text(size = 10), legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1))


p1


#------------------------------------------------------------------------------

#next plot veg fractional coverage 

veg <- ncvar_get(surfdata_multicell, "PCT_NAT_PFT")
veg <- data.frame(veg)

#if using default PFTs
names(veg) <- c("PFT1", "PFT2", "PFT3", "PFT4", "PFT5", "PFT6", "PFT7", "PFT8", "PFT9", "PFT10", "PFT11", "PFT12", 
                "PFT13", "PFT14", "PFT15", "PFT16", "PFT17")

veg$polygon_type <- landcover_types

# Convert from wide to long format
veg <- veg %>%
  pivot_longer(
    cols = PFT1:PFT17,  # Specify columns to pivot
    names_to = "PFT",            
    values_to = "value"   
  )

#if using arctic PFTs
names(veg) <- c("PFT1", "PFT2", "PFT3", "PFT4", "PFT5", "PFT6", "PFT7", "PFT8", "PFT9", "PFT10", "PFT11", "PFT12")

veg$polygon_type <- landcover_types

# Convert from wide to long format
veg <- veg %>%
  pivot_longer(
    cols = PFT1:PFT12,  # Specify columns to pivot
    names_to = "PFT",            
    values_to = "value"   
  )



#fill 0's w/NA so all the PFTs that aren't represented don't get included in the legend
veg <- veg %>%
  mutate(value = na_if(value, 0))

#for default PFTs
pft_colors <- c("PFT1" = "#888888", "PFT10" = "#117733",
                "PFT12" = "#967acc", "PFT13" = "#6699cc")

#for arctic PFTs
#for current BEO setup PFTs 5, 6, 8, 9 aren't represented
pft_colors <- c("PFT1" = "#888888", "PFT2" = "#ddcc77", "PFT3" = "#999933", "PFT4" = "#117733",
                "PFT7" = "#967acc", "PFT10" = "#78b33e", "PFT11" = "#cc6677", "PFT12" = "#6699cc")

ggplot(veg, aes(x = factor(polygon_type), y = value)) +
  geom_bar(stat = "identity", aes(fill = PFT)) +
  scale_fill_manual(values = pft_colors, 
                    labels = c(`PFT1` = "Not Vegetated (PFT1)", `PFT2` = "Arctic Lichen (PFT2)", `PFT3` = "Arctic Bryophyte (PFT3)",
                               `PFT4` = "Arctic Evergreen Shrub Dwarf (PFT4)", `PFT7` = "Arctic Deciduous Shrub Low (PFT7)", 
                               `PFT10` = "Arctic Forb (PFT10)", `PFT11` = "Arctic Dry Graminoid (PFT11)", 
                               `PFT12` = "Arctic Wet Graminoid (PFT12)")) +
  labs(x = NULL, y = "Fractional Coverage", title = "Updated Model Input (BEO, Arctic PFTs)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        legend.title = element_blank(), axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 15), axis.title = element_text(size = 15)) 


pft_colors <- c("PFT1" = "#888888", "PFT10" = "#117733",
                "PFT12" = "#967acc", "PFT13" = "#6699cc")





#next plot inundation level---------------------------------------------

water <- ncvar_get(surfdata_multicell, "ht_above_stream")

water <- data.frame(water, landcover_types)

names(water) <- c("value", "polygon_type")


ggplot(water, aes(x = factor(polygon_type), ymin = -0.1, ymax = value)) +
  geom_linerange(linewidth = 13, color = "tan4") +
  coord_cartesian(ylim = c(-0.10, max(water$value) + 0.10)) + # Adjust y-axis limits
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "steelblue", linewidth = 3) +  # Add horizontal line to show tide height
  labs(x = NULL, y = "Height [m]", title = "Updated Model Input") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        legend.title = element_blank(), axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 15), axis.title = element_text(size = 15))












