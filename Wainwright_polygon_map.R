
#1/28/25

#pulling in polygon data from Haruko Wainwright (Wainwright et al., 2015) to look at high res map of polygons/features for the BEO
#this is the same map that Dengel et al. (2021) used to attribute fluxes to different polygon types
#Haruko shared the data as a .txt file, col 1 is x coords, col 2 is y coords, col 3 is geomorphic type, coordinates are in UTM NAD83, units are m
#there are 9 options for geomorphic type, shown below 
#1 = LCPtrough
#2 = FCPtrough
#3 = HCPtrough
#4 = LCPcenter
#5 = FCPcenter
#6 = ?? there doesn't seem to actually be a 6? ...ask Haruko
#7 = LCPrim 
#8 = FCPrim 
#9 = HCPcenter

#note that 'high and HCP' means HCPcenter here, there is no HCPrim


#----------------------------------------------------------------------------------------------

#----------------
#install/load required packages
#----------------

library(raster) #use to open .tif files and access .tif tags
library(sf) #for using a shapefile to crop raster
library(rcartocolor)#for coming up w/colorblind friendly color palettes
library(terra) #for reprojecting raster data
#install.packages("foreign") #for reading .tif.dbf
library(foreign)
library(dplyr)

#----------------
#Hard-coded variables
#----------------

data.path <- "~/NGEE_modeling/ELM_PFLOTRAN_BEO/2_Analysis"
setwd(data.path)

# Define coordinates for the field site (these are the tower coords)
coords <- data.frame(
  lon = -156.6092, 
  lat = 71.2800    
)

#want to add point for tower location, need the UTM NAD83 version of the coordinates
# Convert the data frame to an sf object
sf_coords <- st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)

# Transform the coordinates 
utm_coords <- st_transform(sf_coords, crs = 26904)
utm_coords <- st_coordinates(utm_coords)

#-------------------------
#Load data & check it out
#-------------------------

polygon_df <- read.table("poly_feat_type_Wainwright.txt", header = FALSE, sep = " ", stringsAsFactors = FALSE)

#reads in as more than 3 cols so have to pull out just the ones w/actual data
polygon_df <- cbind(polygon_df$V4, polygon_df$V7, polygon_df$V10)
polygon_df <- data.frame(polygon_df)
  
colnames(polygon_df) <- c("X", "Y", "land_class")

#want to check out the distribution
hist(polygon_df$land_class,
     main = "Distribution of Polygon Types",
     xlab = "Class", ylab = "Number of Pixels",
     col = "springgreen")

unique(polygon_df$land_class) #shows 8   2   5   3   9  NA NaN   1   7   4

#convert df to raster format for cropping
#need to filter out NA values before converting
polygon_df_omit <- na.omit(polygon_df)
coordinates(polygon_df_omit) <- ~X + Y  # Specify the coordinate columns
crs(polygon_df_omit) <- "+proj=utm +zone=4 +datum=NAD83"  # Set the CRS to UTM NAD83

# Create a raster template
polygon_map <- raster(extent(polygon_df_omit), res = 0.5)  # Define raster resolution (e.g., 0.5m)
crs(polygon_map) <- crs(polygon_df_omit)  # Match the CRS

# Rasterize the points (assign values from 'land_name' column)
raster_layer <- rasterize(polygon_df_omit, polygon_map, field = "land_class")
plot(raster_layer)


#pull in the polygon names
polygon_df <- polygon_df %>%
  mutate(land_name = case_when(land_class == 1 ~ "LCPtrough",
                               land_class == 2 ~ "FCPtrough",
                               land_class == 3 ~ "HCPtrough",
                               land_class == 4 ~ "LCPcenter",
                               land_class == 5 ~ "FCPcenter",
                               land_class == 7 ~ "LCPrim",
                               land_class == 8 ~ "FCPrim",
                               land_class == 9 ~ "HCPcenter",
                               land_class == NaN ~ NA))


#crop map to 400mx400m box around tower, aligning w/Dengel 2021 fig 2
half_side <- 200

# Create the bounding box in projected coordinates around BEO tower
xmin <- utm_coords[1] - half_side
ymin <- utm_coords[2] - half_side
xmax <- utm_coords[1] + half_side
ymax <- utm_coords[2] + half_side

#crop full map to bounding box around flux tower
crop_extent <- extent(xmin, xmax, ymin, ymax)
polygon_map_crop <- crop(raster_layer, crop_extent)
#save the cropped version
writeRaster(polygon_map_crop, "Wainwright_poly_crop.tif", format = "GTiff", overwrite = TRUE)

# Convert raster data to a data frame for plotting...this took like 20Gb, don't do again unless cropped first!!!
polygon_df_crop <- as.data.frame(polygon_map_crop, xy = TRUE)

# Rename columns for clarity
colnames(polygon_df_crop) <- c("X", "Y", "land_class")

#we're calling FCPrim FCPcenter b/c no diff in elevation/hydrology
polygon_df_crop <- polygon_df_crop %>%
  mutate(land_name = case_when(land_class == 1 ~ "LCPtrough",
                               land_class == 2 ~ "FCPtrough",
                               land_class == 3 ~ "HCPtrough",
                               land_class == 4 ~ "LCPcenter",
                               land_class == 5 ~ "FCPcenter",
                               land_class == 7 ~ "LCPrim",
                               land_class == 8 ~ "FCPcenter",
                               land_class == 9 ~ "HCPcenter",
                               land_class == NaN ~ NA))

#calc the fractional coverage of each polygon type/feature combination
length(polygon_df_crop$land_class) #640000 total values

frac_cover <- polygon_df_crop %>%
  group_by(land_name) %>%
  summarise(type_feat_cover = (length(land_name)/640000)*100)

#add a col for categorizing by JUST polygon type instead of polygon type + feature
#so for example call everything w/FCP (FCPtrough, FCPrim, FCPcenter) just 'FCP'
polygon_df_crop <- polygon_df_crop %>%
  mutate(land_name_type = case_when(land_class == 1 ~ "LCP",
                               land_class == 2 ~ "FCP",
                               land_class == 3 ~ "HCP",
                               land_class == 4 ~ "LCP",
                               land_class == 5 ~ "FCP",
                               land_class == 7 ~ "LCP",
                               land_class == 8 ~ "FCP",
                               land_class == 9 ~ "HCP",
                               land_class == NaN ~ NA))

#calc the fractional coverage of each polygon type
frac_cover <- polygon_df_crop %>%
  group_by(land_name_type) %>%
  summarise(type_cover = (length(land_name_type)/640000)*100)

#add a col for categorizing just into wet or dry, call rims dry, troughs wet, FCPcenter and LCPcenter wet, HCPcenter dry
#FCPcenter is technically 'moist', so sometimes would be wet sometimes dry, but Lara 2020 grouped it as 'wet', so going with that
polygon_df_crop <- polygon_df_crop %>%
  mutate(land_name_hydro = case_when(land_class == 1 ~ "wet",
                               land_class == 2 ~ "wet",
                               land_class == 3 ~ "wet",
                               land_class == 4 ~ "wet",
                               land_class == 5 ~ "wet",
                               land_class == 7 ~ "dry",
                               land_class == 8 ~ "dry",
                               land_class == 9 ~ "dry",
                               land_class == NaN ~ NA))

#save the cropped polygon df
write.csv(polygon_df_crop, "Wainwright_poly_crop.csv", row.names = F)


#calc the fractional coverage of wet vs. dry
frac_cover <- polygon_df_crop %>%
  group_by(land_name_hydro) %>%
  summarise(type_cover = (length(land_name_hydro)/640000)*100)


#barplot to see landcover class distribution, filter out NA's, this is for the FULL map
ggplot(polygon_df %>% dplyr::filter(!is.na(land_name)), aes(x= land_name)) +
  #geom_bar(binwidth = 1, fill = PFT, alpha = 0.7) +
  geom_bar(aes(fill = factor(land_name)))+
  #scale_fill_manual(values = color_scheme) +
  #scale_fill_carto_d(name = "Landcover Class", palette = "Safe") +
  labs(title = "BEO Landcover Class Abundance (Wainwright)", x = NULL,
       y = "Pixels per class") +
  #theme_minimal() +
  theme_bw() +
  theme(legend.title = element_blank(), axis.text = element_text(size = 12), axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12), 
        legend.text = element_text(size = 12)) 



#plot map style
ggplot() +
  geom_raster(data = polygon_df %>% dplyr::filter(!is.na(land_name)),
              aes(x = X, y = Y, fill = factor(land_name))) +
  #scale_fill_carto_d(name = "Landcover Class", palette = "Antique") +
  #scale_fill_manual(values = color_scheme) +
  geom_point(data = utm_coords, aes(x = X, y = Y), color = "red", size = 3) + #add point for tower location
  labs(title = "BEO Landcover Classification Map (Wainwright)", x = "Easting (UTM NAD83), m", y = "Northing (UTM NAD83), m") +
  theme_bw() +
  theme(legend.title = element_blank(), axis.text = element_text(size = 12), 
        axis.title.y = element_text(size = 12), 
        legend.text = element_text(size = 12)) 
#coord_sf()



#plot map style, zoomed into 400mx400m around tower
#some of the parts that are shown as white (NA) are shown in Dengel paper as 'drainage', pretty much just the bottom right corner stuff
ggplot() +
  geom_raster(data = polygon_df_crop %>% dplyr::filter(!is.na(land_name_hydro)),
              aes(x = X, y = Y, fill = factor(land_name_hydro))) +
  #scale_fill_carto_d(name = "Landcover Class", palette = "Antique") +
  #scale_fill_manual(values = color_scheme) +
  geom_point(data = utm_coords, aes(x = X, y = Y), color = "red", size = 3) + #add point for tower location
  labs(title = "BEO Landcover Classification Map (Wainwright)", x = "Easting (UTM NAD83), m", y = "Northing (UTM NAD83), m") +
  theme_bw() +
  theme(legend.title = element_blank(), axis.text = element_text(size = 12), 
        axis.title.y = element_text(size = 12), 
        legend.text = element_text(size = 12)) 


#plot barplot again for cropped map
ggplot(polygon_df_crop %>% dplyr::filter(!is.na(land_name_hydro)), aes(x= land_name_hydro)) +
  #geom_bar(binwidth = 1, fill = PFT, alpha = 0.7) +
  geom_bar(aes(fill = factor(land_name_hydro)))+
  #scale_fill_manual(values = color_scheme) +
  #scale_fill_carto_d(name = "Landcover Class", palette = "Safe") +
  labs(title = "BEO Wet vs. Dry Abundance (Wainwright)", x = NULL,
       y = "Pixels per class") +
  #theme_minimal() +
  theme_bw() +
  theme(legend.title = element_blank(), axis.text = element_text(size = 12), axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12), 
        legend.text = element_text(size = 12)) 





