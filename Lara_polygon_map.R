
#1/13/25

#playing around w/Lara et al. (2018) polygon map for BEO

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

#----------------
#Hard-coded variables
#----------------

data.path <- "~/NGEE_modeling/ELM_PFLOTRAN_BEO/2_Analysis/Lara_polygon_map"
setwd(data.path)

# Define coordinates for the field site (these are the tower coords)
field_site_coords <- data.frame(
  lon = -156.6092, 
  lat = 71.2800    
)

# # Convert the data frame to an sf object
# coords_sf <- st_as_sf(field_site_coords, coords = c("lon", "lat"), crs = 4326)  # WGS84 (EPSG:4326)
# 
# # Transform to Albers Equal Area projection (example EPSG: 5070 for North America)
# coords_albers <- st_transform(coords_sf, crs = 5070)
# 
# # View transformed coordinates
# print(coords_albers)
# 
# # Extract coordinates as a data frame
# coords_albers_df <- st_coordinates(coords_albers)
# print(coords_albers_df)



#-------------------------
#Load data & check it out
#-------------------------

# open raster data
polygon_map <- raster(x = "acp_112016.tif") 

#look at the associated database file
dbf_data <- read.dbf(paste0(data.path, "/acp_112016.tif.dbf"))
print(dbf_data)


#look at coordinate reference system, assign to an object so can use later if needed
#says projection is aea (alaska albers equal area), units are in meters
#projection reference: https://www.earthdatascience.org/courses/earth-analytics/lidar-raster-data-r/open-lidar-raster-r/
my_crs <- crs(polygon_map)

#look at x and y resolution, pixels are 30X30m
xres(polygon_map)
yres(polygon_map)

# Input point in latitude and longitude
point_latlon <- st_sfc(st_point(c(-156.6092, 71.28)), crs = 4326) 

# Transform to a projected CRS (e.g., Albers Equal Area)
point_proj <- st_transform(point_latlon, crs = 3338) # NAD83 / Alaska Albers

# Define the bounding box dimensions (10 km = 5 km on each side from the center)
half_side <- 5000 # Half of 10 km in meters
#trying w/400mx400m box around tower, aligning w/Dengel 2021 fig 2
half_side <- 200

# Create the bounding box in projected coordinates that is 100km x 100km around BEO tower lat/lon
bbox_proj <- st_bbox(c(
  xmin = st_coordinates(point_proj)[1] - half_side,
  ymin = st_coordinates(point_proj)[2] - half_side,
  xmax = st_coordinates(point_proj)[1] + half_side,
  ymax = st_coordinates(point_proj)[2] + half_side
), crs = st_crs(point_proj))

#crop full map to bounding box around flux tower
crop_extent <- extent(bbox_proj) 
polygon_map_crop <- crop(polygon_map, crop_extent)

#don't do this yet!! this is what is messing up the categories and making them continuous I think
# # Define the target CRS (latitude and longitude, WGS84)
# target_crs <- CRS("+proj=longlat +datum=WGS84 +no_defs")
# 
# # Reproject the raster to latitude and longitude
# polygon_map <- projectRaster(polygon_map_crop, crs = target_crs)

#spatial extent of file:
st_bbox(polygon_map)

hist(polygon_map,
     main = "Distribution of land cover classes",
     xlab = "Class", ylab = "Number of Pixels",
     col = "springgreen")



# Convert raster data to a data frame for plotting...this took like 20Gb, don't do again unless cropped first!!!
polygon_df <- as.data.frame(polygon_map_crop, xy = TRUE)

# Rename columns for clarity
colnames(polygon_df) <- c("x", "y", "land_class")

# Replace NA values with 0, ocean gets called NA, want to reassign
polygon_df[is.na(polygon_df)] <- 0

#Looks like a bunch of the lakes are labelled as NA, just replace NA with ‘water’, and drop the small/med/large lake
#attach names to the numbers in the df version
polygon_df <- polygon_df %>%
  mutate(land_name = case_when(land_class == 2 ~ "Coalescent Low-center polygon",
                               land_class == 3 ~ "Low-center polygon",
                               land_class == 4 ~ "Drained thaw lake basin",
                               land_class == 7 ~ "Drained slope",
                               land_class == 8 ~ "Flat-center polygon",
                               land_class == 9 ~ "High-center polygon",
                               land_class == 10 ~ "Sandy barren",
                               land_class == 11 ~ "Sand dune",
                               land_class == 12 ~ "Ice",
                               land_class == 13 ~ "Coastal saline water",
                               land_class == 16 ~ "River",
                               land_class == 17 ~ "Urban",
                               land_class == 18 ~ "Riparian corridor",
                               land_class == 20 ~ "Pond",
                               land_class == 21 ~ "Small lake",
                               land_class == 22 ~ "Medium lake",
                               land_class == 23 ~ "Large lake",
                               land_class == 0 ~ "Ocean"))

#hard to tell things apart on map, try replacing ocean with white and using one color for lake/pond
polygon_df <- polygon_df %>%
  mutate(land_name_sub = case_when(land_class == 2 ~ "Coalescent Low-center polygon",
                                   land_class == 3 ~ "Low-center polygon",
                                   land_class == 4 ~ "Drained thaw lake basin",
                                   land_class == 7 ~ "Drained slope",
                                   land_class == 8 ~ "Flat-center polygon",
                                   land_class == 9 ~ "High-center polygon",
                                   land_class == 10 ~ "Sandy barren",
                                   land_class == 11 ~ "Sand dune",
                                   land_class == 12 ~ "Ice",
                                   land_class == 13 ~ NA,
                                   land_class == 16 ~ "River",
                                   land_class == 17 ~ "Urban",
                                   land_class == 18 ~ "Riparian corridor",
                                   land_class == 20 ~ "Lake",
                                   land_class == 21 ~ "Lake",
                                   land_class == 22 ~ "Lake",
                                   land_class == 23 ~ "Lake",
                                   land_class == 0 ~ NA))


#save just so have if need later (to plot in ggplot need it as df)
#write.csv(polygon_df, "polygon_df_100kmx100km.csv", row.names = FALSE)

length(unique(polygon_df$land_class)) #16 unique landcover classes in the 10kmx10km domain

#barplot to see landcover class distribution, filter out NA's
ggplot(polygon_df %>% dplyr::filter(!is.na(land_name_sub)), aes(x= land_name_sub)) +
  #geom_bar(binwidth = 1, fill = PFT, alpha = 0.7) +
  geom_bar(aes(fill = factor(land_name_sub)))+
  scale_fill_manual(values = color_scheme) +
  #scale_fill_carto_d(name = "Landcover Class", palette = "Safe") +
  labs(title = "Landcover Class Abundance (400m X 400m)", x = NULL,
       y = "Pixels per class") +
  #theme_minimal() +
  theme_bw() +
  theme(legend.title = element_blank(), axis.text = element_text(size = 12), axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12), 
        legend.text = element_text(size = 12)) 



#plot map style
#k now there are 13 categories, set up useful colors
#color_scheme <- c("#357d94", "#b7c8b5", "#c5c9bb", "#586f45", "#9a3c35", "#f0d5b7", "#477379", "#a98f5f", "#557067", "#c9bcb5", "#947b6f", "#5d686a")
color_scheme <- c("#855c75", "#d9af6b", "#af6458", "#625377", "#526a83", "#c2d8d8", "#68855c", "#9c9c5e", "#a06177", "#8c785d", "#467378", "#7c7c7c") #these are colorblind friendly

#want to add point for tower location
field_site_coords <- data.frame(
  x = -95812.38, 
  y = 2366607    
)

ggplot() +
  geom_raster(data = polygon_df %>% dplyr::filter(!is.na(land_name_sub)),
              aes(x = x, y = y, fill = factor(land_name_sub))) +
  #scale_fill_carto_d(name = "Landcover Class", palette = "Antique") +
  scale_fill_manual(values = color_scheme) +
  geom_point(data = field_site_coords, aes(x = x, y = y), color = "red", size = 3) + #add point for tower location
  labs(title = "Landcover Classification Map", x = NULL,
       y = NULL) +
  theme_bw() +
  theme(legend.title = element_blank(), axis.text = element_text(size = 12), 
        axis.title.y = element_text(size = 12), 
        legend.text = element_text(size = 12)) 
coord_sf()





















