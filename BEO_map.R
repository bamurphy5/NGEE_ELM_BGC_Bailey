
#1/10/25

#playing around with making a map of AK with an insert showing the location of the BEO study site and another insert 
#showing a picture of the tundra landscape

# Load libraries
library(ggplot2)
#install.packages("ggspatial")
library(ggspatial)
library(cowplot)
#install.packages("ggmap")
library(ggmap)
#install.packages("rnaturalearth")
library(rnaturalearth) #maps() didn't have data for ak so trying this one
library(sf)
#install.packages("jpeg")
library(jpeg)

#downloaded base maps from: https://www.census.gov/cgi-bin/geo/shapefiles/index.php?year=2024&layergroup=States+%28and+equivalent%29 
base_map_path <- "~/NGEE_modeling/ELM_PFLOTRAN_BEO/map_shapefiles/tl_2024_us_state.shp"

# Define coordinates for the field site (these are the tower coords)
field_site_coords <- data.frame(
  lon = -156.6092, 
  lat = 71.2800    
)

# Import base map of Alaska
#alaska_map <- map_data("state") %>% filter(region == "alaska")
#world <- ne_states(country = "United States", returnclass = "sf")
#alaska <- world %>% filter(name == "Alaska")
us_states <- st_read(base_map_path)
# Filter for Alaska
alaska <- subset(us_states, NAME == "Alaska")

# test plot using the geometry column (geometry is 1 of the 15 attributes that come in the shapefile)
#there's a weird island included, will have to set a longitude limit I think
ggplot() +
  geom_sf(data = filtered_alaska, aes(geometry = geometry), fill = "lightblue", color = "black") +
  theme_minimal() +
  labs(title = "Map of Alaska", x = "Longitude", y = "Latitude")

#filter the alaska map to cut off the weird island
# Define the longitude range
longitude_min <- -179
longitude_max <- -120

# Define the bounding box for the longitude range
bbox <- st_bbox(c(
  xmin = longitude_min,
  xmax = longitude_max,
  ymin = st_bbox(alaska)[["ymin"]],  # Keep the full latitude range
  ymax = st_bbox(alaska)[["ymax"]]
), crs = st_crs(alaska))  # Use the same CRS as the Alaska map

# Crop the map data
filtered_alaska <- st_crop(alaska, bbox)


# Define a bounding box for the zoomed-in area around the field site
zoom_bbox <- st_bbox(c(
  xmin = field_site_coords$lon - 0.8, 
  xmax = field_site_coords$lon + 0.8, 
  ymin = field_site_coords$lat - 0.25, 
  ymax = field_site_coords$lat + 0.5
), crs = st_crs(alaska))

zoom_bbox <- st_bbox(c(
  xmin = field_site_coords$lon - 1.5, 
  xmax = field_site_coords$lon + 1.5, 
  ymin = field_site_coords$lat - 1, 
  ymax = field_site_coords$lat + .5
), crs = st_crs(alaska))

# Crop the Alaska map to the zoomed-in area
zoom_area <- st_crop(filtered_alaska, zoom_bbox)

# Convert bbox to an sf polygon (for plottin on top of base map)
zoom_bbox_sf <- st_as_sfc(zoom_bbox)

# Plot the base map of Alaska
base_map <- ggplot() +
  geom_sf(data = filtered_alaska, fill = "lightblue", color = "black") +
  #geom_point(data = field_site_coords, aes(x = lon, y = lat), color = "red", size = 3) +
  geom_sf(data = zoom_bbox_sf, color = "red", fill = NA, linewidth = 1.2) +  # Red outline for zoomed portion of map
  theme_minimal() +
  labs(title = "a)",x = "Longitude", y = "Latitude") +
  annotation_scale(location = "bl", width_hint = 0.5) + # Add a scale bar
  annotation_north_arrow(location = "tl", which_north = "true", 
                         style = north_arrow_fancy_orienteering())
base_map

# Extract longitude range
lon_range <- st_bbox(zoom_area)[c("xmin", "xmax")]

# Create the zoomed-in map
zoomed_map <- ggplot() +
  geom_sf(data = zoom_area, fill = "lightgreen", color = "black") +
  geom_point(data = field_site_coords, aes(x = lon, y = lat), color = "red", size = 3) +
  geom_sf(data = zoom_bbox_sf, color = "red", fill = NA, linewidth = 1.2) +
  theme_minimal() +
  scale_x_continuous(
    breaks = seq(floor(lon_range[1]), ceiling(lon_range[2]), by = 1)) +
  labs(title = "b)", x = NULL, y = NULL) +
  annotation_scale(location = "bl", width_hint = 0.3) 
  #annotation_north_arrow(location = "tl", which_north = "true",
                         #style = north_arrow_fancy_orienteering())

zoomed_map

#Read in the landscape photo
landscape_image_path <- "~/NGEE_modeling/ELM_PFLOTRAN_BEO/map_shapefiles/polygonal_tundra.jpg"
landscape_image <- readJPEG(landscape_image_path)

# Create a ggplot object with the photograph
landscape_photo <- ggplot() +
  annotation_raster(landscape_image, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void() +
  labs(title = "c)")


# Combine the base map with zoomed-in inserts using cowplot
final_plot <- ggdraw() +
  # Add the main Alaska map
  #draw plot uses (thing_to_plot, x location of plot, y location of plot, width of plot, height of plot)
  #0 for x location is left side, 0 for y location is bottom 
  draw_plot(base_map, 0, .5, .5, .5) +
  # Add the zoomed-in field site map (positioned at top-right)
  draw_plot(zoomed_map, 0.5, 0.4, 0.4, 0.4) +
  # Add the landscape photograph (positioned below the zoomed-in map)
  draw_plot(landscape_photo, 0.55, 0.1, 0.3, 0.3)

final_plot

#trying a version w/o the photo so it looks a little cleaner
# Combine the base map with zoomed-in inserts using cowplot
final_plot <- ggdraw() +
  # Add the main Alaska map
  #draw plot uses (thing_to_plot, x location of plot, y location of plot, width of plot, height of plot)
  #0 for x location is left side, 0 for y location is bottom 
  draw_plot(base_map, 0, .3, .6, .6) +
  # Add the zoomed-in field site map (positioned at top-right)
  draw_plot(zoomed_map, 0.35, 0.2, 0.7, 0.7) +
  geom_segment(aes(x = .3, xend = .6, y = .8, yend = .26), color = "red", linewidth = 1.2, alpha = .5) +
  geom_segment(aes(x = .3, xend = .60, y = .83, yend = .78), color = "red", linewidth = 1.2, alpha = .5)
  

final_plot

