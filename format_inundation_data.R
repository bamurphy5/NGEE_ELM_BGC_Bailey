
#9/10/24

#working on formatting the water depth data from Liljedahl and Wilson data for BEO, collected 2012-2014: https://data.ess-dive.lbl.gov/view/doi:10.5440/1183767  
#using this to set up the input files for inundation depth for simulations

#data collected over three summers, measurements taken at 45 locations (polygon centers and troughs), measured manually and w/instruments
#data collected in areas A, B, C, D, need to cross-check with well ID's and see what polygon type/feature each well ID corresponds to

#data is formatted as separate file for each well, need to pull it all together
#inundation level also needs to be calculated, subtract ground surface height from water level to calculate. 
#Water_level and Ground_surface are both in units of meters above sea level (masl)
#missing data represented by -9999

dat.base <- "~/NGEE_modeling/ELM_PFLOTRAN_BEO/2_Analysis/observational_data/inundation/data/water_level_data"
dat.base.main <- "~/NGEE_modeling/ELM_PFLOTRAN_BEO/2_Analysis/observational_data/inundation/data"

# Get a list of all files in the folder
file_list <- list.files(path = dat.base)
#file_list <- file_list[1:3] #trying short version to test

#set up empty dataframe
inundation_data <- data.frame()

for (file in file_list){
  data <- read.csv(file.path(dat.base, file), header = TRUE)
  file_name <- sub("\\.csv$", "", file) #remove extension before splitting
  
  # Split the filename by underscores
  parts <- unlist(strsplit(file_name, "_"))

  # Extract the desired parts (area and well ID)
  area <- parts[6]  
  well_ID <- parts[7]   
  
  data$Site <- "BEO"
  data$Area <- area
  data$Well_ID <- well_ID
  
  inundation_data <- rbind(inundation_data, data)
  
}


#pull in the coordinate data
#compiled coordinates for all well locations in one file, "BEO_Well_Water_Well_Coordinates_ALL.csv"
#lat/lon are in UTM
#combine at end, match by well_ID
coords <- read.csv(file.path(dat.base.main, "BEO_Well_Water_Well_Coordinates_ALL.csv"), header = TRUE)

#drop the elevation col, already have
coords$Field_Surveyed_Elevation <- NULL
inundation_data <- merge(inundation_data, coords, by = c("Area", "Well_ID"))

#replace -9999 with NA
inundation_data$Water_Level <- replace(inundation_data$Water_Level, inundation_data$Water_Level == -9999, NA)
inundation_data$Ground_Surface <- replace(inundation_data$Ground_Surface, inundation_data$Ground_Surface == -9999, NA)

#one random outlier in Water_level still throwing things off, replace w/NA (looks like someone entered -999 instead of -9999...)
inundation_data$Water_Level[inundation_data$Water_Level == min(inundation_data$Water_Level, na.rm = T)] <- NA

#calculate inundation level as the difference between water level and ground elevation
inundation_data$Inundation_level <- inundation_data$Water_Level - inundation_data$Ground_Surface

#there's a random single outlier of -4m for inundation level, replace w/NA
inundation_data$Inundation_level[inundation_data$Inundation_level == min(inundation_data$Inundation_level, na.rm = T)] <- NA


#now need to add in the corresponding polygon type/feature info for each well ID
#polygon type/feature classification for each well ID are annoyingly not included in the dataset, or in the information file for 
#each area (other than random mentions of a couple wells in _information.txt files for each location). 
#Had to look in Liljedahl 2016 supplement…but looks like site/well names are not consistent. SI also only has 29 locations, whereas water dataset has 45 locations. 
#Maybe can just see what matches from the supplement, then look at the map w/elevation shading to guess at the others and check w/Neslihan?

#read in the locations from Liljedahl 2016 supplement
#Liljedahl_loc <- read.csv(file.path(dat.base.main, "Liljedahl_2016_locations.csv"), header = TRUE)

#merge based on location, resulting file will have only the rows where northing/easting are the same in both df's
#matching_rows <- merge(coords, Liljedahl_loc, by = c("Northing_UTM", "Easting_UTM"))

#...k there are no exact matches...the water level files say location info is approximate, so might need to plot everything in google earth and see whats 
#closest to what. Coords in Liljedahl 2016 SI can’t be correct, showing up as off the coast of Russia…
#Looking at the paper there’s only high and low centered polygons, no FCP. Looking at the maps in the water level dataset I’d guess that areas A and D are LCP and areas B and C are HCP. 
#Just going to visually determine whether each location seems to be in a trough or center, and assign it myself.
#adding polygon info to the BEO_Well_Water_Well_Coordinates_ALL.csv file
Well_ID_list <- unique(inundation_data$Well_ID) #45 locations

#read the file in again now that it has the polygon types added, and just run through calculating inundation again
#once complete, go in and remove the well ID's that I couldn't classify (weren't on the map)
inundation_data <- inundation_data[!is.na(inundation_data$Polygon_type), ] #removed 105,315 rows

#save the data before grouping and averaging
write.csv(inundation_data, file.path(dat.base.main, "BEO_inundation_ALL.csv"), row.names = FALSE)
#coming back later to pull in the inundation data...
inundation_data <- read.csv(file.path(dat.base.main, "BEO_inundation_ALL.csv"), header = TRUE)

#find the minimum ground elevation
min_ground <- min(inundation_data$Ground_Surface, na.rm = TRUE)
#calculate an adjusted ground height based on difference from this minimum value, so this minimum value is the 'zero'
inundation_data$Ground_Surface_adj <- inundation_data$Ground_Surface - min_ground
#also calc adjusted water level based on this minimum value
inundation_data$Water_Level_adj <- inundation_data$Water_Level - min_ground

#random -3.31 value, replace w/NA
inundation_data$Water_Level_adj[inundation_data$Water_Level_adj == min(inundation_data$Water_Level_adj, na.rm = T)] <- NA

#datetime col is actually a character, so convert to a date object
inundation_data$Datetime_2 <- as.Date(inundation_data$Datetime) #just date don't really care about time right now

#add a col for month
inundation_data$Month <- month(inundation_data$Datetime_2)


#next group by polygon type and find avg inundation level 
#would it be better to first avg by year? or just throw it all in?
inundation_data_avg <- inundation_data %>%
  group_by(Polygon_type, Year) %>%
  summarise(Inundation_level_avg = mean(Inundation_level, na.rm = TRUE))

#trying avg'ing by year first 
inundation_data_avg <- inundation_data_avg %>%
  group_by(Polygon_type) %>%
  summarise(Inundation_level_avg2 = mean(Inundation_level_avg, na.rm = TRUE))
  
  
  

#tide height is set at 0.10, add this
inundation_data_avg$water_level <- 0.10
inundation_data_avg$soil_height <- inundation_data_avg$water_level - inundation_data_avg$Inundation_level_avg2


#plot to see what it all looks like------------------------------------------------
ggplot(inundation_data_avg, aes(x = factor(Polygon_type), y = soil_height)) +
  geom_bar(stat = "identity", aes(fill = Polygon_type)) +
  #scale_fill_manual(values = pft_colors) +
  #scale_x_discrete(labels=c("Trough", "LCPcenter", "LCPtransition", "HCPcenter", "HCPtransition")) +
  coord_cartesian(ylim = c(-0.10, max(inundation_data_avg$soil_height) + 0.10)) + # Adjust y-axis limits
  labs(x = NULL, y = "Inundation Level (m)", title = "Obs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        legend.title = element_blank(), axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 15), axis.title = element_text(size = 15))

#9/17/24: pull in adjusted soil height data (see slides)
adj_soil_height <- read.csv(file.path(dat.base.main, "adjusted_polygon_soil_height_BEO.csv"), header = TRUE)


#plot with bars that start from -0.1 and extend to the mean for each group
ggplot(inundation_data_avg, aes(x = factor(Polygon_type), ymin = -0.1, ymax = soil_height)) +
  geom_linerange(linewidth = 13, color = "tan4") +
  coord_cartesian(ylim = c(-0.10, max(inundation_data_avg$soil_height) + 0.10)) + # Adjust y-axis limits
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "steelblue", linewidth = 3) +  # Add horizontal line to show tide height
  labs(x = NULL, y = "Inundation Level (m)", title = "Obs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        legend.title = element_blank(), axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 15), axis.title = element_text(size = 15))

#use the outside adjusted soil height file
ggplot(adj_soil_height, aes(x = factor(Polygon_type), ymin = -0.1, ymax = soil_height)) +
  geom_linerange(linewidth = 13, color = "tan4") +
  coord_cartesian(ylim = c(-0.10, max(adj_soil_height$soil_height) + 0.10)) + # Adjust y-axis limits
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "steelblue", linewidth = 3) +  # Add horizontal line to show tide height
  labs(x = NULL, y = "Height [m]", title = "Obs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        legend.title = element_blank(), axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 15), axis.title = element_text(size = 15))





#plot water level over time, grouped by polygon type

# Add a day of year column to plot time continuously
inundation_data$DOY <- yday(inundation_data$Datetime_2)  # Extract day of the year
#there are multiple measurements per day sometimes, so calc daily avg
inundation_data_daily <- inundation_data %>%
  group_by(Year, Polygon_type, DOY) %>%
  summarise(Water_Level = mean(Water_Level, na.rm = TRUE), Ground_Surface = mean(Ground_Surface, na.rm = TRUE),
         Inundation_level = mean(Inundation_level, na.rm = TRUE), Ground_Surface_adj = mean(Ground_Surface_adj, na.rm = TRUE),
         Water_Level_adj = mean(Water_Level_adj, na.rm = TRUE))


inundation_data_daily$Inundation_level_adj <- inundation_data_daily$Water_Level_adj - inundation_data_daily$Ground_Surface_adj

ggplot(inundation_data_daily, aes(x = DOY, y = Inundation_level, color = Polygon_type, group = Polygon_type)) +
  geom_line(linewidth = 1) +              # Add lines, size = 1 for line thickness
  geom_point() +
  facet_grid(Year~.) +
  scale_x_continuous(breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
                     labels = month.abb) +  # Breaks at monthly intervals
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1) +
  # scale_x_date(date_labels = "%b",    # Display month abbreviation (e.g., Jan, Feb)
  #              date_breaks = "1 month") +  # Break the axis at each month
  labs(x = "DOY", y = "Height [m]") +
  theme_minimal() 



#barplot with one bar for ground surface and one for water level grouped by polygon type
#avg'd all together, not grouped by year first (2012 exceptionally dry)
inundation_data_avg <- inundation_data_daily %>%
  group_by(Polygon_type) %>%
  summarise(Water_Level = mean(Water_Level, na.rm = TRUE), Ground_Surface = mean(Ground_Surface, na.rm = TRUE),
            Inundation_level = mean(Inundation_level, na.rm = TRUE), Ground_Surface_adj = mean(Ground_Surface_adj, na.rm = TRUE),
            Water_Level_adj = mean(Water_Level_adj, na.rm = TRUE), Inundation_level_adj = mean(Inundation_level_adj, na.rm = TRUE))



#convert to long format
names(inundation_data_avg)
inundation_data_avg_long <- inundation_data_avg %>%
  pivot_longer(
    cols = Water_Level:Inundation_level_adj, 
    names_to = "variable",
    values_to = "value")

#subset to just soil height and water height for right now
var_list <- c("Ground_Surface", "Water_Level")
inundation_data_avg_long_sub <- subset(inundation_data_avg_long, variable %in% var_list)

ggplot(inundation_data_avg_long_sub, aes(x = factor(Polygon_type), y = value)) +
  geom_bar(stat = "identity", aes(fill = variable), position = "dodge") +
  #scale_fill_manual(values = pft_colors) +
  #scale_x_discrete(labels=c("Trough", "LCPcenter", "LCPtransition", "HCPcenter", "HCPtransition")) +
  #coord_cartesian(ylim = c(-0.10, max(inundation_data_avg$soil_height) + 0.10)) + # Adjust y-axis limits
  labs(x = NULL, y = "Height [m]", title = "Obs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        legend.title = element_blank(), axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 15), axis.title = element_text(size = 15))


  
