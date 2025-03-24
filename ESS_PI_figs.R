
#3/20/25

#making figures for my ESS PI meeting poster ~ELM-PFLOTRAN

library(ggplot2)
library(dplyr)
library(gridExtra)


surfdata_multicell <- nc_open('~/GitHub/NGEE_ELM_BGC_Bailey/BEO_surfdata_multicell_arcticpfts.nc', write = TRUE) 

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

#first make the donut plots for veg fractional cover---------------------------------------------

veg <- ncvar_get(surfdata_multicell, "PCT_NAT_PFT")
veg <- data.frame(veg)

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

#add actual PFT names
veg <- veg %>%
  mutate(PFT_name = case_when(PFT == "PFT1" ~ "Not Vegetated",
                              PFT == "PFT2" ~ "Arctic Lichen",
                              PFT == "PFT3" ~ "Arctic Bryophyte",
                              PFT == "PFT4" ~ "Arctic Evergreen Shrub (dwarf)",
                              PFT == "PFT5" ~ "Arctic Evergreen Shrub (tall)",
                              PFT == "PFT6" ~ "Arctic Deciduous Shrub (dwarf)",
                              PFT == "PFT7" ~ "Arctic Deciduous Shrub (low)",
                              PFT == "PFT8" ~ "Arctic Deciduous Shrub (tall)",
                              PFT == "PFT9" ~ "Arctic Deciduous Shrub (alder)",
                              PFT == "PFT10" ~ "Arctic Forb",
                              PFT == "PFT11" ~ "Arctic Dry Graminoid",
                              PFT == "PFT12" ~ "Arctic Wet Graminoid",))

#fill 0's w/NA so all the PFTs that aren't represented don't get included in the legend
veg <- veg %>%
  mutate(value = na_if(value, 0))


# # Function to create a donut plot
# make_donut_plot <- function(df) {
#   df <- df %>% 
#     arrange(desc(PFT_name)) %>% 
#     mutate(ypos = cumsum(value) - 0.5 * value)
#   
#   ggplot(df, aes(x = 2, y = value, fill = PFT_name)) +
#     geom_bar(stat = "identity", width = 1) +
#     coord_polar("y", start = 0) +
#     geom_text(aes(x = 0, y = 0, label = unique(polygon_type)),
#               size = 8, fontface = "bold") +
#     scale_fill_brewer(palette = "Set3") +
#     theme_void() +
#     xlim(0.5, 2.5)
# }
# 
# 
# # Create donut plots for each landform
# plots <- veg %>% 
#   split(.$polygon_type) %>% 
#   lapply(make_donut_plot)
# 
# # Arrange them together
# grid.arrange(grobs = plots, ncol = 2)


#other approach
library(ggforce)
#k for this version I need to actually remove veg types that don't have coverage so they don't plot

data <- veg %>% 
  filter(value > 0) %>%
  arrange(desc(PFT_name)) %>% 
  mutate(ypos = cumsum(value) - 0.5 * value)

pft_colors <- c("Not Vegetated" = "#e7e6eb", "Arctic Lichen" = "#d8a4a8", "Arctic Bryophyte" = "#8299a5", 
                "Arctic Evergreen Shrub (dwarf)" = "#a6adc0",
                "Arctic Deciduous Shrub (low)" = "#815167", "Arctic Forb" = "#6a7d69", 
                "Arctic Dry Graminoid" = "#382a36", "Arctic Wet Graminoid" = "#2b4e61")

# Donut plot function
donut_plot <- ggplot(data, aes(x = 2, y = value, fill = PFT_name)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  facet_wrap(~polygon_type, ncol = 7) +
  geom_text(aes(x = 0, y = 0, label = polygon_type),
            size = 4, fontface = "bold") +
  theme_void() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = pft_colors) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 15),
        strip.text = element_blank()) 

# Print plot
print(donut_plot)



#kk now plot the OM density---------------------------------------------------

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

#add in the soil depths associated w/the layers
RUNID <- "Alaska_alquimia_arctic_BAM_2_AK-BEO_ICB20TRCNRDCTCBC.elm.h0.1999-01-01-00000.nc"
dat.base <- "~/NGEE_modeling/ELM_PFLOTRAN_BEO/1_Runs/Arctic_PFT/BEO_7cell_RAW_h0"

model_data <- ncdf4::nc_open(file.path(dat.base, RUNID))
soil_layer_depth <- model_data[["dim"]][["levdcmp"]]$vals #depths for each layer are in m, same values for "levdcmp" and "levgrnd"
#subset to top 10 and round them so they're cleaner
soil_layer_depth <- round(soil_layer_depth[1:10], 2)
ncdf4::nc_close(model_data)

#add to OM df
OM <- OM %>%
  mutate(soil_depth = rep(soil_layer_depth, times = 7))


fill_name <- "OM Density \n(kg m-3)"

rng <- range(OM$value)

#need to convert soil depth to factor so it extends the width of the line and fills out the plot
g1 <- ggplot(data = OM, aes(x = polygon_type, y = factor(soil_depth), fill = value))+ 
  #scale_y_continuous(trans = "reverse") + #breaks = seq(1, 10, 1)) +
  scale_y_discrete(limits = rev(unique(factor(OM$soil_depth)))) + #for characters
  labs(x = NULL, y = "Soil Depth (m)", fill = fill_name)

p1 <- g1 + geom_tile(aes(fill = value), color= "white") + #color= "white",size=0.1, alpha = .7, 
  #scale_fill_viridis(option ="C") + #limits=c(rng[1], rng[2])) +
  scale_fill_gradient2(low="#2b4e61", mid="#a6adc0", high="#d8a4a8", 
                       midpoint=mean(rng),
                       limits=c(rng[1], rng[2])) +
  theme_minimal() +
  theme(axis.text = element_text(size = 15), #axis.title.x = element_blank(), axis.title.y = element_blank(), 
        legend.text = element_text(size = 15), legend.title = element_text(size = 15),
        title = element_text(size = 15), legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1))


p1


#k now make the inundated bar plot figure----------------------------------------------------
water <- ncvar_get(surfdata_multicell, "ht_above_stream")

water <- data.frame(water, landcover_types)

names(water) <- c("value", "polygon_type")

water_plot <- ggplot(water, aes(x = factor(polygon_type), ymin = -0.1, ymax = value)) +
  geom_linerange(linewidth = 13, color = "#815167") +
  coord_cartesian(ylim = c(-0.10, max(water$value) + 0.10)) + # Adjust y-axis limits
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "#2b4e61", linewidth = 3) +  # Add horizontal line to show tide height
  # Add shaded water effect below y = 0.1
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.1),
            fill = "#2b4e61", alpha = 0.05) +
  labs(x = NULL, y = "Height Above Water (m)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        legend.title = element_blank(), axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 15), axis.title = element_text(size = 15))

water_plot


#k now to arrange all three plots together....
library(patchwork)
#trying version with patchwork...this works and looks good!
combo_plot <- p1 / p2 / p3 +
  plot_annotation(tag_levels = 'a')

combo_plot <- water_plot / p1 / donut_plot #this stacks all of them on top of each other

combo_plot <- (water_plot / p1) | donut_plot + 
  plot_annotation(tag_levels = 'a')
  #plot_layout(widths = c(1,2))

#need to adjust the margins on the donut plot to remove some white space
donut_plot <- donut_plot +
  theme(
    panel.spacing = unit(0.1, "lines"),  # Reduce space between facets
    plot.margin = margin(2, 2, 2, 2), # Adjust outer plot margins
    aspect.ratio = 1  # Ensures circular appearance, makes donuts look bigger?
  )

combo_plot <- donut_plot / (p1 | water_plot) + 
  plot_annotation(tag_levels = 'a')

#plot_save <- "spatial_AGB_GPP_NEE_NOCONIFERAGB_order.png"
#ggsave(filename = plot_save, plot = combo_plot, dpi = 500)




#----------------------------------------------------------------------------------------------------------------

#now make the polygon scaling figs-------------------------------------------------------------------------------

#pull in the df version of the Wainwright polygon map, cropped to 400X400m domain around the BEO tower
polygon_df_crop <- read.csv("Wainwright_poly_crop.csv", header = T)

#going to plot maps for the three grouping approaches, and bar plots to go underneath each map showing the fractional cover of the different groups
#coordinates are in UTM NAD83easting/northing
#some of the parts that are shown as white (NA) are shown in Dengel paper as 'drainage', pretty much just the bottom right corner stuff

#first plots for polygon type + feature grouping

color_scheme <- c("LCPcenter" = "#d8a4a8", "HCPcenter" = "#8299a5", 
                "FCPcenter" = "#a6adc0",
                "FCPtrough" = "#815167", "HCPtrough" = "#6a7d69", 
                "LCPrim" = "#382a36", "LCPtrough" = "#2b4e61")

type_ft_map <- ggplot() +
  geom_raster(data = polygon_df_crop %>% dplyr::filter(!is.na(land_name)),
              aes(x = X, y = Y, fill = factor(land_name))) +
  #scale_fill_carto_d(name = "Landcover Class", palette = "Antique") +
  scale_fill_manual(values = color_scheme) +
  geom_point(data = utm_coords, aes(x = X, y = Y), color = "white", size = 5) + #add point for tower location
  labs(title = "Polygon Type + Feature Classification", x = "Easting (UTM NAD83), m", y = "Northing (UTM NAD83), m") +
  theme_bw() +
  theme(legend.title = element_blank(), axis.text = element_text(size = 13), 
        axis.title = element_text(size = 13), 
        legend.text = element_text(size = 13), legend.position = "bottom") 


#plot barplot again for cropped map
type_ft_bar <- ggplot(polygon_df_crop %>% dplyr::filter(!is.na(land_name)), aes(x = land_name)) + 
  #geom_bar(binwidth = 1, fill = PFT, alpha = 0.7) +
  geom_bar(aes(fill = factor(land_name)))+
  scale_fill_manual(values = color_scheme) +
  coord_flip() +  # Flip bars horizontally
  scale_y_continuous(labels = scales::scientific) + # Scientific notation
  #scale_fill_carto_d(name = "Landcover Class", palette = "Safe") +
  labs(x = NULL, y = "Pixels per class") +
  #theme_minimal() +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = "none", axis.text = element_text(size = 13), #axis.text.x = element_blank(),
        axis.title.y = element_text(size = 13), 
        legend.text = element_text(size = 13)) 





#next plots for polygon type grouping

color_scheme <- c("FCP" = "#d8a4a8","HCP" = "#6a7d69", 
                  "LCP" = "#2b4e61")

type_map <- ggplot() +
  geom_raster(data = polygon_df_crop %>% dplyr::filter(!is.na(land_name_type)),
              aes(x = X, y = Y, fill = factor(land_name_type))) +
  #scale_fill_carto_d(name = "Landcover Class", palette = "Antique") +
  scale_fill_manual(values = color_scheme) +
  geom_point(data = utm_coords, aes(x = X, y = Y), color = "white", size = 5) + #add point for tower location
  labs(title = "Polygon Type Classification", x = "Easting (UTM NAD83), m", y = "Northing (UTM NAD83), m") +
  theme_bw() +
  theme(legend.title = element_blank(), axis.text = element_text(size = 13), 
        axis.title = element_text(size = 13), 
        legend.text = element_text(size = 13), legend.position = "bottom") 


#plot barplot again for cropped map
type_bar <- ggplot(polygon_df_crop %>% dplyr::filter(!is.na(land_name_type)), aes(x = land_name_type)) + 
  #geom_bar(binwidth = 1, fill = PFT, alpha = 0.7) +
  geom_bar(aes(fill = factor(land_name_type)))+
  scale_fill_manual(values = color_scheme) +
  coord_flip() +  # Flip bars horizontally
  scale_y_continuous(labels = scales::scientific) + # Scientific notation
  #scale_fill_carto_d(name = "Landcover Class", palette = "Safe") +
  labs(x = NULL, y = "Pixels per class") +
  #theme_minimal() +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = "none", axis.text = element_text(size = 13), #axis.text.x = element_blank(),
        axis.title.y = element_text(size = 13), 
        legend.text = element_text(size = 13)) 



#now plots for wet/dry grouping

color_scheme <- c("dry" = "#382a36", 
                  "wet" = "#8299a5")

hydro_map <- ggplot() +
  geom_raster(data = polygon_df_crop %>% dplyr::filter(!is.na(land_name_hydro)),
              aes(x = X, y = Y, fill = factor(land_name_hydro))) +
  #scale_fill_carto_d(name = "Landcover Class", palette = "Antique") +
  scale_fill_manual(values = color_scheme) +
  geom_point(data = utm_coords, aes(x = X, y = Y), color = "white", size = 5) + #add point for tower location
  labs(title = "Polygon Moisture Classification", x = "Easting (UTM NAD83), m", y = "Northing (UTM NAD83), m") +
  theme_bw() +
  theme(legend.title = element_blank(), axis.text = element_text(size = 13), 
        axis.title = element_text(size = 13), 
        legend.text = element_text(size = 13), legend.position = "bottom") 


#plot barplot again for cropped map
hydro_bar <- ggplot(polygon_df_crop %>% dplyr::filter(!is.na(land_name_hydro)), aes(x = land_name_hydro)) + 
  #geom_bar(binwidth = 1, fill = PFT, alpha = 0.7) +
  geom_bar(aes(fill = factor(land_name_hydro)))+
  scale_fill_manual(values = color_scheme) +
  coord_flip() +  # Flip bars horizontally
  scale_y_continuous(labels = scales::scientific) + # Scientific notation
  #scale_fill_carto_d(name = "Landcover Class", palette = "Safe") +
  labs(x = NULL, y = "Pixels per class") +
  #theme_minimal() +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = "none", axis.text = element_text(size = 13), #axis.text.x = element_blank(),
        axis.title.y = element_text(size = 13), 
        legend.text = element_text(size = 13)) 




#combine them all into one figure
combo_plot <- (type_ft_map | type_map | hydro_map) /  (type_ft_bar | type_bar | hydro_bar) + 
  plot_annotation(tag_levels = 'a') +
  plot_layout(heights = c(2, 1)) + # Top row twice as tall as bottom row
  plot_layout(widths = c(2, 1))

combo_plot

#save high res version
plot_save <- "polygon_classification_maps.png"
ggsave(filename = plot_save, plot = combo_plot, dpi = 400)






