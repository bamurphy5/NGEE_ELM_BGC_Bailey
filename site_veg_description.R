
#1/7/25

#pulling in the veg data from intensively sampled plots at BEO to describe the vegetation by polygon type for the paper site description
#want to use the actual species names for the description instead of the PFTs so pulling in the data again

veg_data <- "~/NGEE_modeling/ELM_PFLOTRAN_BEO/2_Analysis/observational_data/vegetation/community_comp/data"
veg_data <- read.csv(file.path(veg_data, "V_1_2_plant_community_composition_Barrow_2012_v1.csv"), header = TRUE)

#see how data is organized
unique(veg_data$polygon_type) #"low-center"   "high-center"  "transitional"

#change polygon names
veg_data <- veg_data %>%
  mutate(polygon_type = case_when(polygon_type == "low-center" ~ "LCP",
                                  polygon_type == "high-center" ~ "HCP",
                                  polygon_type == "transitional" ~ "FCP"))

veg_data$polygon_sub_unit <- ifelse(veg_data$polygon_sub_unit == "edge", "rim", veg_data$polygon_sub_unit)

veg_data$polygon_type_combo <- paste0(veg_data$polygon_type, veg_data$polygon_sub_unit)

#metadata says values are expressed as % cover, but some of them do seem to be >100% for a transect

#convert cols of veg type into long format
veg_data <- veg_data %>%
  pivot_longer(
    cols = Carex_aquatilis:Bare_ground, 
    names_to = "veg_type",
    values_to = "value")

veg_data_summary <- veg_data %>%
  group_by(polygon_type_combo, veg_type) %>%
  summarise(veg_cover = mean(value, na.rm = T))

#need to replace anything less than 1% w/NA otherwise plot is way too busy
veg_data_summary$veg_cover <- ifelse(veg_data_summary$veg_cover > 1, veg_data_summary$veg_cover, NA)

polygon_types <- unique(veg_data_summary$polygon_type_combo)

plot_feature <- polygon_types[9]

veg_data_summary_sub <- subset(veg_data_summary, polygon_type_combo == plot_feature)

# Order categories by value from largest to smallest
veg_data_summary_sub <- veg_data_summary_sub %>%
  mutate(veg_type = factor(veg_type, levels = veg_type[order(-veg_cover)]))

#plot pie charts of coverage, filter out veg types w/NA so they don't get included in the legend
ggplot(veg_data_summary_sub %>% dplyr::filter(!is.na(veg_cover)), aes(x = "", y = veg_cover, fill = veg_type)) +
  #facet_grid(polygon_type_combo) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  labs(title = paste0("% Cover of Vegetation: ", plot_feature)) +
  theme_void() +  # Removes background, grid, and numeric labels
  #scale_fill_brewer(palette = "Set3")  # Optional: Change the color palette
  geom_text(aes(label = paste0(round(veg_cover, 2), "%")), # Display values as labels
            position = position_stack(vjust = 0.5), color = 'white') + # Position the labels
  scale_color_discrete(na.translate = FALSE) + # Exclude NA from the legend
  #scale_fill_manual(values = pft_colors) +
  theme(legend.title = element_blank(),   
        legend.text = element_text(size = 12)) 
