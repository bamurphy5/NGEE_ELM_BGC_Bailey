
#3/18/25

#realized I messed up the polygon type assignments so going back to the hourly output before any averaging was done and fixing the polygon names then re-doing the averaging

site_name <- "BEO"
run_type <- "Arctic"
start_yr <- "1996" 
end_yr <- "1999" 

path.out <- paste0("~/NGEE_modeling/ELM_PFLOTRAN_BEO/1_Runs/", run_type, "_PFT/BEO_7cell_EXTRACTED") 


#don't want 2000 though, moved it to 'BAD' folder, need to extract 1995 later so have full 5 years, just going w/4 years for right now
full_df_2D <- readbulk::read_bulk(directory = path.out, extension = paste0("BEO_hourly_polygon_2D_", run_type, ".*.csv"), header = TRUE)
full_df_2D$File <- NULL
#come back and do 3D, want to make sure everything works first w/2D
#var_df_3D <- readbulk::read_bulk(directory = path.out, extension = paste0("BEO_hourly_polygon_3D_", run_type, ".*.csv"), header = TRUE)

#polygon_type is included in the df's at this point but is incorrect, need to replace
full_df_2D$H2OSFC_TIDE <- as.character(full_df_2D$H2OSFC_TIDE) #convert to char first so rounding doesn't mess up classification
#values have been converted from mm to m already so units are different than before
full_df_2D <- full_df_2D %>%
  mutate(polygon_type = case_when(H2OSFC_TIDE == "0.1" ~ "LCPtrough",
                                  H2OSFC_TIDE == "0.0700899963378906" ~ "LCPcenter",
                                  H2OSFC_TIDE == "0.0910950012207031" ~ "FCPtrough",
                                  H2OSFC_TIDE == "-0.0228581008911133" ~ "FCPcenter",
                                  H2OSFC_TIDE == "-0.0189648399353027" ~ "LCPrim",
                                  H2OSFC_TIDE == "0.0821900024414063" ~ "HCPtrough",
                                  H2OSFC_TIDE == "-0.115806213378906" ~ "HCPcenter"))

full_df_2D$H2OSFC_TIDE <- as.numeric(full_df_2D$H2OSFC_TIDE) #switch back to numeric


#now calculate daily scale for 2D----------------------------------------

#instead of doing daily totals I'm just gonna do daily averages for everything, need units to be comparable to tower
#first convert 2D df into long format
full_df_2D_long <- full_df_2D %>%
  pivot_longer(
    cols = CH4FLUX_ALQUIMIA:FSAT, #MANUALLY CHANGE HERE IF VARIABLE LIST CHANGES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    names_to = "variable",
    values_to = "value")


full_df_2D_long_daily <- full_df_2D_long %>%
  group_by(polygon_type, variable, year, month, day) %>%
  summarise(value_mean = mean(value, na.rm = TRUE))

rm(full_df_2D_long)

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

#CH4FLUX_ALQUIMIA is natively in gC/m^2/s, but FCH4 (default methane model) is in kgC/m2/s
#need to go from kgC to gC
full_df_2D_long_daily$value_mean <- ifelse(full_df_2D_long_daily$variable == "FCH4", full_df_2D_long_daily$value_mean*1000, full_df_2D_long_daily$value_mean)
full_df_2D_long_daily$value_mean <- ifelse(full_df_2D_long_daily$variable == "TSA", full_df_2D_long_daily$value_mean - 273.15, 
                                      full_df_2D_long_daily$value_mean) #convert temp from K to C, came back to add later, forgot


#units note: at this point H2OSFC and H2OSFC_TIDE have been converted from mm to m, TSOI and TSA have been converted from K to C, and 
#FCH4 has been converted to gC/m2/s
#but everything else in the 2D dataset is in its native units, just averaged to daily
#save 2D variable, daily avg values
write.csv(full_df_2D_long_daily, file.path(path.out, paste0(site_name, "_daily_polygon_2D_", run_type, "_", start_yr, "_", end_yr,"_FIXED.csv")), row.names = FALSE)

#want to convert "GPP", "HR", "NEE", "NPP" from gC/m^2/s to umol CO2 m-2 s-1 (daily avg) 
#and "CH4FLUX_ALQUIMIA", "FCH4" from gC/m^2/s to nmol CH4 m-2 s-1
#so things are comparable to flux towers

# Conversion function
convert_to_umol <- function(gC_m2_s) {
  gC_m2_s * (1e6 / 12.01)
}

umol_vars <- c("GPP", "HR", "NEE", "NPP")

full_df_2D_long_daily$value_mean <- ifelse(full_df_2D_long_daily$variable %in% umol_vars, convert_to_umol(full_df_2D_long_daily$value_mean), 
                                           full_df_2D_long_daily$value_mean)


# Conversion function for CH4
convert_to_nmol <- function(gC_m2_s) {
  gC_m2_s * (1e9 / 12.01)
}

nmol_vars <- c("CH4FLUX_ALQUIMIA", "FCH4")

full_df_2D_long_daily$value_mean <- ifelse(full_df_2D_long_daily$variable %in% nmol_vars, convert_to_nmol(full_df_2D_long_daily$value_mean), 
                                           full_df_2D_long_daily$value_mean)



# #want to look at the CH4 flux data real quick to see if outliers are going to be an issue
# test <- subset(full_df_2D_long_daily, variable == "CH4FLUX_ALQUIMIA") #"CH4FLUX_ALQUIMIA"
# 
# test <- test %>%
#   group_by(variable) %>%
#   mutate(mean_value = mean(value_mean, na.rm = TRUE),
#          sd_value = sd(value_mean, na.rm = TRUE),
#          is_outlier = abs(value_mean - mean_value) > 3 * sd_value,
#          value_mean_OR = ifelse(is_outlier, NA, value_mean)) #%>%
#   #select(-mean_value, -sd_value, -is_outlier) # Remove temporary columns used in calculation

#k I'm gonna remove +- 3 SD, since that would align w/the empirical rule (68-95-99.7 rule), which states:
#About 68% of data falls within 1 SD of the mean.
#About 95% falls within 2 SDs.
#About 99.7% falls within 3 SDs.
#also testing 4 vs 3 sd only 2 more datapoints were removed from the alquimia ch4 flux, so I think its pretty reasonable


#filter outliers in everything (only actually should do it for some variables, like shouldn't outlier screen for H2OSFC_TIDE, but adjust after)
full_df_2D_long_daily <- full_df_2D_long_daily %>%
  group_by(variable) %>%
  mutate(mean_value = mean(value_mean, na.rm = TRUE),
         sd_value = sd(value_mean, na.rm = TRUE),
         is_outlier = abs(value_mean - mean_value) > 3 * sd_value,
         value_mean_OR = ifelse(is_outlier, NA, value_mean)) %>%
  select(-mean_value, -sd_value, -is_outlier) # Remove temporary columns used in calculation

#variables that outlier screened version shouldn't be used for: "H2OSFC", "H2OSFC_TIDE", "chem_dt", "RAIN"...I think thats it?


#3D data processing----------------------------------------------------------------------------------------------------------

#k 3D files are so big that I think I need to pull them in one at a time and average instead of doing it all together
#all of the unit conversions and calculating aqeuous forms of things has already been done
#just need to add cols for month, day, hour and avg
var_df_3D_single <- read.csv(file.path(path.out, "BEO_hourly_polygon_3D_Arctic.1996.csv"), header = T)
var_df_3D_single <- read.csv(file.path(path.out, "BEO_hourly_polygon_3D_Arctic.1997.csv"), header = T)
var_df_3D_single <- read.csv(file.path(path.out, "BEO_hourly_polygon_3D_Arctic.1998.csv"), header = T)
var_df_3D_single <- read.csv(file.path(path.out, "BEO_hourly_polygon_3D_Arctic.1999.csv"), header = T)

var_df_3D_single$year <- year(var_df_3D_single$date)
var_df_3D_single$month <- month(var_df_3D_single$date)
var_df_3D_single$day <- day(var_df_3D_single$date)


var_df_3D_single <- var_df_3D_single %>%
  group_by(polygon_type, variable, soil_depth, year, month, day) %>%
  summarise(value_mean = mean(value, na.rm = TRUE))

#rm(full_df_2D_long)

#using Sigrid's season approach
snowmelt <- c(4, 5)
growing <- c(6, 7, 8, 9)
freeze <- c(10, 11)
winter <- c(1, 2, 3, 12)

var_df_3D_single <- var_df_3D_single %>%
  mutate(season = case_when(month %in% snowmelt ~ "Snowmelt",
                            month %in% growing ~ "Growing Season",
                            month %in% freeze ~ "Freeze Up",
                            month %in% winter ~ "Winter Dormant"))

#filter outliers in everything 
var_df_3D_single <- var_df_3D_single %>%
  group_by(variable) %>%
  mutate(mean_value = mean(value_mean, na.rm = TRUE),
         sd_value = sd(value_mean, na.rm = TRUE),
         is_outlier = abs(value_mean - mean_value) > 3 * sd_value,
         value_mean_OR = ifelse(is_outlier, NA, value_mean)) %>%
  select(-mean_value, -sd_value, -is_outlier) # Remove temporary columns used in calculation


full_df_3D_long_daily <- var_df_3D_single #FIRST TIME!!!

full_df_3D_long_daily <- rbind(full_df_3D_long_daily, var_df_3D_single) #SUBSEQUENT TIMES!!!!!!!!!!!!!
dim(full_df_3D_long_daily)

#save to csv
write.csv(full_df_3D_long_daily, file.path(path.out, paste0(site_name, "_daily_polygon_3D_", run_type, "_", start_yr, "_", end_yr,"_FIXED.csv")), row.names = FALSE)

rm(var_df_3D_single)
rm(var_df_3D_full)

#ughghgg need to add back in a date col for plotting
full_df_3D_long_daily$date <- as.Date(with(full_df_3D_long_daily, paste(year, month, day, sep = "-")))

#----------------
#PLOTS 
#----------------

#ridgeplots of 2D variables by season---------------------------------------
season_list <- c("Snowmelt", "Growing Season", "Freeze Up", "Winter Dormant")
length(season_list)
var_list <- unique(full_df_2D_long_daily$variable)
#drop "chem_dt"

for(var_plot in unique(var_list)) {
  plot_list <- list()
  title_plot <- grid::textGrob(var_plot, gp=grid::gpar(fontsize=14, fontface='bold'))
  
  for(season_plot in unique(season_list)) {
    
    full_df_2D_long_daily_sub <- subset(full_df_2D_long_daily, variable == var_plot & season == season_plot)
    
    # Calculate median for each group
    medians <- full_df_2D_long_daily_sub %>%
      group_by(variable, polygon_type, season) %>%
      summarise(median = median(value_mean_OR, na.rm = TRUE))
    
    p <- ggplot(data = full_df_2D_long_daily_sub, aes(x = value_mean_OR, y = season, fill = polygon_type)) +
      geom_density_ridges(alpha = .7, rel_min_height = .01, color = "white") +
      #facet_wrap(~ season, ncol = 2, nrow = 2) + #scales = "free",
      scale_fill_manual(values = c("FCPcenter" = "#abba8f", "LCPcenter" = "#2a7c6e", "LCPtrough" = "#2aaeae", "HCPcenter" = "#781e1f", 
                                   "HCPtrough" = "#c85c35", "LCPrim" = "#b4f0f8", "FCPtrough" = "#cdd0bb")) +
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
variable_list <- c("ALT", "FSAT")

subset_data <- subset(full_df_2D_long_daily, variable %in% variable_list)

p <- ggplot(data = subset_data) + 
  facet_grid(variable ~ ., scales = "free_y") + #, scales = "free_y"allows y axis to vary, but same xaxis for all
  geom_smooth(aes(x = month, y= value_mean_OR, color = polygon_type), linewidth = 0.8) + #linetype = "dotdash") + #add linear fit line
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
#test <- subset(subset_data, year == 1998)

p <- ggplot(data = subset_data) + 
  facet_grid(variable ~ ., scales = "free_y") + #, scales = "free_y"allows y axis to vary, but same xaxis for all
  geom_line(aes(x = date, y= value_mean_OR, color = polygon_type), linewidth = 0.8) + #linetype = "dotdash") + #add linear fit line
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
  scale_x_date(date_labels = "%m", date_breaks = "1 month") + # "%b" Display month names, "%m" is month numbers, expand = c(0, 0) removes extra space on the x-axis
  labs(x = "Month", y = NULL) + 
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

#3D variables:
# [1] "CH4_vr"          "CH4_vr_aq"       "DIC_vr"          "DIC_vr_aq"       "DOC_vr"          "DOC_vr_aq"       "H2OSOI"         
# [8] "LITR1C_vr"       "LITR2C_vr"       "LITR3C_vr"       "SIC_vr"          "SOIL1C_vr"       "SOIL2C_vr"       "SOIL3C_vr"      
# [15] "SOIL4C_vr"       "SOILICE"         "SOILLIQ"         "TSOI"            "VWC"             "frozen_frac"     "soilC_concen"   
# [22] "soil_Fe2"        "soil_Fe2_aq"     "soil_FeOxide"    "soil_FeOxide_aq" "soil_FeS"        "soil_FeS_aq"     "soil_O2"        
# [29] "soil_O2_aq"      "soil_acetate"    "soil_acetate_aq" "soil_pH"         "soil_salinity"   "soil_sulfate"    "soil_sulfate_aq"
# [36] "soil_sulfide"    "soil_sulfide_aq" "watsat"         

subset_data <- subset(full_df_3D_long_daily, variable == "TSOI") #& polygon_type == "HCPcenter")

#fill_name <- expression(atop("Soil Temp.", paste(~gC ~m^{-3})))
fill_name <- "Temp. (C)"

#subset_data$day <- as.numeric(subset_data$day)

subset_data$soil_depth <- factor(subset_data$soil_depth, levels = c("2.86", "1.73", "1.04", "0.62", "0.37", "0.21", "0.12", "0.06", 
                                                                    "0.03", "0.01"))

#need to convert soil depth to factor so it extends the width of the line and fills out the plot
g1 <- ggplot(data = subset_data, aes(x = day, y = soil_depth, z = value_mean_OR))  +
  #scale_y_continuous(breaks = seq(0, 3, 0.5)) + #breaks = seq(0, 3, 0.5)
  #scale_y_reverse() +
  facet_grid(polygon_type~month) +
  #scale_y_continuous(trans = "reverse", breaks = seq(0, 3, 1)) +
  scale_y_discrete(breaks = levels(subset_data$soil_depth)[seq(1, length(levels(subset_data$soil_depth)), by = 2)]) +
  scale_x_continuous(breaks =c(1,10,20,31)) +
  labs(x = "Day", y = "Soil Depth (m)",
       fill = fill_name, title = "Soil Temperature")

p1 <- g1 + geom_tile(aes(fill = value_mean_OR)) + #color= "white",size=0.1, alpha = .7, 
  #scale_fill_viridis_c(option = "mako") +
  scale_fill_viridis(option ="C") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12), #axis.title.x = element_blank(), axis.title.y = element_blank(), 
        legend.text = element_text(size = 13), legend.title = element_text(size = 13),
        title = element_text(size = 13))

p1


#simpler version of the above plot, not broken up by month
#tile plots of variables by depth and time, faceted by polygon type------------------------------------------

# #need a dummy year to set up date object, but avg'd across 1999-2000, just say 2000
# full_df_3D_long_daily$year <- 2000
# full_df_3D_long_daily$month <- sapply(full_df_3D_long_daily$month, add_leading_zero)
# full_df_3D_long_daily$day <- sapply(full_df_3D_long_daily$day, add_leading_zero)
# full_df_3D_long_daily$date <- paste0(full_df_3D_long_daily$year, full_df_3D_long_daily$month, full_df_3D_long_daily$day)
# full_df_3D_long_daily$date <- as.Date(full_df_3D_long_daily$date, format='%Y%m%d')

variable_list <- c("CH4_vr_aq", "DIC_vr_aq", "DOC_vr_aq", "VWC", "TSOI", "frozen_frac", "soilC_concen",
                   "soil_Fe2_aq", "soil_FeOxide_aq", "soil_FeS_aq", "soil_O2_aq", 
                   "soil_pH", "soil_sulfate_aq", "soil_sulfide_aq", "soil_acetate")

variable_plot <- variable_list[15]
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
plot_title <- variable_names[15]

#something weird going on w/CH4 in deepest level...cutting off to see if 0-1m look normal or if I made units error
#subset_data <- subset(subset_data, soil_depth <= 1.04) #1.04

#full depth profile
subset_data$soil_depth <- factor(subset_data$soil_depth, levels = c("2.86", "1.73", "1.04", "0.62", "0.37", "0.21", "0.12", "0.06", 
                                                                    "0.03", "0.01"))

# #this one is w/o deepest 2 soil layers
# subset_data$soil_depth <- factor(subset_data$soil_depth, levels = c("0.62", "0.37", "0.21", "0.12", "0.06", 
#                                                                     "0.03", "0.01"))
# 
# #top meter-ish
# subset_data$soil_depth <- factor(subset_data$soil_depth, levels = c("1.04", "0.62", "0.37", "0.21", "0.12", "0.06", 
#                                                                     "0.03", "0.01"))

#need to convert soil depth to factor so it extends the width of the line and fills out the plot
g1 <- ggplot(data = subset_data, aes(x = date, y = soil_depth, z = value_mean))  +
  #scale_y_continuous(breaks = seq(0, 3, 0.5)) + #breaks = seq(0, 3, 0.5)
  #scale_y_reverse() +
  facet_grid(polygon_type~.) +
  #scale_y_continuous(trans = "reverse", breaks = seq(0, 3, 1)) +
  scale_y_discrete(breaks = levels(subset_data$soil_depth)[seq(1, length(levels(subset_data$soil_depth)), by = 2)]) +
  scale_x_date(date_labels = "%Y-%b") +  # Format to display only month
  labs(x = "Time", y = "Soil Depth (m)",
       fill = fill_name, title = plot_title)

p1 <- g1 + geom_tile(aes(fill = value_mean)) + #color= "white",size=0.1, alpha = .7, 
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
                   "soil_pH", "soil_sulfate_aq", "soil_sulfide_aq", "soil_acetate_aq") #"SIC_vr"

variable_names <- c("Soil Dissolved CH4 (mmol L-1)", "Soil DIC (mmol L-1)", "Soil DOC (mmol L-1)", "VWC",
                    "T Soil (C)", "Frozen Fraction", "Soil C Concentration (gC m-3)", "Soil Fe2 (mmol L-1)",
                    "Soil Fe Oxide (mmol L-1)", "Soil FeS (mmol L-1)", "Soil O2 (mmol L-1)",
                    "Soil pH", "Soil Sulfate (mmol L-1)", "Soil Sulfide (mmol L-1)", "Soil Acetate (mmol L-1)"
                    ) #"Soil Inorganic C (gC m-3)"

#for(i in unique(variable_list)) {

variable_plot <- variable_list[11]
print(variable_plot)

plot_title <- variable_names[11] #"Soil C Concentration (gC m-3)" #"Soil Inorganic C (gC m-3)"

subset_data <- subset(full_df_3D_long_daily, variable == variable_plot)

subset_data <- subset_data %>%
  group_by(variable, polygon_type, soil_depth, season) %>%
  summarise(value_avg = mean(value_mean, na.rm = T))


subset_data <- subset(subset_data, soil_depth <= 1.04) #1.04 #0.5


subset_data$season <- factor(subset_data$season, levels = c("Winter Dormant", "Snowmelt", "Growing Season", "Freeze Up"))



p <- ggplot(subset_data, aes(x = soil_depth, y = value_avg, color = polygon_type)) +
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
p

#}





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

variables_sub <- c("NEE", "HR", "GPP") #both are in gC/m^2/s
variables_sub <- c("FCH4", "CH4FLUX_ALQUIMIA", "ALT")

full_df_2D_long_daily_sub <- subset(full_df_2D_long_daily, variable %in% variables_sub)

#try dropping the last year and replotting, since cH4 and HR are weird for 2000
#full_df_2D_long_daily_sub <- subset(full_df_2D_long_daily_sub, year != 2000)

#want to average all years together first for these figs
full_df_2D_long_daily_sub <- full_df_2D_long_daily_sub %>%
  group_by(polygon_type, variable, season, month, day) %>%
  summarize(value_daily_sd = sd(value_mean, na.rm = TRUE), value_daily = mean(value_mean, na.rm = TRUE))


#now group by polygon type and find the flux associated with each polygon type
match_summary <- full_df_2D_long_daily_sub %>%
  group_by(polygon_type, variable, season) %>%
  summarize(value_daily_avg = mean(value_daily, na.rm = TRUE), value_daily_sd = sd(value_daily, na.rm = TRUE))

match_summary$season <- factor(match_summary$season, levels = c("Winter Dormant", "Snowmelt", "Growing Season", "Freeze Up"))

#polygon_colors <- c("#967acc")
polygon_colors <- c("#abba8f", "#cdd0bb", "#781e1f", "#c85c35", "#2a7c6e", "#b4f0f8", "#2aaeae")

match_summary_CH4 <- subset(match_summary, variable == "CH4FLUX_ALQUIMIA")
match_summary_CH4 <- subset(match_summary, variable == "FCH4")

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
       y = "Default CH4 (nmol m-2 s-1)") +
  geom_hline(yintercept = 0, linetype = 2, color = "black", linewidth = 1) +
  #theme_minimal() +
  theme_bw() +
  theme(legend.title = element_blank(), axis.text = element_text(size = 12),  
        legend.text = element_text(size = 12), legend.position="bottom", strip.background = element_rect(fill = "white")) 



match_summary_ALT <- subset(match_summary, variable == "ALT")

ggplot(match_summary_ALT, aes(x= reorder(polygon_type, value_daily_avg), y = value_daily_avg, fill = factor(polygon_type))) +
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
  labs(title = "Polygon Type Avg. ALT", x = NULL,
       y = "ALT (m)") +
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
       y = "NEE (umol m-2 s-1)") +
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
       y = "HR (umol m-2 s-1)") +
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
       y = "GPP (umol m-2 s-1)") +
  geom_hline(yintercept = 0, linetype = 2, color = "black", linewidth = 1) +
  #theme_minimal() +
  theme_bw() +
  theme(legend.title = element_blank(), axis.text = element_text(size = 12),  
        legend.text = element_text(size = 12), legend.position="bottom", strip.background = element_rect(fill = "white")) 








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








