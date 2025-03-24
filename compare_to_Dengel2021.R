
#1/2/25

#comparing BEO PFLOTRAN simulations to obs data in Dengel et al., 2021. The actual BEO tower data doesn't seem to currently be available on Ameriflux, so just
#looking at Sigrid's avg's in paper

library(ggridges)


site_name <- "BEO" 
run_type <- "Default"

model.data <- paste0("~/NGEE_modeling/ELM_PFLOTRAN_BEO/1_Runs/", run_type, "_PFT/BEO_7cell_EXTRACTED/updated_surface") 

#want the hourly data w/avg per second fluxes over the hour, not the daily totals version
#full_df_2D_daily <- read.csv(file.path(model.data, "BEO_daily_polygon_2D_Default_1998_2000.csv"), header = TRUE)
full_df_2D_hourly <- readbulk::read_bulk(directory = model.data, extension = paste0("BEO_hourly_polygon_2D_", run_type, ".*.csv"), header = TRUE)
full_df_2D_hourly$File <- NULL
#trying w/o 2000 since CH4 is messed up for that year (flat 0 across the board)
full_df_2D_hourly <- subset(full_df_2D_hourly, year != 2000)

#update the season tag to match Sigrid's definitions and obs period:
#April-May = snowmelt period
#June-Sept = growing season
#Oct-Nov = freeze up period
#April-Nov = full obs period 
dim(full_df_2D_hourly)
keep_months <- c(4, 5, 6, 7, 8, 9, 10, 11)
full_df_2D_hourly <- subset(full_df_2D_hourly, month %in% keep_months)
#now data is just 1998-2000, April-Nov
snowmelt <- c(4, 5)
growing <- c(6, 7, 8, 9)
freeze <- c(10, 11)

#convert CO2 fluxes from gC m-2 s-1 to umol CO2 m-2 s-1, and convert CH4 fluxes to nmol CH4 m-2 s-1
#model outputs CH4 and CO2 as gC/m^2/s
#conversion_factor <- (1e6 / 12.01) * (44.01 / 12.01)
conversion_factor <- (1e6 / 12.01)
full_df_2D_hourly$NEE <- conversion_factor * full_df_2D_hourly$NEE
#remove implausible values >50
full_df_2D_hourly$NEE <- ifelse(full_df_2D_hourly$NEE > 50, NA, full_df_2D_hourly$NEE)

# Conversion factor: from gC m-2 s-1 to nmol CH4 m-2 s-1
#conversion_factor <- (1e9 / 12.01) * (16.04 / 12.01)
conversion_factor <- (1e9 / 12.01) #dividing by 12.01 converts from gC to moles C in 1 mole CH4, multiplying by 1e9 converts from moles to nmol
full_df_2D_hourly$CH4FLUX_ALQUIMIA <- conversion_factor * full_df_2D_hourly$CH4FLUX_ALQUIMIA
#remove implausible values >100
full_df_2D_hourly$CH4FLUX_ALQUIMIA <- ifelse(full_df_2D_hourly$CH4FLUX_ALQUIMIA > 100, NA, full_df_2D_hourly$CH4FLUX_ALQUIMIA)

#realized this is probably what's messing up the NEE stuff and making it all look positive...
#this IGR method is way too aggressive, Burba intro to EC book suggests removing values >6 times the SD as implausible, going w/that
#changing to 5SD b/c values seem super high still
#screen for outliers and replace w/NA
# #using inner quartile range method b/c data not normally distributed
# # Function to replace outliers with NA using IQR
# replace_outliers_iqr <- function(data) {
#   Q1 <- quantile(data, 0.25, na.rm = TRUE) # 1st quartile
#   Q3 <- quantile(data, 0.75, na.rm = TRUE) # 3rd quartile
#   IQR_value <- Q3 - Q1                     # Interquartile Range
#   lower_bound <- Q1 - 1.5 * IQR_value      # Lower bound
#   upper_bound <- Q3 + 1.5 * IQR_value      # Upper bound
#   
#   # Replace outliers with NA
#   data[data < lower_bound | data > upper_bound] <- NA
#   return(data)
# }


# #only worrying about replacing outliers in CH4 and CO2 fluxes for now
# #full_df_2D_hourly$CH4FLUX_ALQUIMIA <- replace_outliers_iqr(full_df_2D_hourly$CH4FLUX_ALQUIMIA)
# #full_df_2D_hourly$NEE <- replace_outliers_iqr(full_df_2D_hourly$NEE)
# NEE_sd <- sd(full_df_2D_hourly$NEE, na.rm = T)
# NEE_outlier <- 5 * NEE_sd
# CH4_sd <- sd(full_df_2D_hourly$CH4FLUX_ALQUIMIA, na.rm = T)
# CH4_outlier <- 5 * CH4_sd
# 
# full_df_2D_hourly$NEE <- ifelse(full_df_2D_hourly$NEE > NEE_outlier, NA, full_df_2D_hourly$NEE)
# full_df_2D_hourly$NEE <- ifelse(full_df_2D_hourly$NEE < -NEE_outlier, NA, full_df_2D_hourly$NEE)
# 
# full_df_2D_hourly$CH4FLUX_ALQUIMIA <- ifelse(full_df_2D_hourly$CH4FLUX_ALQUIMIA > CH4_outlier, NA, full_df_2D_hourly$CH4FLUX_ALQUIMIA)
# full_df_2D_hourly$CH4FLUX_ALQUIMIA <- ifelse(full_df_2D_hourly$CH4FLUX_ALQUIMIA < -CH4_outlier, NA, full_df_2D_hourly$CH4FLUX_ALQUIMIA)

#jk might as well do everything, then the distribution plots will be more useful
# Define the function to replace outliers with NA
replace_outliers <- function(data, cols) {
  data[cols] <- lapply(data[cols], function(column) {
    if (is.numeric(column)) {
      sd_col <- sd(column, na.rm = TRUE)
      column[column > 5 * sd_col | column < -5 * sd_col] <- NA
    }
    return(column)
  })
  return(data)
}

# Specify columns to cycle through (e.g., by name or index)
columns_to_adjust <- c("CH4FLUX_ALQUIMIA", "DIC_RUNOFF", "DOC_RUNOFF", "GPP", "H2OSFC", "HR", "LEAFC", "NEE", "NPP","QDRAI", 
                       "QFLX_EVAP_TOT", "QVEGT", "SMINN", "SMINN_TO_PLANT", "SMIN_NO3_RUNOFF", "TOTLITC", "TOTSOMC", "TOTVEGC", "ZWT")

# Apply the function to the specified columns
full_df_2D_hourly <- replace_outliers(full_df_2D_hourly, columns_to_adjust)


#convert to long format
full_df_2D_hourly <- full_df_2D_hourly %>%
  pivot_longer(
    cols = CH4FLUX_ALQUIMIA:chem_dt, 
    names_to = "variable",
    values_to = "value")

full_df_2D_hourly <- full_df_2D_hourly %>%
  mutate(season_sigrid = case_when(month %in% snowmelt ~ "Snowmelt",
                                   month %in% growing ~ "Growing Season",
                                   month %in% freeze ~ "Freeze Up"))


#now group by day and avg to get avg daily fluxes in umol m-2 s-1 and nmol m-2 s-1
full_df_2D_daily <- full_df_2D_hourly %>%
  group_by(variable, polygon_type, year, month, day) %>%
  summarise(value_avg = mean(value, na.rm = T))

full_df_2D_daily <- full_df_2D_daily %>%
  mutate(season_sigrid = case_when(month %in% snowmelt ~ "Snowmelt",
                                   month %in% growing ~ "Growing Season",
                                   month %in% freeze ~ "Freeze Up"))



#calc annual avg's NOT broken down by polygon type
full_df_2D_annual_all <- full_df_2D_daily %>%
  group_by(variable, year) %>%
  summarise(value_mean = mean(value_avg, na.rm = TRUE), value_sd = sd(value_avg, na.rm = TRUE), n = length(value_avg),
            value_se = value_sd/sqrt(n))

#calc annual avg's broken down by polygon type
full_df_2D_annual_poly <- full_df_2D_daily %>%
  group_by(variable, year, polygon_type) %>%
  summarise(value_mean = mean(value_avg, na.rm = TRUE), value_sd = sd(value_avg, na.rm = TRUE), n = length(value_avg),
            value_se = value_sd/sqrt(n))






#PLOTS------------------------------------------------------------------------------------------------

#density ridgeline plots of distributions for the hourly data, grouped by season

# season_list <- c("Snowmelt", "Growing Season", "Freeze Up")
# length(season_list)
# var_list <- c("CH4FLUX_ALQUIMIA", "NEE")
# 
# for(var_plot in unique(var_list)) {
#   plot_list <- list()
#   title_plot <- grid::textGrob(var_plot, gp=grid::gpar(fontsize=14, fontface='bold'))
#   
#   for(season_plot in unique(season_list)) {
#     
#     full_df_2D_daily_sub <- subset(full_df_2D_daily, variable == var_plot)  #& season_sigrid == season_plot)
#     
#     # Calculate mean for each group
#     means <- full_df_2D_daily_sub %>%
#       group_by(variable, polygon_type, season_sigrid) %>%
#       summarise(mean = mean(value_avg, na.rm = TRUE))
#     
#     p <- ggplot(data = full_df_2D_daily_sub, aes(x = value_avg, y = season_sigrid, fill = polygon_type)) +
#       geom_density_ridges(alpha = .7, rel_min_height = .01, color = "white") +
#       #facet_wrap(~ season_sigrid, ncol = 3, nrow = 1) + #scales = "free",
#       scale_fill_manual(values = c("FCPcenter" = "#abba8f", "LCPcenter" = "#2a7c6e", "LCPtrough" = "#2aaeae", "HCPcenter" = "#781e1f", 
#                                    "HCPtrough" = "#c85c35", "LCPrim" = "#b4f0f8", "FCPtrough" = "#cdd0bb")) +
#   
#       #geom_vline(data = means, aes(xintercept = mean, color = polygon_type), linetype = "dashed", linewidth = 1.5) +
#       # scale_color_manual(values = c("FCPcenter" = "#abba8f", "LCPcenter" = "#2a7c6e", "LCPtrough" = "#2aaeae", "HCPcenter" = "#781e1f", 
#       #                               "HCPtrough" = "#c85c35", "LCPrim" = "#b4f0f8", "FCPtrough" = "#cdd0bb"), name = NULL, guide = "none") +
#       theme_ridges(font_size = 13, grid = T) +
#       labs(x = NULL, y = NULL) +
#       theme(legend.position = "top", strip.background = element_blank(), strip.text = element_blank(), legend.text = element_text(size = 15),
#             axis.text.y = element_text(size = 15), axis.title = element_text(size = 15), axis.text.x = element_text(size = 13),
#             legend.title = element_blank()) 
#     
#     plot_list[[season_plot]] <- p
#     
#   }
#   
#   grid_arrange_shared_legend(plot_list[[1]], plot_list[[2]], 
#                              ncol = 1, nrow = 1, position='top', top = title_plot) #bottom = units_plot)
#   
# }
# 


#this is plotting all the data combined, NOT broken up by polygon type
#feel like these should use the hourly data, not the daily avg's from the hourly data
var_list <- c("CH4FLUX_ALQUIMIA")
compare_df_sub <- subset(full_df_2D_hourly, variable %in% var_list)

#also want to plot the distribution for everything combined, so rbind the whole dataset to the end again but w/changed season label
compare_df_sub_all <- compare_df_sub
compare_df_sub_all$season_sigrid <- "All Data"

compare_df_sub <- rbind(compare_df_sub, compare_df_sub_all)

means <- compare_df_sub %>%
  group_by(variable, season_sigrid) %>%
  summarise(mean = mean(value, na.rm = TRUE))

p1 <- ggplot(data = compare_df_sub, aes(x = value, y = season_sigrid, fill = season_sigrid)) +
  geom_density_ridges(alpha = .7, rel_min_height = .01, color = "white") +
  #facet_wrap(~ season_sigrid, ncol = 1) + #scales = "free",
  scale_fill_cyclical(labels = c(`Snowmelt` = "Snowmelt", `Growing Season` = "Growing Season", `Freeze Up` = "Freeze Up", `All Data` = "April-Nov."),
                      values = c("#763c42", "#8b9371", "#4b4456", "#7f632f"), name = NULL, guide = "legend") +
  geom_vline(data = means, aes(xintercept = mean, color = season_sigrid), linetype = "dashed", linewidth = 1.5) +
  scale_color_manual(labels = c(`Snowmelt` = "Snowmelt", `Growing Season` = "Growing Season", `Freeze Up` = "Freeze Up", `All Data` = "April-Nov."),
                     values = c("#763c42", "#8b9371", "#4b4456", "#7f632f"), name = NULL, guide = "none") +
  theme_ridges(font_size = 13, grid = T) +
  labs(x = "CH4 Flux (nmol m-2 s-1)", y = NULL) +
  theme(legend.position = "top", strip.background = element_blank(), strip.text = element_blank(), legend.text = element_text(size = 15),
        axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 12)) 

p1




var_list <- c("NEE")
compare_df_sub <- subset(full_df_2D_hourly, variable %in% var_list)

#also want to plot the distribution for everything combined, so rbind the whole dataset to the end again but w/changed season label
compare_df_sub_all <- compare_df_sub
compare_df_sub_all$season_sigrid <- "All Data"

compare_df_sub <- rbind(compare_df_sub, compare_df_sub_all)

means <- compare_df_sub %>%
  group_by(variable, season_sigrid) %>%
  summarise(mean = mean(value, na.rm = TRUE))

p2 <- ggplot(data = compare_df_sub, aes(x = value, y = season_sigrid, fill = season_sigrid)) +
  geom_density_ridges(alpha = .7, rel_min_height = .01, color = "white") +
  #facet_wrap(~ season_sigrid, ncol = 1) + #scales = "free",
  scale_fill_cyclical(labels = c(`Snowmelt` = "Snowmelt", `Growing Season` = "Growing Season", `Freeze Up` = "Freeze Up", `All Data` = "April-Nov."),
                      values = c("#763c42", "#8b9371", "#4b4456", "#7f632f"), name = NULL, guide = "legend") +
  geom_vline(data = means, aes(xintercept = mean, color = season_sigrid), linetype = "dashed", linewidth = 1.5) +
  scale_color_manual(labels = c(`Snowmelt` = "Snowmelt", `Growing Season` = "Growing Season", `Freeze Up` = "Freeze Up", `All Data` = "April-Nov."),
                     values = c("#763c42", "#8b9371", "#4b4456", "#7f632f"), name = NULL, guide = "none") +
  theme_ridges(font_size = 13, grid = T) +
  labs(x = "CO2 Flux (umol m-2 s-1)", y = NULL) +
  theme(legend.position = "top", strip.background = element_blank(), strip.text = element_blank(), legend.text = element_text(size = 15),
        axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 12)) 

p2


#now making version of above but separated out by polygon type
var_list <- c("CH4FLUX_ALQUIMIA") #"NEE"
compare_df_sub <- subset(full_df_2D_hourly, variable %in% var_list) #& season_sigrid == "Freeze Up")
#testing out dropping the few really NEE high values
compare_df_sub$value <- ifelse(compare_df_sub$value >10, NA, compare_df_sub$value)

means <- compare_df_sub %>%
  group_by(variable, polygon_type) %>%
  summarise(mean = mean(value, na.rm = TRUE))

p1 <- ggplot(data = compare_df_sub, aes(x = value, y = season_sigrid, fill = polygon_type)) +
  geom_density_ridges(alpha = .7, rel_min_height = .01, color = "white") +
  #facet_wrap(~ season_sigrid, ncol = 1) + #scales = "free",
  scale_fill_manual(values = c("FCPcenter" = "#abba8f", "LCPcenter" = "#2a7c6e", "LCPtrough" = "#2aaeae", "HCPcenter" = "#781e1f", 
                               "HCPtrough" = "#c85c35", "LCPrim" = "#b4f0f8", "FCPtrough" = "#cdd0bb"), name = NULL, guide = "legend") +
  geom_vline(data = means, aes(xintercept = mean, color = polygon_type), linetype = "dashed", linewidth = 1.5) +
  scale_color_manual(values = c("FCPcenter" = "#abba8f", "LCPcenter" = "#2a7c6e", "LCPtrough" = "#2aaeae", "HCPcenter" = "#781e1f",
                                "HCPtrough" = "#c85c35", "LCPrim" = "#b4f0f8", "FCPtrough" = "#cdd0bb"), name = NULL, guide = "none") +
  theme_ridges(font_size = 13, grid = T) +
  labs(x = "CH4 Flux (nmol m-2 s-1)", y = NULL) +
  theme(legend.position = "top", strip.background = element_blank(), strip.text = element_blank(), legend.text = element_text(size = 15),
        axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 12)) 

p1


plot_list <- list(p1, p2, p3)

plot_final <- grid_arrange_shared_legend(plot_list[[1]], plot_list[[2]], plot_list[[3]],
                                         ncol = 1, nrow = 3, position='top')





#create the rotated bar plots for CO2 and CH4-----------------------------------------------
#units are diff and outliers have been removed so will look different than other ones I plotted previously
variables_sub <- c("NEE", "CH4FLUX_ALQUIMIA") 

full_df_2D_hourly_sub <- subset(full_df_2D_hourly, variable %in% variables_sub)

#now group by polygon type and find the flux associated with each polygon type
match_summary <- full_df_2D_hourly_sub %>%
  group_by(polygon_type, variable, season_sigrid) %>%
  summarize(value_avg = mean(value, na.rm = TRUE), value_sd = sd(value, na.rm = TRUE), n = length(value), value_se = value_sd/sqrt(n))

match_summary$season_sigrid <- factor(match_summary$season_sigrid, levels = c("Snowmelt", "Growing Season", "Freeze Up"))

#polygon_colors <- c("#967acc")
polygon_colors <- c("#abba8f", "#cdd0bb", "#781e1f", "#c85c35", "#2a7c6e", "#b4f0f8", "#2aaeae")

match_summary_CH4 <- subset(match_summary, variable == "CH4FLUX_ALQUIMIA")

ggplot(match_summary_CH4, aes(x= reorder(polygon_type, value_avg), y = value_avg, fill = factor(polygon_type))) +
  #geom_bar(binwidth = 1, fill = PFT, alpha = 0.7) +
  facet_grid(.~season_sigrid) +
  geom_bar(stat="identity")+
  #geom_bar(aes(fill = factor(PFT_name))+
  #scale_fill_manual(values = pft_colors, labels = names) +
  #scale_fill_manual(values = pft_colors, labels = c("CH4 (gC m-2 s-1)")) +
  scale_fill_manual(values = polygon_colors) +
  geom_errorbar(aes(ymin = value_avg-value_se, ymax = value_avg+value_se), width=.2,
                position=position_dodge(.9)) +
  coord_flip() +
  labs(title = "Polygon Type Avg. CH4 Flux", x = NULL,
       y = "CH4 (nmol m-2 s-1)") +
  geom_hline(yintercept = 0, linetype = 2, color = "black", linewidth = 1) +
  #theme_minimal() +
  theme_bw() +
  theme(legend.title = element_blank(), axis.text = element_text(size = 12),  
        legend.text = element_text(size = 12), legend.position="bottom", strip.background = element_rect(fill = "white")) 

#LCP mean CH4 flux 2013-2019 from Sigrids paper = 19.22
#FCP CH4 = 19.76

#this is for a version NOT broken up by season-------------------
match_summary <- full_df_2D_hourly_sub %>%
  group_by(polygon_type, variable) %>%
  summarize(value_avg = mean(value, na.rm = TRUE), value_sd = sd(value, na.rm = TRUE), n = length(value), value_se = value_sd/sqrt(n))


polygon_colors <- c("#abba8f", "#cdd0bb", "#781e1f", "#c85c35", "#2a7c6e", "#b4f0f8", "#2aaeae")

match_summary_CH4 <- subset(match_summary, variable == "CH4FLUX_ALQUIMIA")

obs_groups <- c("FCPcenter", "LCPcenter")
obs_avg <- c(19.76, 19.22)
obs_avg <- data.frame(obs_groups, obs_avg)

ggplot(match_summary_CH4, aes(x= reorder(polygon_type, value_avg), y = value_avg, fill = factor(polygon_type))) +
  #facet_grid(.~season_sigrid) +
  geom_bar(stat="identity")+
  scale_fill_manual(values = polygon_colors) +
  geom_errorbar(aes(ymin = value_avg-value_se, ymax = value_avg+value_se), width=.2,
                position=position_dodge(.9)) +
  geom_point(data = obs_avg, aes(x = obs_groups, y = obs_avg), color = "red", size = 4, inherit.aes = FALSE) + # Add points for Sigrid's avg's
  coord_flip() +
  labs(title = "Polygon Type Avg. CH4 Flux", x = NULL,
       y = "CH4 (nmol m-2 s-1)") +
  geom_hline(yintercept = 0, linetype = 2, color = "black", linewidth = 1) +
  #theme_minimal() +
  theme_bw() +
  theme(legend.title = element_blank(), axis.text = element_text(size = 12),  
        legend.text = element_text(size = 12), legend.position="bottom", strip.background = element_rect(fill = "white")) 




#NEE version of barplot broken up by season
match_summary <- full_df_2D_hourly_sub %>%
  group_by(polygon_type, variable, season_sigrid) %>%
  summarize(value_avg = mean(value, na.rm = TRUE), value_sd = sd(value, na.rm = TRUE), n = length(value), value_se = value_sd/sqrt(n))

match_summary$season_sigrid <- factor(match_summary$season_sigrid, levels = c("Snowmelt", "Growing Season", "Freeze Up"))

match_summary_NEE <- subset(match_summary, variable == "NEE")

ggplot(match_summary_NEE, aes(x= reorder(polygon_type, value_avg), y = value_avg, fill = factor(polygon_type))) +
  facet_grid(.~season_sigrid) +
  geom_bar(stat="identity")+
  scale_fill_manual(values = polygon_colors) +
  geom_errorbar(aes(ymin = value_avg-value_se, ymax = value_avg+value_se), width=.2,
                position=position_dodge(.9)) +
  coord_flip() +
  labs(title = "Polygon Type Avg. CO2 Flux", x = NULL,
       y = "CO2 (umol m-2 s-1)") +
  geom_hline(yintercept = 0, linetype = 2, color = "black", linewidth = 1) +
  #theme_minimal() +
  theme_bw() +
  theme(legend.title = element_blank(), axis.text = element_text(size = 12),  
        legend.text = element_text(size = 12), legend.position="bottom", strip.background = element_rect(fill = "white")) 





#NEE version NOT broken up by season-------------------
match_summary <- full_df_2D_hourly_sub %>%
  group_by(polygon_type, variable) %>%
  summarize(value_avg = mean(value, na.rm = TRUE), value_sd = sd(value, na.rm = TRUE), n = length(value), value_se = value_sd/sqrt(n))

match_summary_NEE <- subset(match_summary, variable == "NEE")

obs_groups <- c("FCPcenter", "LCPcenter")
obs_avg <- c(-1.04, -1.15)
obs_avg <- data.frame(obs_groups, obs_avg)

ggplot(match_summary_NEE, aes(x= reorder(polygon_type, value_avg), y = value_avg, fill = factor(polygon_type))) +
  geom_bar(stat="identity")+
  scale_fill_manual(values = polygon_colors) +
  geom_errorbar(aes(ymin = value_avg-value_se, ymax = value_avg+value_se), width=.2,
                position=position_dodge(.9)) +
  geom_point(data = obs_avg, aes(x = obs_groups, y = obs_avg), color = "red", size = 4, inherit.aes = FALSE) + # Add points for Sigrid's avg's
  coord_flip() +
  labs(title = "Polygon Type Avg. CO2 Flux", x = NULL,
       y = "CO2 (umol m-2 s-1)") +
  geom_hline(yintercept = 0, linetype = 2, color = "black", linewidth = 1) +
  #theme_minimal() +
  theme_bw() +
  theme(legend.title = element_blank(), axis.text = element_text(size = 12),  
        legend.text = element_text(size = 12), legend.position="bottom", strip.background = element_rect(fill = "white")) 



#set up function for facet tags
tag_facet <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf, 
                      hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
  
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
  p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE) 
}




#full_df_2D_hourly_sub$date <- as.Date(full_df_2D_hourly_sub$date, format='%Y-%m-%d')
full_df_2D_hourly_sub$month <- as.numeric(full_df_2D_hourly_sub$month)
# Create a combined date-like column
full_df_2D_hourly_sub$plot_date <- make_date(month = full_df_2D_hourly_sub$month, day = full_df_2D_hourly_sub$day,
                                             year = full_df_2D_hourly_sub$year)

# Convert Month to a factor with ordered levels
#test$month <- factor(test$month, ordered = TRUE)
test$plot_date <- as.character(test$plot_date)
test$plot_date <- as.Date(test$plot_date)

p <- ggplot(data = test) + #full_df_2D_hourly_sub
  facet_grid(variable ~ year, scales = "free_y") + #, scales = "free_y"allows y axis to vary, but same xaxis for all
  geom_smooth(aes(x = month, y= value, color = polygon_type)) + #linetype = "dotdash") + #add linear fit line
  scale_color_manual(values = c("FCPcenter" = "#abba8f", "LCPcenter" = "#2a7c6e", "LCPtrough" = "#2aaeae", "HCPcenter" = "#781e1f", 
                                "HCPtrough" = "#c85c35", "LCPrim" = "#b4f0f8", "FCPtrough" = "#cdd0bb")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #scale_x_continuous(breaks=seq(4,11,1)) +
  labs(x = "Month", y = NULL) + 
  #guides(fill = "none") + #don't display fill legend
  # scale_x_date(
  #   limits = c(as.Date("1998-04-01"), as.Date("1998-11-30")), # Set date range
  #   date_labels = "%b", # Show abbreviated month names on the x-axis
  #   date_breaks = "1 month" # Breaks for each month
  # ) +
  scale_x_continuous(
    limits = c(4, 11), # April (4) to November (11)
    breaks = 4:11, # Show ticks for April to November
    labels = month.abb[4:11] # Abbreviated month names
  ) +
  theme_bw() +
  theme(legend.title = element_blank(), axis.text = element_text(size = 10), axis.title.x = element_text(size = 12), 
        legend.position = "top", legend.text = element_text(size = 12), strip.text.x = element_text(size = 12), 
        strip.text.y = element_text(size = 12), strip.background = element_blank()) 

#tag_facet(p) 








p <- ggplot(data = full_df_2D_hourly_sub) + 
  facet_grid(variable ~ year(plot_date), scales = "free_y") + #, scales = "free_y"allows y axis to vary, but same xaxis for all
  geom_smooth(aes(x = plot_date, y= value, color = polygon_type), linewidth = 0.8) + #linetype = "dotdash") + #add linear fit line
  scale_color_manual(values = c("FCPcenter" = "#abba8f", "LCPcenter" = "#2a7c6e", "LCPtrough" = "#2aaeae", "HCPcenter" = "#781e1f", 
                                "HCPtrough" = "#c85c35", "LCPrim" = "#b4f0f8", "FCPtrough" = "#cdd0bb")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  labs(x = "Month", y = NULL) + 
  #guides(fill = "none") + #don't display fill legend
  theme_bw() +
  theme(legend.title = element_blank(), axis.text = element_text(size = 10), axis.title.x = element_text(size = 12), 
        legend.position = "top", legend.text = element_text(size = 12), strip.text.x = element_text(size = 12), 
        strip.text.y = element_text(size = 12), strip.background = element_blank()) 
















