# WISE Switzerland data analysis script
# Version: {0.2}
# Author: {Hisham Ben Hamidane}
# Changes:
# 1. Addition of script metadata and version
# 2. Improvement of script comments
# 3. Removal of unnecessary packages: openxlsx, readxl
# 4. Coercing "below_LOQ" variable to logical (previously character) and grouping variables (AT_name, AT_code, WB_type, site_name, etc...) to factors
# 5. Addition of preprocessing for k-means
# 6. Elbow plot for cluster number estimation according to the within sum square method


# Loading the required packages
library(tidyverse)
library(lubridate)
library(janitor)
library(visdat)
library(ggplot2)
library(stats)
library(caret)
library(factoextra)

# 0. Parameters for the analysis
##Defining variables for dataselection
# meas_count_threshold <- 5000

# 1. Loading the data and metadata, merging the 2 together and manipulating variables (selecting/renaming column names, adding new columns) 
## EU dataset for water quality; filtered available database for Switzerland entries (country code = CH)
setwd("C:/Users/shami/Documents/CAS 2021/")
wd <- getwd()
##Reading the data and metadata
data <- read.csv("DataExtract_Switzerland.csv", header = T)
metadata <- read.csv("WISE_spatialobject_CH.csv", header = T)
## Selecting the usefull columns of the metadata to add to data
meta_cols_to_keep <- c("monitoringSiteIdentifier","monitoringSiteName", "waterBodyIdentifier", "waterBodyIdentifierScheme",
                       "waterBodyName", "specialisedZoneType", "subUnitName",
                       "rbdIdentifier", "rbdName", "lon", "lat")
metadata <- metadata %>% select(all_of(meta_cols_to_keep))
##Dropping the useless columns in the data to make it lighter
data_cols_to_drop <- c("parameterSedimentDepthSampled","parameterSpecies",
                       "resultMoisture", "resultFat", "resultExtractableLipid",
                       "resultLipid")
##And left joining data to metadata
data <- data %>% select(-all_of(data_cols_to_drop)) %>% left_join(metadata, by = "monitoringSiteIdentifier")
rm(meta_cols_to_keep, data_cols_to_drop, metadata)
## Selecting and renaming columns
data <- data %>% transmute(site = as.factor(monitoringSiteIdentifier),
                  site = as.factor(site),
                  site_name = as.factor(monitoringSiteName),
                  site_lat =lat,
                  site_lon =lon,
                  WB_name =as.factor(waterBodyName),
                  WB_system_name =as.factor(rbdName),
                  WB_type =as.factor(parameterWaterBodyCategory), 
                  AT_code =as.factor(observedPropertyDeterminandCode),
                  AT_name = as.factor(observedPropertyDeterminandLabel),
                  meas_date =as.Date(data$phenomenonTimeSamplingDate, format="%d/%m/%Y %H:%M:%S"),
                  meas_year =phenomenonTimeSamplingDate_year,
                  measured_value = resultObservedValue,
                  measured_value_unit = as.factor(resultUom),
                  below_LOQ = as.logical(resultQualityObservedValueBelowLOQ),
                  LOQ_value = as.numeric(procedureLOQValue))
## Converting meas_year to a factor vector instead of character
data$meas_year <- as.factor(as.numeric(data$meas_year))
## Creating a month column to look for yearly periodicity in the data
data$meas_month <- as.factor(as.numeric(str_extract(data$meas_date, pattern = "(?<=^.{5})[[:digit:]]{2}")))

# 2. Creating data overview and selecting observations to keep
## Graphical representation of dataframe after edits in 1 using Visdat
# vis_dat(data, warn_large_data = F)
## Percentage of missing values (total)
percent_missing_values <- sum(is.na(data$measured_value))/nrow(data)*100
print(paste("The overall percentage of missing values is", round(percent_missing_values,2), "%", sep = " "))
## Creating a recap table summarizing the mMissing values and values at or below LOQ per analytical target, number of sites where the AT is measured as well as corresponding percentages
data_completeness <- data %>% group_by(AT_code) %>% 
                              summarize(total_values = n(),
                                        missing_values = sum(is.na(measured_value)),
                                        not_measured_values = sum(below_LOQ==T),
                                        number_of_sites = n_distinct(site),
                                        number_of_years = n_distinct(meas_year)) %>% 
                              mutate(missing_perc = missing_values/total_values,
                                     not_measured_perc = not_measured_values/total_values) %>%
                              arrange(desc(total_values), not_measured_perc, missing_perc)  
## Additional information combined between data and data_completeness
data_completeness <- data %>% select(AT_name, AT_code) %>% filter(!duplicated(AT_code)) %>% right_join(data_completeness, by = "AT_code") %>% 
                              mutate(avg_meas_per_site = total_values/number_of_sites,
                                     meas_site_ratio = number_of_sites/length(unique(data$site))*100) %>%
                              arrange(desc(total_values), number_of_sites, desc(avg_meas_per_site), not_measured_values, missing_values)
#Filtering data to retain only the 11 most sampled ATs after removing missing values
data_f <- data %>% filter(!is.na(measured_value) & AT_name %in% data_completeness$AT_name[1:11])
# Creating a conservative data frame where AT values below LOQ are dropped
data_f_conservative <- data_f %>% filter(below_LOQ == F)

# 3. Preparing the data for k-means clustering
## For a first clustering attempt, values per site and AT will be averaged (removing the time dimension for now)
## Step 1: calculate the median and average value for each combination of site and AT
data_f <- data_f %>% group_by(site, AT_code) %>% 
                   mutate(median_measured_value = median(measured_value),
                          avg_measured_value = mean(measured_value)) %>%
                   select(site, WB_name, WB_system_name, WB_type, AT_code, avg_measured_value) %>% 
                   arrange(site, AT_code)
## Remove duplicated entries (each row is still a single measurements but we are only interested in the average measured value per site and AT)
data_f <- data_f[!duplicated(data_f),]
## Further filtering to retain only sites that have a value for all 11 analytical targets
data_f_recap <- data_f %>% mutate(WB_site =str_c(WB_system_name, site, sep="-")) %>% group_by(WB_site) %>% summarize(num_AT_measured = n()) %>% filter(num_AT_measured == 11)
data_f <- data_f %>% mutate(WB_site =str_c(WB_system_name, site, sep="-")) %>% filter(WB_site %in% data_f_recap$WB_site)
## Conversion to a dataframe where each observation is a site and each variable is an AT
data_f <- data_f %>% select(site, AT_code, avg_measured_value) %>% pivot_wider(id_cols = site, names_from = AT_code, values_from = avg_measured_value)
## Extracting site column to use as row name
data_f.index <- data_f$site
data_f <- data_f[-1]
row.names(data_f) <- data_f.index
## Scaling (normalizing the data and converting to a matrix)
data_f.scaled <- scale(data_f)
## Creating the distance matrix
data_f.dist <- dist(data_f)
## Determining the number of clusters using an elbow plot
### synthax from factoextra
fviz_nbclust(data_f, kmeans, method = "wss")





## Overview of the distribution of data between the higher level categories (Water Body system)
table(data_f$WB_system_name, data_f$WB_type)
### uneven distribution of measurements between water bassins; will have to be accounted for when sampling 





# ## Recap tables for the different dataset dimensions: AT, location and date (year)
# # Recap table per AT_name
# recap_table_AT <- data %>% filter(below_LOQ == "False") %>%
#   group_by(AT_name) %>%
#   summarize(min = min(measured_value, na.rm=T), 
#             max = max(measured_value, na.rm=T),
#             mean = mean(measured_value, na.rm=T),
#             sd = sd(measured_value, na.rm=T),
#             unit = unique(measured_value_unit)[1],
#             sample_size = n()) %>%
#   arrange(desc(sample_size))
# # Vector of analytical targets (AT) with a number of measurements above meas_count_threshold
# valid_ATs <- recap_table_AT %>% filter(sample_size >= meas_count_threshold) %>% select(AT_name)
# valid_ATs <- as.vector(valid_ATs$AT_name)
# 
# # Recap table per site_name and top 1 AT
# recap_table_SN <- data %>% filter(below_LOQ == "False") %>%
#   group_by(site_name) %>%
#   summarize(min = min(measured_value, na.rm=T), 
#             max = max(measured_value, na.rm=T),
#             mean = mean(measured_value, na.rm=T),
#             sd = sd(measured_value, na.rm=T),
#             sample_size = n()) %>%
#   arrange(desc(sample_size))
# 
# # Recap table per WB_type
# recap_table_WB <- data %>% filter(below_LOQ == "False") %>%
#   group_by(WB_type) %>%
#   summarize(min = min(measured_value, na.rm=T), 
#             max = max(measured_value, na.rm=T),
#             mean = mean(measured_value, na.rm=T),
#             sd = sd(measured_value, na.rm=T),
#             sample_size = n()) %>%
#   arrange(desc(sample_size))
# 
# # Recap table per meas_year
# recap_table_YR <- data %>% filter(below_LOQ == "False") %>%
#   group_by(meas_year) %>%
#   summarize(min = min(measured_value, na.rm=T), 
#             max = max(measured_value, na.rm=T),
#             mean = mean(measured_value, na.rm=T),
#             sd = sd(measured_value, na.rm=T),
#             sample_size = n()) %>%
#   arrange(desc(sample_size))
# 
# ## Defining different data subsets for different statistical analysis.
# 
# # 1. Measurements of the same AT at the same site over different dates -> time series
# 
# # Creating a recap table for time series
# data_top3AT <- data %>% filter(AT_name %in% valid_ATs[1:3])
# 
# recap_table_top3AT_site_year <- data_top3AT %>% group_by(site, meas_year) %>% 
#   summarize(min = min(measured_value, na.rm=T), 
#             max = max(measured_value, na.rm=T),
#             mean = mean(measured_value, na.rm=T),
#             sd = sd(measured_value, na.rm=T),
#             sample_size = n()) %>%
#   arrange(desc(sample_size), desc(meas_year))
# 
# # CHNTG15 has 48 measurements for the top 3 AT per year for a 5 year period from 2003 to 2008
# data_top3AT_CHNTG15 <- data_top3AT %>% filter(site == "CHNTG15" & meas_year %in% c("2003","2004", "2005", "2006", "2007", "2008"))
# data_top3AT_CHNTG15 <- data_top3AT_CHNTG15 %>% group_by(meas_month, AT_name) %>% mutate(sampled_per_month = n())
# 
# ggplot(data_top3AT_CHNTG15, aes(x=meas_month, y=measured_value)) + geom_point(size=2) + 
#    geom_smooth(aes(group=meas_year), method="glm", se=F) +
#   facet_wrap(AT_name ~ meas_year, nrow=3, scales = "free_y") +
#   theme_light() +
#   labs(title = paste(paste("Measurements of electrical conductivity, nitrate concentration and pH at", data_top3AT_CHNTG15$site_name, sep=" "), "between 2003 and 2008", sep=" "),
#        x = "Measurement month (January to December)",
#        y ="Measured value (resp. uS/cm, mg/L and pH units)") 
# 
# ggplot(data_top3AT_CHNTG15, aes(x=as.factor(meas_month), y=measured_value)) + geom_boxplot(size=1.2, outlier.colour = "red", outlier.size = 2, varwidth = TRUE) + facet_wrap(~AT_name, nrow=1, scales = "free") +
# theme_light() +
#   labs(title = paste(paste("Distributions of measurements of electrical conductivity, nitrate concentration and pH at", data_top3AT_CHNTG15$site_name, sep=" "), "between 2003 and 2008", sep=" "),
#        x = "Measurement month (January to December)",
#        y ="Measured value (resp. uS/cm, mg/L and pH units)")
# 
# # Is the data normally distributed? And at which scale/level?
# ## Per Body of Water type?
# # Plot of distribution of measured valutes for top 3 AT per WB type
# # Facet wrap with free scale
# ggplot(data_top3AT, aes(x=measured_value)) + geom_density(stat="density", na.rm = T) + facet_wrap(AT_name~WB_type, scales="free_y")
# # Due to poor readability, split by AT and faceted only per WB
# # Conductivity
# ggplot(data_top3AT[data_top3AT$AT_name == "Electrical conductivity",], aes(x=measured_value)) + geom_density(stat="density", na.rm = T) + facet_wrap(~WB_type, scales="free_y")
# # pH
# ggplot(data_top3AT[data_top3AT$AT_name == "pH",], aes(x=measured_value)) + geom_density(stat="density", na.rm = T) + facet_wrap(~WB_type, scales="free_y")
# # Nitrate concentration
# ggplot(data_top3AT[data_top3AT$AT_name == "Nitrate",], aes(x=measured_value)) + geom_density(stat="density", na.rm = T) + facet_wrap(~WB_type, scales="free_y")
# 
# # Subsetting to keep only the pH as it appears to follow a pseudo normal distribution
# data_pH <- data %>% filter(AT_name == "pH" & WB_type != "LW")
# 
# # Plotting of density distribution for pH values in either groundwater (GW) or river water (RW)
# ggplot(data_pH, aes(x=measured_value)) + geom_density(stat="density", na.rm = T) + facet_wrap(~WB_type) + coord_cartesian(xlim = c(5, 10)) + theme_light() +
#   labs(title = "Distributions of pH values across all sites and dates for groundwater (GW) and river water (RW)",
#        x = "Measured pH value",
#        y ="Probability density")
# 
# # Performing the Shapiro-Wilk normality test for values from both bodies of water type (GW and RW)
# shapiro.test(x = data_pH[data_pH$WB_type == "GW",]$measured_value)
# ## Note: shapiro wilk test is limited to dataset of max 5000 values (in R), therefore a random sampling of 5000 values has been used to perform the test on RW values 
# shapiro.test(x = sample(data_pH[data_pH$WB_type == "RW",]$measured_value, 5000))
# ## Both test give a very low p-value (10e-16), thus rejecting the null hypothesis (that the data is normally distributed)
# 
# # Visual alternative using a qqplot
# ggplot(data_pH, aes(sample=measured_value)) + stat_qq() + stat_qq_line() + facet_wrap(~WB_type) + theme_light() +
#   labs(title = "QQ plots of pH values across all sites and dates for groundwater (GW) and river water (RW)",
#        x = "Theoretical quantiles",
#        y ="Sample quantiles")
# 
# # Is it any different if we only look at one site (paired independent values)
# data_pH %>% group_by(site) %>% count() %>% arrange(desc(n))
# # 1 CHNTG15     153; site with the most pH measurements overall
# data_pH_CHNTG15 <- data_pH[data_pH$site == "CHNTG15",]
# shapiro.test(x = data_pH_CHNTG15$measured_value)
# # No, the data is still not normal at site level
# 
# # Comparing values based on the 3rd variable: date
# View(data_top3AT_CHNTG15)
# # Creating a new group for date: semester; 
# # winter semester months 10, 11, 12, 1, 2 and 3
# # summer semester months 4, 5, 6, 7, 8 and 9
# data_top3AT_CHNTG15$semester <- if_else(data_top3AT_CHNTG15$meas_month %in% c("1", "2", "3", "10", "11", "12"), "winter", "summer")
# 
# # Checking the data normality per subset
# shapiro.test(x = data_top3AT_CHNTG15[data_top3AT_CHNTG15$AT_name == "pH",]$measured_value)
# data_top3AT_CHNTG15_winter <- data_top3AT_CHNTG15 %>% filter(semester == "winter")
# data_top3AT_CHNTG15_summer <- data_top3AT_CHNTG15 %>% filter(semester == "summer")
# shapiro.test(x = data_top3AT_CHNTG15_winter[data_top3AT_CHNTG15_winter$AT_name == "pH",]$measured_value)
# shapiro.test(x = data_top3AT_CHNTG15_winter[data_top3AT_CHNTG15_summer$AT_name == "pH",]$measured_value)
# 
# # Performing a Wilcoxon test to compare values at Kappelen (CHNTG15) between summer and winter semester
# # for pH values
# wilcox.test(data_top3AT_CHNTG15_summer[data_top3AT_CHNTG15_summer$AT_name == "pH",]$measured_value,
#             data_top3AT_CHNTG15_winter[data_top3AT_CHNTG15_winter$AT_name == "pH",]$measured_value,
#             paired=TRUE,
#             correct = TRUE,
#             conf.int = TRUE,
#             conf.level = 0.95,
#             )
# 
# # for conductivity values
# wilcox.test(data_top3AT_CHNTG15_summer[data_top3AT_CHNTG15_summer$AT_name == "Electrical conductivity",]$measured_value,
#             data_top3AT_CHNTG15_winter[data_top3AT_CHNTG15_winter$AT_name == "Electrical conductivity",]$measured_value,
#             paired=TRUE,
#             correct = TRUE,
#             conf.int = TRUE,
#             conf.level = 0.95,
# )
# # for nitrate concentration values
# wilcox.test(data_top3AT_CHNTG15_summer[data_top3AT_CHNTG15_summer$AT_name == "Nitrate",]$measured_value,
#             data_top3AT_CHNTG15_winter[data_top3AT_CHNTG15_winter$AT_name == "Nitrate",]$measured_value,
#             paired=TRUE,
#             correct = TRUE,
#             conf.int = TRUE,
#             conf.level = 0.95,
# )
# 
# 
# 
# # 2. Comparing measurements of the same AT at the same site  between 2 years -> two groups paired
# ## Comparing pH values between GW and RW using the Mann-Whitney U
# wilcox.test(data_pH$measured_value ~ data_pH$WB_type) 
# 
# ## Comparing the Mann-Whitney U test result with an unpaired t test with Welch correction
# t.test(data_pH$measured_value ~ data_pH$WB_type,
#        paired = FALSE,
#        var.equal = FALSE,
#        conf = 0.95)
# 
# 
# 
# 
# 
# # 3. Comparing measurements of the same AT across n different sites -> n groups unpaired
# 
# # Comparing different ATs across 1 or more sites -> meaningless
# 
# 
# 
# 


       