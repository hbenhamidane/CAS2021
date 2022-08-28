# Loading libraries
library(tidyverse)
library(lubridate)
library(janitor)
library(visdat)
library(ggplot2)
library(stats)
library(caret)
library(factoextra)
library(caret)
library(mclust)
library(gridExtra)

########################################################################################################################################################################################
########################################################################################################################################################################################
########################################################################################################################################################################################

# Defining work directory and general variables
setwd("C:/Users/shami/Documents/CAS 2021/CAS2021")
wd <- getwd()
na_string <- c("", " ", NA, "NA", "N/A")

# Reading and cleaning up data; manually aggregated dataset for year 2020, split by site, AT and semester
data_ludo <- read.csv("Data_EU_aggregated_perSiteTargetSeason_2020.csv", header = T, na.strings = na_string)

# Reading the metadata
metadata <- read.csv("Waterbase_v2021_1_S_WISE6_SpatialObject_DerivedData.csv", header = T, na.strings = na_string)
metadata <- metadata[!duplicated(metadata),]

# Adding the metadata info
data_agg <- data_ludo %>% left_join(metadata, by = "monitoringSiteIdentifier")
# removing unused df
rm(data_ludo, metadata)
# Standardizing column names
data_agg <- clean_names(data_agg)
# Removing unused columns
data_agg <- data_agg %>% select(1:9, 17, 24, 26, 28, 29)
# Removing all the missing values
data_agg <- data_agg %>% filter(!is.na(lat) & !is.na(lon) & !is.na(rbd_name))
# Removing all duplicated values
data_agg <- data_agg[!duplicated(data_agg),]
# Converting season to a factor column
data_agg[str_detect(data_agg$season, pattern="winter"),]$season <- "W"
data_agg[str_detect(data_agg$season, pattern="summer"),]$season <- "S"
data_agg$season <- as.factor(data_agg$season)
# Visualizing the extent of missing values in agg
vis_dat(data_agg, warn_large_data = F)
# Creating a reference dataset for the different analysis
data_agg_seed <- data_agg

########################################################################################################################################################################################
########################################################################################################################################################################################
########################################################################################################################################################################################

# Creating a summary of the data 
agg <- data_agg %>% group_by(country_code, rbd_name, monitoring_site_identifier, season) %>% 
                    summarize(AT_measured=n_distinct(observed_property_determinand_label),
                              num_measurements = sum(result_observed_value_count),
                              ATs = list(observed_property_determinand_label))
# Trying to identify comment ATs measured across all groups i.e. combination of country_code/rbd_name/monitoring_site_identifer/season
Reduce(intersect, agg$ATs)
# character(0)
# No common AT measured accross all groups

# Create a targeted subset of the data retaining only the most commonly measured ATs or specific ATs of interest
agg_top10 <- data_agg %>% group_by(observed_property_determinand_label) %>% 
  summarize(sites = n_distinct(monitoring_site_identifier)) %>% arrange(desc(sites))

# Taking the 10 most measured ATs
targets <- agg_top10$observed_property_determinand_label[1:10]
data_agg <- data_agg %>% filter(observed_property_determinand_label %in% targets)
rm(agg_top10)

# Some entries are "duplicated" due to some decimal increment in the gps coordinates
data_agg$sorter <- str_c(data_agg$monitoring_site_identifier, data_agg$observed_property_determinand_label, data_agg$season, sep="_")
data_agg <- data_agg[!duplicated(data_agg$sorter),]
data_agg <- data_agg %>% select(-sorter)

# Comparing which of the 10 top measured ATs are measured accross the most countries
data_agg <- data_agg %>% ungroup()
data_agg_countries <- data_agg %>% group_by(observed_property_determinand_label, season) %>% summarise(n_country = n_distinct(country_code), countries = list(country_code))
countries <- Reduce(intersect, data_agg_countries$countries)
# Using countries as a subset for data_agg
data_agg <- data_agg %>% filter(country_code %in% countries)
rm(data_agg_countries)

# Recap of filtered again dataset
data_agg_recap <-data_agg %>% group_by(observed_property_determinand_label, country_code, season) %>% summarize(n_sites=n_distinct(monitoring_site_identifier), 
                                                                                                              n_measurements = sum(result_observed_value_count),
                                                                                                              med_measurements = median(result_observed_value_count),
                                                                                                              min_measurements = min(result_observed_value_count),
                                                                                                              max_measurements = max(result_observed_value_count)) %>%
  arrange(observed_property_determinand_label, desc(n_measurements), desc(n_sites))

# Visualizing the measurement count per AT and country
ggplot(data=data_agg_recap, aes(x=as.factor(observed_property_determinand_label), y=n_measurements, fill=country_code)) + 
  geom_bar(position="stack", stat="identity") + theme(axis.text.x = element_text(angle = 90))

# Plotting distribution of measurement counts
ggplot(data=data_agg, aes(x=result_observed_value_count, fill=season)) +geom_histogram(binwidth=1) +xlim(0,25) + xlab("Number of distinct analytical targets measured") + ylab("Number of samples measured")

# Pivoting data_agg with 0 inputation of missing values
data_agg_w <- data_agg %>% 
  pivot_wider(id_cols = c(monitoring_site_identifier, season, lat, lon, country_code, specialised_zone_type, rbd_identifier, rbd_name),
              names_from = observed_property_determinand_label,
              values_from = c(result_observed_value_mean, result_observed_value_std, result_observed_value_count, result_quality_observed_value_below_loq_sum, result_quality_observed_value_below_loq_perc),
              values_fill = 0)

# Quick visualization of data_agg_2020_w
# vis_dat(data_agg_w, warn_large_data = F)
# Lot of missing values again
# NAs still present in std columns, removing them due to incompleteness and to reduce data dimension for PCA
data_agg_w <- data_agg_w %>% select(-starts_with("resultObservedValue_std"))

# Extracting starts_with("resultObservedValue_mean") columns and converting to a matrix, passing monitoringSiteIdentifier as row.names
matrix <- as.matrix(data_agg_w %>% select(starts_with("result_observed_value_mean")))
row.names(matrix) <- str_c(data_agg_w$monitoring_site_identifier, data_agg_w$season, sep="_")

# Running a PCA on matrix
pca <- prcomp(matrix, center = T, scale. = T)
fviz_pca_biplot(pca, axes = c(1,2), repel = T)

# K-means
## Step 1: determining the optimal number of clusters using elbow/silhouette plots
fviz_nbclust(matrix, kmeans, method = "silhouette")
## Step 2: kmeans clustering with the optimal number of clusters determined in Step 1 and 100 random cluster center starting positions
km_out <- kmeans(matrix, centers=2, nstart=200)
km_clust <- km_out$cluster
## Step 3: Creating a visualization to evaluate the clustering along the first 2 PCA dimensions
fviz_cluster(list(data=matrix, cluster = km_clust))

rm(matrix, data_agg, data_agg_recap, data_agg_w, countries, targets, km_out, pca, km_clust)

# No conclusive groups or underlying structure in the dataset

########################################################################################################################################################################################
########################################################################################################################################################################################
########################################################################################################################################################################################

# Revisiting the pivot, not using fill_values and dropping all NAs in resultObservedValue_mean columns
data_agg <- data_agg_seed
# Create a targeted subset of the data retaining only the most commonly measured ATs or specific ATs of interest
agg_top10 <- data_agg %>% group_by(observed_property_determinand_label) %>% 
  summarize(sites = n_distinct(monitoring_site_identifier)) %>% arrange(desc(sites))

# Taking the 10 most measured ATs
targets <- agg_top10$observed_property_determinand_label[1:10]
data_agg <- data_agg %>% filter(observed_property_determinand_label %in% targets)
rm(agg_top10)

# Some entries are "duplicated" due to some decimal increment in the gps coordinates
data_agg$sorter <- str_c(data_agg$monitoring_site_identifier, data_agg$observed_property_determinand_label, data_agg$season, sep="_")
data_agg <- data_agg[!duplicated(data_agg$sorter),]
data_agg <- data_agg %>% select(-sorter)

# Comparing which of the 10 top measured ATs are measured accross the most countries
data_agg <- data_agg %>% ungroup()
data_agg_countries <- data_agg %>% group_by(observed_property_determinand_label, season) %>% summarise(n_country = n_distinct(country_code), countries = list(country_code))
countries <- Reduce(intersect, data_agg_countries$countries)
# Using countries as a subset for data_agg
data_agg <- data_agg %>% filter(country_code %in% countries)
rm(data_agg_countries)

# Pivoting without any inputation of missing values
data_agg_w <- data_agg %>% 
  pivot_wider(id_cols = c(monitoring_site_identifier, season, lat, lon, country_code, specialised_zone_type, rbd_identifier, rbd_name),
              names_from = observed_property_determinand_label,
              values_from = c(result_observed_value_mean, result_observed_value_std, result_observed_value_count, result_quality_observed_value_below_loq_sum, result_quality_observed_value_below_loq_perc))
# Cleaning column names for easier manipulation
data_agg_w <- clean_names(data_agg_w)
# removing all entries where an AT is not measured (missing mean value)
data_agg_w <- data_agg_w %>% filter(!is.na(result_observed_value_mean_nitrate) &
                                    !is.na(result_observed_value_mean_nitrite) &
                                    !is.na(result_observed_value_mean_phosphate) &
                                    !is.na(result_observed_value_mean_ammonium) &
                                    !is.na(result_observed_value_mean_dissolved_oxygen) &
                                    !is.na(result_observed_value_mean_water_temperature) &
                                    !is.na(result_observed_value_mean_electrical_conductivity) &
                                    !is.na(result_observed_value_mean_p_h) &
                                    !is.na(result_observed_value_mean_chloride) &
                                    !is.na(result_observed_value_mean_sulphate))

# Extracting starts_with("resultObservedValue_mean") columns and converting to a matrix, passing monitoringSiteIdentifier as row.names
matrix <- as.matrix(data_agg_w %>% select(starts_with("result_observed_value_mean")))
row.names(matrix) <- str_c(data_agg_w$monitoring_site_identifier, data_agg_w$season, sep="_")

# Running a PCA on matrix
pca <- prcomp(matrix, center = T, scale. = T)
fviz_pca_biplot(pca, axes = c(1,2), repel = T)

# K-means
## Step 1: determining the optimal number of clusters using elbow/silhouette plots
fviz_nbclust(matrix, kmeans, method = "silhouette")
## Step 2: kmeans clustering with the optimal number of clusters determined in Step 1 and 100 random cluster center starting positions
km_out <- kmeans(matrix, centers=2, nstart=1000)
km_clust <- km_out$cluster
## Step 3: Creating a visualization to evaluate the clustering along the first 2 PCA dimensions
fviz_cluster(list(data=matrix, cluster = km_clust))

rm(matrix, data_agg, data_agg_w, countries, targets, pca, km_out, matrix, km_clust)

# Excluding imputation of missing values and more stringent filtering did not allow for better clustering of the data

########################################################################################################################################################################################
########################################################################################################################################################################################
########################################################################################################################################################################################

# Create a targeted subset of the data retaining only specific ATs most likely to exhibit seasonal variations (targeted analysis): electrical conductivity, pH and water temperature
targets <- c("Electrical conductivity", "pH", "Water temperature")
data_agg <- data_agg_seed
data_agg <- data_agg %>% filter(observed_property_determinand_label %in% targets)

# Some entries are "duplicated" due to some decimal increment in the gps coordinates
data_agg$sorter <- str_c(data_agg$monitoring_site_identifier, data_agg$observed_property_determinand_label, data_agg$season, sep="_")
data_agg <- data_agg[!duplicated(data_agg$sorter),]
data_agg <- data_agg %>% select(-sorter)

# Recap of filtered dataset
data_agg_recap <-data_agg %>% group_by(observed_property_determinand_label, country_code, season) %>% summarize(n_sites=n_distinct(monitoring_site_identifier), 
                                                                                                                n_measurements = sum(result_observed_value_count),
                                                                                                                med_measurements = median(result_observed_value_count),
                                                                                                                min_measurements = min(result_observed_value_count),
                                                                                                                max_measurements = max(result_observed_value_count))

# Comparing which of the 10 top measured ATs are measured accross the most countries
data_agg_recap <- data_agg_recap %>% ungroup()
data_agg_recap_country <- data_agg_recap %>% group_by(observed_property_determinand_label, season) %>% summarise(n_country = n_distinct(country_code), countries = list(country_code))
countries <- Reduce(intersect, data_agg_recap_country$countries)
# Using countries as a subset for data_agg
data_agg <- data_agg %>% filter(country_code %in% countries)
rm(data_agg_recap, data_agg_recap_country)

# Recap of filtered again dataset
data_agg_recap <-data_agg %>% group_by(observed_property_determinand_label, country_code, season) %>% summarize(n_sites=n_distinct(monitoring_site_identifier), 
                                                                                                                n_measurements = sum(result_observed_value_count),
                                                                                                                med_measurements = median(result_observed_value_count),
                                                                                                                min_measurements = min(result_observed_value_count),
                                                                                                                max_measurements = max(result_observed_value_count)) %>%
  arrange(observed_property_determinand_label, desc(n_measurements), desc(n_sites))

# Visualizing the measurement count per AT and country
ggplot(data=data_agg_recap, aes(x=as.factor(observed_property_determinand_label), y=n_measurements, fill=country_code)) + 
  geom_bar(position="stack", stat="identity") + theme(axis.text.x = element_text(angle = 90))

# Pivoting data_agg without inputation of missing values
data_agg_w <- data_agg %>% 
  pivot_wider(id_cols = c(monitoring_site_identifier, season, lat, lon, country_code, specialised_zone_type, rbd_identifier, rbd_name),
              names_from = observed_property_determinand_label,
              values_from = c(result_observed_value_mean, result_observed_value_std, result_observed_value_count, result_quality_observed_value_below_loq_sum, result_quality_observed_value_below_loq_perc))

# Quick visualization of data_agg_2020_w
vis_dat(data_agg_w, warn_large_data = F)
# Lot of missing values again; Removing all NAs accross mean and std columns
data_agg_w <- data_agg_w %>% filter_at(vars(starts_with("result_observed_value")), all_vars(!is.na(.)))
# Extent of missing data removed: 58%

# Extracting starts_with("resultObservedValue_mean") columns and converting to a matrix, passing monitoringSiteIdentifier as row.names
matrix <- as.matrix(data_agg_w %>% select(starts_with("result_observed_value_mean"),  starts_with("result_observed_value_std")))
row.names(matrix) <- str_c(data_agg_w$monitoring_site_identifier, data_agg_w$season, sep="_")

# Running a PCA on matrix
pca <- prcomp(matrix, center = T, scale. = T)
fviz_pca_biplot(pca, axes = c(1,2), repel = T)

# K-means
## Step 1: determining the optimal number of clusters using elbow/silhouette plots
fviz_nbclust(matrix, kmeans, method = "silhouette")
## Step 2: kmeans clustering with the optimal number of clusters determined in Step 1 and 100 random cluster center starting positions
km_out <- kmeans(matrix, centers=4, nstart=100)
km_clust <- km_out$cluster
## Step 3: Creating a visualization to evaluate the clustering along the first 2 PCA dimensions
fviz_cluster(list(data=matrix, cluster = km_clust))

fviz_cluster(list(data=matrix, cluster = km_clust),axes = c(2, 3))
