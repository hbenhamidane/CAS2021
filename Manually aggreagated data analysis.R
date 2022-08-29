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

# Reading and cleaning up data
data_ludo <- read.csv("Data_EU_aggregated_custom_perYear_from_disaggregated.csv", header = T, na.strings = na_string)

# Reading the metadata
metadata <- read.csv("Waterbase_v2021_1_S_WISE6_SpatialObject_DerivedData.csv", header = T, na.strings = na_string)
metadata <- metadata[!duplicated(metadata),]

# Adding the metadata info
data_agg <- data_ludo %>% left_join(metadata, by = "monitoringSiteIdentifier")
# removing unused df
rm(data_ludo, metadata)
# Removing unused columns
data_agg <- data_agg %>% select(1:9, 17, 24, 26, 28, 29)
# Removing all the missing values
data_agg <- data_agg %>% filter(!is.na(lat) & !is.na(lon) & !is.na(rbdName))
# Removing all duplicated values
data_agg <- data_agg[!duplicated(data_agg),]
# Plotting distribution of measurements across years
agg <- data_agg %>% group_by(year, countryCode) %>% summarize(count=n())
ggplot(data=agg, aes(x=year, y=count, fill=countryCode)) + geom_bar(position="stack", stat="identity") + xlim(1990, 2022)


# Arbitrary filtering: year = 2020
data_agg_2020 <- data_agg %>% filter(year == 2020)
# Taking the 10 most measured ATs
data_agg_2020_recap <- data_agg_2020 %>% group_by(observedPropertyDeterminandLabel) %>% 
  summarize(sites = n_distinct(monitoringSiteIdentifier)) %>% arrange(desc(sites))
targets <- data_agg_2020_recap$observedPropertyDeterminandLabel[1:10]
data_agg_2020 <- data_agg_2020 %>% filter(observedPropertyDeterminandLabel %in% targets)
rm(data_agg_2020_recap)

# Some entries are "duplicated" due to some decimal increment in the gps coordinates -_-
data_agg_2020$sorter <- str_c(data_agg_2020$monitoringSiteIdentifier, data_agg_2020$observedPropertyDeterminandLabel, sep="_")
data_agg_2020 <- data_agg_2020[!duplicated(data_agg_2020$sorter),]
data_agg_2020 <- data_agg_2020 %>% select(-sorter)

# Recap of filtered dataset
data_agg_2020_recap <-data_agg_2020 %>% group_by(observedPropertyDeterminandLabel, countryCode) %>% summarize(n_sites=n_distinct(monitoringSiteIdentifier), 
                                                                                             n_measurements = sum(resultObservedValue_count),
                                                                                             med_measurements = median(resultObservedValue_count),
                                                                                             min_measurements = min(resultObservedValue_count),
                                                                                             max_measurements = max(resultObservedValue_count))

# Comparing which of the 10 top measured ATs are measured accross the most countries
data_agg_2020_recap_country <- data_agg_2020_recap %>% summarise(n_country = n_distinct(countryCode), countries = list(countryCode))
count_per_country <- table(unlist(data_agg_2020_recap_country$countries))
countries <- row.names(count_per_country[count_per_country == 10])
rm(data_agg_2020_recap_country, count_per_country)

# Using countries as a subset for data_agg_2020
data_agg_2020 <- data_agg_2020 %>% filter(countryCode %in% countries)
# Recap of filtered dataset
data_agg_2020_recap <-data_agg_2020 %>% group_by(observedPropertyDeterminandLabel, countryCode) %>% summarize(n_sites=n_distinct(monitoringSiteIdentifier), 
                                                                                                              n_measurements = sum(resultObservedValue_count),
                                                                                                              med_measurements = median(resultObservedValue_count),
                                                                                                              min_measurements = min(resultObservedValue_count),
                                                                                                              max_measurements = max(resultObservedValue_count)) %>%
  arrange(observedPropertyDeterminandLabel, desc(n_measurements), desc(n_sites))


# Visualizing the measurement count per AT and country
ggplot(data=data_agg_2020_recap, aes(x=as.factor(observedPropertyDeterminandLabel), y=n_measurements, fill=countryCode)) + 
  geom_bar(position="stack", stat="identity") + theme(axis.text.x = element_text(angle = 90))
# Plotting distribution of measurement counts
ggplot(data=data_agg_2020, aes(x=resultObservedValue_count)) +geom_histogram(binwidth=1) +xlim(0,50)

# Pivoting data_agg_2020
data_agg_2020_w <- data_agg_2020 %>% 
  pivot_wider(id_cols = c(monitoringSiteIdentifier, year, lat, lon, countryCode, specialisedZoneType, rbdIdentifier, rbdName),
              names_from = observedPropertyDeterminandLabel,
              values_from = c(resultObservedValue_mean, resultObservedValue_std, resultObservedValue_count, resultQualityObservedValueBelowLOQ_sum, resultQualityObservedValueBelowLOQ_perc),
              values_fill = 0)

# Quick visualization of data_agg_2020_w
vis_dat(data_agg_2020_w, warn_large_data = F)
# Lot of missing values again
# imputation of missing values to 0 --> rerun pivot_wider with values_fill = 0
# NAs still present in std columns, removing them due to incompleteness and to reduce data dimensionality for PCA
data_agg_2020_w <- data_agg_2020_w %>% select(-starts_with("resultObservedValue_std"))

# Extracting starts_with("resultObservedValue_mean") columns and converting to a matrix, passing monitoringSiteIdentifier as row.names
matrix <- as.matrix(data_agg_2020_w %>% select(starts_with("resultObservedValue_mean")))
row.names(matrix) <- data_agg_2020_w$monitoringSiteIdentifier

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

# No conclusive separation of data points --> revisiting the pivot, not using fill_values and dropping all NAs in resultObservedValue_mean columns
data_agg_2020_w2 <- data_agg_2020 %>% 
  pivot_wider(id_cols = c(monitoringSiteIdentifier, year, lat, lon, countryCode, specialisedZoneType, rbdIdentifier, rbdName),
              names_from = observedPropertyDeterminandLabel,
              values_from = c(resultObservedValue_mean, resultObservedValue_std, resultObservedValue_count, resultQualityObservedValueBelowLOQ_sum, resultQualityObservedValueBelowLOQ_perc))
# Cleaning column names for easier manipulation
data_agg_2020_w2 <- clean_names(data_agg_2020_w2)
# removing all entries where an AT is not measured (missing mean value)
data_agg_2020_w2 <- data_agg_2020_w2 %>% filter(!is.na(result_observed_value_mean_nitrate) &
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
matrix2 <- as.matrix(data_agg_2020_w2 %>% select(starts_with("result_observed_value_mean")))
row.names(matrix2) <- data_agg_2020_w2$monitoring_site_identifier

# Running a PCA on matrix
pca2 <- prcomp(matrix2, center = T, scale. = T)
fviz_pca_biplot(pca2, axes = c(1,2), repel = T)

# K-means
## Step 1: determining the optimal number of clusters using elbow/silhouette plots
fviz_nbclust(matrix2, kmeans, method = "silhouette")
## Step 2: kmeans clustering with the optimal number of clusters determined in Step 1 and 100 random cluster center starting positions
km_out2 <- kmeans(matrix2, centers=2, nstart=1000)
km_clust2 <- km_out2$cluster
## Step 3: Creating a visualization to evaluate the clustering along the first 2 PCA dimensions
fviz_cluster(list(data=matrix2, cluster = km_clust2))



