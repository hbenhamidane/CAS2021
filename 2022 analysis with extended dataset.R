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

<<<<<<< HEAD
########################################################################################################################################################################################
########################################################################################################################################################################################
########################################################################################################################################################################################

# Defining work directory and general variables
setwd("C:/Users/shami/Documents/CAS 2021")
wd <- getwd()
na_string <- c("", " ", NA, "NA", "N/A")

=======

# Defining work directory and general variables
setwd("C:/Users/hisham.benhamidane/OneDrive - Thermo Fisher Scientific/Documents/R/projects/CAS2021/Raw data/2021")
wd <- getwd()
na_string <- c("", " ", NA, "NA", "N/A")


>>>>>>> 16d5e90d9eac54fc080fdfd83b3ab428bc088d73
# Reading and cleaning up data
data <- read.csv("Waterbase_v2021_1_T_WISE6_AggregatedData.csv", header = T, na.strings = na_string)
## renaming the 1st column (special char issue)
colnames(data)[1] <- "monitoringSiteIdentifier"
## dropping the useless columns in the data to make it lighter
data_cols_to_drop <- c( "remarks","metadata_versionId","metadata_beginLifeSpanVersion","metadata_statusCode","metadata_observationStatus",
<<<<<<< HEAD
                        "metadata_statements", "UID")
## and removing rows without site id
data <- data %>% select(-all_of(data_cols_to_drop)) %>% filter(!is.na(monitoringSiteIdentifier))
=======
                        "metadata_statements", "UID"   )
## and removing rows without site id
data <- data %>% select(-all_of(data_cols_to_drop)) %>% filter(!is.na(monitoringSiteIdentifier))


>>>>>>> 16d5e90d9eac54fc080fdfd83b3ab428bc088d73
# Reading and cleaning up metadata
metadata <- read.csv("Waterbase_v2021_1_S_WISE6_SpatialObject_DerivedData.csv", header = T, na.strings = na_string)
## renaming the 1st column (special char issue)
colnames(metadata)[1] <- "countrycode"
## selecting the useful columns of the metadata to add to data
meta_cols_to_keep <- c("countrycode", "monitoringSiteIdentifier","monitoringSiteName", "waterBodyIdentifier",
                       "waterBodyName", "specialisedZoneType", "subUnitName",
                       "rbdIdentifier", "rbdName", "lon", "lat")
metadata <- metadata %>% select(all_of(meta_cols_to_keep)) %>% filter(!is.na(monitoringSiteIdentifier))
## Remove duplicated rows from metadata
metadata <- metadata[!duplicated(metadata),]

# Removing ambiguous entries from metadata
## Some monitoringSiteIdentifier in metadata are repeated 2 times; Solution was to exclude them from analysis
metadata_dup <- metadata %>% group_by(monitoringSiteIdentifier) %>% summarize(site_count = n()) %>% filter(site_count>1) %>% select(monitoringSiteIdentifier)
metadata <- metadata[!metadata$monitoringSiteIdentifier %in% metadata_dup$monitoringSiteIdentifier,]
## Dropping the corresponding sites in data
data <- data %>% filter(monitoringSiteIdentifier %in% unique(metadata$monitoringSiteIdentifier))
<<<<<<< HEAD
=======


>>>>>>> 16d5e90d9eac54fc080fdfd83b3ab428bc088d73
# Left joining data to metadata
data <- data %>% left_join(metadata, by = "monitoringSiteIdentifier" )

# Quality checks post merge
## Checking for the % of unique(monitoringSiteIdentifier) from data in metadata
sum(unique(data$monitoringSiteIdentifier) %in% unique(metadata$monitoringSiteIdentifier))/length(unique(data$monitoringSiteIdentifier))*100
# [1] 100

<<<<<<< HEAD
# Removing all non-necessary variables
rm(meta_cols_to_keep, data_cols_to_drop, metadata, metadata_dup)

=======

# removing all non-necessary variables
rm(meta_cols_to_keep, data_cols_to_drop, metadata, metadata_dup)


>>>>>>> 16d5e90d9eac54fc080fdfd83b3ab428bc088d73
# Selecting and renaming columns in the merged data 
data <- data %>% transmute(site = as.factor(monitoringSiteIdentifier),
                           site_name = as.factor(monitoringSiteName),
                           site_lat =lat,
                           site_lon =lon,
                           country = as.factor(countrycode),
                           WB_name =as.factor(waterBodyName),
                           WB_system_name =as.factor(rbdName),
                           WB_type =as.factor(parameterWaterBodyCategory), 
                           AT_code =as.factor(observedPropertyDeterminandCode),
                           AT_name = as.factor(observedPropertyDeterminandLabel),
                           meas_year =phenomenonTimeReferenceYear,
                           number_of_samples = as.numeric(resultNumberOfSamples),
                           measured_value_min = as.numeric(resultMinimumValue),
                           measured_value_mean = as.numeric(resultMeanValue),
                           measured_value_median = as.numeric(resultMedianValue),
                           measured_value_max = as.numeric(resultMaximumValue),
                           measured_value_sd = as.numeric(resultStandardDeviationValue),
                           measured_value_unit = as.factor(resultUom),
                           sampling_depth = parameterSampleDepth,
<<<<<<< HEAD
                           meas_period = parameterSamplingPeriod,
                           number_of_samples_below_LOQ = as.numeric(resultQualityNumberOfSamplesBelowLOQ),
                           LOQ_value = as.numeric(procedureLOQValue))

# Reassigning NAs in data$number_of_samples_below_LOQ to 0
data$number_of_samples_below_LOQ <- coalesce(data$number_of_samples_below_LOQ, 0) 

# Adding a column describing the percentage of measurements below LOQ
data$below_LOQ_perc <- data$number_of_samples_below_LOQ/data$number_of_samples

# Changing the format of meas_period to factor
data$meas_period <- as.factor(data$meas_period)

# Removing missing values
data <- data %>% filter(!is.na(measured_value_mean) & !is.na(number_of_samples))

# Plotting the distribution of measurements below LOQ percentage for the whole dataset
ggplot(data[!is.na(data$below_LOQ_perc),], aes(x=below_LOQ_perc)) +geom_histogram(stat="bin")

# Exporting data as .csv file
# write.csv(data, file="WISE EU data after merging with metadata and clean up.csv")


########################################################################################################################################################################################
########################################################################################################################################################################################
########################################################################################################################################################################################

# Filtering to keep only the 3 most measured ATs
filter <- data %>% group_by(AT_name) %>% summarize(meas_per_AT = n()) %>% arrange(desc(meas_per_AT)) %>% select(AT_name) %>% head(3)
data <- data %>% filter(AT_name %in% filter$AT_name)
rm(filter)
# Identifying and removing sites to drop based on multiplicity after grouping by site, year and AT
sites_to_drop <- data %>% group_by(site, meas_year, AT_code) %>% summarize(repetition = n(), tot_number_of_samples = sum(number_of_samples)) %>% arrange(desc(repetition)) %>% filter(repetition > 1)
data <- data %>% filter(!site %in% sites_to_drop$site)
rm(sites_to_drop)
# Pivoting the dataset to have a column for each combination of AT and measured value type (min, mean, max, sd, unit, number of samples, etc...)
data_w <- data %>% group_by(site, meas_year) %>% select(-c(AT_code, sampling_depth, meas_period,  measured_value_median,  below_LOQ_perc)) %>% 
  pivot_wider(id_cols = c(site, site_name, site_lat, site_lon, country, WB_name, WB_system_name, WB_type, meas_year) ,
              names_from = AT_name,
              values_from = c(measured_value_min, measured_value_mean, measured_value_max, measured_value_sd, measured_value_unit, number_of_samples, number_of_samples_below_LOQ, LOQ_value))

# Removing entries for which not all 3 ATs are measured
data_w <- data_w %>% filter(!is.na(number_of_samples_Phosphate) & !is.na(number_of_samples_Nitrate) & !is.na(number_of_samples_Ammonium))

# Removing rows that have an NA value in the sd columns
data_w <- data_w %>% filter(!is.na(measured_value_sd_Nitrate) & !is.na(measured_value_sd_Phosphate) & !is.na(measured_value_sd_Ammonium))

# Removing rows that have an NA value in the min and/or max columns
data_w <- data_w %>% filter(!is.na(measured_value_min_Nitrate) & !is.na(measured_value_min_Phosphate) & !is.na(measured_value_min_Ammonium))
data_w <- data_w %>% filter(!is.na(measured_value_max_Nitrate) & !is.na(measured_value_max_Phosphate) & !is.na(measured_value_max_Ammonium))

# Creating a unique identifier
data_w$identifier <- str_c(data_w$site, data_w$WB_type, data_w$meas_year)
row_names <- data_w$identifier

# Ungrouping data_w
data_w <- ungroup(data_w)

# Creating a column subsetter to split the data into a numerical matrix and a metadata df
# col_sub <- colnames(data_w %>%  select(starts_with("measured_value_min") |starts_with("measured_value_mean") | starts_with("measured_value_max") | starts_with("measured_value_sd"))) 
# col_sub <- colnames(data_w %>%  select(starts_with("measured_value_mean") | starts_with("measured_value_sd"))) 
col_sub <- colnames(data_w %>%  select(starts_with("measured_value_min") | starts_with("measured_value_max")))

## splitting into a matrix for k means by retaining only the mean and sd columns
metadata_w <- data_w %>%  select(-col_sub)
data_w <- data_w %>% select(col_sub)

# Converting to a matrix and passing identifier as row names
data_w <- as.matrix(data_w)
rownames(data_w) <- row_names

# Performing PCA
data_w_pca <- prcomp(data_w, center = T, scale. = T)
fviz_pca_biplot(data_w_pca, axes = c(1,2), repel = T)





# K means on matrix_s1wcf
## Step 1: determining the optimal number of clusters using elbow/silhouette plots
fviz_nbclust(data_w, kmeans, method = "silhouette")
# ERROR MESSAGE
# Error in silhouette.default(cluster, d) : 
#   long vectors (argument 1) are not supported in .C

## Step 2: kmeans clustering with the optimal number of clusters determined in Step 1 and 100 random cluster center starting positions
km_out <- kmeans(data_w, centers=2, nstart=5000)
km_clust <- km_out$cluster

## Step 3: Creating a visualization to evaluate the clustering along the first 2 PCA dimensions
### Approach 1: site average; optimal number of clusters k =
fviz_cluster(list(data=matrix_1w, cluster = km_clust))
# fviz_cluster(list(data=matrix_s1wcf, cluster = km_clust), repel = T, show.clust.cent = T, ellipse.alpha = 0.15, ggtheme = theme_classic())


















# Analytical scenario 1
## Group 1: usual analysis: pH, T, dissolved Oxygen, Phosphate, Nitrate and Ammonium (alternative to phosphate & nitrate would be total P/N respectively)  
AT_group1 <- c("pH", "Water Temperature", "Oxygen saturation", "Electrical conductivity", "Ammonium", "Phosphate", "Nitrate", "BOD5")

## FILTER:  AT_group 1, below_LOQ_perc < 1, meas_year >= 2000
data_1 <- data %>% filter(AT_name %in% AT_group1 & meas_year >= 2000 & below_LOQ_perc < 1)

## Creation of a unique identifier column based on the combination of site, AT_code and meas_year for unambiguous identification of entries
data_1$identifier <- str_c(data_1$site, data_1$AT_code, data_1$meas_year, sep="_")

# PROBLEM 1:  FIX THE PROBLEM OF MULTIPLE ENTRIES PER SITE + MEAS_YEAR + AT_NAME !!!
# Number of unique idenfiers in data_S1
# length(unique(data_s1$identifier))
# [1] 115029
# After investigation it was identified that the cause for replicated identifier is multiple sampling depth per identifier
# Summarizing the replication count per identifier and providing the number of rows (i.e. number of unique identifiers) with a replication > 1
# data_s1 %>% group_by(identifier) %>% summarize(sampling_depths = n_distinct(sampling_depth)) %>% arrange(desc(sampling_depths)) %>%
#   filter(sampling_depths > 1) %>% summarize(row_total = sum(sampling_depths), number_of_rows = nrow(.))
# A tibble: 1 × 2
# row_total number_of_rows
# <int>          <int>
#   1       275             67
# Verification: 275 (total replication count) - 67 (sum of each unique identifier with replication >1) = 208
# Subtracting the result (208) from the nrow(data_S1) = length(unique(data_S1$identifier)) thus accounting for all identifier replications
# SOLUTION 1: remove those identifiers alltogether
## Creating a df of identifiers to drop
to_drop <- data_1 %>% group_by(identifier) %>% summarize(sampling_depths = n_distinct(sampling_depth)) %>% arrange(desc(sampling_depths)) %>%
  filter(sampling_depths > 1) %>% select(identifier)
## Removing those identifiers from data_S1
data_1 <- data_1 %>% filter(!identifier %in% to_drop$identifier)
rm(to_drop)
# SOLUTION 1 SUCCESSFULLY IMPLEMENTED

# Pivoting data_S1 into a wider daframe with a column for each AT mean and sd
data_1w <- data_1 %>% group_by(site, meas_year) %>% select(-c(AT_code, sampling_depth, meas_period, identifier, measured_value_min, measured_value_median, measured_value_max, below_LOQ_perc)) %>% 
    pivot_wider(id_cols = c(site, site_name, site_lat, site_lon, country, WB_name, WB_system_name, WB_type, meas_year) ,
                names_from = AT_name,
                values_from = c(measured_value_mean, measured_value_sd, measured_value_unit, number_of_samples, number_of_samples_below_LOQ, LOQ_value))

# Ensuring column names are tidy compliant post pivot
data_1w <- clean_names(data_1w)
## correcting the unecessary space introduced in pH by clean_names
colnames(data_1w) <- str_replace(colnames(data_1w), "p_h", "ph")

# Replacing NAs in number_of_samples and number of samples_below_LOQ
subset <- which(str_detect(string = colnames(data_1w), pattern = "^number_of_samples"))
cols <- colnames(data_1w)[subset]
for(i in 1:length(cols)){
  data_1w[cols[i]] <- replace_na(data = unlist(data_1w[cols[i]]), replace = 0)
}

## COMMENT: this implementation is not elegant and was difficult due to input data format restriction (thus the use of 'unlist' as tibble columns are embedded as lists); 
## nonetheless this works and all the number_of_samples columns (14 in total) have NAs converted to 0

# Adding a column for total number of samples and total number of samples below LOQ per site & year combination as well as the percentage of below LOQ values
# Adding measurement ratio column per AT (number of measurement of AT over number of measurements of the 7 ATs for each unique combination of site and year)
data_1w <- data_1w  %>% ungroup() %>% mutate(total_number_of_samples = select(., cols[1]:cols[7]) %>% rowSums(na.rm = TRUE),
                                             total_number_of_samples_below_LOQ = select(., cols[8]:cols[14]) %>% rowSums(na.rm = TRUE),
                                             below_LOQ_ratio = total_number_of_samples_below_LOQ/total_number_of_samples,
                                             bod5_ratio = number_of_samples_bod5/total_number_of_samples,
                                             cond_ratio = number_of_samples_electrical_conductivity/total_number_of_samples,
                                             ph_ratio = number_of_samples_ph/total_number_of_samples,
                                             nitrate_ratio = number_of_samples_nitrate/total_number_of_samples,
                                             ammonium_ratio = number_of_samples_ammonium/total_number_of_samples,
                                             o2sat_ratio = number_of_samples_oxygen_saturation/total_number_of_samples,
                                             phosphate_ratio = number_of_samples_phosphate/total_number_of_samples,
                                             identifier = as.factor(str_c(site, meas_year, sep="_")))
# Removing the unnecessary variables
rm(subset, cols, i)
# Arranging the dataframe by desc(total_number_of_samples) and below_LOQ_ratio
data_1w <- data_1w %>% arrange(desc(total_number_of_samples), below_LOQ_ratio)
# Visdat of data_S1_wide
vis_dat(data_1w, warn_large_data = F)
# Visualization of the distribution of below LOQ ratio
ggplot(data=data_1w, aes(x = below_LOQ_ratio)) + geom_histogram(stat="bin", binwidth = 0.02)
# Exporting data_s1w as csv
# write.csv(data_s1w, file="WISE EU data subset for group 1 AT after pivot.csv")
# Filtering to retain only the entries for which a non 0 zero ratio for each AT is observed
data_1w <- data_1w %>% filter(bod5_ratio !=0 & cond_ratio !=0 & ph_ratio !=0 & nitrate_ratio !=0 & ammonium_ratio !=0 &o2sat_ratio !=0 & phosphate_ratio !=0)
# Visdat of data_S1_wide
vis_dat(data_1w, warn_large_data = F)
# Exporting data_s1w as csv
# write.csv(data_s1wc, file="WISE EU data subset for group 1 all measured AT after pivot.csv")

# Creating a summary table for visualization support
data_1w_summary <- data_1w %>% group_by(country, wb_type, meas_year) %>% summarize(number_of_sites=n())
# Visualizing the number of sites per year and country for which all 7 ATs in group 1 are measured
ggplot(data = data_1w_summary, aes(x=as.factor(meas_year), y=number_of_sites, fill = wb_type)) +geom_histogram(stat="identity") +facet_wrap(~country)

# Filtering data_s1wc to retain only years 2013 and 2014 for BE, EE and ES
data_1w <- data_1w %>% filter(meas_year %in% c(2013, 2014) & country %in% c("BE", "EE", "ES") & wb_type == "RW")

# Reshaping the data for k means clustering
## Creating row names based on the unique identifier obtained by aggregation of site and meas_year
# rownames(data_s1w) <- str_c(data_s1w$site, data_s1w$meas_year, sep="_")
row_names <- data_1w$identifier
## splitting into a matrix for k means by retaining only the mean and sd columns
matrix_1w <- data_1w %>%  select(starts_with("measured_value_mean") | starts_with("measured_value_sd"))
## and a dataframe for the metadata
metadata_1w <- data_1w %>%  select(-starts_with("measured_value_mean") & -starts_with("measured_value_sd"))


## NZV not needed
# Applying filters to improve matrix_s1wcf data quality
# 1. Near zero variance filtering
# freqCut <- 95/5
# uniqueCut <- 5
# matrix_1w_nzv_subsetter <- nearZeroVar(matrix_1w, freqCut = freqCut, uniqueCut = uniqueCut ,saveMetrics = F)
# matrix_1w <- matrix_1w %>% select(-all_of(matrix_1w_nzv_subsetter))

# This results in all sd columns being dropped except for nitrate
# I will manually remove measured_value_sd_nitrate and convert to a matrix for a first evaluation of the silhouette plot
# Do I scale the data?
# For the 1st analysis, the following line of code was run
# data_a1 <- scale(data_a1)
matrix_1w <- matrix_1w %>% select(-starts_with("measured_value_sd"))
# Replacing NAs by 0
# matrix_1w <- matrix_1w %>% replace_na(list(measured_value_mean_bod5 = 0,
#                                            measured_value_mean_electrical_conductivity = 0,
#                                            measured_value_mean_ph  = 0,
#                                            measured_value_mean_nitrate = 0,
#                                            measured_value_mean_ammonium = 0,
#                                            measured_value_mean_oxygen_saturation = 0,
#                                            measured_value_mean_phosphate = 0))
# Converting to a matrix and adding the row names
matrix_1w <- as.matrix(matrix_1w)
rownames(matrix_1w) <- row_names

# K means on matrix_s1wcf
## Step 1: determining the optimal number of clusters using elbow/silhouette plots
fviz_nbclust(matrix_1w, kmeans, method = "silhouette")

## Step 2: kmeans clustering with the optimal number of clusters determined in Step 1 and 100 random cluster center starting positions
km_out <- kmeans(matrix_1w, centers=2, nstart=5000)
km_clust <- km_out$cluster

## Step 3: Creating a visualization to evaluate the clustering along the first 2 PCA dimensions
### Approach 1: site average; optimal number of clusters k =
fviz_cluster(list(data=matrix_1w, cluster = km_clust))
# fviz_cluster(list(data=matrix_s1wcf, cluster = km_clust), repel = T, show.clust.cent = T, ellipse.alpha = 0.15, ggtheme = theme_classic())


matrix_1w_pca <- prcomp(matrix_1w, center = T, scale. = T)
fviz_pca_biplot(matrix_1w_pca, axes = c(1,2), repel = T)


########################################################################################################################################################################################
########################################################################################################################################################################################
########################################################################################################################################################################################






=======
                           number_of_samples_below_LOQ = as.numeric(resultQualityNumberOfSamplesBelowLOQ),
                           LOQ_value = as.numeric(procedureLOQValue))




# Quick overview of the new dataset

### 1. by country
View(data %>% group_by(country) %>% summarize(meas_per_country = n()) %>% arrange(desc(meas_per_country)))

### 2. by Water system
view(data %>% group_by(WB_system_name) %>% summarize(meas_per_WB = n()) %>% arrange(desc(meas_per_WB)))

### 3. by AT
view(data %>% group_by(AT_name) %>% summarize(meas_per_AT = n()) %>% arrange(desc(meas_per_AT)))
>>>>>>> 16d5e90d9eac54fc080fdfd83b3ab428bc088d73



#### Idea: defining AT groups for more meaning full comparisons between sites/countries/WB/years
<<<<<<< HEAD
=======
#### Group 1: usual analysis: pH, T, dissolved Oxygen, Phosphate, Nitrate and Ammonium (alternative to phosphate & nitrate would be total P/N respectively)  
AT_group1 <- c("pH", "Water Temperature", "Oxygen saturation", "Electrical conductivity", "Ammonium", "Phosphate", "Nitrate", "BOD5")
data_g1 <- data %>% filter(AT_name %in% AT_group1)
>>>>>>> 16d5e90d9eac54fc080fdfd83b3ab428bc088d73

#### Group 2: metals
AT_group2 <- c("Lead and its compounds", "Copper and its compounds", "Cadmium and its compounds", "Zinc and its compounds",
               "Nickel and its compounds", "Arsenic and its compounds", "Mercury and its compounds", "Chromium and its compounds",
               "Iron and its compounds")
data_g2 <- data %>% filter(AT_name %in% AT_group2)

<<<<<<< HEAD

# Plotting the distribution of measurements below LOQ percentage for the whole dataset
ggplot(data_g2, aes(x=below_LOQ_perc)) +geom_histogram(stat="bin")

# Creating 3 groups based on LOQ percentage distribution
## Groupe A: measurements are always Above the LOQ limit
data_g2_A  <- data_g2 %>% filter(below_LOQ_perc == 1)
## Group B: measurements are always Below the LOQ limit
data_g2_B <- data_g2 %>% filter(below_LOQ_perc == 0)
## Groupe P: measurements are Partially above/below the LOQ limit
data_g2_P <- data_g2 %>% filter(below_LOQ_perc > 0 & below_LOQ_perc < 1)

## Group NA: measurements for which the LOQ_value is NA
data_g2_NA <- data_g2 %>% filter(is.na(LOQ_value))




View(data_g2_A %>% group_by(AT_name) %>% summarize(mean = fivenum(measured_value_mean)))


## Extra trim for visualization purposes
data_g2_A <- data_g2_A %>% filter(meas_year > 1999)

ggplot(data_g2_A[data_g2_A$country == "ES",], aes(x=AT_name, y=measured_value_mean)) + geom_violin(position="dodge") + facet_wrap(~meas_year) + theme(axis.text.x = element_text(angle = 90))







## FILTER:  AT_group 1, below_LOQ_perc < 1, meas_year >= 2000
data_s1 <- data %>% filter(AT_name %in% AT_group1 & meas_year >= 2000 & below_LOQ_perc < 1)

## Creation of a unique identifier column based on the combination of site, AT_code and meas_year for unambiguous identification of entries
data_s1$identifier <- str_c(data_s1$site, data_s1$AT_code, data_s1$meas_year, sep="_")

# PROBLEM 1:  FIX THE PROBLEM OF MULTIPLE ENTRIES PER SITE + MEAS_YEAR + AT_NAME !!!
# Number of unique idenfiers in data_S1
length(unique(data_s1$identifier))
# [1] 115029

# After investigation it was identified that the cause for replicated identifier is multiple sampling depth per identifier
# Summarizing the replication count per identifier and providing the number of rows (i.e. number of unique identifiers) with a replication > 1
data_s1 %>% group_by(identifier) %>% summarize(sampling_depths = n_distinct(sampling_depth)) %>% arrange(desc(sampling_depths)) %>%
  filter(sampling_depths > 1) %>% summarize(row_total = sum(sampling_depths), number_of_rows = nrow(.))
# A tibble: 1 × 2
# row_total number_of_rows
# <int>          <int>
#   1       275             67























ggplot(data_g2_A, aes(x=measured_value_mean)) +
  geom_histogram(stat="bin", binwidth = 0.25)+ facet_wrap(~AT_name) + theme(axis.text.x = element_text(angle = 90), xlim = 5)







=======
>>>>>>> 16d5e90d9eac54fc080fdfd83b3ab428bc088d73
# TEST: retaining only Spain and Austria data (most number of measurements)
data_g2_f <- data_g2 %>% filter(country == "ES" | country == "AT")
common_years <- unique(data_g2_f[data_g2_f$country == "ES",]$meas_year)
data_g2_f <- data_g2_f %>% filter(meas_year %in% common_years)


ggplot(data = data_g2_f[data_g2_f$AT_name == "Iron and its compounds",],
       aes(x=measured_value_mean,
           fill = country)) +
  geom_histogram(position="dodge", binwidth = 5)

ggplot(data = data_g2_f[data_g2_f$AT_name == "Zinc and its compounds",],
       aes(x=measured_value_mean, col=country)) + geom_density()



## Needs rework: either reshaping of data or reworking of ggplot expr
# ggplot(data = data_g2_f, 
#        aes(x = as.factor(meas_year),
#            y = measured_value_mean,
#            color = country)) +
#   geom_boxplot(position = "dodge") +
#   facet_wrap(~AT_name, scales = "free")


data_g2_f_recap <- data_g2_f %>% group_by(country, meas_year, AT_name) %>% 
  summarize(mean = mean(measured_value_mean),
            sd = sd(measured_value_mean),
            tot_meas_number = sum(number_of_samples),
            tot_meas_number_below_LOQ = sum(number_of_samples_below_LOQ))

ggplot(data = data_g2_f_recap,
       aes(x=as.factor(meas_year), y = mean, col = country, fill = country)) + 
  geom_violin(position = "dodge") +
  facet_wrap(~AT_name, scales = "free")

#### Group 3: halogen containing compounds
AT_group3 <- c("")


data_g1_recap <- data_g1 %>% group_by(country, AT_name) %>% summarise(count = n()) %>% arrange(desc(count))
ggplot(data = data_g1_recap, aes(x = AT_name, y=count, fill = country)) + geom_histogram(stat="identity", position="stack") + theme(element_text(angle = 90))


data_g1_recap <- data_g1 %>% group_by(country, meas_year, AT_name) %>% summarise(count = n()) %>% arrange(desc(count))
ggplot(data = data_g1_recap, aes(x = AT_name, y=count, fill = country)) + geom_histogram(stat="identity", position="stack") + facet_wrap(~meas_year) + theme(axis.text.x = element_text(angle = 90))



data_g1 <- data_g1 %>% filter(meas_year > 2004 & meas_year < 2013) 

