# WISE Switzerland data analysis script
# Version: {0.8}
# Author: {Hisham Ben Hamidane}
#____________________________________________________________________________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________________________________________________________________________
#
# Changes:
# 1. Addition of script metadata and version
# 2. Improvement of script comments
# 3. Removal of unnecessary packages: openxlsx, readxl
# 4. Coercing "below_LOQ" variable to logical (previously character) and grouping variables (AT_name, AT_code, WB_type, site_name, etc...) to factors
# 5. Addition of preprocessing for k-means
# 6. Elbow plot for cluster number estimation according to the within sum square method
# 7. Expansion to 44 dimensions adding the min, max and sd values to the average for each AT and site and performing of K-means
# 8.Adding a semester column (winter/summer) to account for seasonality in the measurements and try to cluster along
# 9. Reworked the overall Kmeans workflow distinguishing 4 approaches and renaming variables according to improve readability
# 10. Addition of new approach, retaining all the data, pivoting to a wider dataframe after grouping according to the time variable. 
# 11. Use of near zero variance filter (caret) to drop the variables with zero or near zero variance
# 12. Rework of non zero variance filter in the main script (rows 252 to 282) for option 1; to be used as template for options 2 to 4 and possibly write up a function
# 13. Restructuring of code and renaming of data frames pertaining to the different data grouping approaches
# 14. Addition of standard deviation per AT as an additional dimension
# 15. Implementation of supervised clustering using k nearest neighbours
# 16. Evaluation of differences between 0 and mean inputation
# 17. Reorganization of code and addition of comments
# 18. Addition of To DO, to Investigate and Problems comment sections to track analysis progress
# 19. Reorganization of general parameters for the analysis in 0.
# 20. Re-implementing knn modeling on the raw data (after pivot_wider transformation) rather than PCA transformed data
# 21. Determination of the optimal number of k nearest neighbours: k_opt
# 22. KNN model using k_opt and predicting group values: test_pred
# 23. Confusion matrix for test_y vs test_pred
# 24. Creation of dedicated df for PCA using both 0 and avg inputation
# 25. Implementation of average imputation using colMeans after reshaping but prior to nZv filtering
#
#____________________________________________________________________________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________________________________________________________________________
#
# To DO:
# - Compare the 0 inputation with the mean inputation per AT (mean across all sites, or mean across all sites of a specific category (WB_type or WB_name))
# - Recap the parameter evaluation for the nearZeroVar filtering (optimization)
# - Replicated the knn clustering analysis using WB_system_name instead of WB_type
#
# To Investigate:
# - bootstraping vs cross-validation for knn model train control
# - impact of including sd vs clustering efficiency
# - alternatives to knn for clustering
#
# # # Problems # # #
## Overview of the distribution of data between the higher level categories (Water Body system)
# table(data$WB_system_name, data$WB_type)
### uneven distribution of measurements between water basins; will have to be accounted for when sampling 
#           GW    LW    RW
# DANUBE   353     0   688
# PO      1384    26  3202
# RHINE  16147    48 74313
# RHONE   2210   477 12009
#
#____________________________________________________________________________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________________________________________________________________________
#
# Loading the required packages
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
#
#____________________________________________________________________________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________________________________________________________________________
#
# 0. Parameters for the analysis
##Defining variables for dataselection
# meas_count_threshold <- 5000
## Selected analytical targets 
curated_ATs <- c("EEA_3152-01-0", "EEA_3142-01-6", "CAS_16887-00-6", "CAS_14797-55-8")
## setting the random seed
set.seed(1234)
## Defining the train:test data splitting for supervised and unsupervised model training (0 < split_p < 1)
split_p <- 0.6
## parameters for the near zero variance filtering:
freqCut <- 80/20
uniqueCut <- 33
#
#____________________________________________________________________________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________________________________________________________________________
#
# 1. Loading the data and metadata, merging the 2 together and manipulating variables (selecting/renaming column names, adding new columns)
## EU dataset for water quality; filtered available database for Switzerland entries (country code = CH)
### path 1 (hbh private cpu)
setwd("C:/Users/shami/Documents/CAS 2021/")
### path 1 (hbh professional cpu)
# setwd("C:/Users/hisham.benhamidane/OneDrive - Thermo Fisher Scientific/Documents/R/projects/CAS2021")
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
## Creating a new factor column grouping by semester 
### winter semester months 10, 11, 12, 1, 2 and 3 --> W
### summer semester months 4, 5, 6, 7, 8 and 9 --> S
data$semester <- as.factor(if_else(data$meas_month %in% c("1", "2", "3", "10", "11", "12"), "W", "S"))
site_num <- str_extract(data$site, pattern = "[[:digit:]]{2,4}[[:alpha:]]*$")
# Creating unique and compact obs_id identifyers for rows
data$obs_id <- str_c(str_extract(data$WB_type, pattern = "^[[:alpha:]]{1}"), site_num, sep="-")
# Adding semester information to the unique obs_id identifyer
data$obs_id_S <- str_c(data$semester, str_extract(data$WB_type, pattern = "^[[:alpha:]]{1}"), site_num, sep="-")  
# Adding month information to the unique obs_id identifyer
data$obs_id_M <- str_c(data$meas_month, str_extract(data$WB_type, pattern = "^[[:alpha:]]{1}"), site_num, sep="-")
rm(site_num)
#
#____________________________________________________________________________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________________________________________________________________________
#
# 2. Creating data overview, assessing data completeness and missing values 
## Graphical representation of dataframe after edits in 1 using Visdat
# vis_dat(data, warn_large_data = F)
## Percentage of missing values (total)
percent_missing_values <- sum(is.na(data$measured_value))/nrow(data)*100
print(paste("The overall percentage of missing values is", round(percent_missing_values,1), "%", sep = " "))
## Percentage of measured values at or below the LOQ of their respective analytical methods
percent_NQ_values <- sum(data$below_LOQ)/nrow(data)*100
print(paste("The overall percentage of values at or below LOQ is", round(percent_NQ_values,1), "%", sep = " "))
## removing percent_missing_values and percent_NQ_values
rm(percent_missing_values, percent_NQ_values)
## Creating a recap table summarizing the mMissing values and values at or below LOQ per analytical target, number of sites where the AT is measured as well as corresponding percentages
data_completeness <- data %>% group_by(AT_code, semester) %>% 
                              summarize(total_values = n(),
                                        missing_values = sum(is.na(measured_value)),
                                        not_measured_values = sum(below_LOQ==T),
                                        number_of_sites = n_distinct(site),
                                        number_of_years = n_distinct(meas_year)) %>% 
                              mutate(missing_perc = missing_values/total_values*100,
                                     not_measured_perc = not_measured_values/total_values*100) %>%
                              arrange(desc(total_values), not_measured_perc, missing_perc)  
## Additional information combined between data and data_completeness
data_completeness <- data %>% select(AT_name, AT_code) %>% filter(!duplicated(AT_code)) %>% right_join(data_completeness, by = "AT_code") %>% 
                              mutate(avg_meas_per_site = total_values/number_of_sites,
                                     meas_site_ratio = number_of_sites/length(unique(data$site))*100) %>%
                              arrange(desc(total_values), number_of_sites, desc(avg_meas_per_site), not_measured_perc, missing_perc)
#Filtering data to remove missing values
data_f <- data %>% filter(!is.na(measured_value))
#
#____________________________________________________________________________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________________________________________________________________________
#
# 3. Preparing the data for k-means clustering based on the different grouping conditions we wish to use in the analysis:
## Approach 1: 1 measurement per site (all measurements for a given site ID will be averaged); missing values will be 0 filled
data_a1 <- data  %>% group_by(obs_id, AT_code) %>% mutate(avg = mean(measured_value), sd = sd(measured_value)) %>%  ungroup() %>%   
  select(obs_id, AT_code, avg, sd) %>% 
  distinct() %>%
  mutate(avg = replace_na(data = avg, 0),
         sd = replace_na(data = sd, 0)) %>% 
  pivot_wider(id_cols = obs_id, names_from = AT_code, values_from = c(avg, sd), values_fill = 0)
### checking the test df structure and content
vis_dat(data_a1)
### creating a unique row index and passing it as rownames and retaining only value columns (to retain a num matrix)
data_a1_index <- data_a1$obs_id
data_a1 <- data_a1 %>% select(-obs_id)
### Filtering and removing variables based on near zero variance
data_a1_nzv_subsetter <- nearZeroVar(data_a1, freqCut = freqCut, uniqueCut = uniqueCut ,saveMetrics = F)
data_a1 <- data_a1 %>% select(-all_of(data_a1_nzv_subsetter))
rm(data_a1_nzv_subsetter)
### Scaling (normalizing the data and converting to a matrix)
data_a1 <- scale(data_a1)
### Adding a unique row identifier
row.names(data_a1) <- data_a1_index

## Approach 2: 1 measurement per site AND per semester (to account for the seasonal variability but still limit the PCA projection space); missing values will be 0 filled
data_a2 <- data  %>% group_by(obs_id_S, AT_code) %>% mutate(avg = mean(measured_value), sd = sd(measured_value)) %>%  ungroup() %>%   
  select(obs_id_S, AT_code, avg, sd) %>% 
  distinct() %>%
  mutate(avg = replace_na(data = avg, 0), sd = replace_na(data = sd, 0)) %>% 
  pivot_wider(id_cols = obs_id_S, names_from = AT_code, values_from = c(avg, sd), values_fill = 0)
### checking the test df structure and content
vis_dat(data_a2)
### creating a unique row index and passing it as rownames and retaining only value columns (to retain a num matrix)
data_a2_index <- data_a2$obs_id_S
data_a2 <- data_a2 %>% select(-obs_id_S)
### Filtering and removing variables based on near zero variance
data_a2_nzv_subsetter <- nearZeroVar(data_a2, freqCut = freqCut, uniqueCut = uniqueCut, saveMetrics = F)
data_a2 <- data_a2 %>% select(-all_of(data_a2_nzv_subsetter))
rm(data_a2_nzv_subsetter)
### Scaling (normalizing the data and converting to a matrix)
data_a2 <- scale(data_a2)
### Adding a unique row identifier
row.names(data_a2) <- data_a2_index

## Approach 3: 1 measurement per site AND per month (to provide the maximum meaningful time resolution per site, might resulte in crowded PCA projection space); missing values will be 0 filled
data_a3 <- data  %>% group_by(obs_id_M, AT_code) %>% mutate(avg = mean(measured_value), sd = sd(measured_value)) %>%  ungroup() %>%   
  select(obs_id_M, AT_code, avg, sd) %>% 
  distinct() %>%
  mutate(avg = replace_na(data = avg, 0), sd = replace_na(data = sd, 0)) %>% 
  pivot_wider(id_cols = obs_id_M, names_from = AT_code, values_from = c(avg, sd), values_fill = 0)
### checking the test df structure and content
vis_dat(data_a3)
### creating a unique row index and passing it as rownames and retaining only value columns (to retain a num matrix)
data_a3_index <- data_a3$obs_id_M
data_a3 <- data_a3 %>% select(-obs_id_M)
### Filtering and removing variables based on near zero variance
data_a3_nzv_subsetter <- nearZeroVar(data_a3, freqCut = freqCut, uniqueCut = uniqueCut, saveMetrics = F)
data_a3 <- data_a3 %>% select(-all_of(data_a3_nzv_subsetter))
rm(data_a3_nzv_subsetter)
### Scaling (normalizing the data and converting to a matrix)
data_a3 <- scale(data_a3)
### Adding a unique row identifier
row.names(data_a3) <- data_a3_index
#
#____________________________________________________________________________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________________________________________________________________________
#
# 4. Performing K-means clustering (unsupervised learning approach)
## Step 1: determining the optimal number of clusters using elbow/silhouette plots
### Approach 1: site average
fviz_nbclust(data_a1, kmeans, method = "silhouette")
### Approach 2: site and semester average
fviz_nbclust(data_a2, kmeans, method = "silhouette")
### Approach 3: site and month average
fviz_nbclust(data_a3, kmeans, method = "silhouette")

## Step 2: kmeans clustering with the optimal number of clusters determined in Step 1 and 100 random cluster center starting positions
### Approach 1: site average; optimal number of clusters k = 2
km_out_a1 <- kmeans(data_a1, centers=2, nstart=100)
km_clust_a1 <- km_out_a1$cluster
### Approach 2: site average; optimal number of clusters k = 2
km_out_a2 <- kmeans(data_a2, centers=2, nstart=100)
km_clust_a2 <- km_out_a2$cluster
### Approach 3: site average; optimal number of clusters k = 2
km_out_a3 <- kmeans(data_a3, centers=2, nstart=100)
km_clust_a3 <- km_out_a3$cluster

## Step 3: Creating a visualization to evaluate the clustering along the first 2 PCA dimensions
### Approach 1: site average; optimal number of clusters k =
fviz_cluster(list(data=data_a1, cluster = km_clust_a1))
### Approach 1: site average; optimal number of clusters k =
fviz_cluster(list(data=data_a2, cluster = km_clust_a2))
### Approach 1: site average; optimal number of clusters k =
fviz_cluster(list(data=data_a3, cluster = km_clust_a3))
#
#____________________________________________________________________________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________________________________________________________________________
#
# 5. Performing and visualize a principal component analysis and attempting to cluster using mclust (multiple models)

# 5.1 Creating 2 distinct input dataframe with respectively 0 and avg inputation for avg(measured_value) and sd(measured value)
## 0 inputation dataset
data_0  <- data  %>% group_by(obs_id, AT_code) %>% mutate(avg = mean(measured_value, na.rm = T), sd = sd(measured_value, na.rm = T)) %>%  ungroup() %>%   
  select(obs_id, AT_code, avg, sd) %>% 
  distinct() %>%
  mutate(avg = replace_na(data = avg, 0),
         sd = replace_na(data = sd, 0)) %>% 
  pivot_wider(id_cols = obs_id, names_from = AT_code, values_from = c(avg, sd), values_fill = 0)
### checking the test df structure and content
vis_dat(data_0)
### creating a unique row index and passing it as rownames and retaining only value columns (to retain a num matrix)
data_0_index <- data_0$obs_id
data_0 <- data_0 %>% select(-obs_id)
### Filtering and removing variables based on near zero variance
data_0_nzv_subsetter <- nearZeroVar(data_0, freqCut = freqCut, uniqueCut = uniqueCut ,saveMetrics = F)
data_0 <- data_0 %>% select(-all_of(data_0_nzv_subsetter))
### Converting to a matrix class object
data_0 <- as.matrix(data_0)
### Scaling (normalizing the data and converting to a matrix)
data_0 <- scale(data_0)
### Adding a unique row identifier
row.names(data_0) <- data_0_index

## Avg inputation dataset
data_avg  <- data  %>% group_by(obs_id, AT_code) %>% mutate(avg = mean(measured_value, na.rm=T), sd = sd(measured_value, na.rm=T)) %>%  ungroup() %>%   
  select(obs_id, AT_code, avg, sd) %>% 
  distinct() %>% 
  pivot_wider(id_cols = obs_id, names_from = AT_code, values_from = c(avg, sd), values_fill = NA)
### Visualizing the extent of missing (NA) values
vis_dat(data_avg)
### creating a unique row index and passing it as rownames and retaining only value columns (to retain a num matrix)
data_avg_index <- data_avg$obs_id
data_avg <- data_avg %>% select(-obs_id)
### Creating a vector of column means for inputation
col_means <- colMeans(data_avg, na.rm=T)
### Passing the colMeans value to the NAs for each variable (replacement needs to be a list if input is a df)
data_avg <- data_avg %>% replace_na(as.list(col_means))
### Visualizing the extent of missing (NA) values AFTER average inputation but BEFORE nzV filtering
vis_dat(data_avg)
### Filtering and removing variables based on near zero variance
data_avg_nzv_subsetter <- nearZeroVar(data_avg, freqCut = freqCut, uniqueCut = uniqueCut ,saveMetrics = F)
data_avg <- data_avg %>% select(-all_of(data_avg_nzv_subsetter))
### Visualizing the extent of missing (NA) values AFTER average inputation and AFTER nzV filtering
vis_dat(data_avg)
### Converting to a matrix class object
data_avg <- as.matrix(data_avg)
### Scaling (normalizing the data and converting to a matrix)
data_avg <- scale(data_avg)
### Adding a unique row identifier
row.names(data_avg) <- data_avg_index

# 5.2 Performing the PCA on the site aggregated data for the 0 and average inputated matrices
data_0_PCA <- prcomp(data_0, center = T, scale. = T)
data_avg_PCA <- prcomp(data_avg, center = T, scale. = T)

# 5.3 Visualizing the Variance by eigenvector (Visualizationwith package factoextra
fviz_eig(X = data_0_PCA, addlabels = T, main = "Variance explained by principal component dimension for the 0 inputed data", xlab = "Principal component dimension")
fviz_eig(X = data_avg_PCA, addlabels = T, main = "Variance explained by principal component dimension for the Avg inputed data", xlab = "Principal component dimension")

# 5.4 Creating PCA biplots for 0 and Avg inputed data, projecting on PC1 and PC2
## 0 inputed PCA biplot
cat_WB_0 <- as.factor(str_extract(string = rownames(data_0_PCA$x), pattern = "^[:alpha:]{1}"))
fviz_pca_biplot(data_0_PCA, axes = c(1,2), repel = T, col.ind = cat_WB_0, title = "PCA projection along PC1 and PC2 for the 0-inputed dataset", geom = c("point"))
## Avg inputed PCA biplot
cat_WB_avg <- as.factor(str_extract(string = rownames(data_avg_PCA$x), pattern = "^[:alpha:]{1}"))
fviz_pca_biplot(data_avg_PCA, axes = c(1,2), repel = T, col.ind = cat_WB_avg, title = "PCA projection along PC1 and PC2 for the Avg-inputed dataset")
#
#____________________________________________________________________________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________________________________________________________________________
#
# 6. Supervised learning using k-nearest neighbours
## Create a dataframe containing all obs x var and additional factor columns for the WB type and WB name 
data_knn <- data  %>% group_by(obs_id, AT_code) %>% mutate(avg = mean(measured_value, na.rm=T), sd = sd(measured_value, na.rm=T)) %>%  ungroup() %>%   
  select(obs_id, AT_code, avg, sd) %>% 
  distinct() %>%
  mutate(avg = replace_na(data = avg, 0), sd = replace_na(data = sd, 0)) %>% 
  pivot_wider(id_cols = obs_id, names_from = AT_code, values_from = c(avg, sd), values_fill = 0)
### Adding the WB type and WB name columns
data_knn <- data %>% select(obs_id, WB_type, WB_system_name) %>% distinct() %>% right_join(data_knn, by="obs_id")
### Filtering and removing variables based on near zero variance
data_knn_nzv_subsetter <- nearZeroVar(data_knn, freqCut = freqCut, uniqueCut = uniqueCut, saveMetrics = F)
data_knn <- data_knn %>% select(-all_of(data_knn_nzv_subsetter))
rm(data_knn_nzv_subsetter)

## 6.1 Creating a train/test data set based on the PC1 from the PCA analysis
### Creating a subsetter to split the data between training and testing subsets using createDataPartition with a split p defined by split_p (defautl used 75:25 train:test)
train_subsetter <- createDataPartition(y=data_knn$WB_type, times = 1,  p=split_p, list=F)
## 6.2 Defining the training subset input (avg and sd of measured values for each AT) and output (factor column used to group the data: WB_type or WB_system_name)
### For train_x (training model input), keeping all the variables 
train_x <- data_knn[train_subsetter, 3:ncol(data_knn)]
### For train_y (training model correct output), the WB_type value is extracted from the row.names of the PCA matrix
train_y <- data_knn[train_subsetter, "WB_type"]
### Alternate train_y (training model correct output), the WB_type value is extracted from the row.names of the PCA matrix
# train_y <- data_knn[train_subsetter, "WB_system_name"]
## 6.3 Defining the testing subset input (measured values transformed after PCA) and correct output (factor column used to group the data: WB_type or WB_system_name)
### For test_x (training model input), keeping all the variables from the PCA analysis
test_x <- data_knn[-train_subsetter, 3:ncol(data_knn)]
### For test_y (testing model correct output), the WB_type value is extracted from the row.names of the PCA matrix
test_y <- data_knn[-train_subsetter, "WB_type"]
## 6.4 Defining the train control argument for the model. Here using the cross validation method; more investigation into bootstrapping ("boot") and cross-validation ("cv") is required
train_control <- trainControl(method="cv", number = 10)
## 6.5 Defining a range for parameter k (number of nearest neighbours taken into consideration when evaluating an unknown observation)  to test for (model optimization)
###  Note: A suggested starting point for k is sqrt(nrow(data_knn)); k should preferably be an odd number for better decision power (since knn is vote based)
k_range <- data.frame(k=seq(from=1, to=2*floor(sqrt(nrow(data_knn))), by=1))
## 6.6 Running a first model to determine the optimal number of k
knn_model_Kopt <- train(x = train_x,
                        y = train_y,
                        method="knn",
                        tuneGrid = k_range,
                        trControl = train_control)
## Optimal number of k determined by the model:
k_opt <- data.frame(k=knn_model_Kopt$bestTune$k)
## Visualizing the model performance as a function of k
ggplot(knn_model_Kopt$results, aes(x=k, y=Accuracy, color = "Red")) + geom_point() +geom_line() +ggtitle("KNN model accuracy as a function of the number of k-nearest neighbors considered")

## 6.7 Training a model based on the optimal number of k
knn_model <- train(x = train_x,
                   y = train_y,
                   method="knn",
                   tuneGrid = k_opt,
                   trControl = train_control)

## 6.8 Predicting values from test_x and comparing them with the expected output test_y
### Creating a vector of predicted values
test_pred <- predict(object = knn_model, newdata = test_x)
### Creating a dataframe model_values combining the actual and predicted values for the test subset
model_values <- data.frame(obs=test_y, pred=test_pred)
### Visualizing the confusion matrix and summary statistics
confusionMatrix(test_pred, test_y)


#
#____________________________________________________________________________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________________________________________________________________________
# Performing the same workflow but with mean instead of 0 inputed missing values
data_knn_2 <- data  %>% group_by(obs_id, AT_code) %>% mutate(avg = mean(measured_value), sd = sd(measured_value)) %>%  ungroup() %>%   
  select(obs_id, AT_code, avg, sd) %>% 
  distinct() %>%
  mutate(avg = replace_na(data = avg, mean(avg, na.rm=T)), sd = replace_na(data = sd, mean(sd))) %>% 
  pivot_wider(id_cols = obs_id, names_from = AT_code, values_from = c(avg, sd), values_fill = 0)
### Adding the WB type and WB name columns
data_knn_2 <- data %>% select(obs_id, WB_type, WB_system_name) %>% distinct() %>% right_join(data_knn_2, by="obs_id")
### Filtering and removing variables based on near zero variance
data_knn_nzv_subsetter_2 <- nearZeroVar(data_knn_2, freqCut = freqCut, uniqueCut = uniqueCut, saveMetrics = F)
data_knn_2 <- data_knn_2 %>% select(-all_of(data_knn_nzv_subsetter_2))
rm(data_knn_nzv_subsetter_2)
## 6.1 Creating a train/test data set based on the PC1 from the PCA analysis
### Creating a subsetter to split the data between training and testing subsets using createDataPartition with a split p defined by split_p (defautl used 75:25 train:test)
train_subsetter <- createDataPartition(y=data_knn_2$WB_type, times = 1,  p=split_p, list=F)

## 6.2 Defining the training subset input (avg and sd of measured values for each AT) and output (factor column used to group the data: WB_type or WB_system_name)
### For train_x (training model input), keeping all the variables 
train_x <- data_knn_2[train_subsetter, 4:ncol(data_knn_2)]
### For train_y (training model correct output), the WB_type value is extracted from the row.names of the PCA matrix
train_y <- data_knn_2[train_subsetter, "WB_type"]
### Alternate train_y (training model correct output), the WB_type value is extracted from the row.names of the PCA matrix
# train_y <- data_knn_2[train_subsetter, "WB_system_name"]

## 6.3 Defining the testing subset input (measured values transformed after PCA) and correct output (factor column used to group the data: WB_type or WB_system_name)
### For test_x (training model input), keeping all the variables from the PCA analysis
test_x <- data_knn_2[-train_subsetter, 4:ncol(data_knn_2)]
### For test_y (testing model correct output), the WB_type value is extracted from the row.names of the PCA matrix
test_y <- data_knn_2[-train_subsetter, "WB_type"]

## 6.6 Running a first model to determine the optimal number of k
knn_model_Kopt_2 <- train(x = train_x,
                        y = train_y,
                        method="knn",
                        tuneGrid = k_range,
                        trControl = train_control)
## Optimal number of k determined by the model:
k_opt_2 <- data.frame(k=knn_model_Kopt_2$bestTune$k)
## Visualizing the model performance as a function of k
ggplot(knn_model_Kopt$results, aes(x=k, y=Accuracy, color = "Red")) + geom_point() +geom_line()

## 6.7 Training a model based on the optimal number of k
knn_model <- train(x = train_x,
                   y = train_y,
                   method="knn",
                   tuneGrid = k_opt,
                   trControl = train_control)

## 6.8 Predicting values from test_x and comparing them with the expected output test_y
### Creating a vector of predicted values
test_pred <- predict(object = knn_model, newdata = test_x)
### Creating a dataframe model_values combining the actual and predicted values for the test subset
model_values <- data.frame(obs=test_y, pred=test_pred)
### Visualizing the confusion matrix and summary statistics
confusionMatrix(test_pred, test_y)
















#
#____________________________________________________________________________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________________________________________________________________________
#____________________________________________________________________________________________________________________________________________________________________________________________________________________
## REWORK AND REORGANISATION NEEDED BEYOND THIS LINE


## 1st attempt at model based clustering
data_a1_mclust <- Mclust(data_a1,G = 1:12)
fviz_mclust_bic(data_a1_mclust)
fviz_mclust(data_a1_mclust, "classification", "point")
## optimal number of clusters found to be 1; further look into BIC required to understand this result



## Creating the training set: X and Y 
data_knn_train_X <- data_knn %>% select(all_of(starts_with("avg"))) %>% slice(train_subsetter) %>% data.frame()
# TRy to rework the code to also select sd columns & all_of(starts_with("sd"))
data_knn_train_Y <- data_knn %>% select(WB_type) %>% slice(train_subsetter) %>% data.frame()
## Creating the testing set: X and Y 
data_knn_test_X <- data_knn %>% select(all_of(starts_with("avg"))) %>% slice(-train_subsetter) %>% data.frame()
data_knn_test_Y <- data_knn %>% select(WB_type) %>% slice(-train_subsetter) %>% data.frame()

## Defining the train_control argument
train_control <- trainControl(method="cv", number = 10)

k_lots <- data.frame(k=seq(from=1, to=500, by=5))

set.seed(1234)

knn_reg_cv_10 <- train(x = data_knn_train_X,
                       y = data_knn_train_Y$WB_type,
                       method="knn",
                       tuneGrid = k_lots,
                       trControl = train_control)

ggplot(knn_reg_cv_10$results, aes(x=k, y=Accuracy, color = "Red")) + geom_point()


test <- prcomp(data_knn, center = T, scale = T)



#_________________________________________________________
# Supervised learning using k nearest neighbours using knn3 from the caret package
# Defining training and testing data subsets 
# Clustering attempt on water body

# Data preparation according to the previous steps defined in 3
data_knn <- data  %>% group_by(obs_id, AT_code) %>% mutate(avg = mean(measured_value), 
                                                          sd = sd(measured_value),
                                                          min = min(measured_value),
                                                          max = max(measured_value)) %>%  ungroup() %>%   
  select(obs_id, AT_code, avg, sd, min, max) %>% 
  distinct() %>%
  mutate(avg = replace_na(data = avg, 0),
         sd = replace_na(data = sd, 0),
         min = replace_na(data = min, 0),
         max = replace_na(data = sd, 0)) %>% 
  pivot_wider(id_cols = obs_id, names_from = AT_code, values_from = c(avg, sd, min, max), values_fill = 0)

# Default near zero var parameters
data_knn_nzv_subsetter <- nearZeroVar(data_knn, freqCut = 95/5, uniqueCut = 10 ,saveMetrics = F)
# Creating a table with the group information: water body type and name
data_knn_clusterID <- data %>% select(obs_id, WB_system_name, WB_type)
# applying nZv filter and adding WB identifiers
data_knn <- data_knn %>% select(-all_of(data_knn_nzv_subsetter)) %>% left_join(data_knn_clusterID, by = "obs_id")
# scaling the data
data_knn <- scale(data_knn)
data_knn_index <- data_knn$obs_id
data_knn <- data_knn %>% select(-obs_id)
row.names(data_knn) <- data_knn_index
# Creating a seed for the random number generator
set.seed(123)
# Splitting the data set between train and test with a 80:20 ratio
size <- floor(0.8*nrow(data_knn))
sample(x = seq_len(nrow(data)), size = size)



#______________________________________________________________________________________________________________________________________________________________________________________________________________________
#______________________________________________________________________________________________________________________________________________________________________________________________________________________
#______________________________________________________________________________________________________________________________________________________________________________________________________________________

# 5. Creating aggregated statistical values beyond the mean for the different grouping conditions
## For a first clustering attempt, values per site and AT will be averaged (removing the time dimension for now)
## Step 1: calculate the min, max, median, average and standard deviation value as well as the total count of measurement for each combination of site and AT
data_f <- data_f %>% group_by(AT_code, obs_id) %>% 
                   mutate(median_measured_value = median(measured_value),
                          avg_measured_value = mean(measured_value),
                          sd_measured_value = sd(measured_value),
                          min_measured_value = min(measured_value),
                          max_measured_value = max(measured_value),
                          count_measured_value = n_distinct(measured_value),
                          alt_count_measured_value = n()) %>%
                   select(obs_id, site, WB_name, WB_system_name, WB_type, AT_code, semester, min_measured_value, avg_measured_value, max_measured_value, sd_measured_value, count_measured_value, alt_count_measured_value) %>% 
                   arrange(obs_id, AT_code)
## Remove duplicated entries (each row is still a single measurements but we are only interested in the average measured value per site and AT)
data_f <- data_f[!duplicated(data_f),]
## Getting a better feel for the data distribution
data_f_10meas <- data_f %>% filter(alt_count_measured_value >= 10)


# 4. Reshaping and selecting the data for Kmeans

# Approach 1: Further filtering to retain only sites that have a value for all 11 top analytical targets
data_f_recap <- data_f %>% group_by(obs_id) %>% summarize(num_AT_measured = n()) %>% filter(num_AT_measured == 11)

data_f_recap_2 <- data_f %>% group_by(obs_id) %>% filter(alt_count_measured_value >= 10) %>% summarize(num_AT_measured = n()) 
data_f_recap_3 <- data_f_conservative %>% group_by(obs_id) %>% filter(alt_count_measured_value >= 10) %>% summarize(num_AT_measured = n())

## Filtering to retain only the sites with the top 11 AT measured (results in a biais with only RW sites)
data_f <- data_f %>% filter(obs_id %in% data_f_recap$obs_id)
## Method 1: Conversion to a dataframe where each observation is a site and each variable is an AT retaining only the average measured value per site
### Data dimension: 11 dim, 1 for each AT measured, returning the average measured value
data_1 <- data_f %>% select(obs_id, AT_code, avg_measured_value) %>% pivot_wider(id_cols = obs_id, names_from = AT_code, values_from = avg_measured_value)
### Passing the site code as row name and removing the site column to retain only a numerical matrix
data_1_rownames <- as.character(data_1$obs_id)
data_1 <- data_1 %>% drop_na() %>% select(-obs_id)
row.names(data_1) <- data_1_rownames
## Scaling (normalizing the data and converting to a matrix)
data_1_scaled <- scale(data_1)
## Creating the distance matrix
data_1_dist <- dist(data_1_scaled)
## Determining the number of clusters using an elbow plot
### synthax from factoextra
fviz_nbclust(data_1_scaled, kmeans, method = "wss")
# 4.Performing kmeans clustering using k=2 and nstart=100
km_out_1 <- kmeans(data_1_scaled, centers=2, nstart=100)
print(km_out_1)
## Visualizing
km_clusters_1 <- km_out_1$cluster
fviz_cluster(list(data=data_1_scaled, cluster = km_clusters_1))

# Approach 2: Conversion to a dataframe where each observation is either the min, avg, max or sd of a specific AT measured at a given  site
### Data dimension: Expanding from 11 dimensions from approach 1 to 44 dimensions (4 per AT: min, avg, max and sd) 
### Converting to a long dataframe to have all the computed values (avg, min, max, sd) in 2 columns, 1 for the values the other for the key 
data_2 <- data_f %>% select(obs_id, AT_code, min_measured_value, avg_measured_value, max_measured_value, sd_measured_value) %>%
  pivot_longer(cols = c(min_measured_value, avg_measured_value, max_measured_value, sd_measured_value), values_to = "value", names_to = "type")
### Converting back to a wide dataframe, combining the AT code and the type (which of the aggregated measured value: min, avg, max and sd) to have a unique observation per case
data_2 <- data_2 %>% pivot_wider(id_cols = obs_id, names_from = c(AT_code, type), names_sep = "-")
### Passing the site code as row name and removing the site column to retain only a numerical matrix
data_2_rownames <- as.character(data_2$obs_id)
data_2 <- data_2 %>% drop_na() %>% select(-obs_id)
rownames(data_2) <- data_2_rownames
## Scaling (normalizing the data and converting to a matrix)
data_2_scaled <- scale(data_2)
## Creating the distance matrix
data_2_dist <- dist(data_2_scaled)
## Determining the number of clusters using an elbow plot
### synthax from factoextra
fviz_nbclust(data_2_scaled, kmeans, method = "wss")
# 4.Performing kmeans clustering using k=3 and nstart=100
km_out_2 <- kmeans(data_2_scaled, centers=3, nstart=100)
print(km_out_2)
## Visualizing
km_clusters_2 <- km_out_2$cluster
fviz_cluster(list(data=data_2_scaled, cluster = km_clusters_2))

## Approach 3: curating top ATs, dropping all ATs with significant (>2%) missing OR measured at/below LOQ values
## Retaining only pH (EEA_3152-01-0), conductivity (EEA_3142-01-6), chloride (CAS_16887-00-6), nitrate (CAS_14797-55-8) defined in step 0.
data_3 <- data_f %>% filter(AT_code %in% curated_ATs)
# Converting to a wider df with AT as columns
data_3 <- data_3 %>% select(obs_id, AT_code, avg_measured_value) %>% pivot_wider(id_cols = obs_id, names_from = AT_code, values_from = avg_measured_value)
### Passing the site code as row name and removing the site column to retain only a numerical matrix
data_3_rownames <- as.character(data_3$obs_id)
data_3 <- data_3 %>% drop_na() %>% select(-obs_id)
row.names(data_3) <- data_3_rownames
## Scaling (normalizing the data and converting to a matrix)
data_3_scaled <- scale(data_3)
## Creating the distance matrix
data_3_dist <- dist(data_3_scaled)
## Determining the number of clusters using an elbow plot
### synthax from factoextra
fviz_nbclust(data_3_scaled, kmeans, method = "wss")
# 4.Performing kmeans clustering using k=2 and nstart=100
km_out_3 <- kmeans(data_3_scaled, centers=2, nstart=100)
print(km_out_3)
## Visualizing
km_clusters_3 <- km_out_3$cluster
fviz_cluster(list(data=data_3_scaled, cluster = km_clusters_3))

## Approach 4: curating top ATs, dropping all ATs with significant (>2%) missing OR measured at/below LOQ values
## Retaining only pH (EEA_3152-01-0), conductivity (EEA_3142-01-6), chloride (CAS_16887-00-6), nitrate (CAS_14797-55-8)
data_4 <- data_f %>% filter(AT_code %in% curated_ATs)
# Converting to a wider df with AT as columns
data_4 <- data_4 %>% select(obs_id, AT_code, min_measured_value, avg_measured_value, max_measured_value, sd_measured_value) %>%
  pivot_longer(cols = c(min_measured_value, avg_measured_value, max_measured_value, sd_measured_value), values_to = "value", names_to = "type")
### Converting back to a wide dataframe, combining the AT code and the type (which of the aggregated measured value: min, avg, max and sd) to have a unique observation per case
data_4 <- data_4 %>% pivot_wider(id_cols = obs_id, names_from = c(AT_code, type), names_sep = "-")
### Passing the site code as row name and removing the site column to retain only a numerical matrix
data_4_rownames <- as.character(data_4$obs_id)
data_4 <- data_4 %>% drop_na() %>% select(-obs_id)
row.names(data_4) <- data_4_rownames
## Scaling (normalizing the data and converting to a matrix)
data_4_scaled <- scale(data_4)
## Creating the distance matrix
data_4_dist <- dist(data_4_scaled)
## Determining the number of clusters using an elbow plot
### synthax from factoextra
fviz_nbclust(data_4_scaled, kmeans, method = "wss")
# 4.Performing kmeans clustering using k=2 and nstart=100
km_out_4 <- kmeans(data_4_scaled, centers=2, nstart=100)
print(km_out_4)
## Visualizing
km_clusters_4 <- km_out_4$cluster
fviz_cluster(list(data=data_4_scaled, cluster = km_clusters_4))


#_________________________________________________________________________________________________________________________________
#_________________________________________________________________________________________________________________________________









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
