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
setwd("C:/Users/hisham.benhamidane/OneDrive - Thermo Fisher Scientific/Documents/R/projects/CAS2021/Raw data/20220622")
wd <- getwd()
na_string <- c("", " ", NA, "NA", "N/A")

##Reading the data and metadata

data <- read.csv("Waterbase_v2020_1_T_WISE6_AggregatedData.csv", header = T, na.strings = na_string)
colnames(data)[1] <- "monitoringSiteIdentifier"
metadata <- read.csv("Waterbase_v2020_1_S_WISE6_SpatialObject_DerivedData.csv", header = T, na.strings = na_string)
colnames(metadata)[1] <- "countrycode"

## Selecting the usefull columns of the metadata to add to data
meta_cols_to_keep <- c("countrycode", "monitoringSiteIdentifier","monitoringSiteName", "waterBodyIdentifier", "waterBodyIdentifierScheme",
                       "waterBodyName", "specialisedZoneType", "subUnitName",
                       "rbdIdentifier", "rbdName", "lon", "lat")
metadata <- metadata %>% select(all_of(meta_cols_to_keep)) %>% filter(!is.na(monitoringSiteIdentifier))

##Dropping the useless columns in the data to make it lighter
data_cols_to_drop <- c( "remarks","metadata_versionId","metadata_beginLifeSpanVersion","metadata_statusCode","metadata_observationStatus",
                        "metadata_statements", "UID"   )
##And left joining data to metadata
data <- data %>% select(-all_of(data_cols_to_drop)) %>% filter(!is.na(monitoringSiteIdentifier))
# Creating a vector of unique monitoringSiteIdentifier in data
dSID <- unique(data$monitoringSiteIdentifier)
# Checking for the % of unique(monitoringSiteIdentifier) from data in metadata
sum(dSID %in% unique(metadata$monitoringSiteIdentifier))/length(dSID)*100
# [1] 61.56156
# Creating a corresponding index for the site identifiers included in metadata
dSID_index <- which(dSID %in% unique(metadata$monitoringSiteIdentifier))
dSID2 <- dSID[dSID_index]
# length(dSID)
# length(dSID2)
# Reducing data to keep only sites id found in metadata as well
data <- data %>% filter(monitoringSiteIdentifier %in% dSID2)

# Merging data and metadata
data <- data %>% left_join(metadata, by = "monitoringSiteIdentifier" )
# removing all non-necessary variables
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
#Filtering data to remove missing values
data_f <- data %>% filter(!is.na(measured_value))