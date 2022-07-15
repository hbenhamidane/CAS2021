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


# Defining work directory and general variables
setwd("C:/Users/hisham.benhamidane/OneDrive - Thermo Fisher Scientific/Documents/R/projects/CAS2021/Raw data/2021")
wd <- getwd()
na_string <- c("", " ", NA, "NA", "N/A")


# Reading and cleaning up data
data <- read.csv("Waterbase_v2021_1_T_WISE6_AggregatedData.csv", header = T, na.strings = na_string)
## renaming the 1st column (special char issue)
colnames(data)[1] <- "monitoringSiteIdentifier"
## dropping the useless columns in the data to make it lighter
data_cols_to_drop <- c( "remarks","metadata_versionId","metadata_beginLifeSpanVersion","metadata_statusCode","metadata_observationStatus",
                        "metadata_statements", "UID"   )
## and removing rows without site id
data <- data %>% select(-all_of(data_cols_to_drop)) %>% filter(!is.na(monitoringSiteIdentifier))


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


# Left joining data to metadata
data <- data %>% left_join(metadata, by = "monitoringSiteIdentifier" )

# Quality checks post merge
## Checking for the % of unique(monitoringSiteIdentifier) from data in metadata
sum(unique(data$monitoringSiteIdentifier) %in% unique(metadata$monitoringSiteIdentifier))/length(unique(data$monitoringSiteIdentifier))*100
# [1] 100


# removing all non-necessary variables
rm(meta_cols_to_keep, data_cols_to_drop, metadata, metadata_dup)


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
                           number_of_samples_below_LOQ = as.numeric(resultQualityNumberOfSamplesBelowLOQ),
                           LOQ_value = as.numeric(procedureLOQValue))




# Quick overview of the new dataset

### 1. by country
View(data %>% group_by(country) %>% summarize(meas_per_country = n()) %>% arrange(desc(meas_per_country)))

### 2. by Water system
view(data %>% group_by(WB_system_name) %>% summarize(meas_per_WB = n()) %>% arrange(desc(meas_per_WB)))

### 3. by AT
view(data %>% group_by(AT_name) %>% summarize(meas_per_AT = n()) %>% arrange(desc(meas_per_AT)))



#### Idea: defining AT groups for more meaning full comparisons between sites/countries/WB/years
#### Group 1: usual analysis: pH, T, dissolved Oxygen, Phosphate, Nitrate and Ammonium (alternative to phosphate & nitrate would be total P/N respectively)  
AT_group1 <- c("pH", "Water Temperature", "Oxygen saturation", "Electrical conductivity", "Ammonium", "Phosphate", "Nitrate", "BOD5")
data_g1 <- data %>% filter(AT_name %in% AT_group1)

#### Group 2: metals
AT_group2 <- c("Lead and its compounds", "Copper and its compounds", "Cadmium and its compounds", "Zinc and its compounds",
               "Nickel and its compounds", "Arsenic and its compounds", "Mercury and its compounds", "Chromium and its compounds",
               "Iron and its compounds")
data_g2 <- data %>% filter(AT_name %in% AT_group2)

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

