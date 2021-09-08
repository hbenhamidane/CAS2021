library(janitor)
library(plyr)
library(tidyverse)
library(stringr)
library(readxl)
library(openxlsx)
library(visdat)
library(ggplot2)
library(stats)
library(lubridate)

# EU dataset for water quality; filtered available database for Switzerland entries (country code = CH)
setwd("C:/R/CAS/EU water data/")
wd <- getwd()
#Reading the data and metadata
data <- read.csv("DataExtract_Switzerland.csv", header = T)
metadata <- read.csv("WISE_spatialobject_CH.csv", header = T)

# Selecting the usefull columns of the metadata to add to data
meta_cols_to_keep <- c("monitoringSiteIdentifier","monitoringSiteName", "waterBodyIdentifier", "waterBodyIdentifierScheme",
                       "waterBodyName", "specialisedZoneType", "subUnitName",
                       "rbdIdentifier", "rbdName", "lon", "lat")
metadata <- metadata %>% select(all_of(meta_cols_to_keep))

# Dropping the useless columns in the data to make it lighter
data_cols_to_drop <- c("parameterSedimentDepthSampled","parameterSpecies",
                       "resultMoisture", "resultFat", "resultExtractableLipid",
                       "resultLipid")
# And left joining data to metadata
data <- data %>% select(-all_of(data_cols_to_drop)) %>% left_join(metadata, by = "monitoringSiteIdentifier")
rm(meta_cols_to_keep, data_cols_to_drop, metadata)

## Preliminary questions about the data
#How many uniquely monitored locations
length(unique(data$monitoringSiteIdentifier))
# [1] 186
# How many observed properties
length(unique(data$observedPropertyDeterminandCode))
# [1] 130
# What is the distribution of the different types of water bodies?
table(data$parameterWaterBodyCategory)
# GW     LW     RW 
# 113274   1380 119547 

# How complete is the dataset (missing/incomplete data)?
vis_dat(data, warn_large_data = F)

# Target columns for analysis: 
# "resultUom", "resultObservedValue", "observedPropertyDeterminandCode","observedPropertyDeterminandLabel","procedureAnalysedMatrix"
# How many observations per monitoring site
data %>% group_by(monitoringSiteIdentifier) %>% count() %>% arrange(desc(n))
# How many observation per monitoring site and year:
data %>% group_by(monitoringSiteIdentifier, phenomenonTimeSamplingDate_year) %>% count()%>% arrange(desc(n))
# How many observation per monitoring site, year and target compound (observedPropertyDeterminantLabel):
data %>% group_by(monitoringSiteIdentifier, phenomenonTimeSamplingDate_year, observedPropertyDeterminandLabel) %>% count() %>% arrange(desc(n))

# Alternatively for comparison, how many measurements per year across all BoW:
meas_per_year <- data %>% group_by(phenomenonTimeSamplingDate_year) %>% count() %>% arrange(desc(n))
ggplot(meas_per_year, aes(x=phenomenonTimeSamplingDate_year, y=n)) + geom_bar(stat="identity")
# For each target compound, how many measurements
meas_per_target <- data %>% group_by(observedPropertyDeterminandLabel) %>% count() %>% arrange(desc(n))
ggplot(meas_per_target[1:10,], aes(x=as.factor(observedPropertyDeterminandLabel), y=n)) + 
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90))

# For each target, how many quantifiable measurements (observed value > LOQ value)
meas_per_quant_target <- data %>% filter(resultQualityObservedValueBelowLOQ == "False") %>%
  group_by(observedPropertyDeterminandLabel) %>%
  count() %>% arrange(desc(n))
ggplot(meas_per_quant_target[1:10,], aes(x=as.factor(observedPropertyDeterminandLabel), y=n)) + 
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90))

# Statistical overview of the measured values per analytical target
recap_table <- data %>% filter(resultQualityObservedValueBelowLOQ == "False") %>%
         group_by(observedPropertyDeterminandLabel) %>%
         summarize(min = min(resultObservedValue, na.rm=T), 
                   max = max(resultObservedValue, na.rm=T),
                   mean = mean(resultObservedValue, na.rm=T),
                   sd = sd(resultObservedValue, na.rm=T))

# Creating a vector of names for the most sampled analytical targets
analytical_targets <- meas_per_quant_target[1:10,]$observedPropertyDeterminandLabel
# Subsetting the dataset for those analytical target and filtering out the values at or below LOQ
data_top10_analyticaltargets <- data %>% filter(observedPropertyDeterminandLabel %in% analytical_targets & resultQualityObservedValueBelowLOQ == "False")
# Parsing the full date column (phenomenonTimeSamplingDate) and expressing as POSIXCT object (lubridate package)
data_top10_analyticaltargets$date <- as_date(data_top10_analyticaltargets$phenomenonTimeSamplingDate)
# won't work, dmy format needs to be explicitely specified, look into lubridate vignette






    
# Taking a peak at a water body of interest
View(data %>% filter(monitoringSiteIdentifier == "CHRW-2078") %>% group_by(phenomenonTimeSamplingDate_year) %>% count(observedPropertyDeterminandLabel))
CHRW2078 <- data %>% filter(monitoringSiteIdentifier == "CHRW-2078")
# What kind of water body is CHRW2078
unique(CHRW2078$parameterWaterBodyCategory)
#RW : river
# Keeping only the quantifiable targets and counting unique targets
CHRW2078_q <- CHRW2078 %>% filter(!as.logical(CHRW2078$resultQualityObservedValueBelowLOQ))
unique(CHRW2078_q$observedPropertyDeterminandLabel)

target_sample <- c("pH", "Dissolved oxygen", "Chlorine Cl-","Dissolved organic carbon (DOC)","Total phosphorus","Total nitrogen")
CHRW2078_q_t <- CHRW2078_q %>% filter(observedPropertyDeterminandLabel %in% target_sample)
ggplot(CHRW2078_q_t, aes(x=phenomenonTimeSamplingDate_year, y=resultObservedValue)) + geom_line() + facet_wrap(~observedPropertyDeterminandLabel)

pH <- data %>% filter(observedPropertyDeterminandLabel == "pH" & !as.logical(resultQualityObservedValueBelowLOQ))
sites_pH <- pH %>% group_by(monitoringSiteIdentifier) %>% count() %>% arrange(desc(n))
sites <- sites_pH$monitoringSiteIdentifier[1:6]

ggplot(data=pH[pH$monitoringSiteIdentifier %in% sites,], aes(x=as.factor(phenomenonTimeSamplingDate_year), y=resultObservedValue)) +geom_boxplot() + facet_wrap(~monitoringSiteIdentifier)

