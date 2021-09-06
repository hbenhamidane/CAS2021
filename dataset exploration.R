library(janitor)
library(plyr)
library(tidyverse)
library(stringr)
library(readxl)
library(openxlsx)
library(visdat)
library(ggplot2)
library(stats)

# EU dataset for water quality; filtered available database for Switzerland entries (country code = CH)
setwd("C:/R/Datasets/EU water data/")
wd <- getwd()
#Reading the data
data <- read.csv("DataExtract_Switzerland.csv", header = T)

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

# Empty Columns 
# "parameterSedimentDepthSampled", "parameterSpecies", "resultMoisture","resultFat","resultExtractableLipid","resultLipid"

# Target columns for analysis: 
# "resultUom", "resultObservedValue", "observedPropertyDeterminandCode","observedPropertyDeterminandLabel","procedureAnalysedMatrix"

# How many observations per body of water
data %>% group_by(monitoringSiteIdentifier) %>% count()
# How many observation per body of water and year:
data %>% group_by(monitoringSiteIdentifier, phenomenonTimeSamplingDate_year) %>% count()
# How many observation per body of water, year and target compound (observedPropertyDeterminantLabel):
data %>% group_by(monitoringSiteIdentifier, phenomenonTimeSamplingDate_year, observedPropertyDeterminandLabel) %>% count() %>% arrange(desc(n))


# How many measurements did not allow for quantification (result below LOQ)
sum(as.logical(data$resultQualityObservedValueBelowLOQ))
# [1] 87213
# Percentage of unquantifiable data (below LOQ)
uq <- sum(as.logical(data$resultQualityObservedValueBelowLOQ))/nrow(data)*100
# [1] 37.23853
# Percentage of quantifiable data 
q <- sum(!as.logical(data$resultQualityObservedValueBelowLOQ))/nrow(data)*100
# [1] 62.76147
q+uq
# [1] 100
# Confirming 100% of the data falls in either quantifiable or unquantifiable


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

