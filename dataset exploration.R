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
# Visdat on dataframe prior to any modification
vis_dat(data, warn_large_data = F)
# Converting the full date column (phenomenonTimeSamplingDate) from character into a date object
data$date <- as.Date(x = data$phenomenonTimeSamplingDate, format ="%d/%m/%Y %H:%M:%S") 
# Creating a subset of data for quantified analytical targets
data_q <- data %>% filter(resultQualityObservedValueBelowLOQ == "False")

# Statistical overview of the measured values per analytical target
recap_table <- data %>% filter(resultQualityObservedValueBelowLOQ == "False") %>%
  group_by(observedPropertyDeterminandLabel) %>%
  summarize(min = min(resultObservedValue, na.rm=T), 
            max = max(resultObservedValue, na.rm=T),
            mean = mean(resultObservedValue, na.rm=T),
            sd = sd(resultObservedValue, na.rm=T),
            sample_size = n()) %>%
  arrange(desc(sample_size))

# Statistical overview of the measured values per analytical target
recap_table_WB <- data %>% filter(resultQualityObservedValueBelowLOQ == "False") %>%
  group_by(parameterWaterBodyCategory) %>%
  summarize(min = min(resultObservedValue, na.rm=T), 
            max = max(resultObservedValue, na.rm=T),
            mean = mean(resultObservedValue, na.rm=T),
            sd = sd(resultObservedValue, na.rm=T),
            sample_size = n()) %>%
  arrange(desc(sample_size))


# Creating a vector of names for the most sampled analytical targets
sampling_threshold <- 5000
analytical_targets <- recap_table %>% filter(sample_size > sampling_threshold)  %>% select(observedPropertyDeterminandLabel)
analytical_targets <- as.vector(analytical_targets$observedPropertyDeterminandLabel)
# Subsetting the dataset for those analytical target and filtering out the values at or below LOQ
data_topAT <- data_q %>% filter(observedPropertyDeterminandLabel %in% analytical_targets)
# Assessing the normality of data per analytical target across all measurements (site & year)
# ggplot(data=data_topAT
#        aes())

# Special case of Nitrate
data_no3 <- data_q[data_q$observedPropertyDeterminandLabel == "Nitrate",] 
data_no3$sim_norm_distr <- rnorm(n = recap_table[recap_table$observedPropertyDeterminandLabel=="Nitrate",]$sample_size,
                                 mean = recap_table[recap_table$observedPropertyDeterminandLabel=="Nitrate",]$mean,
                                 sd = recap_table[recap_table$observedPropertyDeterminandLabel=="Nitrate",]$sd) 

ggplot(data=data_no3, aes(x=resultObservedValue)) + geom_histogram(stat="density")
ggplot(data=data_no3, aes(x=resultObservedValue)) + geom_histogram(stat="density", position ="dodge") +facet_wrap(~parameterWaterBodyCategory) +xlim(0, 50)

qqplot(x=data_no3$resultObservedValue, y=data_no3$sim_norm_distr)
# Overall distribution of NO3 measurement across all years and site is not normally distributed

# Looking at the distribution of NO3 measurements per site and year
data_no3 %>% group_by(monitoringSiteIdentifier,phenomenonTimeSamplingDate_year) %>% summarise(meas_count=n()) %>% arrange(desc(meas_count))
# Site CHNTG15 has between 4 and 16 measurements per year over a 15-year period (total number of observations = 159)
data_no3_CHNTG15 <- data_no3 %>% filter(monitoringSiteIdentifier == "CHNTG15")
data_no3_CHNTG15 %>% group_by(phenomenonTimeSamplingDate_year) %>% summarise(meas_count=n()) %>% arrange(desc(meas_count))
ggplot(data=data_no3_CHNTG15, aes(x=as.factor(phenomenonTimeSamplingDate_year), y=resultObservedValue)) + geom_boxplot()

# Let's check the normality of pH measurements distribution over all sites and years (as it is the 2nd most sampled analytical target after Nitrate)
data_pH <- data_q[data_q$observedPropertyDeterminandLabel == "pH",] 
data_pH$sim_norm_distr <- rnorm(n = recap_table[recap_table$observedPropertyDeterminandLabel=="pH",]$sample_size,
                                 mean = recap_table[recap_table$observedPropertyDeterminandLabel=="pH",]$mean,
                                 sd = recap_table[recap_table$observedPropertyDeterminandLabel=="pH",]$sd) 
# Plotting pH value distribution and qqplot for all sites and years
ggplot(data=data_pH, aes(x=resultObservedValue)) + geom_histogram(stat="density") + facet_wrap(~parameterWaterBodyCategory)
qqplot(x=data_pH$resultObservedValue, y=data_pH$sim_norm_distr)
# Plotting the pH value distributions per type of Water Body
ggplot(data=data_pH, aes(x=resultObservedValue)) + geom_histogram(stat="density") + facet_wrap(~parameterWaterBodyCategory)
# Creating a separate theoretical normal distribution for each group (i.e. each body of water category)

recap_pH_per_BoWtype <- data_pH %>% group_by(parameterWaterBodyCategory) %>%
  summarize(min = min(resultObservedValue, na.rm=T), 
            max = max(resultObservedValue, na.rm=T),
            mean = mean(resultObservedValue, na.rm=T),
            sd = sd(resultObservedValue, na.rm=T),
            sample_size = n()) %>%
  arrange(desc(sample_size))

# QQ Plot for Ground Water pH measurements
qqplot(x = data_pH[data_pH$parameterWaterBodyCategory == "GW",]$resultObservedValue,
       y = rnorm(n = recap_pH_per_BoWtype[recap_pH_per_BoWtype$parameterWaterBodyCategory=="GW",]$sample_size,
                 mean = recap_pH_per_BoWtype[recap_pH_per_BoWtype$parameterWaterBodyCategory=="GW",]$mean,
                 sd = recap_pH_per_BoWtype[recap_pH_per_BoWtype$parameterWaterBodyCategory=="GW",]$sd)
)
## QQ Plot for River Water pH measurements
qqplot(x = data_pH[data_pH$parameterWaterBodyCategory == "RW",]$resultObservedValue,
       y = rnorm(n = recap_pH_per_BoWtype[recap_pH_per_BoWtype$parameterWaterBodyCategory=="RW",]$sample_size,
                 mean = recap_pH_per_BoWtype[recap_pH_per_BoWtype$parameterWaterBodyCategory=="RW",]$mean,
                 sd = recap_pH_per_BoWtype[recap_pH_per_BoWtype$parameterWaterBodyCategory=="RW",]$sd)
)


obs_val <- "Electrical conductivity"
data_cond <- data_q[data_q$observedPropertyDeterminandLabel == obs_val,] 
data_cond$sim_norm_distr <- rnorm(n = recap_table[recap_table$observedPropertyDeterminandLabel==obs_val,]$sample_size,
                                 mean = recap_table[recap_table$observedPropertyDeterminandLabel==obs_val,]$mean,
                                 sd = recap_table[recap_table$observedPropertyDeterminandLabel==obs_val,]$sd) 

ggplot(data=data_cond, aes(x=resultObservedValue)) + geom_histogram(stat="density")
ggplot(data=data_cond, aes(x=resultObservedValue)) + geom_histogram(stat="density", position ="dodge") +facet_grid(parameterWaterBodyCategory~rbdName)

qqplot(x=data_cond$resultObservedValue, y=data_no3$sim_norm_distr)




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

