library(janitor)
library(plyr)
library(tidyverse)
library(stringr)
library(readxl)
library(openxlsx)
library(visdat)

setwd("C:/R/Datasets/")
wd <- getwd()

water_data <- read_csv("Waterbase_v2020_1_T_WISE6_AggregatedData") 