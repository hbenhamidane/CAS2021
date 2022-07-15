"""
Created on 28 June 2022

@author = "Ludovic Le Reste"
@credits = ["Hisham Ben Hamidane", "Ludovic Le Reste"]
@status = "Prototype"

This scripts loads and analyses the Waterbase WISE 6 database from the European Environment Agency.
Analysis is for the whole dataset (not limited to Switzerland as previously)
Data sources:
2021 WISE dataset (Measurements, disaggregated data, spatial data):
    https://www.eea.europa.eu/data-and-maps/data/waterbase-water-quality-icm-2

TO DO:
- equivalent of visdat
    + looks like pandas-profiling is a good candidate (run first on a small dataset)
    + perhaps in conjunction with great_expectations?
- number of point per site per year per molecule 
- check which parameters are common to all sites or to all waterbody types
- merge df with spatial (only country and river basin info) for grouping
    + missing and duplicate entries to be managed
    
TO DO (great to have/explore):    
    - rewrite script using functions and the workflow will follow from it
    - use disaggregated data to spot best candidates for time trace analysis
    - plot time traces
    - find interesting questions to answer
        + is there a pattern in terms of metal pollution?
        [10:27] Ben Hamidane, Hisham
Groupe 1: cibles analytiques générales

AT_group1 <- c("pH", "Water Temperature", "Oxygen saturation", "Electrical conductivity", "Ammonium", "Phosphate", "Nitrate", "BOD5")

[10:27] Ben Hamidane, Hisham
Groupe 2: métaux

AT_group2 <- c("Lead and its compounds", "Copper and its compounds", "Cadmium and its compounds", "Zinc and its compounds",
               "Nickel and its compounds", "Arsenic and its compounds", "Mercury and its compounds", "Chromium and its compounds",
               "Iron and its compounds")



"""
__author__ = "Ludovic Le Reste"
__credits__ = ["Hisham Ben Hamidane", "Ludovic Le Reste"]
__status__ = "Prototype"

  
# %% PACKAGES
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"
import pandas as pd
import os
import sys
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
# from mpl_toolkits.basemap import Basemap
import geopandas as gpd
# import plot_style
import statsmodels.api as sm
from scipy import stats
from sklearn import metrics
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.mixture import GaussianMixture
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn import svm
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import train_test_split
import umap
from imblearn.over_sampling import RandomOverSampler
from imblearn.under_sampling import RandomUnderSampler
from collections import Counter

# %% FUNCTONS
def load_csv_disaggregated_data(save=False) -> pd.DataFrame:
    """
    - load csv disaggregated data
    - filter columns and compress size
    - Optional: saves as Pickle file for later reload

    Parameters
    ----------
    save : boolean, optional.
        If True, saves output as pickle file in current directory.
        The default is False.

    Returns
    -------
    df : pd.DataFrame
        DataFrame with filtered columns and specified data types to reduce file
        size from 15 GB to 1.1 GB (with float32; only marginal number of results with >7 significant digits).
        
    """
    # WISE tabular data
    # filter columns and specify data types to reduce file size from 15.3 GB to 1.1 GB
    data_types = {"monitoringSiteIdentifier": "category",
                  "parameterWaterBodyCategory": "category",
                  "observedPropertyDeterminandLabel": "category",
                  "resultUom": "category",
                  "phenomenonTimeSamplingDate": "int32",
                  "resultObservedValue": "float32",
                  "resultQualityObservedValueBelowLOQ": "boolean",
                  "parameterSpecies": "category",
                  "resultObservationStatus": "category"}
        
    df = pd.read_csv("WISE/Waterbase_v2021_1_T_WISE6_DisaggregatedData.csv", usecols=list(data_types.keys()), dtype=data_types)
    if save == True:
         df.to_pickle("WISE/Data_EU_disaggregated_colFiltered.pkl")
    return df


def load_csv_aggregated_data(save=False) -> pd.DataFrame:
    """
    - load csv aggregated data
    - filter columns and compress size
    - Optional: saves as Pickle file for later reload

    Parameters
    ----------
    save : boolean, optional.
        If True, saves output as pickle file in current directory.
        The default is False.

    Returns
    -------
    df : pd.DataFrame
        DataFrame with filtered columns and specified data types to reduce file
        size from 1.1 GB to 210 MB (with float32) or 314 MB (eith float64).
        
    """
    float_type = "float64"
    data_types = {"monitoringSiteIdentifier": "category",
                  "parameterWaterBodyCategory": "category",
                  "observedPropertyDeterminandLabel": "category",
                  "resultUom": "category",
                  "phenomenonTimeReferenceYear": "category", 
                  "parameterSamplingPeriod": "category",
                  "procedureLOQValue": float_type,
                  "resultNumberOfSamples": "Int16",
                  "resultQualityNumberOfSamplesBelowLOQ": "Int16",
                  "resultQualityMinimumBelowLOQ": "boolean",
                  "resultMinimumValue": float_type,                   
                  "resultQualityMeanBelowLOQ": "boolean",
                  "resultMeanValue": float_type,
                  "resultQualityMaximumBelowLOQ": "boolean",
                  "resultMaximumValue": float_type,
                  "resultQualityMedianBelowLOQ": "boolean",
                  "resultMedianValue": float_type,
                  "resultStandardDeviationValue": float_type,
                  "metadata_observationStatus": "category"}     
    df = pd.read_csv("WISE/Waterbase_v2021_1_T_WISE6_AggregatedData.csv", usecols=list(data_types.keys()), dtype=data_types)
    if save == True:
         df.to_pickle("WISE/Data_EU_aggregated_colFiltered.pkl")
    return df


def investigate_data():
    # ------------------
    # General info
    # ------------------
    df.dtypes
    df.memory_usage()
    
    site_counts = df.monitoringSiteIdentifier.value_counts()
    target_counts = df.observedPropertyDeterminandLabel.value_counts()
    
    # ------------------
    # df <-> spatial match
    # ------------------
    # Are all sites in df listed in spatial? => No
    df.monitoringSiteIdentifier.isin(spatial.monitoringSiteIdentifier).all()
    # Are there NaNs in df sites? => No
    df.monitoringSiteIdentifier.isnull().sum()
    # Which sites in df are not listed in spatial? => 244 sites
    mask = site_counts.index.isin(spatial.monitoringSiteIdentifier)
    sites_unmatch = site_counts.index[~mask]
    
    # Are there NaNs in spatial? => Yes, 1631
    spatial.monitoringSiteIdentifier.isnull().sum()
    # Are there duplicates in spatial sites? => Yes
    # due to different monitoringSiteIdentifierScheme (the euMonitoringSiteCode has lat and lon, not the other) and NAs
    sum(spatial.monitoringSiteIdentifier.duplicated(), )
    mask = spatial.monitoringSiteIdentifier.duplicated(keep=False)
    sites_duplicated = spatial[mask].reset_index(drop=True)
    sites_idscheme = spatial.monitoringSiteIdentifierScheme.value_counts()
    sites_duplicated_idscheme = sites_duplicated.monitoringSiteIdentifierScheme.value_counts()


def prep_data():
    # ------------------
    # Purge df and spatial from unclear data and merge
    # ------------------
    """spatial:
        - remove all NAs site IDs
        - remove duplicates (keep in priority sites with euMonitoringSiteCode scheme)"""
    spatial_trim = spatial.dropna(subset = ["monitoringSiteIdentifier"])
    spatial_trim = spatial_trim.astype({'monitoringSiteIdentifier':'string', 'countryCode':'category', 'rbdName':'category'})
    mask = spatial_trim.monitoringSiteIdentifier.duplicated(keep=False)
    spatial_nondup = spatial_trim[~mask]
    spatial_dup = spatial_trim[mask]
    spatial_dup.sort_values(by='monitoringSiteIdentifierScheme', ascending= False, inplace=True)
    spatial_dup.drop_duplicates(subset=["monitoringSiteIdentifier"], inplace=True)
    spatial_trim = pd.concat([spatial_nondup, spatial_dup])[['monitoringSiteIdentifier','countryCode', 'rbdName']]
    
    df = df.loc[df.monitoringSiteIdentifier.isin(spatial.monitoringSiteIdentifier)].reset_index(drop=True)
    
    dfm = pd.merge(df, spatial_trim, how='left', on='monitoringSiteIdentifier').reset_index(drop=True)
    dfm.to_pickle("WISE/Data_EU_disaggregated_mergedSpatial.pkl")


def explore_data():
    # ------------------
    # Specific questions
    # ------------------
    
    # check counts for results above/below LOQ against
    # metadata_observation status (A: Normal record;U: Record with lower reliability;V: Unvalidated record)
    # metatdata_statusCode (experimental, stable, valid)
    test = df.groupby(["resultQualityObservedValueBelowLOQ", "metadata_statusCode",
                      "metadata_observationStatus"], as_index=False).size()
    
    # No of measurement per site, water body, river basin district (with names)
    wbody_counts = df.parameterWaterBodyCategory.value_counts()
    rbdIdent_counts = df.rbdIdentifier.value_counts()
    site_name_counts = df.monitoringSiteName.value_counts()
    wbody_name_counts = df.waterBodyName.value_counts()
    rbdIdent_name_counts = df.rbdName.value_counts()
    print("There are {} monitoring sites".format(site_counts.shape[0]))
    print("There are {} monitoring site names".format(site_name_counts.shape[0]))
    print("There are {} water bodies".format(wbody_counts.shape[0]))
    print("There are {} water bodies names".format(wbody_name_counts.shape[0]))
    print("There are {} river basin districts".format(rbdIdent_counts.shape[0]))
    print("There are {} river basin districts names".format(rbdIdent_name_counts.shape[0]))
    
    # Is there any measurement with no associated river basin disctict? => No
    df.rbdIdentifier.isnull().any()
    
    # check if targets have more than one associated measuring unit. => Problem with 2 targets
    # Chloride	['mg/L' 'mmol/L']	2
    # Dissolved oxygen	['mg/L' 'mg{O2}/L']	2
    # Phosphate	['mg{P}/L' 'mg{PO4}/L']	2
    target_uom = df.groupby('observedPropertyDeterminandLabel')['resultUom'].unique()
    target_uom_count = df.groupby('observedPropertyDeterminandLabel')['resultUom'].nunique()
    target_uom_mult = pd.merge(target_uom, target_uom_count, how='left', on='observedPropertyDeterminandLabel')
    target_uom_mult = target_uom_mult[target_uom_mult.resultUom_y > 1]
    pd.unique(df.observedPropertyDeterminandLabel).shape
    pd.unique(df.resultUom).shape
    
    # date of first and last measurement
    min(df.phenomenonTimeSamplingDate)
    max(df.phenomenonTimeSamplingDate)
    
    # how many NAs? 4887 NAs in result value (2%)
    print("There are {} empty result values, i.e. {}% of all results"
          .format(np.sum(df.resultObservedValue.isna()), np.sum(df.resultObservedValue.isna()) / df.shape[0] *100))
    
    """
    potential interesting targets
    Oxygen saturation: for fauna
    MCPA is a herbicide
    What is MTBE, AOX, NTA?
    """


if __name__ == "__main__":
    # %% LOAD FILES
    # path = "D:\Ludo\Docs\programming\CAS_applied_data_science\project_Water\Datasets".replace(
    #     "\\", "/")
    path = r"C:\Users\ludovic.lereste\Documents\CAS_applied_data_science\project_Water\Datasets" \
            .replace("\\", "/")
    os.chdir(path)
    
    
    # spatial = pd.read_csv("WISE/Waterbase_v2021_1_S_WISE6_SpatialObject_DerivedData.csv")
    # df_agg = load_csv_aggregated_data(save=True)

    # df = pd.read_pickle("WISE/Data_EU_disaggregated_colFiltered.pkl")
    # dfm = pd.read_pickle("WISE/Data_EU_disaggregated_mergedSpatial.pkl")
    df_agg = pd.read_pickle("WISE/Data_EU_aggregated_colFiltered.pkl")
    
    
    # %% IDENTIFY BEST CANDIDATES FOR TIME TRACE ANALYSIS
    """
    - filter out:
        + most results at LOQ
        + records with low reliability
        + sites not listed in spatial
    - get max number of results per site, per AT, per year
    - identify most represented AT for comparison between sites
    """
    meas_hist = df_agg.resultNumberOfSamples.hist(bins=100)
    df_agg.sort_values(by='resultNumberOfSamples', axis='index', ascending=False, inplace=True, na_position='first')
