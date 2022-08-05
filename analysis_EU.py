"""
Created on 28 June 2022

@author = "Ludovic Le Reste"
@credits = ["Hisham Ben Hamidane", "Ludovic Le Reste"]
@status = "Prototype"

This scripts loads and analyses the Waterbase WISE 6 database from the European Environment Agency.
Analysis is for the whole dataset (not limited to Switzerland as previously)
Data sources:
2021 WISE dataset (disaggregated data, aggregated data, spatial data):
    https://www.eea.europa.eu/data-and-maps/data/waterbase-water-quality-icm-2

TO DO:
    - check completion of data

TO DO for classification analysis:
    - check which parameters are common to all sites or to all waterbody types

    
TO DO for time traces analysis:    
    - read dates from original csv file as strings and convert to datetime before saving pickle file
    - plot time traces
    - find interesting questions to answer
        + is there a pattern in terms of metal pollution?
        [10:27] Ben Hamidane, Hisham

Targets of interest from Hisham:
    - Groupe 1: cibles analytiques générales
    AT_group1 <- c("pH", "Water Temperature", "Oxygen saturation", "Electrical conductivity", "Ammonium", "Phosphate", "Nitrate", "BOD5")
    - Groupe 2: métaux
    AT_group2 <- c("Lead and its compounds", "Copper and its compounds", "Cadmium and its compounds", "Zinc and its compounds",
               "Nickel and its compounds", "Arsenic and its compounds", "Mercury and its compounds", "Chromium and its compounds",
               "Iron and its compounds")

"""

__author__ = "Ludovic Le Reste"
__credits__ = ["Hisham Ben Hamidane", "Ludovic Le Reste"]
__status__ = "Prototype"

# %% PACKAGES
# from collections import Counter
# from imblearn.under_sampling import RandomUnderSampler
# from imblearn.over_sampling import RandomOverSampler
# import umap
# from sklearn.model_selection import train_test_split
# from sklearn.feature_selection import VarianceThreshold
# from sklearn.preprocessing import StandardScaler
# from sklearn import svm
# from sklearn.manifold import TSNE
# from sklearn.decomposition import PCA
# from sklearn.mixture import GaussianMixture
# from sklearn.metrics import silhouette_score
# from sklearn.cluster import KMeans
# from sklearn import metrics
# from scipy import stats
# import statsmodels.api as sm
# import geopandas as gpd
import matplotlib.pyplot as plt
# import seaborn as sns
# import numpy as np
import os
import pandas as pd



# %% PACKAGES
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"
# from mpl_toolkits.basemap import Basemap
# import plot_style

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
                  "phenomenonTimeSamplingDate": "str",
                  "resultObservedValue": "float32",
                  "resultQualityObservedValueBelowLOQ": "boolean",
                  "parameterSpecies": "category",
                  "resultObservationStatus": "category"}
    df1 = pd.read_csv("WISE/Waterbase_v2021_1_T_WISE6_DisaggregatedData.csv",
                     usecols=list(data_types.keys()), dtype=data_types)
    df2 = pd.read_csv("WISE/Waterbase_v2021_2_T_WISE6_DisaggregatedData.csv",
                     usecols=list(data_types.keys()), dtype=data_types)
    # df = 
    
    df["phenomenonTimeSamplingDate"] = pd.to_datetime(df["phenomenonTimeSamplingDate"], format='%Y%m%d')
    
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
        size from 1.1 GB to 210 MB (with float32) or 317 MB (with float64).

    """
    float_type = "float64"
    data_types = {"monitoringSiteIdentifier": "category",
                  "parameterWaterBodyCategory": "category",
                  "observedPropertyDeterminandLabel": "category",
                  "resultUom": "category",
                  "phenomenonTimeReferenceYear": "int16",
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
                  "parameterSampleDepth": float_type,
                  "metadata_observationStatus": "category"}
    df = pd.read_csv("WISE/Waterbase_v2021_1_T_WISE6_AggregatedData.csv",
                     usecols=list(data_types.keys()), dtype=data_types)
    if save == True:
        df.to_pickle("WISE/Data_EU_aggregated_colFiltered.pkl")
    return df


def investigate_data(df):
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


def prep_data(df, spatial):
    """
    Purge df (dsaggregated data) and spatial from unclear data and merge
    spatial:
        - remove all NAs site IDs in spatial
        - remove duplicated sites in spatial (keep in priority sites with euMonitoringSiteCode scheme)
    result:
        - 60'171 rows are removed from df to create dfm (0.1% of rows in df)
    """
    spatial_trim = spatial.dropna(subset=["monitoringSiteIdentifier"])
    spatial_trim = spatial_trim.astype(
        {'monitoringSiteIdentifier': 'string', 'countryCode': 'category', 'rbdName': 'category'})
    mask = spatial_trim.monitoringSiteIdentifier.duplicated(keep=False)
    spatial_nondup = spatial_trim[~mask]
    spatial_dup = spatial_trim[mask]
    spatial_dup.sort_values(
        by='monitoringSiteIdentifierScheme', ascending=False, inplace=True)
    spatial_dup.drop_duplicates(
        subset=["monitoringSiteIdentifier"], inplace=True)
    spatial_trim = pd.concat([spatial_nondup, spatial_dup])[
        ['monitoringSiteIdentifier', 'countryCode', 'rbdName']]

    df = df.loc[df.monitoringSiteIdentifier.isin(
        spatial.monitoringSiteIdentifier)].reset_index(drop=True)

    dfm = pd.merge(df, spatial_trim, how='left',
                   on='monitoringSiteIdentifier').reset_index(drop=True)
    dfm.to_pickle("WISE/Data_EU_disaggregated_mergedSpatial.pkl")
    
    return dfm


def find_time_traces(df_agg, thresh_samples_per_year=100, thresh_sites_per_target=10) -> pd.DataFrame:
    """
    Computes a summary of intersting sites, targets, and years for time-trace analysis.
    Returns a pivot table from filtered aggregated data.
    
    - filter out:
        + records with low reliability
        + sites not listed in spatial
        + most results at LOQ
        + time traces with few points
    - get max number of results per site, per AT, per year
    - identify most represented AT for comparison between sites
        + identifiers: AT, site
        + filter out AT with too few sites meeting criteria

    Parameters
    ----------
    df_agg : pd.DataFrame
        aggregated data
    thresh_samples_per_year : int, optional
        The default is 100.
    thresh_sites_per_target : int, optional
        The default is 10.

    Returns
    -------
    tt_id : pd.DataFrame
        DataFrame with sites as index and targets as columns.
    tt_year_id : pd.DataFrame
        DataFrame with sites and year as index and targets as columns.
    """
    
    # filter out results with low reliability
    df_agg = df_agg.loc[df_agg['metadata_observationStatus'] == 'A']

    # filter out sites that are not in spatial => There are none!!
    df_agg = df_agg.loc[df_agg.monitoringSiteIdentifier.isin(spatial.monitoringSiteIdentifier)]

    # filter out results too close to LOQ
    df_agg = df_agg.loc[(df_agg['resultQualityMedianBelowLOQ']==False) | (df_agg['resultQualityMedianBelowLOQ'].isnull())]
    df_agg = df_agg.loc[(df_agg['resultQualityMaximumBelowLOQ']==False) | (df_agg['resultQualityMaximumBelowLOQ'].isnull())]
    df_agg = df_agg.loc[(df_agg['resultQualityMeanBelowLOQ']==False) | (df_agg['resultQualityMeanBelowLOQ'].isnull())]

    # filter out time traces with few points
    df_agg = df_agg.loc[df_agg.resultNumberOfSamples>thresh_samples_per_year]
    
    # convert some categorical data into relevant dtypes for further analysis
    """ This is done at this stage as df_agg has now a much smaller size"""
    df_agg = df_agg.astype({'phenomenonTimeReferenceYear': 'int16',
                   'monitoringSiteIdentifier': 'str',
                   'observedPropertyDeterminandLabel': 'str'})
    
    # compute table with AT and site identifiers
    df_agg_colFil = df_agg[['monitoringSiteIdentifier', 
                      'observedPropertyDeterminandLabel',
                      'phenomenonTimeReferenceYear',
                      'resultNumberOfSamples']]
    tt_id_raw = pd.pivot_table(df_agg_colFil,
                          values='resultNumberOfSamples',
                          index='monitoringSiteIdentifier',
                          columns='observedPropertyDeterminandLabel',
                          aggfunc='count')
    tt_id = tt_id_raw.dropna(axis=1, thresh=thresh_sites_per_target)
    df_agg_colFil_2 = df_agg_colFil[df_agg_colFil.observedPropertyDeterminandLabel.isin(tt_id.columns)]
    tt_year_id = pd.pivot_table(df_agg_colFil_2,
                          values='resultNumberOfSamples',
                          index=['monitoringSiteIdentifier', 'phenomenonTimeReferenceYear'],
                          columns='observedPropertyDeterminandLabel',
                          aggfunc='count')
    
    return tt_id, tt_year_id

def select_time_trace(dfm, tt_year_id, site: str, target: str) -> pd.DataFrame:
    """
    WEIRD, there is no match in years between the filtered aggregated data and the raw disaggregated data

    Parameters
    ----------
    dfm : TYPE
        DESCRIPTION.
    tt_year_id : TYPE
        DESCRIPTION.
    site : str
        DESCRIPTION.
    target : str
        DESCRIPTION.

    Returns
    -------
    tt : TYPE
        DESCRIPTION.

    """
    tt_year_id_fil = tt_year_id.loc[site, target].dropna()
    mask = (dfm["monitoringSiteIdentifier"]==site) #& \
            # (dfm["observedPropertyDeterminandLabel"]==target) #& \
            # (dfm["phenomenonTimeSamplingDate"].dt.year.isin(tt_year_id_fil.index))
    tt = dfm[mask]
    
    return tt, tt_year_id_fil

def dump():
    # status and visual aids
    targets = df_agg.observedPropertyDeterminandLabel.value_counts()
    sites = df_agg.monitoringSiteIdentifier.value_counts()
    n_sites_per_target = df_agg \
        .groupby(['observedPropertyDeterminandLabel'], observed=True) \
        .monitoringSiteIdentifier \
        .nunique() \
        .sort_values(ascending=False)
        
    df_agg.phenomenonTimeReferenceYear.value_counts().sort_index().plot(kind='bar')
        
    targets_id = targets[targets>100]
    
    df_agg.loc[df_agg.observedPropertyDeterminandLabel==n_sites_per_target.index[0]] \
        .monitoringSiteIdentifier \
        .unique() \
        .sort_values(ascending=False)
        
    df_agg.resultNumberOfSamples.hist(bins=50)
    df_agg.sort_values(by='resultNumberOfSamples', axis='index',
                            ascending=False, inplace=True, na_position='last')
    
    

if __name__ == "__main__":
    # %% LOAD FILES
    # path = "D:\Ludo\Docs\programming\CAS_applied_data_science\project_Water\Datasets".replace(
    #     "\\", "/")
    path = r"C:\Users\ludovic.lereste\Documents\CAS_applied_data_science\project_Water\Datasets" \
        .replace("\\", "/")
    os.chdir(path)

    # FROM CSV
    spatial = pd.read_csv("WISE/Waterbase_v2021_1_S_WISE6_SpatialObject_DerivedData.csv")
    # df = load_csv_disaggregated_data(save=True)
    # df_agg = load_csv_aggregated_data(save=True)
    # dfm = prep_data(df, spatial)
    
    # FROM PICKLE
    df = pd.read_pickle("WISE/Data_EU_disaggregated_colFiltered.pkl")
    dfm = pd.read_pickle("WISE/Data_EU_disaggregated_mergedSpatial.pkl")
    df_agg = pd.read_pickle("WISE/Data_EU_aggregated_colFiltered.pkl")
    
    # %% IDENTIFY BEST CANDIDATES FOR TIME TRACE ANALYSIS
    tt_id, tt_year_id = find_time_traces(df_agg)
    
    # %% PLOT TIME TRACES
    tt, tt_year_id_fil = select_time_trace(dfm, tt_year_id, site='ES020ESPF004300171', target='Dissolved oxygen')
    dfm.phenomenonTimeSamplingDate.dt.year.value_counts().plot(kind='bar')
    plt.xlabel(xlabel, kwargs)
    
    

