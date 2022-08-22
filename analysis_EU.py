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


TO DO for classification analysis:
    - check which parameters are common to all sites or to all waterbody types
   
TO DO for time traces analysis:    
    - read dates from original csv file as strings and convert to datetime before saving pickle file
    - plot time traces
    - find interesting questions to answer
        + maybe for Electrical conductivity, but results are more interesting to show:
            - annual variation
            - depth variation

TO DO for mapping:
    - create map of all sites including
        + label color with water type
        + label details with site identifier, river body name, etc.
        
TO DO for documentation:
    - get detailed doc on functions using Sphinx
    - troubleshoot doc update on ReadTheDoc (or just import at the end)

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
import numpy as np
import os
import pandas as pd
from matplotlib.widgets import Slider
import matplotlib.dates as mdates
import matplotlib.colors as mcol



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
        size from 15 GB to 1.3 GB (with float32; only marginal number of results with >7 significant digits).

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
                  "parameterSampleDepth": "float32",
                  "resultObservationStatus": "category"}
    df = pd.read_csv("WISE/Waterbase_v2021_1_T_WISE6_DisaggregatedData.csv",
                     usecols=list(data_types.keys()), dtype=data_types)
    
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


def investigate_data(df, spatial):
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


def prep_data(df, spatial, save=False):
    """
    Purge df (dsaggregated data) and spatial from unclear data and merge
    spatial:
        - remove all NAs site IDs in spatial
        - remove duplicated sites in spatial (keep in priority sites with euMonitoringSiteCode scheme)
    result:
        - 60'171 rows are removed from df due to no site match in spatial (0.1% of rows in df)
        - 8'041 rows are dropped because result value is NA
        - Added column 'year' (for subsequent aggregation criteria)
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
    dfm.dropna(subset=['resultObservedValue'], inplace=True)
    
    dfm['year'] = dfm.phenomenonTimeSamplingDate .dt.year
    
    if save == True:
        dfm.to_pickle("WISE/Data_EU_disaggregated_mergedSpatial.pkl")
    
    return dfm


def find_time_traces(df_agg, spatial, thresh_samples_per_year=100, thresh_sites_per_target=10) -> pd.DataFrame:
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
    df_agg_fil = df_agg.loc[df_agg['metadata_observationStatus'] == 'A']

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
                      'parameterSamplingPeriod',
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
    
    return tt_id, tt_year_id, df_agg_colFil_2

def select_time_trace(dfm, tt_id_year, site: str, target: str) -> pd.DataFrame:
    """
    *** function unused since WISE aggregated data does not match disaggregated data :-( !!! ***
    simply filters disaggregated data based on site, target, and year

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
    tt_id_year_fil = tt_id_year.loc[site, target].dropna()
    mask = (dfm["monitoringSiteIdentifier"]==site) & \
            (dfm["observedPropertyDeterminandLabel"]==target) & \
            (dfm["phenomenonTimeSamplingDate"].dt.year.isin(tt_id_year_fil.index))
    tt = dfm[mask]
    
    return tt, tt_id_year_fil


def aggregate(df, groups: list, save=False) -> pd.DataFrame:
    """
    "C'est marqué dessus comme le port-salut"
    
    Aggregates raw data based on groups

    Parameters
    ----------
    df : pd.DataFrame
        Disaggregated data.
    groups : list
        list of columns to group.
    save : Boolean, optional
        Saves the output as a pickle file. The default is False.

    Returns
    -------
    df_agg : TYPE
        DESCRIPTION.

    """

    df_agg = df.groupby(by=groups, observed=True) \
        .agg({'resultObservedValue': ['count', 'mean', 'std'],
              'resultQualityObservedValueBelowLOQ': 'sum'})
        
    df_agg['resultQualityObservedValueBelowLOQ', 'perc'] = df_agg['resultQualityObservedValueBelowLOQ', 'sum'] / df_agg['resultObservedValue', 'count'] * 100
    df_agg.sort_values(by=('resultObservedValue', 'count'), ascending=False, inplace=True)
    df_agg.columns = ["_".join(a) for a in df_agg.columns.to_flat_index()]
             
    if save == True:
        df_agg.to_pickle("WISE/Data_EU_aggregated_custom_from_disaggregated.pkl")
        
    return df_agg


def find_time_traces_ca(dfm_agg,
                        dfm_agg_year,
                        thresh_LOQ=30,
                        thresh_samples=300,
                        thresh_samples_per_year=100,
                        thresh_sites_per_target=10,
                        thresh_siteYears_per_target=50):
    """
    Computes a summary of intersting sites, targets, and years for time-trace analysis.
    Takes custom aggregated (hence "_ca" suffix) data as input
    Returns pivot tables (actually unstacked tables) from filtered aggregated data.
    
    - filter out:
        + results too close to LOQ
        + time traces with too few points
    - computes number of results per site, per AT, per year
    - identify most represented AT for comparison between sites
        + identifiers: AT, site (,and year)
        + filter out AT with too few sites meeting criteria
        + filter out AT='Other chemical parameter' because wtf is that.

    Parameters
    ----------
    dfm_agg : TYPE
        DESCRIPTION.
    dfm_agg_year : TYPE
        DESCRIPTION.
    thresh_LOQ : TYPE, optional
        DESCRIPTION. The default is 30.
    thresh_samples : TYPE, optional
        DESCRIPTION. The default is 300.
    thresh_samples_per_year : TYPE, optional
        DESCRIPTION. The default is 100.
    thresh_sites_per_target : TYPE, optional
        DESCRIPTION. The default is 10.
    thresh_sites_per_target : TYPE, optional
        DESCRIPTION. The default is 10.

    Returns
    -------
    tt_id_ca : pd.DataFrame
        DataFrame with sites as index and targets as columns.
    tt_id_year_ca : pd.DataFrame
        DataFrame with sites and year as index and targets as columns.

    """
    # filter out results too close to LOQ
    # 11'049'525 ($071%) rows removed for 30% LOQ for dfm_agg_year
    dfm_agg = dfm_agg[dfm_agg['resultQualityObservedValueBelowLOQ_perc'] < thresh_LOQ]
    dfm_agg_year = dfm_agg_year[dfm_agg_year['resultQualityObservedValueBelowLOQ_perc'] < thresh_LOQ]
    
    # Filter out time traces with too few points
    dfm_agg = dfm_agg[dfm_agg['resultObservedValue_count'] > thresh_samples]
    dfm_agg_year = dfm_agg_year[dfm_agg_year['resultObservedValue_count'] > thresh_samples_per_year]
    
    # computes the "pivot table"
    tt_id_ca = dfm_agg['resultObservedValue_count'].unstack(level=1)
    tt_id_year_ca = dfm_agg_year['resultObservedValue_count'].unstack(level=1)
    
    # filter out targets with too few sites or site-year 
    tt_id_ca = tt_id_ca.dropna(axis=1, thresh=thresh_sites_per_target)
    tt_id_year_ca = tt_id_year_ca.dropna(axis=1, thresh=thresh_siteYears_per_target)
    
    tt_id_ca = tt_id_ca.drop(columns='Other chemical parameter')
    tt_id_year_ca = tt_id_year_ca.drop(columns='Other chemical parameter')
    
    # remove empty rows due to previous column (i.e. target) drops
    tt_id_ca = tt_id_ca.dropna(axis=0, how='all')
    tt_id_year_ca = tt_id_year_ca.dropna(axis=0, how='all')
    
    return tt_id_ca, tt_id_year_ca

    
def select_time_trace_ca(dfm, tt_id_year, site: str, target: str) -> pd.DataFrame:
    """
    simply filters disaggregated data based on site, target, and year

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
    
    tt_id_year_fil = tt_id_year.loc[site, target].dropna()
    mask = (dfm["monitoringSiteIdentifier"]==site) & \
            (dfm["observedPropertyDeterminandLabel"]==target) & \
            (dfm["year"].isin(tt_id_year_fil.index))
                
    tt = dfm[mask]
    
    return tt

def prep_plot(dfm, tt_id_year, target):
    """N.B.: generating a MultiIndex (site, year) is much faster than .apply(isin())"""
    site_year_filter = tt_id_year.loc[:, target].dropna().index
    tts = dfm[dfm["observedPropertyDeterminandLabel"]==target]
    site_year = pd.MultiIndex.from_frame(tts[['monitoringSiteIdentifier', 'year']])
    tts = tts[site_year.isin(site_year_filter)]
    
    sites = site_year_filter.get_level_values(0).unique()
    tts_per_site = []
    for i, el in enumerate(sites):
        tts_per_site.append(tts[tts.monitoringSiteIdentifier == sites[i]]
                            .sort_values(by='phenomenonTimeSamplingDate'))
                                 
    return tts, tts_per_site

def dump():
    """miscalleneous commands"""
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
    
    # compare year distributions between aggregated and disaggregated data
    years_agg = df_agg.groupby(by=['phenomenonTimeReferenceYear']).sum()
    years_disagg = df.phenomenonTimeSamplingDate.dt.year.value_counts().sort_index()
    years = pd.merge(left=years_agg.resultNumberOfSamples,
                     right=years_disagg,
                     how='outer',
                     left_on='phenomenonTimeReferenceYear',
                     right_index=True)
    years.plot(x="phenomenonTimeReferenceYear", y=["resultNumberOfSamples", "phenomenonTimeSamplingDate"], kind="bar")
    

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
    # dfm = prep_data(df, spatial, save=True)
    
    # FROM PICKLE
    # df = pd.read_pickle("WISE/Data_EU_disaggregated_colFiltered.pkl")
    dfm = pd.read_pickle("WISE/Data_EU_disaggregated_mergedSpatial.pkl")
    df_agg = pd.read_pickle("WISE/Data_EU_aggregated_colFiltered.pkl")
    dfm_agg = pd.read_pickle("WISE/Data_EU_aggregated_custom_from_disaggregated.pkl")
    dfm_agg_year = pd.read_pickle("WISE/Data_EU_aggregated_custom_perYear_from_disaggregated.pkl")
    
    
    
    # %% IDENTIFY BEST CANDIDATES FOR TIME TRACE ANALYSIS
    
    # # From WISE aggregated data
    # """This is useless since WISE aggregated data do not match WISE disaggregated data"""
    # tt_id, tt_year_id = find_time_traces(df_agg, spatial)
    
    # # From custom-aggregated data
    # dfm_agg = aggregate(dfm, groups=['monitoringSiteIdentifier', 'observedPropertyDeterminandLabel'], 
    #                     save=True)
    # dfm_agg_year = aggregate(dfm, groups=['monitoringSiteIdentifier',
    #                                       'observedPropertyDeterminandLabel',
    #                                       'year'])
    # dfm_agg_year.to_pickle("WISE/Data_EU_aggregated_custom_perYear_from_disaggregated.pkl")
    tt_id, tt_id_year = find_time_traces_ca(dfm_agg, dfm_agg_year)
    
    
    # %% PREP FOR PLOTS
    """Nota Bene:
        - Electrical conductivity has many different sampling depths
            + color data points with sample depth value
        - Conductivity seem to go throuh yearly cycles. INTERESTING !!
        - tt_pH[24] is taken in 1973-1974 and is an exception
        - tt_pH[0] there is A LOT of results taken the same day (30 per day on average...)"""
        
    target = 'pH'
    units = dfm.loc[dfm["observedPropertyDeterminandLabel"]==target, "resultUom"].unique()
    tts, tts_per_site = prep_plot(dfm, tt_id_year, target=target)
    
    
    
    # %% PLOT - line plots
    # """to be done with a slide"""
    # for i in range(3):
    #     tts_per_site[i].plot(x='phenomenonTimeSamplingDate', y='resultObservedValue')
    
    # # test_tt = select_time_trace_ca(dfm,
    # #                           tt_id_year_ca,
    # #                           site='PT19E02',
    # #                           target='Oxygen saturation')
                                                                           
    # # plt.plot(tts.phenomenonTimeSamplingDate, tts.resultObservedValue)
    
    
    # n_rows = 4
    # n_cols = 1
    # n_plots = n_rows * n_cols
    # fig, axs = plt.subplots(nrows=n_rows, ncols=n_cols) #, sharey=True
    # plt.get_current_fig_manager().window.state('zoomed')
    # plt.subplots_adjust(top=0.90)

    # # inititate plot lines (ls_xxx) that will later be updated by the slider
    # ls_data = []

    # for i, ax in enumerate(axs.flat):
    #     ls_data.append(ax.plot(tts_per_site[i]['phenomenonTimeSamplingDate'], 
    #                            tts_per_site[i]['resultObservedValue']))

    # ax_slider = plt.axes([0.25, 0.95, 0.65, 0.03])
    # slider_packet = Slider(ax=ax_slider,
    #                        label='Test packet',
    #                        valmin=0,
    #                        valmax=len(tts_per_site)-n_plots,
    #                        valstep=n_plots,
    #                        valinit=0)

    # def update_slider(val):
    #     for i in np.arange(n_plots):
    #         ls_data[i][0].set_data(tts_per_site[i+val]['phenomenonTimeSamplingDate'],
    #                                 tts_per_site[i+val]['resultObservedValue'])

            
    # slider_packet.on_changed(update_slider)
    
    # %% PLOT 2 - scatter plots
    """I think it is a bit buggy since there are time traces with 100+ points that I do not see"""
  
   # fig, axs = plt.subplots(2,2)  
   # for i in range(4):
   #       # plt.figure()
   #       # fig, ax = plt.subplots(2,2)
   #       axs.flat[i].scatter(tts_per_site[i]['phenomenonTimeSamplingDate'], 
   #                  tts_per_site[i]['resultObservedValue'])
   #       # plt.show()
  
    n_rows = 4
    n_cols = 2
    n_plots = n_rows * n_cols
    fig, axs = plt.subplots(nrows=n_rows, ncols=n_cols) #, sharey=True
    # plt.get_current_fig_manager().window.state('zoomed')
    plt.get_current_fig_manager().window.showMaximized()
    plt.subplots_adjust(top=0.90, bottom=0.05)
    fig.suptitle(target)
    
    # colormap = mcol.ListedColormap('cool')#.reversed()
    # colormap.get_bad()
    # colormap.set_bad("black")
    # colormap = mcol.LinearSegmentedColormap.from_list("Blue-Red-Colormap", ["b", "r"])
    # colormap.set_bad("black")
    # colormap = mcol.LinearSegmentedColormap.from_list("Blues", "b")
    # colormap.set_bad("yellow")
    
    ax_slider = plt.axes([0.25, 0.92, 0.65, 0.03])
    slider_tt = Slider(ax=ax_slider,
                            label='Time traces',
                            valmin=0,
                            valmax=len(tts_per_site)-n_plots,
                            valstep=n_plots,
                            valinit=0)

    for i, ax in enumerate(axs.flat):
        axs.flat[i].scatter(tts_per_site[i]['phenomenonTimeSamplingDate'], 
                   tts_per_site[i]['resultObservedValue'],
                    c=tts_per_site[i]['parameterSampleDepth'],
                    cmap='Blues',
                    plotnonfinite=True,
                    edgecolors='black',
                    linewidths=0.1,
                    s=10)    
        ax.xaxis.set_major_locator(mdates.YearLocator())
        unit = tts_per_site[i].loc[tts_per_site[i]["observedPropertyDeterminandLabel"]==target, "resultUom"].unique().to_numpy()
        ax.set_ylabel(f"{unit}")
        ax.text(0.01, 0.04, 
                f"{tts_per_site[i].monitoringSiteIdentifier.unique()[0]}, {tts_per_site[i].parameterWaterBodyCategory.unique()[0]}",
                transform=ax.transAxes)
        

    def update_slider(val):
        for i, ax in enumerate(axs.flat):
            ax.clear()
            ax.scatter(tts_per_site[i+val]['phenomenonTimeSamplingDate'],
                        tts_per_site[i+val]['resultObservedValue'],
                        c=tts_per_site[i+val]['parameterSampleDepth'],
                        cmap='Blues',
                        plotnonfinite=True,
                        edgecolors='black',
                        linewidths=0.1,
                        s=10)       
            ax.xaxis.set_major_locator(mdates.YearLocator())
            unit = tts_per_site[i+val].loc[tts_per_site[i+val]["observedPropertyDeterminandLabel"]==target, "resultUom"].unique().to_numpy()
            ax.set_ylabel(f"{unit}")
            ax.text(0.01, 0.04, 
                    f"{tts_per_site[i+val].monitoringSiteIdentifier.unique()[0]}, {tts_per_site[i+val].parameterWaterBodyCategory.unique()[0]}",
                    transform=ax.transAxes)
            
    slider_tt.on_changed(update_slider)
    
    
    

