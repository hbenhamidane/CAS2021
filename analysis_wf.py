"""
Created on 21 December 2021

@author = "Ludovic Le Reste"
@credits = ["Hisham Ben Hamidane", "Ludovic Le Reste"]
@status = "Prototype"

This module loads and analyses the Waterbase WISE 6 database from the European Environment Agency.
Analysis was limited to data from Switzerland.
This module is based on the first Analysis.py script. This module restructure 
the previous script into funcntions to facilitate various workflows.

Data sources
Measurements, disaggregated data:
    https://discomap.eea.europa.eu/App/DiscodataViewer/?fqn=[WISE_SOE].[v1r1].[Waterbase_T_WISE6_DisaggregatedData]
Spatial Metadata:
    https://discomap.eea.europa.eu/App/DiscodataViewer/?fqn=[WISE_SOE].[v1r1].[Waterbase_S_WISE_SpatialObject_DerivedData]
    

TO DO:
- equivalent of visdat
    + looks like pandas-profiling is a good candidate (run first on a small dataset)
    + perhaps in conjunction with great_expectations?

- check map coordinate conversion with function to_crs() instead of manual transformation
- use ecomorphology map (continuous rivers) instead of typology map
    + troubleshoot LV03 Vs LV95
"""
__author__ = "Ludovic Le Reste"
__credits__ = ["Hisham Ben Hamidane", "Ludovic Le Reste"]
__status__ = "Prototype"

 
# %% PACKAGES
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"
import pandas as pd
import os
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
# from imblearn.over_sampling import RandomOverSampler
# from imblearn.under_sampling import RandomUnderSampler
# from collections import Counter

# Use my plotting custom style
# plot_style.custom_style()  


# %% LOAD FILES
def load_files(path):
    """
    Sets "path" as root directory and load files from there.

    Parameters
    ----------
    path : 'string'
        Path of the root folder where datasets are located.

    Returns
    -------
    None.

    """
    os.chdir(path)
    
    # WISE tabular data
    df = pd.read_csv("WISE/waterbase_t_wise6_disagregateddata_CH.csv")
    spatial = pd.read_csv("WISE/waterbase_s_wise_spatialobject_deriveddata_CH.csv")
    df_def = pd.read_csv(
        "WISE/Waterbase_v2020_1_dataset_definition_T_WISE6_DisaggregatedData.csv")
    spatial_def = pd.read_csv(
        "WISE/Waterbase_v2020_1_dataset_definition_S_WISE6_SpatialObject_DerivedData.csv")
    
    return(df, spatial, df_def, spatial_def)


# %% DATA PREP
def merge_data(df, spatial):
    """
    Merges the main data (measurment values) with the spatial metadata.

    Parameters
    ----------
    df : pandas DataFrame
        Main tabular data from WISE.
    spatial : pandas DataFrame
        Spatial metadata from WISE.

    Returns
    -------
    dfm : merged data)
    df : modified version of df
    spatial : modified version of spatial
    """
    
    # ------------------
    # Clean up main data
    # ------------------
    # drop empty columns
    empty_cols = [col for col in df.columns if df[col].isnull().all()]
    df.drop(empty_cols, axis=1, inplace=True)
    
    # drop uninteresting columns
    """
    col 2 monitoringSiteIdentifierScheme is 'eionetMonitoringSiteCode'
    col 6 procedureAnalysedMatrix is W-DIS or W
    col 10 sampleIdentifier is 0 or na
    """
    test = df.iloc[:, 16].value_counts(dropna=False)
    test = df.loc[:, "resultQualityObservedValueBelowLOQ"].value_counts(
        dropna=False)
    
    drops = ["monitoringSiteIdentifierScheme",
             "procedureAnalysedMatrix", "sampleIdentifier"]
    df.drop(drops, axis=1, inplace=True)
    
    # convert data types
    df.phenomenonTimeSamplingDate = pd.to_datetime(df.phenomenonTimeSamplingDate)
    
    # ------------------
    # Clean up spatial data
    # ------------------
    # drop empty columns
    empty_cols_spatial = [
        col for col in spatial.columns if spatial[col].isnull().all()]
    spatial.drop(empty_cols_spatial, axis=1, inplace=True)
    
    # drop other uninteresting columns (with IsentifierScheme number and duplicated columns with df)
    drops_spatial = spatial.filter(
        regex='Scheme').columns.append(pd.Index(['countryCode']))
    spatial.drop(drops_spatial, axis=1, inplace=True)
    # N.B.: can also be done using spatial.columns.str.contains('Scheme')
    
    # convert latitude and longitude (from WGS84 coordinate system) to LV95 swiss system
    # More info: https://www.swisstopo.admin.ch/en/maps-data-online/calculation-services.html
    """
    Apparently I could use function to_crs() with crs='EPSG = 2056'
    """
    spatial['lon_arcsec'] = (spatial.lon * 3600 - 26782.5) / 10000
    spatial['lat_arcsec'] = (spatial.lat * 3600 - 169028.66) / 10000
    spatial['swiss_E'] = 2600072.37 \
        + 211455.93 * spatial.lon_arcsec \
        - 10938.51 * spatial.lon_arcsec * spatial.lat_arcsec \
        - 0.36 * spatial.lon_arcsec * spatial.lat_arcsec**2 \
        - 44.54 * spatial.lon_arcsec**3
    spatial['swiss_N'] = 1200147.07 \
        + 308807.95 * spatial.lat_arcsec \
        + 3745.25 * spatial.lon_arcsec**2 \
        + 76.63 * spatial.lat_arcsec**2 \
        - 194.56 * spatial.lon_arcsec**2 * spatial.lat_arcsec \
        + 119.79 * spatial.lat_arcsec**3
    
    # monitoring sites Vs water bodies
    
    
    # ------------------
    # merge datasets
    # ------------------
    # Are all sites in df listed in spatial? => Yes
    df.monitoringSiteIdentifier.isin(spatial.monitoringSiteIdentifier).all()
    
    # Are there NaNs in df sites? => No
    df.monitoringSiteIdentifier.isnull().sum()
    
    # Are there duplicates in spatial sites? => No, only NaNs are duplicates
    sum(spatial.monitoringSiteIdentifier.dropna().duplicated(), )
    
    # Merge
    dfm = pd.merge(df, spatial, how='left', on='monitoringSiteIdentifier')
    
    # Rename river body types
    dfm.loc[df.parameterWaterBodyCategory == 'GW', 'parameterWaterBodyCategory'] = 'ground'
    dfm.loc[df.parameterWaterBodyCategory == 'RW', 'parameterWaterBodyCategory'] = 'river'
    dfm.loc[df.parameterWaterBodyCategory == 'LW', 'parameterWaterBodyCategory'] = 'lake'
    
    return dfm, df, spatial


# %% DATA EXPLORATION
def exolore_data():
    """
    Creates arrays/series and print some characteristics of the datasets.
    
    TBD, potential interesting targets:
    Oxygen saturation for fauna
    MCPA is a herbicide
    What is MTBE, AOX, NTA?

    Returns
    -------
    target_uom: pandas Series
        Unit of measures of all analytical targets (used in explore_per_target())

    """
    # ------------------
    # General info
    # ------------------
    df.dtypes
    df.memory_usage()
    
    meas_persite_count = df.monitoringSiteIdentifier.value_counts()
    target_counts = df.observedPropertyDeterminandLabel.value_counts()
    
    # ------------------
    # Specific questions
    # ------------------
    # sites with target Ibuprofene
    site_ibu = df.loc[df.observedPropertyDeterminandLabel ==
                      "Ibuprofen", "monitoringSiteIdentifier"].value_counts()
    
    # check counts for results above/below LOQ against
    # metadata_observation status (A: Normal record;U: Record with lower reliability;V: Unvalidated record)
    # metatdata_statusCode (experimental, stable, valid)
    test = df.groupby(["resultQualityObservedValueBelowLOQ", "metadata_statusCode",
                      "metadata_observationStatus"], as_index=False).size()
    
    # No of measurement per site, water body, river basin district (with names)
    site_counts = dfm.monitoringSiteIdentifier.value_counts()
    wbody_counts = dfm.waterBodyIdentifier.value_counts()
    rbdIdent_counts = dfm.rbdIdentifier.value_counts()
    site_name_counts = dfm.monitoringSiteName.value_counts()
    wbody_name_counts = dfm.waterBodyName.value_counts()
    rbdIdent_name_counts = dfm.rbdName.value_counts()
    print("There are {} monitoring sites".format(site_counts.shape[0]))
    print("There are {} monitoring site names".format(site_name_counts.shape[0]))
    print("There are {} water bodies".format(wbody_counts.shape[0]))
    print("There are {} water bodies names".format(wbody_name_counts.shape[0]))
    print("There are {} river basin districts".format(rbdIdent_counts.shape[0]))
    print("There are {} river basin districts names".format(rbdIdent_name_counts.shape[0]))
    
    # Is there any measurement with no associated river basin disctict? => No
    dfm.rbdIdentifier.isnull().any()
    
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
    min(dfm.phenomenonTimeSamplingDate)
    max(dfm.phenomenonTimeSamplingDate)
    
    # how many NAs? 4887 NAs in result value (2%)
    print("There are {} empty result values, i.e. {}% of all results"
          .format(np.sum(df.resultObservedValue.isna()), np.sum(df.resultObservedValue.isna()) / df.shape[0] *100))
    
    return target_uom

def explore_per_target(dfm, target_uom):
    """
    Analyse and plots various characteristics per analytical target.

    Parameters
    ----------
    dfm : pandas DataFrame
        merged data.

    Returns
    -------
    None.

    """
    target = "pH"
    bins = {"Nitrate": np.linspace(0, 60, 100),
            "pH": np.linspace(5, 10, 50),
            "Electrical conductivity": np.linspace(0, 1200, 100)}
    
    # Histograms: value distibution + time distribution
    axs = dfm.loc[(dfm.observedPropertyDeterminandLabel == target) 
                  & (dfm.resultQualityObservedValueBelowLOQ == 0), [
                      'resultObservedValue', 'phenomenonTimeSamplingDate']
                  ].hist(bins=100)
    axs[0][0].set_xlabel('{} in {}'.format(target, target_uom.loc[target][0]))
    axs[0][0].set_ylabel('counts')
    axs[0][1].set_xlabel('year')
    axs[0][1].set_ylabel('counts')
    plt.suptitle(target)
    
    # Histogram: value distibution
    plt.figure()
    ax = dfm.loc[(dfm.observedPropertyDeterminandLabel == target) 
                  & (dfm.resultQualityObservedValueBelowLOQ == 0), [
                      'resultObservedValue']
                  ].hist(bins=bins[target],
                         grid=False)
    ax[0][0].set_xlabel('{} in {}'.format(target, target_uom.loc[target][0]))
    ax[0][0].set_ylabel('counts')
    plt.suptitle(target)
    
    # Histogram: value distribution, colour-coded with water body
    plt.figure()
    ax = plt.subplot()
    par = [{"wbody": "ground", "color": "brown", "alpha": 0.5, "label": "Ground Water"},
           {"wbody": "lake", "color": "green", "alpha": 0.8, "label": "Lake Water"},
           {"wbody": "river", "color": "blue", "alpha": 0.3, "label": "River Water"}]
    for el in par:
        plt.hist(dfm.loc[(dfm.observedPropertyDeterminandLabel == target)
                  & (dfm.resultQualityObservedValueBelowLOQ == 0)
                  & (dfm.parameterWaterBodyCategory == el["wbody"]),
                  'resultObservedValue'],
                 color=el["color"],
                 alpha=el["alpha"],
                 histtype='stepfilled',
                 edgecolor='black',
                 label=el["label"],
                 bins=bins[target])
    plt.legend(loc='upper right')
    ax.set_ylabel('counts')
    ax.set_xlabel('{} in {}'.format(target, target_uom.loc[target][0]))
    plt.title(target)
    
    # Box plot: single target by type of water body
    t1 = dfm.loc[dfm.observedPropertyDeterminandLabel == target, ['resultObservedValue', 'specialisedZoneType']]
    ax2 = t1.boxplot(by='specialisedZoneType')
    ax2.set_title("{}".format(target))
    # plt.ylim(6, 10)
    
    # Time series per monitoring site
    # select site with most measurements
    site_count_per_target = dfm[dfm.observedPropertyDeterminandLabel == target].monitoringSiteIdentifier.value_counts()
    site_ID = site_count_per_target.index[0]
    site_name = dfm.loc[dfm.monitoringSiteIdentifier == site_ID].monitoringSiteName.iloc[0]
    
    plt.figure()
    dfm_t = dfm.set_index('phenomenonTimeSamplingDate')
    ax = dfm_t.loc[(dfm_t.monitoringSiteIdentifier == site_ID) &
              (dfm_t.observedPropertyDeterminandLabel == target),
              'resultObservedValue'].plot()
    ax.set_ylabel('{} in {}'.format(target, target_uom.loc[target][0]))
    plt.title("Time series for {}".format(site_name))
    
    # QQ plots
    fig = sm.qqplot(dfm.loc[(dfm.observedPropertyDeterminandLabel == target) 
                  & (dfm.resultQualityObservedValueBelowLOQ == 0),
                      'resultObservedValue'
                  ], line='q')
    plt.title("QQ plot for {}".format(target))
 
    
# %% PLOTS
def plot_data(dfm, df):
    """
    Plots various type of data for exploration purposes.
    
    Remark: Use preferably the following iPython magic commands
    - %matplotlib auto (not qt as it struggles with large figure size) for
    seperate figures.
    - %matplotlib inline for plots in interactive window
    
    Parameters
    ----------
    dfm : pandas DataFrame
        merged data.
    df : pandas DataFrame
        Main measurement data from WISE (might have been modified by 
                                         explore_data()).

    Returns
    -------
    None.

    """
    
    # Sites: histrogram
    plt.figure()
    ax = dfm.monitoringSiteName.value_counts().plot(kind='barh', figsize=(10, 30))
    ax.invert_yaxis()
    ax.set_xlabel('number of measurements')
    ax.set_title('Number of measurements per site')
    
    # Time: distribution of measurements through time
    plt.figure()
    dfm.phenomenonTimeSamplingDate.hist(bins=100)
    
    # Year: histogram ordered by counts
    plt.figure()
    ax = dfm.phenomenonTimeSamplingDate_year.value_counts().plot(kind='barh')
    ax.invert_yaxis()
    ax.set_xlabel('number of measurements')
    ax.set_title('Number of measurements per year')
    
    # Year: histogram ordered chronologically
    plt.figure()
    ax = dfm.phenomenonTimeSamplingDate_year.value_counts().sort_index().plot(kind='bar')
    ax.set_xlabel('number of measurements')
    ax.set_title('Number of measurements per year')
    
    # Water bodies: histogram
    plt.figure()
    ax = dfm.specialisedZoneType.value_counts().plot(kind='barh')
    ax.invert_yaxis()
    ax.set_xlabel('number of measurements')
    ax.set_title('Number of measurements per type of water body')
    plt.tight_layout()
    
    # River basin districts: histogram
    plt.figure()
    ax = dfm.rbdName.value_counts().plot(kind='barh')
    ax.invert_yaxis()
    ax.set_xlabel('number of measurements')
    ax.set_title('Number of measurements per river basin disctrict')
    plt.tight_layout()
    
    # Targets
    plt.figure()
    ax = df.observedPropertyDeterminandLabel.value_counts().plot(kind='barh',figsize=(10, 30))
    ax.invert_yaxis()
    ax.set_xlabel('number of measurements')
    ax.set_title('Number of measurements per target')
    
    plt.figure()
    ax = df.observedPropertyDeterminandLabel.value_counts()[0:10].plot(kind='barh')
    ax.invert_yaxis()
    ax.set_xlabel('number of measurements')
    ax.set_title('Number of measurements for top-10 targets')
    
    # Box plots: multiple targets by type of water body
    # TO DO:
    # - add counts to each boxplot
    # - shorten xticklabels to "ground", "lake", "river" (might do that on dfm directly)
    target_list = df.observedPropertyDeterminandLabel.value_counts().index.tolist()
    fig, axs = plt.subplots(nrows=2, ncols=5)
    axs = axs.ravel()
    t = [None] * 10
    counts = [None] * 10
    for i in range(10):
        t[i] = dfm.loc[dfm.observedPropertyDeterminandLabel == target_list[i], ['resultObservedValue', 'parameterWaterBodyCategory']]
        counts[i] = t[i].parameterWaterBodyCategory.value_counts().sort_index()
        axs[i] = t[i].boxplot(by='parameterWaterBodyCategory', ax=axs[i])
        axs[i].set_title("{}".format(target_list[i]))


# %% MAPS
def maps():
    """
    open map files and plot multiple layers
    
    maps investigated (not used at the moment):   
    # open FOEN river typology map data
    # This map has data on discharge volume, slope, etc.
    geo_df = gpd.read_file("maps/FOEN_river_typology/Typisierung_LV95/FGT.shp")
        
    swiss_rivers_simple = gpd.read_file("maps/FOEN_swiss_water_simple/Gewaesser_500_simplify_80.shp")
    
    # open swiss map (weird)
    swiss_map_10k = gpd.read_file("maps/EEA_switzerland_shapefile/ch_10km.shp") 
    
    # open swiss topo maps
    # from https://www.swisstopo.admin.ch/fr/geodata/landscape/boundaries3d.html
    swiss_canton = gpd.read_file("maps/swissboundaries3d_2021-07_2056_5728.shp_SHAPEFILE_LV95_LN02/swissBOUNDARIES3D_1_3_TLM_KANTONSGEBIET.shp")
    
    # open topo map with water (not bound to swiss border)
    # from https://www.swisstopo.admin.ch/en/geodata/maps/smv/smv1000.html#photo_video
    swiss_topo_rivers = gpd.read_file("maps/SwissTopo_20161001_SMV1000_SHAPE_CHLV95/Shapefiles/21_DKM1M_GEWAESSER_LIN.shp")
    swiss_topo_lakes = gpd.read_file("maps/SwissTopo_20161001_SMV1000_SHAPE_CHLV95/Shapefiles/20_DKM1M_GEWAESSER_PLY.shp")

    Returns
    -------
    site_map : geopandas DataFrame
        locations of all monitoring sites.

    """
    # ------------------
    # Open files
    # ------------------
    
    # open FOEN simple river and lake map data
    swiss_lakes_simple = gpd.read_file("maps/FOEN_swiss_water_simple/See_500.shp")
    swiss_lakes_simple_LV95 = swiss_lakes_simple.geometry.translate(xoff=2000000, yoff=1000000)
    
    # open ecomorphology map
    # This map shows continuous rivers Vs discontinued rivers for the typology map
    # has info on river width, depth, altitude?
    # https://map.geo.admin.ch/?lang=en&topic=bafu&zoom=8.18474730471186&bgLayer=ch.swisstopo.pixelkarte-farbe&catalogNodes=2771,2772,2780,2818&layers=ch.bafu.hydroweb-messstationen_zustand,ch.bafu.hydroweb-messstationen_temperatur,ch.bafu.oekomorphologie-f_abstuerze,ch.bafu.oekomorphologie-f_bauwerke,ch.bafu.oekomorphologie-f_abschnitte&layers_visibility=false,false,false,false,true&E=2608337.25&N=1125282.67
    geo_df = gpd.read_file("maps/FOEN_Ecomorphology/Ökomorphologie_Stufe_F_mit_Geometrie (feature classes)/omch_vec2507_2013.gdb")
    
    swiss_topo_lakes = gpd.read_file("maps/SwissTopo_20161001_SMV1000_SHAPE_CHLV95/Shapefiles/20_DKM1M_GEWAESSER_PLY.shp")
    
   
    # ------------------
    # Plots
    # ------------------
    # TO DO: site map should be created differently (from spatial rather than dfm) to reduce its size significantly
    site_map = gpd.GeoDataFrame(dfm, geometry=gpd.points_from_xy(dfm.swiss_E, dfm.swiss_N))
    site_map_2 = gpd.GeoDataFrame(dfm, geometry=gpd.points_from_xy(dfm.lon, dfm.lat))
    fig, ax = plt.subplots()
    ax = geo_df.plot(zorder=1)
    swiss_lakes_simple.plot(ax=ax, color='green', zorder=2)
    site_map_2.plot(ax=ax,
                  marker='p', 
                  edgecolor='black',
                  c='yellow', 
                  markersize=20,
                  zorder=3)
    
    swiss_topo_lakes.plot()
    
    return site_map


# %% HYPOTESIS TESTING
def stat_test(dfm, target, wbody):
    """
    Runs statistical tests on specified target
    Currently runs D'agostino Pearson normality test

    Parameters
    ----------
    dfm : pandas DataFrame
        Merged data.
    target : string
        Analytical target used for testing
    wbody: 'string'
        Water Body Type

    Returns
    -------
    None.

    """
    # Normality test: D'agostino Pearson
    x = dfm.loc[(dfm.observedPropertyDeterminandLabel == target)  & 
                (dfm.resultQualityObservedValueBelowLOQ == 0) &
                (dfm.parameterWaterBodyCategory == wbody), [
                      'resultObservedValue']
                  ].dropna()
    k2, p = stats.normaltest(x) # D Agostino-Pearson. The method returns the test statistic value and the p-value
    alpha = 0.001 # Rejection criterion defined by you
    print("Normality test (D'agostino Pearson) for {} in {}".format(target, wbody))
    print('Alpha = ',alpha)
    print('p = ',p)
    if p < alpha:  # null hypothesis: x comes from a normal distribution
         print("The null hypothesis can be rejected")
    else:
      print("The null hypothesis cannot be rejected")
    print("For info, the tiniest float64 number is {}".format(np.finfo(np.float64).tiny))


# %% MACHINE LEARNING
'''
Workflow:
    Drop targets that are not common to all categories
    Group by monitoring sites, average values
    Scale data
    Replace NAs by zeros
    TBD: rebalance dataset
    Split training and test datasets
        50/50 split
        Keep proportionnality of labelled data in the train and test data
        TBD: cross-validation
    Run classification model
        Linear SVC
        TBD: K-neighbors
    Results
        Confusion matrix
        Scores
'''

def prep_data(dfm):
    """
    Prepares data for machine learning.
    
    TBD:
        Parameters will be added to configure the variations of data prep 
        workflows

    Parameters
    -------
    dfm : pandas DataFrame
        merged data.

    Returns
    -------
    data : pandas DataFrame
        numerical tidy data ready for machine learning (sites x targets)
    train_perc : numeric
        percentage of train vs test data for supervised learning.
    x_train : pandas DataFrame
        split data to be used for model training.
    x_test : pandas DataFrame
        split data to be used for model testing.
    y_train : pandas Index
        output (usually labels) for training the model.
    y_test : pandas Index
        output (usually labels) for testing the model.

    """      
    # ------------------------------
    # preliminary investigation
    # ------------------------------
    dfa = dfm.copy()
    
    # drop empty results 
    dfa.dropna(subset = ["resultObservedValue"], inplace=True)
    target_perWB_counts = (dfa
                        .groupby(by="parameterWaterBodyCategory")
                        .observedPropertyDeterminandLabel
                        .value_counts()
                        .unstack(0)
                        .sort_values(by=['lake', 'ground', 'river'], ascending=False)
                        .reindex(columns=['lake', 'ground', 'river'])
                        )
    target_perRBD_counts = (dfa
                        .groupby(by="rbdName")
                        .observedPropertyDeterminandLabel
                        .value_counts()
                        .unstack(0)
                        .sort_values(by=['DANUBE', 'PO', 'RHONE', 'RHINE'], ascending=False)
                        .reindex(columns=['DANUBE', 'PO', 'RHONE', 'RHINE'])
                        )
    
    # ------------------------------
    # Data prep
    # ------------------------------
    # # keep only essential columns
    # dfa['month'] = dfm.phenomenonTimeSamplingDate.dt.month
    # cols_to_keep = ['monitoringSiteIdentifier', 'parameterWaterBodyCategory',
    #        'observedPropertyDeterminandCode', 'observedPropertyDeterminandLabel',
    #        'resultUom', 
    #        'phenomenonTimeSamplingDate_year', 'resultObservedValue',
    #        'resultQualityObservedValueBelowLOQ',
    #        'monitoringSiteName', 'waterBodyName',
    #        'specialisedZoneType', 'rbdName', 'month']
    # dfa = dfa[cols_to_keep]
    
    # Selection of categories
    # dfa = dfa[dfa.parameterWaterBodyCategory != "lake"]
    
    # Selection of analytical targets
    '''
    There are several options:
        - Option A: select targets that have data (measured values) on all
        categories to classify => keeps ~max number of targets
            + drop NAs or not? => comparison to be done
            At the moment, I am more in favour of removing targets that are not
            at all measured for one or more of the categories. It makes it more
            fair when I susbsequently impute NAs with 0.
        - Option B: select specific (e.g. most promising) targets based on
        analysis/exploration
    '''
    # Option A - select top targets (many sub-options...)
    # targets = target_counts.iloc[:19].index
    # targets = target_perWB_counts.dropna().index
    targets = target_perRBD_counts.dropna().index
    # Option B - select specific targets
    # targets = ["Dissolved oxygen", "Electrical conductivity", "Total phosphorus", "pH"]
    
    target_nb = len(targets)
    dfa = dfa[dfa.observedPropertyDeterminandLabel.isin(targets)]
    
    # Aggregate data
    dfa_g = dfa.groupby(by= ["monitoringSiteIdentifier", 
                             "observedPropertyDeterminandLabel", 
                             "parameterWaterBodyCategory", 
                             "rbdName"]).mean()
    dfa_g = dfa_g["resultObservedValue"]
    dfa_g = dfa_g.unstack(1)
    
    # Scale data
    '''
    - Necessary for training supervised models.
    Without scaling, the LinearSVC method fails to converge.
    - Scaling is done before filling NAs.
    '''
    data = StandardScaler().fit_transform(dfa_g.values)
    data = pd.DataFrame(data, index=dfa_g.index, columns=dfa_g.columns)
       
    # Manage Na values / Imputation
    '''
    Keeping all targets leads to a lot of NAs:
        -Option A: regouping per WB types (with large target selection), filtering
        out all sites with at least one Na means means retaining only 33 sites 
        (32 GW and 1 RW) => useless.
        - One alternative is replacing Na with 0
            + and then optionnaly, filter targets using various feature
            selection methods (see https://scikit-learn.org/stable/modules/feature_selection.html)
    '''
    # Option A: Drop NA
    # dfa_g = dfa_g.dropna()
    # Option B: impute and select features 
    data.fillna(value=0, inplace=True)
    
    # Optional: filter targets using a feature selection method
    '''
    filtering based on variance introduces can filter out targets...:
        - with lots of NAs
        - with lots of measurements at LOQ
    '''
    # sel = VarianceThreshold(threshold=(0))
    # var = dfa_g.var(axis='rows')
    # test = sel.fit_transform(dfa_g)
    
    # # Optional plot
    # dfa_g_plot = dfa_g.copy()
    # dfa_g_plot['WBtype'] = dfa_g_plot.index.get_level_values(1)
    # dfa_g_plot['RBdistrict'] = dfa_g_plot.index.get_level_values(2)
    # # sns.relplot(data=dfa_g, x='pH', y='Electrical conductivity', hue="parameterWaterBodyCategory")
    # sns.pairplot(data=dfa_g_plot, hue="WBtype")
    # sns.pairplot(data=dfa_g_plot, hue="RBdistrict")
    
    # Manage imbalanced Datasets
    '''
    rebalanced datasets to even (or reduce imbalance) representation of categories
    '''
    # how balanced are the categories?
    dist_WB = data.index.get_level_values('parameterWaterBodyCategory').value_counts()
    dist_RBD = data.index.get_level_values('rbdName').value_counts()
    
    # print(Counter(dist_WB))
    # print(Counter(dist_RBD))
    
    plt.figure()
    ax = dist_WB.plot(kind='barh')
    ax.invert_yaxis()
    ax.set_xlabel('number of monitoring sites')
    ax.set_title('Number of sites per category')
    
    plt.figure()
    ax = dist_RBD.plot(kind='barh')
    ax.invert_yaxis()
    ax.set_xlabel('number of monitoring sites')
    ax.set_title('Number of sites per category')
    
    # define categories to build model on
    labels = data.index.get_level_values('parameterWaterBodyCategory')#.to_numpy()
    
    # # mix of over- and under-sampling
    # '''!! try converting labels to integers to troubleshoot fit_resample error on random_state!!
    # '''
    # data_orig = data.copy()
    # over = RandomOverSampler(sampling_strategy='auto')
    # data, labels = over.fit_resample(data, labels)
    # # TBD: reconstruct dataframe using over.sample_indices_
    
    # Split train and test datasets
    train_perc = 0.5
    x_train, x_test, y_train, y_test = train_test_split(data, labels, train_size=train_perc, stratify=labels)
    
    return data, train_perc, x_train, x_test, y_train, y_test

def classify_data(data, train_perc: float, x_train, x_test, y_train, y_test):
    """
    Runs classifier method (supervised learning). At the moment the method
    used is LinearSVC
    
    TBD:
        - options for different classification methods (LinearSVC, K-neighbors,
          etc.)
        - plot resembling a pair plot with seperation vectors from svc
        - run cross validation (might have to do it in workflow and not 
          directly in function)

    Parameters
    ----------
    data : pandas DataFrame
        numerical tidy data ready for machine learning (sites x targets)
    train_perc : float
        percentage of train vs test data for supervised learning.
    x_train : pandas DataFrame
        split data to be used for model training.
    x_test : pandas DataFrame
        split data to be used for model testing.
    y_train : pandas Index
        output (usually labels) for training the model.
    y_test : pandas Index
        output (usually labels) for testing the model.

    Returns
    -------
    y_pred : numpy.ndarray
        vetor of predicted labels.
    """
    classifier = svm.LinearSVC()
    classifier.fit(x_train, y_train)
    y_pred = classifier.predict(x_test)
    
    cm = metrics.confusion_matrix(y_test, y_pred, labels=classifier.classes_)
    cm_disp = metrics.ConfusionMatrixDisplay(cm, display_labels=classifier.classes_)
    cm_disp.plot(cmap='Greens')
    plt.title('Confusion Matrix, train/test split = {:n}/{:n}'.format(train_perc*100, (1-train_perc)*100))
    
    scores = metrics.classification_report(y_test, y_pred, target_names=classifier.classes_)
    print(scores)
    
    return y_pred

def cluster_data(data, method: 'str'):
    """
    Run cluster analysis (unsupervised learning). The methods ''gaussian' and 'k-means' run 
    Gaussian_mixtures and KMeans for various number of cluster (default is 2 
    to 15 clusters) and compute scores (silhouette, aic and bic (gaussian
    mixtures only).

    Parameters
    ----------
    data : pandas DataFrame
        numerical data to be clustered. Should be the output of prep_data().
    method : 'str'
        Clustering method to be used ('custom', 'gaussian', 'k-means')
        - 'custom' applies a method with specificied number of clusters.
        - 'gaussian' and 'k-means' run analysis for and compute scores for max 
        15 clusters.

    Returns
    -------
    cluster_labels : 'numpy.ndarray of int32'
        1D numpy array of integer corresponding to the cluster number for reach data input.

    """
    # initialisation of variables
    aic_list = []
    bic_list = []
    sil_k = []
    sil_g =[]
    max_cl = 15
    
    # Kmeans
    if method == 'k-means':
        for i in range(2, max_cl):
            clusterer = KMeans(n_clusters=i)
            clusterer.fit(data)
            cluster_centers = clusterer.cluster_centers_
            cluster_labels = clusterer.predict(data)
            score = metrics.silhouette_score(data, cluster_labels)
            sil_k.append(score)
            
        plt.plot(range(2, max_cl), sil_k, '-o')
        plt.xlabel('number of clusters')
        plt.ylabel('score')
        plt.title('silhouette - Kmeans')
        plt.show()
    
    # Gaussian mixtures
    elif method == 'gaussian':
        for i in range(2, max_cl):
            clusterer = GaussianMixture(n_components=i)
            clusterer.fit(data)
            cluster_labels = clusterer.predict(data)
            score = metrics.silhouette_score(data, cluster_labels)
            aic_list.append(clusterer.aic(data))
            bic_list.append(clusterer.bic(data))
            sil_g.append(score)
        
        plt.plot(range(2, max_cl), sil_g, '-o')
        plt.xlabel('number of clusters')
        plt.ylabel('score')
        plt.title('silhouette - Gaussian Mixture')
        plt.show()
        
        plt.plot(range(2, max_cl), aic_list, '-o')
        plt.xlabel('number of clusters')
        plt.ylabel('score')
        plt.title('aic - Gaussian Mixture')
        plt.show()
        
        plt.plot(range(2, max_cl), bic_list, '-o')
        plt.xlabel('number of clusters')
        plt.ylabel('score')
        plt.title('bic - Gaussian Mixture')
        plt.show()
    
    # Clustering to plot
    elif method == 'custom':
        clusterer = KMeans(n_clusters=5)
        clusterer.fit(data)
        cluster_labels = clusterer.predict(data)
        
        return cluster_labels

def proj_data(data: 'pd.DataFrame or np.array', target_nb, cluster_labels: 'np.array'):
    """
    Run dimension reduction methods (UMAP and PCA) and plot results.
    
    TBD:
        - Make use of cluster_labels optional, i.e. only when available.

    Parameters
    ----------
    data : 'pd.DataFrame or np.array'
        Data to be analysed.
    cluster_labels : 'np.array'
        labels coming from classifier analysis.

    Returns
    -------
    None.

    """
    # UMAP projection
    '''
    UMAP can randomly output very different projections
        + use random generator to see if results are indeed consistent
    '''
    umap_model = umap.UMAP()
    embedding = umap_model.fit_transform(data)
    embedding.shape   
    
    plt.figure()
    plt.scatter(
        embedding[:, 0],
        embedding[:, 1],
        c=cluster_labels)
    plt.title('UMAP projection of the top-{} targets, K-means cluster labels'.format(target_nb), fontsize=12)
    
    plt.figure()
    sns.relplot(
        embedding[:, 0],
        embedding[:, 1],
        hue=data.index.get_level_values('parameterWaterBodyCategory'),
        style=data.index.get_level_values('rbdName'),
        s=40,
        alpha=0.8,
        palette='muted')
    plt.title('UMAP projection of the top-{} targets, real label'.format(target_nb), fontsize=12)
    
    # PCA
    pca = PCA()
    pca.fit(data)
    data_pcaed = pca.transform(data)
    
    plt.plot(pca.explained_variance_ratio_, '-o')
    plt.ylabel('percentage of explained variance')
    plt.xlabel('component number')
    plt.title('Scree plot')
    plt.show()
    
    sns.scatterplot(data_pcaed[:,0],data_pcaed[:,1], hue=cluster_labels, palette='bright')
    plt.xlabel('First component')
    plt.ylabel('Second component')
    plt.title('PCA - K-means clusters')
    plt.show()
    
    sns.relplot(data_pcaed[:,0],
                data_pcaed[:,1],
                hue=data.index.get_level_values('parameterWaterBodyCategory'),
                style=data.index.get_level_values('rbdName'),
                s=40,
                alpha=0.8,
                palette='muted')
    plt.xlabel('First component')
    plt.ylabel('Second component')
    plt.title('PCA- real labels')
    plt.show()
    
if __name__ == "__main__":
    path = "D:\Ludo\Docs\programming\CAS_applied_data_science\project_Water\Datasets".replace(
        "\\", "/")
    # path = r"C:\Users\ludovic.lereste\Documents\CAS_applied_data_science\project_Water\Datasets".replace("\\", "/")
    df, spatial, df_def, spatial_def = load_files(path)
    dfm, df, spatial = merge_data(df, spatial)
    data, train_perc, x_train, x_test, y_train, y_test = prep_data(dfm)
    y_pred = classify_data(data, train_perc, x_train, x_test, y_train, y_test)
    cluster_labels = cluster_data(data, method='custom')
    
