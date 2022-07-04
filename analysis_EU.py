"""
Created on 07 Sep 2021

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
"""
__author__ = "Ludovic Le Reste"
__credits__ = ["Hisham Ben Hamidane", "Ludovic Le Reste"]
__status__ = "Prototype"

if __name__ == "__main__":   
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
    
    
    # %% LOAD FILES - read csv
    # path = "D:\Ludo\Docs\programming\CAS_applied_data_science\project_Water\Datasets".replace(
    #     "\\", "/")
    path = r"C:\Users\ludovic.lereste\Documents\CAS_applied_data_science\project_Water\Datasets".replace("\\", "/")
    os.chdir(path)
    
    # WISE tabular data
    # filter columns and specify data types to reduce file size from 15.3 GB to 1.1 GB
    cols_kept = ["monitoringSiteIdentifier",
                  "parameterWaterBodyCategory",
                  "observedPropertyDeterminandLabel",
                  "resultUom",
                  "phenomenonTimeSamplingDate",
                  "resultObservedValue",
                  "resultQualityObservedValueBelowLOQ",
                  "parameterSpecies",
                  "resultObservationStatus"]
    data_types = {"monitoringSiteIdentifier": "category",
                  "parameterWaterBodyCategory": "category",
                  "observedPropertyDeterminandLabel": "category",
                  "resultUom": "category",
                  "phenomenonTimeSamplingDate": "int32",
                  "resultObservedValue": "float32",
                  "resultQualityObservedValueBelowLOQ": "boolean",
                  "parameterSpecies": "category",
                  "resultObservationStatus": "category"}
    
    # df = pd.read_csv("WISE/Waterbase_v2021_1_T_WISE6_DisaggregatedData.csv", usecols=cols_kept, dtype=data_types)
    # df.to_pickle("WISE/Data_EU_disaggregated_colFiltered.pkl")

    spatial = pd.read_csv("WISE/Waterbase_v2021_1_S_WISE6_SpatialObject_DerivedData.csv")

    # %% LOAD FILES - read pickle
    # path = "D:\Ludo\Docs\programming\CAS_applied_data_science\project_Water\Datasets".replace(
    #     "\\", "/")
    path = r"C:\Users\ludovic.lereste\Documents\CAS_applied_data_science\project_Water\Datasets".replace("\\", "/")
    os.chdir(path)
    
    df = pd.read_pickle("WISE/Data_EU_disaggregated_colFiltered.pkl")
    
    
    # %% DATA PREP
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
    
    # Are there NaNs in spatial? => Yes
    spatial.monitoringSiteIdentifier.isnull().sum()
    # Are there duplicates in spatial sites? => Yes
    # due to different monitoringSiteIdentifierScheme (the euMonitoringSiteCode has lat and lon, not the other)
    sum(spatial.monitoringSiteIdentifier.duplicated(), )
    mask = spatial.monitoringSiteIdentifier.duplicated()
    sites_duplicated = spatial.monitoringSiteIdentifier.dropna()[mask].reset_index(drop=True)
    
    
    # %% DATA EXPLORATION
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
    
    
    # %% PLOTS
    # Use iPython magic §matplotlib auto (not qt as it struggles with large figure size) for seperate figures
    # %matplotlib inline for plots in interactive window
    
    # Use my plotting custom style
    # plot_style.custom_style()
    
    # Sites: histrogram
    plt.figure()
    ax = site_name_counts.plot(kind='barh', figsize=(10, 30))
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
    target_list = target_counts.index.tolist()
    fig, axs = plt.subplots(nrows=2, ncols=5)
    axs = axs.ravel()
    t = [None] * 10
    counts = [None] * 10
    for i in range(10):
        t[i] = dfm.loc[dfm.observedPropertyDeterminandLabel == target_list[i], ['resultObservedValue', 'parameterWaterBodyCategory']]
        counts[i] = t[i].parameterWaterBodyCategory.value_counts().sort_index()
        axs[i] = t[i].boxplot(by='parameterWaterBodyCategory', ax=axs[i])
        axs[i].set_title("{}".format(target_list[i]))
    
    # %% ANALYSIS PER TARGET
    
    target = "pH"
    bins = {"Nitrate": np.linspace(0, 60, 100),
            "pH": np.linspace(5, 10, 50),
            "Electrical conductivity": np.linspace(0, 1200, 100)}
    # # Histograms: value distibution + time distribution
    # axs = dfm.loc[(dfm.observedPropertyDeterminandLabel == target) 
    #               & (dfm.resultQualityObservedValueBelowLOQ == 0), [
    #                   'resultObservedValue', 'phenomenonTimeSamplingDate']
    #               ].hist(bins=100)
    # axs[0][0].set_xlabel('{} in {}'.format(target, target_uom.loc[target][0]))
    # axs[0][0].set_ylabel('counts')
    # axs[0][1].set_xlabel('year')
    # axs[0][1].set_ylabel('counts')
    # plt.suptitle(target)
    
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
    par = [{"wbody": "GW", "color": "brown", "alpha": 0.5, "label": "Ground Water (GW)"},
           {"wbody": "LW", "color": "green", "alpha": 0.8, "label": "Lake Water (LW)"},
           {"wbody": "RW", "color": "blue", "alpha": 0.3, "label": "River Water (RW)"}]
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
    
    # %% HYPOTESIS TESTING
    # Normality test: D'agostino Pearson
    wbody = "LW"
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
    # ------------------------------
    # preliminary investigation
    # ------------------------------
    dfa = dfm.copy()
    
    # Rename river body types
    dfa.loc[df.parameterWaterBodyCategory == 'GW', 'parameterWaterBodyCategory'] = 'ground'
    dfa.loc[df.parameterWaterBodyCategory == 'RW', 'parameterWaterBodyCategory'] = 'river'
    dfa.loc[df.parameterWaterBodyCategory == 'LW', 'parameterWaterBodyCategory'] = 'lake'
    
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
    # dfa = dfa[dfa.parameterWaterBodyCategory != "LW"]
    
    
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
    targets = target_perWB_counts.dropna().index
    # targets = target_perRBD_counts.dropna().index
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
    
    print(Counter(dist_WB))
    print(Counter(dist_RBD))
    
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
    
    # mix of over- and under-sampling
    '''!! try converting labels to integers to troubleshoot fit_resample error on random_state!!
    '''
    data_orig = data.copy()
    over = RandomOverSampler(sampling_strategy='auto')
    data, labels = over.fit_resample(data, labels)
    # TBD: reconstruct dataframe using over.sample_indices_
    
    # Split train and test datasets
    train_perc = 0.5
    x_train, x_test, y_train, y_test = train_test_split(data, labels, train_size=train_perc, stratify=labels)
    
    # ------------------------------
    # Classification (supervised)
    # ------------------------------
    '''
    TO DO:
        - get scores of fit/confusion matrix
        - plot like pair plot with seperation vectors from svc
        - repeat exercise with more/less targets (depending on feature selection)
        - run cross validation
    '''
    # x_train = data.sample(frac=0.8).sort_index(level=1)
    # y_train = x_train.index.get_level_values(1)
    # x_train = StandardScaler().fit_transform(x)
    
    # classifier = svm.LinearSVC()
    # classifier.fit(x, y)
    # x_test = StandardScaler().fit_transform(data)
    # y_test = classifier.predict(x_test)
    
    classifier = svm.LinearSVC()
    classifier.fit(x_train, y_train)
    y_pred = classifier.predict(x_test)
    
    cm = metrics.confusion_matrix(y_test, y_pred, labels=classifier.classes_)
    cm_disp = metrics.ConfusionMatrixDisplay(cm, display_labels=classifier.classes_)
    cm_disp.plot(cmap='Greens')
    plt.title('Confusion Matrix, train/test split = {:n}/{:n}'.format(train_perc*100, (1-train_perc)*100))
    
    scores = metrics.classification_report(y_test, y_pred, target_names=classifier.classes_)
    print(scores)
    
    # ------------------------------
    # Clustering (unsupervised)
    # ------------------------------
    
    
    aic_list = []
    bic_list = []
    sil_k = []
    sil_g =[]
    max_cl = 15
    
    # Kmeans
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
    clusterer = KMeans(n_clusters=5)
    clusterer.fit(data)
    cluster_labels = clusterer.predict(data)
    
    # ------------------------------
    # 2D projections
    # ------------------------------
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
