# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 13:33:52 2022

@author: ludovic.lereste

This scripts plots interactive map(s) of waters

TO DO:
- check map coordinate conversion with function to_crs() instead of manual transformation
- use ecomorphology map (continuous rivers) instead of typology map
    + troubleshoot LV03 Vs LV95
"""
__author__ = "Ludovic Le Reste"
__credits__ = ["Hisham Ben Hamidane", "Ludovic Le Reste"]
__status__ = "Prototype"

import os
import pandas as pd
import geopandas as gpd
import sys
import matplotlib.pyplot as plt
import folium
    
if __name__ == "__main__":
    # path = "D:\Ludo\Docs\programming\CAS_applied_data_science\project_Water\Datasets".replace(
    #     "\\", "/")
    path = r"C:\Users\ludovic.lereste\Documents\CAS_applied_data_science\project_Water\Datasets".replace("\\", "/")
    os.chdir(path)
    
    spatial = pd.read_csv("WISE/Waterbase_v2021_1_S_WISE6_SpatialObject_DerivedData.csv")
    sites = gpd.GeoDataFrame(spatial, geometry=gpd.points_from_xy(spatial.lon, spatial.lat), crs='EPSG:4326')
    sites = sites.to_crs(crs=21781)
    
    fig, ax = plt.subplots()
    sites.plot(ax=ax,
                  marker='p', 
                  edgecolor='black',
                  c='yellow', 
                  markersize=20,
                  zorder=3)
    center = sites[['lon', 'lat']].mean()
    fmap = folium.Map(location=center)
    fmap.save("./maps/map.html")
