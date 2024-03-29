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
# import geopandas as gpd
# import sys
import matplotlib.pyplot as plt
# import folium

class Human:
    """Machine builders :-)."""

    race = "Homo sapiens sapiens"

    def __init__(self, name, age, energy):
        # super().__init__(**kwargs)
        self.name = name
        self.age = age
        self.energy = energy

    # Methods
    def grow(self, years):
        """
        Add years to age of a human.

        Parameters
        ----------
        years : years alive

            That would be the age.

        Returns
        -------
        age

            This the age the human instance (i.e. a person).
        """
        self.age = self.age + years
        return self.age
    
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
    center = sites[['lon', 'lat']].median()
    fmap = folium.Map(location=center)
    
    # site_
    fmap.save("./maps/map.html")
