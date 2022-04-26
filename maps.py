# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 13:33:52 2022

@author: ludovic.lereste

This scripts plots interactive map(s) of waters
"""
__author__ = "Ludovic Le Reste"
__credits__ = ["Hisham Ben Hamidane", "Ludovic Le Reste"]
__status__ = "Prototype"

import os
import pandas as pd
import geopandas as gpd
import sys
import matplotlib.pyplot as plt
    
if __name__ == "__main__":
    # path = "D:\Ludo\Docs\programming\CAS_applied_data_science\project_Water\Datasets".replace(
    #     "\\", "/")
    path = r"C:\Users\ludovic.lereste\Documents\CAS_applied_data_science\project_Water\Datasets".replace("\\", "/")
    os.chdir(path)
    
    spatial = pd.read_csv("WISE/waterbase_s_wise_spatialobject_deriveddata_CH.csv")
    # open ecomorphology map
    # This map shows continuous rivers Vs discontinued rivers for the typology map
    # has info on river width, depth, altitude?
    # https://map.geo.admin.ch/?lang=en&topic=bafu&zoom=8.18474730471186&bgLayer=ch.swisstopo.pixelkarte-farbe&catalogNodes=2771,2772,2780,2818&layers=ch.bafu.hydroweb-messstationen_zustand,ch.bafu.hydroweb-messstationen_temperatur,ch.bafu.oekomorphologie-f_abstuerze,ch.bafu.oekomorphologie-f_bauwerke,ch.bafu.oekomorphologie-f_abschnitte&layers_visibility=false,false,false,false,true&E=2608337.25&N=1125282.67
    rivers = gpd.read_file("maps/FOEN_Ecomorphology/Ökomorphologie_Stufe_F_mit_Geometrie (feature classes)/omch_vec2507_2013.gdb")
    lakes = gpd.read_file("maps/FOEN_swiss_water_simple/See_500.shp")
    sites = gpd.GeoDataFrame(spatial, geometry=gpd.points_from_xy(spatial.lon, spatial.lat), crs=4326)
    sites = sites.to_crs(crs=21781)
    
    fig, ax = plt.subplots()
    rivers.plot(ax=ax, zorder=1)
    lakes.plot(ax=ax, color='green', zorder=2)
    sites.plot(ax=ax,
                  marker='p', 
                  edgecolor='black',
                  c='yellow', 
                  markersize=20,
                  zorder=3)
