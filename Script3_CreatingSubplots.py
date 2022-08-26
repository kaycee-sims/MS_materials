#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 12:20:19 2022

@author: kayceeschaper
"""

#Goal is to make a subroutine for plotting clusters/lineations as subplots


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import datetime
from obspy import read_events
from mpl_toolkits.basemap import Basemap
from easyQuake import simple_cat_df                                                            
from obspy.core import UTCDateTime
from sklearn.cluster import DBSCAN
from sklearn import metrics

from sklearn import linear_model, datasets
from obspy.clients.fdsn import Client

#Loading Catalogs#
    # easyQuake Catalog_2010 
cat_eQuake= read_events('/Users/kayceeschaper/Earthquake_Catalogs/catalog_2010_mag.xml')     
catdf_eQuake=simple_cat_df(cat_eQuake )

#catdf_eQuake['origintime'] = catdf_eQuake.index

catdf_eQuake['origintime'] = pd.to_datetime(catdf_eQuake['origintime'])
catdf_eQuake['magnitude'] = catdf_eQuake['magnitude'].astype(float)

catdf_eQuake_Reloc = pd.read_csv('/Users/kayceeschaper/HYPODD/hypodd.relocpanfixed_KSBest', header=None, delimiter=r"\s+") #delimiter=r"\s+" means working with txt file/using spaces instead of commas
catdf_eQuake_Reloc['latitude']=catdf_eQuake_Reloc[1]
catdf_eQuake_Reloc['longitude']=catdf_eQuake_Reloc[2]


    # OGS Catalog 
catdf_OGS = pd.read_csv('http://wichita.ogs.ou.edu/eq/catalog/complete/complete.csv')
catdf_OGS['origintime'] = pd.to_datetime(catdf_OGS['origintime'])
catdf_OGS['magnitude'] = catdf_OGS['magnitude'].str.replace("None", "0").astype(float)
catdf_OGS['magnitude'] = catdf_OGS['magnitude'].astype(float)
catdf_OGS['depth'] = catdf_OGS['depth_km'].astype(float)

#to limit the data to only 2010
catdf_OGS2 = catdf_OGS[(catdf_OGS['origintime']>'2010-01-01') & 
                         (catdf_OGS['origintime']<'2011-01-01')]

catdf_OGS_Reloc= pd.read_csv('/Users/kayceeschaper/HYPODD/OK2010/hypoDD.relocpanfixed_2010', header=None, delimiter=r"\s+")
catdf_OGS_Reloc['latitude']=catdf_OGS_Reloc[1]
catdf_OGS_Reloc['longitude']=catdf_OGS_Reloc[2]

#Fetching the desired catalog from USGS via Obspy Client
starttime= UTCDateTime("2010-01-01")
endtime = UTCDateTime("2011-01-01")
client = Client("USGS")
cat_USGS2010 = client.get_events(starttime=starttime, endtime=endtime, minlatitude=33.,
                         maxlatitude=37.5, minlongitude=-104, maxlongitude=-94, 
                         minmagnitude=1)
catdf_USGS2010=simple_cat_df(cat_USGS2010)
catdf_USGS2010['origintime'] = catdf_USGS2010.index

#Keranen Catalog
catdf_Ker=pd.read_excel('/Users/kayceeschaper/Earthquake_Catalogs/Keranen_catalog-3.xlsx')

date1 = []
for date in catdf_Ker['Date']:
    date1.append((datetime.datetime(int(str(date)[0:2])+2000,int(str(date)[2:4]),
                                    int(str(date)[4:6]),0,0,0)))
catdf_Ker['origintime']=pd.to_datetime(date1)
catdf_Ker['latitude']=catdf_Ker.Latitude
catdf_Ker['longitude']=catdf_Ker.Longitude

catdf_Ker2 = catdf_Ker[(catdf_Ker['origintime']>'2010-01-01') & (catdf_Ker['origintime']<'2011-01-01')]


#Limiting to Jones
#Injection site Coord 35.43,-97.46
#Keranen Jones zone 35.4-35.8, -97.5- -96.75

#easyQuake submodule; creates new dataframe geographic limitations
def catdf_narrowbounds(catdf=None,lat_a=None,lat_b=None,lon_a=None,lon_b=None):
    catdf = catdf[(catdf['latitude']>lat_a) & (catdf['latitude']<lat_b) & (catdf['longitude']>lon_a) & (catdf['longitude']<lon_b)]
    return catdf

catdf_eQuakeNarrow=catdf_narrowbounds(catdf=catdf_eQuake, lat_a=35.4, lat_b=35.8, lon_a=-97.5, lon_b=-96.75)
catdf_eQuake_RelocNarrow=catdf_narrowbounds(catdf=catdf_eQuake_Reloc, lat_a=35.4, lat_b=35.8, lon_a=-97.5, lon_b=-96.75)

catdf_OGS2Narrow=catdf_narrowbounds(catdf=catdf_OGS2, lat_a=35.4, lat_b=35.8, lon_a=-97.5, lon_b=-96.75)
catdf_USGSNarrow=catdf_narrowbounds(catdf=catdf_USGS2010, lat_a=35.4, lat_b=35.8, lon_a=-97.5, lon_b=-96.75)
catdf_Ker2Narrow=catdf_narrowbounds(catdf=catdf_Ker2, lat_a=35.4, lat_b=35.8, lon_a=-97.5, lon_b=-96.75)
catdf_OGS_RelocNarrow=catdf_narrowbounds(catdf=catdf_OGS_Reloc, lat_a=35.4, lat_b=35.8, lon_a=-97.5, lon_b=-96.75)


def lat_lon_clusters(catdf=None,sub1=None,sub2=None):
        
    #easyQuake
    ##with geopgraphic bounds
    #adding injection site location and stations to map
    #Station locations for USGS catalog 
    client = Client("IRIS")
    starttime= UTCDateTime("2010-01-01")
    endtime = UTCDateTime("2010-07-01")
    lat1=35.4
    lat2=35.8
    lon1=-96.75
    lon2=-97.5
    inv = client.get_stations(minlatitude=lat1, minlongitude=lon2 ,maxlatitude=lat2, maxlongitude=lon1,
                                    starttime=starttime,
                                    endtime=endtime)
    latitudes = []
    longitudes = []
    for net in inv:
        for sta in net:
            latitudes.append(sta.latitude)
            longitudes.append(sta.longitude)
    
    map = Basemap(llcrnrlon=-104,llcrnrlat=33,urcrnrlon=-93,urcrnrlat=38, resolution='l', projection='tmerc',lat_0=35,lon_0=-98)
    x,y = map(np.array((catdf.longitude)),np.array((catdf.latitude)))
    injlon,injlat=-97.5,35.4
    xinj,yinj = map(np.array(injlon),np.array(injlat)) #injection site coordinates: 35.4N -97.5W
    xinv,yinv = map(np.array(longitudes),np.array(latitudes))
    X=np.transpose(np.array((x,y)))
    #depth = -np.array(catdf.depth)*1000
    
    #plt.figure()
    #np.concatenate((a, b.T), axis=1)
    # Compute DBSCAN
    db = DBSCAN(eps=800, min_samples=5).fit(X)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_
    
    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)
    print('Estimated number of clusters: %d' % n_clusters_)
    print('Estimated number of noise points: %d' % n_noise_)
    
    # Black removed and is used for noise instead.
    unique_labels = set(labels)
    colors = [plt.cm.Spectral(each)
              for each in np.linspace(0, 1, len(unique_labels))]
    
    for k, col in zip(unique_labels, colors):
        #print(k)
        if k == -1:
            # Black used for noise.
            col = [0, 0, 0, 1]
    
        class_member_mask = (labels == k)
        xy = X[class_member_mask & core_samples_mask]
    
    
       # plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
        #         markeredgecolor='k', markersize=14) 
    
        xy = X[class_member_mask & ~core_samples_mask]
        #plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
         #        markeredgecolor='k', markersize=6)
    
    #plt.scatter(xinj, yinj, marker='^', color='r', label='Injection Site 35.43N, -97.46W')
    #plt.scatter(xinv, yinv, marker='v', color='g', label='Stations')
    
    #plt.legend()
    #plt.title('eps=800; min_samples=5 - Estimated number of clusters: %d' % n_clusters_)
    
    
    #Injection site Coord 35.43,-97.46
    #Keranen Jones zone 35.4-35.8, -97.5- -96.75
    
    #Back to lat/long
    #plt.figure()
    
    starttime= UTCDateTime("2010-01-01")
    endtime = UTCDateTime("2010-07-01")
    lat1=35.4
    lat2=35.8
    lon1=-96.75
    lon2=-97.5
    inv = client.get_stations(minlatitude=lat1, minlongitude=lon2 ,maxlatitude=lat2, maxlongitude=lon1,
                                    starttime=starttime,
                                    endtime=endtime)
    
    latitudes = []
    longitudes = []
    for net in inv:
        for sta in net:
            latitudes.append(sta.latitude)
            longitudes.append(sta.longitude)
    
    map = Basemap(llcrnrlon=lon1,llcrnrlat=lat1,urcrnrlon=lon2,urcrnrlat=lat2, resolution='l', projection='tmerc',lat_0=35.5,lon_0=-97)
    x,y = (np.array((catdf['longitude'])),np.array((catdf['latitude'])))
    injlon,injlat=-97.5,35.4
    #xinj,yinj = map(np.array(injlon),np.array(injlat), inverse=True) #injection site coordinates: 35.4N -97.5W
    xinv,yinv = (np.array(longitudes),np.array(latitudes))
    X=np.transpose(np.array((x,y)))
    #depth = -np.array(catdf.depth)*1000
    #plt.plot(lonx,latx)
    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)
    print('Estimated number of clusters: %d' % n_clusters_)
    print('Estimated number of noise points: %d' % n_noise_)
    
    # Black removed and is used for noise instead.
    unique_labels = set(labels)
    colors = [plt.cm.Spectral(each)
              for each in np.linspace(0, 1, len(unique_labels))]
    
    for k, col in zip(unique_labels, colors):
        #print(k)
        if k == -1:
            # Black used for noise.
            col = [0, 0, 0, 1]
    
        class_member_mask = (labels == k)
        xy = X[class_member_mask & core_samples_mask]
    
    
        axs[sub1,sub2].plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
                 markeredgecolor='k', markersize=14) 
    
        xy = X[class_member_mask & ~core_samples_mask]
        axs[sub1,sub2].plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
                 markeredgecolor='k', markersize=6)
    
    axs[sub1,sub2].scatter(injlon, injlat, marker='^', color='r', label='Injection Site 35.43N, -97.46W')
    axs[sub1,sub2].scatter(xinv, yinv, marker='v', color='g', label='Stations')
    
    axs[sub1,sub2].legend()





fig, axs = plt.subplots(2, 2, figsize=(10,10))
fig.suptitle('NonRelocated Catalog Clusters')
lat_lon_clusters(catdf=catdf_eQuakeNarrow,sub1=0,sub2=0)
axs[0,0].title.set_text('easyQuake 2010')

lat_lon_clusters(catdf_OGS2Narrow,0,1)
axs[0,1].title.set_text('OGS 2010')

lat_lon_clusters(catdf_USGSNarrow,1,0)
axs[1,0].title.set_text('USGS 2010')

lat_lon_clusters(catdf_Ker2Narrow,1,1)
axs[1,1].title.set_text('Keranen 2010')





fig, axs = plt.subplots(3,2, figsize=(10,10))
fig.suptitle('Non-Relocated vs Relocated Catalog Clusters')

lat_lon_clusters(catdf=catdf_eQuakeNarrow,sub1=0,sub2=0)
axs[0,0].title.set_text('easyQuake 2010')

lat_lon_clusters(catdf=catdf_eQuake_RelocNarrow,sub1=0,sub2=1)
axs[0,1].title.set_text('easyQuake, Relocated 2010')

lat_lon_clusters(catdf_OGS2Narrow,1,0)
axs[1,0].title.set_text('OGS 2010')

lat_lon_clusters(catdf_OGS_RelocNarrow,1,1)
axs[1,1].title.set_text('OGS, Relocated 2010')

lat_lon_clusters(catdf_Ker2Narrow,2,1)
axs[2,1].title.set_text('Keranen14 (Relocated) 2010')

#NOTE: Only relocated Keranen14 data available
