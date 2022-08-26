#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 14:17:01 2022

@author: kayceeschaper
"""


#copy of "Clusters_and_Regressions..."scrip to play with/ show lineations lat/lon conversion issues without messing up original



#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 14:35:52 2022

@author: kayceeschaper
"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime
from obspy import read_events
from mpl_toolkits.basemap import Basemap
from easyQuake import simple_cat_df                                                            
from obspy.core import UTCDateTime
from sklearn.cluster import DBSCAN

from sklearn import linear_model, datasets #used in function
from obspy.clients.fdsn import Client



#Injection site Coord 35.43,-97.46
#Keranen Jones zone 35.4-35.8, (-97.5)-(-96.75)


#Loading Catalogs#
    # easyQuake Catalog_2010 
cat_eQuake= read_events('/Users/kayceeschaper/Earthquake_Catalogs/catalog_2010_mag.xml')     
catdf_eQuake=simple_cat_df(cat_eQuake )

catdf_eQuake['origintime'] = pd.to_datetime(catdf_eQuake['origintime'])
catdf_eQuake['magnitude'] = catdf_eQuake['magnitude'].astype(float)


    #easyQuake Big catalog (Relocated)
BigCat = pd.read_csv('/Users/kayceeschaper/HYPODD/bigRelocInputs/bigrelocations/eQuakeBig.csv', header= None) 
BigCat['longitude']=BigCat[2]
BigCat['latitude']=BigCat[1]
BigCat['magnitude']=BigCat[17] 
BigCat['origintime']=pd.to_datetime(BigCat[10])

BigCat2010=BigCat[(BigCat['origintime']>'2010-01-01') & 
                         (BigCat['origintime']<'2011-01-01')]

catdf_eQuake_Reloc = BigCat2010


    # OGS Catalog 
catdf_OGS = pd.read_csv('http://wichita.ogs.ou.edu/eq/catalog/complete/complete.csv')
catdf_OGS['origintime'] = pd.to_datetime(catdf_OGS['origintime'])
catdf_OGS['magnitude'] = catdf_OGS['magnitude'].str.replace("None", "0").astype(float)
catdf_OGS['magnitude'] = catdf_OGS['magnitude'].astype(float)

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
catdf_OGS_RelocNarrow=catdf_narrowbounds(catdf=catdf_OGS_Reloc, lat_a=35.4, lat_b=35.8, lon_a=-97.5, lon_b=-96.75)
catdf_USGSNarrow=catdf_narrowbounds(catdf=catdf_USGS2010, lat_a=35.4, lat_b=35.8, lon_a=-97.5, lon_b=-96.75)
catdf_Ker2Narrow=catdf_narrowbounds(catdf=catdf_Ker2, lat_a=35.4, lat_b=35.8, lon_a=-97.5, lon_b=-96.75)


#Constant Variables
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



#Functions
def compute_clusters(catdf=None, Title=None):
        
    map = Basemap(llcrnrlon=-104,llcrnrlat=33,urcrnrlon=-93,urcrnrlat=38, resolution='l', projection='tmerc',lat_0=35,lon_0=-98)
    x,y = map(np.array((catdf.longitude)),np.array((catdf.latitude)))
    injlon,injlat=-97.5,35.4
    xinj,yinj = map(np.array(injlon),np.array(injlat)) #injection site coordinates: 35.4N -97.5W
    xinv,yinv = map(np.array(longitudes),np.array(latitudes))
    X=np.transpose(np.array((x,y)))
    #depth = -np.array(catdf1.depth)*1000
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
    
        plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
             markeredgecolor='k', markersize=14) 

        xy = X[class_member_mask & ~core_samples_mask]
        plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
                 markeredgecolor='k', markersize=6)
   
    plt.scatter(xinj, yinj, marker='^', color='r', label='Injection Site 35.43N, -97.46W')
    plt.scatter(xinv, yinv, marker='v', color='g', label='Stations')
    plt.title(Title)    
    
def clusters_to_latlon(catdf=None, Title=None):
    map = Basemap(llcrnrlon=-104,llcrnrlat=33,urcrnrlon=-93,urcrnrlat=38, resolution='l', projection='tmerc',lat_0=35,lon_0=-98)
    x,y = map(np.array((catdf.longitude)),np.array((catdf.latitude)))
    injlon,injlat=-97.5,35.4
    xinj,yinj = map(np.array(injlon),np.array(injlat)) #injection site coordinates: 35.4N -97.5W
    xinv,yinv = map(np.array(longitudes),np.array(latitudes))
    X=np.transpose(np.array((x,y)))

    db = DBSCAN(eps=800, min_samples=5).fit(X)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_
        
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
    
        xy = X[class_member_mask & ~core_samples_mask]
      

    x,y = (np.array((catdf['longitude'])),np.array((catdf['latitude'])))
    injlon,injlat=-97.5,35.4
    xinv,yinv = (np.array(longitudes),np.array(latitudes))
    X=np.transpose(np.array((x,y)))

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
    
        plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
                 markeredgecolor='k', markersize=14) 
    
        xy = X[class_member_mask & ~core_samples_mask]
        plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
                 markeredgecolor='k', markersize=6)
        
    plt.scatter(injlon, injlat, marker='^', color='r', label='Injection Site 35.43N, -97.46W')
    plt.scatter(xinv, yinv, marker='v', color='g', label='Stations')
    plt.title(Title)    


        
    

def clusters_with_lineations(catdf=None, Title=None):    
    map = Basemap(llcrnrlon=-104,llcrnrlat=33,urcrnrlon=-93,urcrnrlat=38, resolution='l', projection='tmerc',lat_0=35,lon_0=-98)
    x,y = map(np.array((catdf.longitude)),np.array((catdf.latitude)))
    injlon,injlat=-97.5,35.4
    xinj,yinj = map(np.array(injlon),np.array(injlat)) #injection site coordinates: 35.4N -97.5W
    xinv,yinv = map(np.array(longitudes),np.array(latitudes))
    X=np.transpose(np.array((x,y)))

    #np.concatenate((a, b.T), axis=1)
    # Compute DBSCAN
    db = DBSCAN(eps=800, min_samples=5).fit(X)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_

    fault_lengths=[]
    
    Xmax=[]
    Xmin=[]
    Ymax=[]
    Ymin=[]
    
                
    plt.plot(X[:,0],X[:,1],'ko')
    plt.xlabel('(m)')
    plt.ylabel('(m)')

    def cluster_regression(xy):
        X = xy[:,0].reshape(-1,1)
        y = xy[:,1]
        # Add outlier data
        # Fit line using all data
        lr = linear_model.LinearRegression()
        lr.fit(X, y)
        # Robustly fit linear model with RANSAC algorithm
        ransac = linear_model.RANSACRegressor()
        ransac.fit(X, y)
        # Predict data of estimated models
        line_X = np.arange(X.min(), X.max())[:, np.newaxis]
        line_y = lr.predict(line_X)
        line_y_ransac = ransac.predict(line_X)
        return line_X, line_y_ransac, line_y
    
    # Black removed and is used for noise instead.
    unique_labels = set(labels)
    colors = [plt.cm.Spectral(each)
              for each in np.linspace(0, 1, len(unique_labels))]
    #for k in range(10):
    #tab 146-175 to use above loop
    for k, col in zip(unique_labels, colors):
        #print(k, col)
        if k == -1:
            # Black used for noise.
            col = [0, 0, 0, 1]
            
        class_member_mask = (labels == k) 
        xy = X[class_member_mask & ~core_samples_mask]
        xy = X[class_member_mask & core_samples_mask]
        try:
            line_X, line_y_ransac, line_y = cluster_regression(xy)
            #print("yes")
            plt.plot(xy[:,0],xy[:,1],'ro',color='cornflowerblue')
            plt.plot(line_X, line_y_ransac, color='r', linewidth=2) #, label='RANSAC regressor')
            fault_lengths.append(np.sqrt((line_X.max()-line_X.min())**2+(line_y_ransac.max()-line_y_ransac.min())**2))
            Xmax.append(line_X[-1][0])
            Xmin.append(line_X[0][0])
            Ymax.append(line_y_ransac[-1])
            Ymin.append(line_y_ransac[0])
        except:
            pass
    plt.scatter(xinj, yinj, marker='^', color='r', label='Injection Site')
    plt.scatter(xinv, yinv, marker='v', color='g', label='Station')
    plt.title(Title)
  
def cluster_regression(xy):
    X = xy[:,0].reshape(-1,1)
    y = xy[:,1]
    # Add outlier data
    # Fit line using all data
    lr = linear_model.LinearRegression()
    lr.fit(X, y)
    # Robustly fit linear model with RANSAC algorithm
    ransac = linear_model.RANSACRegressor()
    ransac.fit(X, y)
    # Predict data of estimated models
    line_X = np.arange(X.min(), X.max())[:, np.newaxis]
    line_y = lr.predict(line_X)
    line_y_ransac = ransac.predict(line_X)
    return line_X, line_y_ransac, line_y
    
def clusters_with_lineations_latlon(catdf=None, Title=None):    
    map = Basemap(llcrnrlon=-104,llcrnrlat=33,urcrnrlon=-93,urcrnrlat=38, resolution='l', projection='tmerc',lat_0=35,lon_0=-98)
    x,y = map(np.array((catdf.longitude)),np.array((catdf.latitude)))
    #injlon,injlat=-97.5,35.4
    #xinj,yinj = map(np.array(injlon),np.array(injlat)) #injection site coordinates: 35.4N -97.5W
    #xinv,yinv = map(np.array(longitudes),np.array(latitudes))
    X=np.transpose(np.array((x,y)))
    #np.concatenate((a, b.T), axis=1)
    # Compute DBSCAN
    db = DBSCAN(eps=800, min_samples=5).fit(X)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_
  
    #x,y = (np.array((catdf.longitude)),np.array((catdf.latitude)))
    injlon,injlat=-97.5,35.4
    xinv,yinv = (np.array(longitudes),np.array(latitudes))
    #X=np.transpose(np.array((x,y)))
 
    plt.plot(np.array((catdf.longitude)),np.array((catdf.latitude)),'ko')
    fault_lengths=[]
    
    Xmax=[]
    Xmin=[]
    Ymax=[]
    Ymin=[]

    
    # Black removed and is used for noise instead.
    unique_labels = set(labels)
    colors = [plt.cm.Spectral(each)
              for each in np.linspace(0, 1, len(unique_labels))]


    for k, col in zip(unique_labels, colors):
        #print(k, col)
        if k == -1:
            # Black used for noise.
            col = [0, 0, 0, 1]
        class_member_mask = (labels == k) 
        xy = X[class_member_mask & ~core_samples_mask]
        xy = X[class_member_mask & core_samples_mask]
        print(xy)

        try:
            line_X, line_y_ransac, line_y = cluster_regression(xy)
            #print("yes")
            lonev,latev=map(xy[:,0],xy[:,1],inverse=True)
            plt.plot(lonev,latev,'o',color='cornflowerblue')

            fault_lengths.append(np.sqrt((line_X.max()-line_X.min())**2+(line_y_ransac.max()-line_y_ransac.min())**2))
            print(line_X)
            Xmax.append(line_X[-1][0])
            Xmin.append(line_X[0][0])
            Ymax.append(line_y_ransac[-1])
            Ymin.append(line_y_ransac[0]) 
        except:
            pass
        #plt.plot(line_X,line_y_ransac,color='r',linewidth=2)

    plt.scatter(injlon, injlat, marker='^', color='r', label='Injection Site')
    plt.scatter(xinv, yinv, marker='v', color='g', label='Station')
    #plt.title(Title)
    lonmin, latmin = map(Xmin,Ymin,inverse=True)

    lonmax,latmax = map(Xmax,Ymax,inverse=True)
    #Below is what was used to plot the lineations separately/wihout the clusters;
    #was in cartesian though, not lat/lon
    for i in range(len(Xmax)):
        plt.plot([lonmax[i],lonmin[i]], [latmax[i],latmin[i]],'r')
    plt.title(Title)
     
     
#data: easyQuake with geopgraphic bounds
plt.figure()
compute_clusters(catdf_eQuakeNarrow,'easyQuake, Non-Relocated 2010')
plt.legend()

#Back to lat/long
plt.figure()
clusters_to_latlon(catdf_eQuakeNarrow,'easyQuake, Non-Relocated 2010')
plt.legend()

#Now with lineations
plt.figure()
clusters_with_lineations(catdf_eQuakeNarrow,'easyQuake, Non-Relocated 2010')
plt.legend()

#Now Lineations in lat/lon
plt.figure()
clusters_with_lineations_latlon(catdf_eQuakeNarrow,'easyQuake, Non-Relocated 2010')
plt.legend()

#Relocated easyQuake Data
plt.figure()
clusters_with_lineations_latlon(catdf_eQuake_RelocNarrow, 'easyQuake, Relocated 2010')
plt.legend()

#Relocated OGS Data
plt.figure()
clusters_with_lineations_latlon(catdf_OGS_RelocNarrow, 'OGS, Relocated 2010')
plt.legend()

