#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 12:58:17 2022

@author: kayceeschaper
"""
'''
Script creates plots of lineations in Lat/Lon, Histograms of length and azimuth estimations,
and Rose diagrams for both OGS and easyQuake data
'''
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import datetime
from mpl_toolkits.basemap import Basemap

from obspy import read_events
from easyQuake import simple_cat_df                                                            
from obspy.core import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth, degrees2kilometers

from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn import linear_model, datasets #used in function


#catalog_okla_hyp_all
catdfml= pd.read_csv('/Users/kayceeschaper/Earthquake_Catalogs/catalog_okla_hyp_all.csv')# header=None delimiter=r"\s+")

catdfml = catdfml.dropna(subset=['horizontal_error'])
catdfml = catdfml.reset_index(drop=True)
#can play around with dropping more events based on errors:
catdfml = catdfml[catdfml['rms']<5]
catdfml = catdfml.reset_index(drop=True)

eQuakeall = catdfml

eQuake2010= catdfml[(catdfml.origintime>'2010-01-01') & (catdfml.origintime<'2011-01-01')]

OGSall = pd.read_csv('http://wichita.ogs.ou.edu/eq/catalog/complete/complete.csv')

OGSall['origintime'] = pd.to_datetime(OGSall['origintime'])
OGSall['magnitude'] = OGSall['magnitude'].str.replace("None", "0").astype(float)
OGSall['magnitude'] = OGSall['magnitude'].astype(float)
OGSall['depth'] = OGSall['depth_km'].astype(float)
#to limit the data to only 2010
OGS2010 = OGSall[(OGSall['origintime']>'2010-01-01') & 
                         (OGSall['origintime']<'2011-01-01')]
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

def clusters_with_lineations_latlon(catdf=None, Title=None, RoseTitle=None):    
    map = Basemap(llcrnrlon=-105,llcrnrlat=30.5,urcrnrlon=-90,urcrnrlat=40, resolution='l', projection='tmerc',lat_0=35,lon_0=-98)
    x,y = map(np.array((catdf.longitude)),np.array((catdf.latitude)))
    map.drawstates(linewidth=0.5, color='gray')
    map.drawcounties(linewidth=0.25, color='gray')
    #adding lat/lon markers
    parallels = np.arange(-90.,91.,.5)
    meridians = np.arange(-180.,181., 1)
    map.drawmeridians(meridians,labels=[True,False,False,True])
    map.drawparallels(parallels,labels=[False,True,True,False])
    
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
 
    map.plot(np.array((catdf.longitude)),np.array((catdf.latitude)),'ko',latlon=True)
    
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
        #print(xy)

        try:
            line_X, line_y_ransac, line_y = cluster_regression(xy)
            #print("yes")
            lonev,latev=map(xy[:,0],xy[:,1],inverse=True)
            map.plot(lonev,latev,'o',color='cornflowerblue',latlon=True)

            fault_lengths.append(np.sqrt((line_X.max()-line_X.min())**2+(line_y_ransac.max()-line_y_ransac.min())**2))
            #print(line_X)
            Xmax.append(line_X[-1][0])
            Xmin.append(line_X[0][0])
            Ymax.append(line_y_ransac[-1])
            Ymin.append(line_y_ransac[0]) 
        except:
            pass
        #plt.plot(line_X,line_y_ransac,color='r',linewidth=2)
        
   
    map.scatter(injlon, injlat, marker='^', color='r', label='Injection Site',latlon=True)
    map.scatter(xinv, yinv, marker='v', color='g', label='Station',latlon=True)
    #plt.title(Title)
    lonmin, latmin = map(Xmin,Ymin,inverse=True)

    lonmax,latmax = map(Xmax,Ymax,inverse=True)
    #Below is what was used to plot the lineations separately/wihout the clusters;
    #was in cartesian though, not lat/lon
    for i in range(len(Xmax)):
        map.plot([lonmax[i],lonmin[i]], [latmax[i],latmin[i]],'r',latlon=True)
    plt.title(Title)
    
        
    plt.figure()
    plt.hist(fault_lengths, bins= 50)
    plt.title('Fault Lengths (km)')
    plt.xlabel('Length (km)')
    plt.ylabel('# of Occurrences')

        
    deltaX=[]
    deltaY=[] 
    for i in range(len(Xmax)):
        #print(i)
        differenceX=(Xmax[i]-Xmin[i])
        deltaX.append(differenceX)
        differenceY=(Ymax[i]-Ymin[i])
        deltaY.append(differenceY)
        
    array_azimuths=(np.degrees(np.arctan(np.array(deltaX)/np.array(deltaY))))
    
    finAz=[]
    for j in range(len(array_azimuths)):
        #print(j)
        if array_azimuths[j] < 0:
            finAz.append((180+array_azimuths[j]))
        else:
            finAz.append(array_azimuths[j])
    #print(finAz)
    
    plt.figure()
    plt.hist(finAz, bins=50) 
    plt.title('Fault Azimuths')
    plt.xlabel('Azimuth (degrees)')
    plt.ylabel('# of Occurrences')
            
    
    #Calculate the number of directions (strikes) every 10° using numpy.histogram.
    bin_edges = np.arange(-5, 366, 10)
    number_of_strikes, bin_edges = np.histogram(finAz, bin_edges)
    
    #Sum the last value with the first value.
    number_of_strikes[0] += number_of_strikes[-1]
    
    #Sum the first half 0-180° with the second half 180-360° to achieve the "mirrored behavior" of Rose Diagrams.
    half = np.sum(np.split(number_of_strikes[:-1], 2), 0)
    two_halves = np.concatenate([half, half])
    
    fig = plt.figure(figsize=(8,8))
    
    ax = fig.add_subplot(projection='polar')
    
    ax.bar(np.deg2rad(np.arange(0, 360, 10)), two_halves, 
           width=np.deg2rad(10), bottom=0.0, color='.8', edgecolor='k')
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_thetagrids(np.arange(0, 360, 10), labels=np.arange(0, 360, 10))
    ax.set_rgrids(np.arange(1, two_halves.max() + 1, 2), angle=0, weight= 'black')
    ax.set_title( RoseTitle, y=1.10, fontsize=15)




plt.figure()
clusters_to_latlon(OGSall, 'OGS')

plt.figure()
clusters_to_latlon(eQuakeall,'eQuake')

plt.figure()
clusters_with_lineations_latlon(OGSall,Title='OGS', RoseTitle='OGS Orientations') 

plt.figure()
clusters_with_lineations_latlon(OGSall, Title='eQuake', RoseTitle='easyQuake Orientations') 

