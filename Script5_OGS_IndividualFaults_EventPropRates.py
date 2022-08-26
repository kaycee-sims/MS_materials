#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 10:12:44 2022

@author: kayceeschaper
"""

'''
This script looks at event propagation rates for the OGS 2010 clusters and
provides some statistics

'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import datetime
from mpl_toolkits.basemap import Basemap
from statistics import mean, stdev

from obspy import read_events
from easyQuake import simple_cat_df                                                            
from obspy.core import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth, degrees2kilometers

from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn import linear_model, datasets #used in function

from scipy import stats
from scipy.stats import linregress, t



#Injection site Coord 35.43,-97.46
#Keranen Jones zone 35.4-35.8, (-97.5)-(-96.75)

###Loading Catalogs
catdf_OGS_Reloc= pd.read_csv('/Users/kayceeschaper/HYPODD/OK2010/hypoDD.relocpanfixed_2010', header=None, delimiter=r"\s+")
catdf_OGS_Reloc['latitude']=catdf_OGS_Reloc[1]
catdf_OGS_Reloc['longitude']=catdf_OGS_Reloc[2]


#easyQuake submodule; creates new dataframe geographic limitations
def catdf_narrowbounds(catdf=None,lat_a=None,lat_b=None,lon_a=None,lon_b=None):
    catdf = catdf[(catdf['latitude']>lat_a) & (catdf['latitude']<lat_b) & (catdf['longitude']>lon_a) & (catdf['longitude']<lon_b)]
    return catdf

catdf_OGS_RelocNarrow=catdf_narrowbounds(catdf=catdf_OGS_Reloc, lat_a=35.4, lat_b=35.8, lon_a=-97.5, lon_b=-96.75)


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



map = Basemap(llcrnrlon=-104,llcrnrlat=33,urcrnrlon=-93,urcrnrlat=38, resolution='l', projection='tmerc',lat_0=35,lon_0=-98)


#condenses the columns with datetime info into one column
catdfr = catdf_OGS_RelocNarrow
#catdfr = catdfr.dropna()
#catdfr = catdfr.reset_index(drop=True)
#rutc = np.zeros((len(catdfr.index),1))
rutc = []
for i in range(0,len(catdfr.index)):
    rutc.append(UTCDateTime(int(catdfr.iloc[i,10]),int(catdfr.iloc[i,11]),
                            int(catdfr.iloc[i,12]),int(catdfr.iloc[i,13]),int(catdfr.iloc[i,14]),catdfr.iloc[i,15]).datetime)
catdfr['rutc'] = rutc
catdfr.sort_values(by=['rutc'], inplace=True)
catdfr = catdfr.reset_index(drop=True)
x,y = map(np.array((catdfr.longitude)),np.array((catdfr.latitude)))
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
#get lables out,
#labels = cluster_labelsout(catdf_eQuake_RelocNarrow)
unique_labels = set(labels)
for z in range(len(list(unique_labels)[:-1])):
    df = catdfr[labels==z]
    #map = Basemap(llcrnrlon=-104,llcrnrlat=33,urcrnrlon=-93,urcrnrlat=38, resolution='l', projection='tmerc',lat_0=35,lon_0=-98)
    #x,y = map(np.array((df.longitude)),np.array((df.latitude)))
    X=np.transpose(np.array((x,y)))
    class_member_mask = (labels == z) 
    xy = X[class_member_mask & ~core_samples_mask]
    xy = X[class_member_mask & core_samples_mask]
    try:
        line_X, line_y_ransac, line_y = cluster_regression(xy)
        lonmin, latmin = map(line_X[0][0],line_y_ransac[0],inverse=True)
        lonmax,latmax = map(line_X[-1][0],line_y_ransac[-1],inverse=True)
        # plt.figure()
        # plt.scatter(df.longitude, df.latitude,s=df.iloc[:,16].values**3*8,c=df.index,marker='o',alpha=0.5)
        # plt.plot([lonmax,lonmin], [latmax,latmin],'r')
        # plt.show()
        # cbar = plt.colorbar()
        # N_TICKS=8
        # indexes = [catdfr['rutc'].iloc[i].strftime('%Y-%m-%d') for i in np.linspace(0,catdfr.shape[0]-1,N_TICKS).astype(int)]
           
        # #indexes = [catdfr.index[i].strftime('%Y-%m-%d') for i in np.linspace(0,catdfr.shape[0]-1,N_TICKS).astype(int)]
        # cbar.ax.set_yticklabels(indexes)
        # #plt.savefig('hypoDDmap.png')
        # #plt.show()
    except:
        pass
        
    plt.figure()
    ##needs to be the lat/lon of the endpoints of the fault line in question
    #line = [(-98.21,35.725), (-98.21, 36.16)] # kingfisher
    line = [(lonmin,latmin),(lonmax,latmax)]
    
    p1 = np.array(line[0])
    p2 = np.array(line[1])
    
    total_dist = (gps2dist_azimuth(p1[1],p1[0],p2[1],p2[0]))[0]/1000
    
    p3=np.transpose(np.array([df[2],df[1]]))
    
    d = np.cross(p2-p1,p3-p1)/np.linalg.norm(p2-p1)
    dkm = degrees2kilometers(d)
    
    d_along = []
    dtime = [] 
    for idx2, val in enumerate(df[1]):
        #if np.abs(dkm[idx2])<.1:
        epi_dist1, az1, baz1 = gps2dist_azimuth(p1[1],p1[0],df[1].iloc[idx2],df[2].iloc[idx2])
        d_along.append(epi_dist1)
        dtime1 = df['rutc'].iloc[idx2]-df['rutc'].iloc[0]
        dtime.append(dtime1.total_seconds()/60)
    
    #plt.plot(d_along,dtime,'.')
    #d_along = np.array(d_along)
    #df['rutc'] = pd.to_datetime(df['rutc'])

    plt.xlabel('meters')
    plt.ylabel('hours')
    #catdfr['d_along'] = d_along
    #catdfr['d'] = d
    #total_dist = (gps2dist_azimuth(p1[1],p1[0],p2[1],p2[0]))[0]/1000
    #depth = -(catdfalong[3])*1000
    #plt.plot(catdfalong['d_along'],depth,'.')
    
    #Printing name/label of cluster to make ouput readable
    print(f"z={z}")
    
    res = stats.linregress(d_along, dtime)
    print(f"R-squared: {res.rvalue**2:.6f}")
    plt.plot(d_along, dtime, '.', label='original data')
    plt.plot(d_along, res.intercept + res.slope*(np.array(d_along)), 'r', label='fitted line')
    plt.legend()
    plt.show()
    
    tinv = lambda p, df: abs(t.ppf(p/2, df))
    ts = tinv(0.05, len(x)-2)

    
    print(f"slope (95%): {res.slope:.6f} +/- {ts*res.stderr:.6f}")
    print(f"intercept (95%): {res.intercept:.6f}"
          f" +/- {ts*res.intercept_stderr:.6f}")

    
    #Standard Deviation for distance ALONG fault lineation
    print(f"StanDev d_along:{np.std(d_along)}")
    #Standard Deviation for distance FROM fault lineation
    print(f"StanDev d:{np.std(d)}")
  

    # #trying to make a list of outliers to be able to remove them
    # HighOutlier_d = []
    # HighOutlier_d_along = []
    # HighOutlier_dtime=[]
    # LowOutlier_d = []
    # LowOutlier_d_along = []    
    # LowOutlier_dtime=[]
    
    # HighOutlier_d.append(np.where(d > (np.mean(d) + (2*np.std(d)))))
    # LowOutlier_d.append(np.where(d < (np.mean(d) + (2*np.std(d)))))
        
    # HighOutlier_d_along.append(np.where(d_along > (np.mean(d_along) + (2*np.std(d_along)))))
    # LowOutlier_d_along.append(np.where(d_along < (np.mean(d_along) + (2*np.std(d_along)))))
    
    
    # HighOutlier_dtime.append(np.where(dtime > (np.mean(dtime) + (2*np.std(dtime)))))
    # LowOutlier_dtime.append(np.where(dtime < (np.mean(dtime) + (2*np.std(dtime)))))
    
 
    
    
