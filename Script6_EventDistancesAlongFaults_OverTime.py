#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 11:50:29 2022

@author: kayceeschaper
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 14:51:04 2022

@author: kayceeschaper
"""
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

from scipy import stats
from scipy.stats import linregress, t
from scipy.optimize import curve_fit as CF

#catalog_okla_hyp_all
catdfml= pd.read_csv('/Users/kayceeschaper/Earthquake_Catalogs/catalog_okla_hyp_all.csv')# header=None delimiter=r"\s+")

catdfml = catdfml.dropna(subset=['horizontal_error'])
catdfml = catdfml.reset_index(drop=True)
#can play around with dropping more events based on errors:
catdfml = catdfml[catdfml['rms']<7]
catdfml = catdfml.reset_index(drop=True)


#Injection site Coord 35.43,-97.46
#Keranen Jones zone 35.4-35.8, (-97.5)-(-96.75)


#easyQuake submodule; creates new dataframe geographic limitations
def catdf_narrowbounds(catdf=None,lat_a=None,lat_b=None,lon_a=None,lon_b=None):
    catdf = catdf[(catdf['latitude']>lat_a) & (catdf['latitude']<lat_b) & (catdf['longitude']>lon_a) & (catdf['longitude']<lon_b)]
    return catdf

#catdf_okla_hyp_all_Narrow=catdf_narrowbounds(catdf=catdf_okla_hyp_all, lat_a=35.4, lat_b=35.8, lon_a=-97.5, lon_b=-96.75)



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


map = Basemap(llcrnrlon=-105,llcrnrlat=30.5,urcrnrlon=-90,urcrnrlat=40, resolution='l', projection='tmerc',lat_0=35,lon_0=-98)
#catdf=catdf_okla_hyp_all
#condenses the columns with datetime info into one column
catdfr = catdfml
catdfr = catdfr.dropna()
catdfr = catdfr.reset_index(drop=True)
#rutc = np.zeros((len(catdfr.index),1))
rutc = []
for i in range(0,len(catdfr.index)):
    rutc.append(UTCDateTime(catdfr.origintime[i]).datetime)
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

    except:
        pass  
    
   #plt.figure()
   
    ##needs to be the lat/lon of the endpoints of the fault line in question
    #line = [(-98.21,35.725), (-98.21, 36.16)] # kingfisher
    line = [(lonmin,latmin),(lonmax,latmax)]
    
    p1 = np.array(line[0])
    p2 = np.array(line[1])
    
    total_dist = (gps2dist_azimuth(p1[1],p1[0],p2[1],p2[0]))[0]/1000
    
    p3=np.transpose(np.array([df.longitude,df.latitude]))
    
    d = np.cross(p2-p1,p3-p1)/np.linalg.norm(p2-p1)
    dkm = degrees2kilometers(d)
    
    d_along = []
    dtime = [] 
    for idx2, val in enumerate(df.latitude):
        #if np.abs(dkm[idx2])<.1:
        epi_dist1, az1, baz1 = gps2dist_azimuth(p1[1],p1[0],df.latitude.iloc[idx2],df.longitude.iloc[idx2])
        d_along.append(epi_dist1)
        dtime1 = df['rutc'].iloc[idx2]-df['rutc'].iloc[0]
        dtime.append(dtime1.total_seconds()/3600/24/365.25)
    
    #plt.plot(d_along,dtime,'.')
    #d_along = np.array(d_along)
    #df['rutc'] = pd.to_datetime(df['rutc'])

    plt.xlabel('meters along fault')
    plt.ylabel('years')
    #Printing name/label of cluster to make ouput readable
    #print(f"z={z}")
    
    res = stats.linregress(d_along, dtime)
    #print(f"R-squared: {res.rvalue**2:.6f}")
    plt.plot(d_along, dtime, '.', label='original data')
   #plt.plot(d_along, res.intercept + res.slope*(np.array(d_along)), 'r', label='fitted line')
    #plt.legend()
    plt.title('Event Distances Along Faults Over Time')
    plt.show()
    
    tinv = lambda p, df: abs(t.ppf(p/2, df))
    ts = tinv(0.05, len(x)-2)

    
    #print(f"slope (95%): {res.slope:.6f} +/- {ts*res.stderr:.6f}")
    #print(f"intercept (95%): {res.intercept:.6f}"
         # f" +/- {ts*res.intercept_stderr:.6f}")

    
    #Standard Deviation for distance ALONG fault lineation
    #print(f"StanDev d_along:{np.std(d_along)}")
    #Standard Deviation for distance FROM fault lineation
    #print(f"StanDev d:{np.std(d)}")
   
    
#plt.legend()