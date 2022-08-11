#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 15:25:46 2022

@author: kayceeschaper
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Original Pieces Created on Mon Jun 28 16:52:56 2021

@author: kayceeschaper
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

from obspy import read_events
from mpl_toolkits.basemap import Basemap
from easyQuake import simple_cat_df                                                            
from obspy.core import UTCDateTime
from sklearn.cluster import DBSCAN

from obspy.clients.fdsn import Client

from easyQuake import quakeml_to_hypodd


#Load Catalog
 # easyQuake Catalog_2010 
cat_eQuake= read_events('/Users/kayceeschaper/Earthquake_Catalogs/catalog_2010_mag.xml')     
catdf_eQuake=simple_cat_df(cat_eQuake )

#catdf_eQuake['origintime'] = catdf_eQuake.index

catdf_eQuake2 = catdf_eQuake[(catdf_eQuake.origintime>'2010-03-01') & (catdf_eQuake.origintime<'2010-10-01')]


#Getting HypoDD file from the easyQuake QuakeML file

#quakeml_to_hypodd(cat=cat, project_folder='write_out_current_full_path', project_code='ok', download_station_metadata=True)
quakeml_to_hypodd(cat=cat_eQuake, project_folder='/Users/kayceeschaper/HYPODD', project_code='ok')


#Plotting catalog data
def plot_hypodd_catalog(file=None):    
    catdfr = pd.read_csv(file,delimiter=r"\s+")
    catdfr = catdfr.dropna()
    catdfr = catdfr.reset_index(drop=True)
    rutc = []
    for i in range(0,len(catdfr.index)):
        rutc.append(UTCDateTime(int(catdfr.iloc[i,10]),int(catdfr.iloc[i,11]),int(catdfr.iloc[i,12]),int(catdfr.iloc[i,13]),int(catdfr.iloc[i,14]),catdfr.iloc[i,15]))
   
    catdfr['rutc'] = rutc
    catdfr.sort_values(by=['rutc'], inplace=True)
    catdfr = catdfr.reset_index(drop=True)
   

    # 1. Draw the map background
    lat0 = 35.55 #np.median(catdfr.iloc[:,1].values)
    lon0 = -97.4 #np.median(catdfr.iloc[:,2].values)
    m = Basemap(projection='lcc', resolution='h',
                lat_0=lat0, lon_0=lon0,
                width=.8E5, height=.55E5)
    m.drawcountries(color='gray')
    m.drawcounties(color='gray')
    m.drawstates(color='gray')
   
    #adding lat/lon markers
    parallels = np.arange(-90.,91.,.05)
    meridians = np.arange(-180.,181., .05)
    m.drawmeridians(meridians,labels=[True,False,False,True], rotation=45)
    m.drawparallels(parallels,labels=[False,True,True,False])

    # 2. scatter city data, with color reflecting population and size reflecting area
    m.scatter(catdfr.iloc[:,2].values,catdfr.iloc[:,1].values,s=catdfr.iloc[:,16].values**3*8,c=catdfr.index,marker='o',alpha=0.5,latlon=True)
   

    cbar = plt.colorbar(pad=0.15)
    N_TICKS=8
    indexes = [catdfr['rutc'].iloc[i].strftime('%Y-%m-%d') for i in np.linspace(0,catdfr.shape[0]-1,N_TICKS).astype(int)]
    
    cbar.ax.set_yticklabels(indexes)
    #plt.savefig('hypoDDmap.png')
    plt.show()
  
#Unrelocated/input data plot
plt.figure(1)
plot_hypodd_catalog('/Users/kayceeschaper/HYPODD/hypodd.loc')
plt.title('Non-relocated Data')



#Relocated data plot #Best relocation
plt.figure(2)
plot_hypodd_catalog(file='/Users/kayceeschaper/HYPODD/hypoDD.relocpanfixed_KSBest')
plt.title('HypoDD Relocated Events; KS Best')



# Clustering

def augment(xyzs):
    axyz = np.ones((len(xyzs), 4))
    axyz[:, :3] = xyzs
    return axyz

def estimate(xyzs):
    axyz = augment(xyzs[:3])
    return np.linalg.svd(axyz)[-1][-1, :]

def is_inlier(coeffs, xyz, threshold):
    return np.abs(coeffs.dot(augment([xyz]).T)) < threshold

def plot_plane(a, b, c, d):
    xx, yy = np.meshgrid(np.arange(np.min(xyzs[:,0]),np.max(xyzs[:,0])),np.arange(np.min(xyzs[:,1]),np.max(xyzs[:,1])))
    return xx, yy, (-d - a * xx - b * yy) / c
def calc_strikedip(pts):
    ptA, ptB, ptC = pts[0], pts[1], pts[2]
    x1, y1, z1 = float(ptA[0]), float(ptA[1]), float(ptA[2])
    x2, y2, z2 = float(ptB[0]), float(ptB[1]), float(ptB[2])
    x3, y3, z3 = float(ptC[0]), float(ptC[1]), float(ptC[2])


    u1 = float(((y1 - y2) * (z3 - z2) - (y3 - y2) * (z1 - z2)))
    u2 = float((-((x1 - x2) * (z3 - z2) - (x3 - x2) * (z1 - z2))))
    u3 = float(((x1 - x2) * (y3 - y2) - (x3 - x2) * (y1 - y2)))

    '''
    Calculate pseudo eastings and northings from origin
    these are actually coordinates of a new point that represents
    the normal from the plane's origin defined as (0,0,0).

    If the z value (u3) is above the plane we first reverse the easting 
    then we check if the z value (u3) is below the plane, if so
    we reverse the northing.

    This is to satisfy the right hand rule in geology where dip is always
    to the right if looking down strike.
    '''
    if u3 < 0:
        easting = u2
    else:
        easting = -u2

    if u3 > 0:
        northing = u1
    else:
        northing = -u1

    if easting >= 0:
        partA_strike = math.pow(easting, 2) + math.pow(northing, 2)
        strike = math.degrees(math.acos(northing / math.sqrt(partA_strike)))
    else:
        partA_strike = northing / math.sqrt(math.pow(easting, 2) + math.pow(northing, 2))
        strike = math.degrees(2 * math.pi - math.acos(partA_strike))

    # determine dip
    print(strike, 'strike')
    part1_dip = math.sqrt(math.pow(u2, 2) + math.pow(u1, 2))
    part2_dip = math.sqrt(math.pow(u1,2) + math.pow(u2,2) + math.pow(u3,2))
    dip = math.degrees(math.asin(part1_dip / part2_dip))

    return strike, dip




#Relocated and clustered data

catdf1 = pd.read_csv('/Users/kayceeschaper/HYPODD/hypodd.relocpanfixed_KSBest', header=None, delimiter=r"\s+") #delimiter=r"\s+" means working with txt file/using spaces instead of commas
catdf = catdf1
map = Basemap(llcrnrlon=-104,llcrnrlat=33,urcrnrlon=-93,urcrnrlat=38, resolution='l', projection='tmerc',lat_0=35,lon_0=-98)
x,y = map(np.array((catdf[2])),np.array((catdf[1])))
X=np.transpose(np.array((x,y)))
depth = -np.array(catdf[3])*1000

#yeilds 18 clusters; 159 noise points 
plt.figure(3)
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


    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
             markeredgecolor='k', markersize=14) #need to make each color its own label, but this is where I can add labels

    xy = X[class_member_mask & ~core_samples_mask]
    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
             markeredgecolor='k', markersize=6)

plt.legend()
plt.title('eps=800; min_samples=5 - Estimated number of clusters: %d' % n_clusters_)
plt.show()





#Initial dataset; not relocated events. Clustered

catdf2 = pd.read_csv('/Users/kayceeschaper/HYPODD/hypodd.loc', header=None, delimiter=r"\s+") #delimiter=r"\s+" means working with txt file/using spaces instead of commas
catdf = catdf2
map = Basemap(llcrnrlon=-104,llcrnrlat=33,urcrnrlon=-93,urcrnrlat=38, resolution='l', projection='tmerc',lat_0=35,lon_0=-98)
x,y = map(np.array((catdf[2])),np.array((catdf[1])))

X=np.transpose(np.array((x,y)))
depth = -np.array(catdf[3])*1000


plt.figure(4)
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


    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
             markeredgecolor='k', markersize=14)

    xy = X[class_member_mask & ~core_samples_mask]
    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
             markeredgecolor='k', markersize=6)

plt.title('eps=800; min_samples=5 - Estimated number of clusters: %d' % n_clusters_)
plt.show()



#Relocated and clustered data; with injection site and seismic stations plotted

#adding injection site location and stations to map
#Station locations for USGS catalog 
client = Client("IRIS")
starttime= UTCDateTime("2010-01-01")
endtime = UTCDateTime("2010-07-01")
lat1=34
lat2=39
lon1=-95
lon2=-100
inv = client.get_stations(minlatitude=lat1, minlongitude=lon2 ,maxlatitude=lat2, maxlongitude=lon1,
                                starttime=starttime,
                                endtime=endtime)

latitudes = []
longitudes = []
for net in inv:
    for sta in net:
        latitudes.append(sta.latitude)
        longitudes.append(sta.longitude)

catdf1 = pd.read_csv('/Users/kayceeschaper/HYPODD/hypodd.relocpanfixed_KSBest', header=None, delimiter=r"\s+") #delimiter=r"\s+" means working with txt file/using spaces instead of commas
catdf = catdf1
map = Basemap(llcrnrlon=-104,llcrnrlat=33,urcrnrlon=-93,urcrnrlat=38, resolution='l', projection='tmerc',lat_0=35,lon_0=-98)
x,y = map(np.array((catdf[2])),np.array((catdf[1])))
xinj,yinj = map(-97.5,35.4) #injection site coordinates: 35.4N -97.5W
xinv,yinv = map(np.array(longitudes),np.array(latitudes))
X=np.transpose(np.array((x,y)))
depth = -np.array(catdf[3])*1000


plt.figure(5)
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


    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
             markeredgecolor='k', markersize=14) #need to make each color its own label, but this is where I can add labels

    xy = X[class_member_mask & ~core_samples_mask]
    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
             markeredgecolor='k', markersize=6)

plt.scatter(xinj, yinj, marker='^', color='r', label='Injection Site')
plt.scatter(xinv, yinv, marker='v', color='g', label='Station')
plt.legend()
plt.title('eps=800; min_samples=5 - Estimated number of clusters: %d' % n_clusters_)
plt.show()

