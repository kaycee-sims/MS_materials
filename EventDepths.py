#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 28 18:10:10 2021

@author: kayceeschaper
"""

#import obspy
import numpy as np
import matplotlib.pyplot as plt
#import easyQuake
import pandas as pd
#import datetime


#from obspy import read_events
#from mpl_toolkits.basemap import Basemap
from easyQuake import simple_cat_df                                                            
from obspy.core import UTCDateTime

from obspy.clients.fdsn import Client

    # OGS Catalog 
catdf_ogs = pd.read_csv('http://wichita.ogs.ou.edu/eq/catalog/complete/complete.csv')
catdf_ogs['origintime'] = pd.to_datetime(catdf_ogs['origintime'])
catdf_ogs['magnitude'] = catdf_ogs['magnitude'].str.replace("None", "0").astype(float)
catdf_ogs['magnitude'] = catdf_ogs['magnitude'].astype(float)

catdf_ogs2 = catdf_ogs[(catdf_ogs['origintime']>'2010-01-01') & 
                         (catdf_ogs['origintime']<'2021-01-01')]

#Fetching the desired catalog from USGS via Obspy Client

#Tried 2021-01-01 but given "FDSNException"/"Bad Request" 20016 events exceeds search limit of 20000
#Error said "Modify the search to mathc fewer events" so tried removing one day at a time until no longer getting this error

starttime= UTCDateTime("2010-01-01")
endtime = UTCDateTime("2020-12-29")
client = Client("USGS")
cat_USGS = client.get_events(starttime=starttime, endtime=endtime, minlatitude=33.,
                         maxlatitude=37.5, minlongitude=-104, maxlongitude=-94, 
                         minmagnitude=1)
catdf_USGS=simple_cat_df(cat_USGS)
catdf_USGS['origintime'] = catdf_USGS.index[::-1]

catdf_USGS['depth_km'] = catdf_USGS['depth']/1000 



#OGS
plt.figure(1)
plt.title('OGS 2010-2020 Event Depths Over Time')
plt.ylabel('Depth(km)')
plt.xlabel('Time')
plt.scatter(catdf_ogs2.origintime,catdf_ogs2.depth_km, marker='.',c='r')

#USGS
plt.figure(2)
plt.title('USGS 2010-2020 Oklahoma Event Depths Over Time')
plt.ylabel('Depth(km)')
plt.xlabel('Time')
plt.scatter(catdf_USGS.origintime,catdf_USGS.depth_km, marker='.',c='b')

#USGS & OGS
plt.figure(3)
plt.title('USGS and OGS 2010-2020 Oklahoma Event Depths Over Time')
plt.ylabel('Depth(km)')
plt.xlabel('Time')
plt.scatter(catdf_ogs2.origintime,catdf_ogs2.depth_km, marker='.',c='r', label='OGS')
plt.scatter(catdf_USGS.origintime,catdf_USGS.depth_km, marker='.',c='b', label='USGS')
#Checked which order to list scatters; USGS needs to be on top
# ls = linestyle
plt.grid(which='major', ls='dotted')
plt.minorticks_on()
plt.grid(which='minor', ls='dotted')
plt.legend(loc='best')

#Looking at only events 40km or shallower
plt.figure(4)
plt.title('USGS and OGS 2010-2020 Oklahoma Event Depths Over Time')
plt.ylabel('Depth(km)')
plt.xlabel('Time')
plt.grid(which='major', ls='dotted')
plt.minorticks_on()
plt.grid(which='minor', ls='dotted')
plt.scatter(catdf_ogs2.origintime,catdf_ogs2.depth_km, marker='.',c='r', label='OGS')
plt.scatter(catdf_USGS.origintime,catdf_USGS.depth_km, marker='.',c='b', label='USGS')
plt.ylim(bottom=-2,top=40)
plt.legend(loc='best')

#USGS & OGS w/magnitudes
plt.figure(5)
plt.title('USGS and OGS 2010-2020 Oklahoma Event Depths Over Time')
plt.ylabel('Depth(km)')
plt.xlabel('Time')
plt.scatter(catdf_ogs2.origintime,catdf_ogs2.depth_km, marker='.',c='r',
            s=catdf_ogs2.magnitude.values**3*10, label='OGS')
plt.scatter(catdf_USGS.origintime,catdf_USGS.depth_km, marker='.',c='b', 
            s=catdf_USGS.magnitude.values**3*10, label='USGS')
plt.legend(loc='best')

plt.figure(6)
plt.title('USGS and OGS 2010-2020 Oklahoma Event Depths Over Time')
plt.ylabel('Depth(km)')
plt.xlabel('Time')
plt.scatter(catdf_ogs2.origintime,catdf_ogs2.depth_km, marker='.',c='r',
            s=catdf_ogs2.magnitude.values*3*10, label='ogs')
plt.scatter(catdf_USGS.origintime,catdf_USGS.depth_km, marker='x',c='b', 
            s=catdf_USGS.magnitude.values**3*10, label='USGS')
plt.ylim(bottom=-2,top=40)
plt.legend(loc='best')

    #Depths vs Magnitudes
plt.figure(7)
plt.title('USGS and OGS 2010-2020 Oklahoma Event Depths vs Magnitude')
plt.ylabel('Depth(km)')
plt.xlabel('Magnitude')
plt.scatter(catdf_ogs2.magnitude,catdf_ogs2.depth_km, marker='.',c='r',label='OGS')
plt.scatter(catdf_USGS.magnitude,catdf_USGS.depth_km, marker='.',c='b',label='USGS')
plt.legend(loc='best')

plt.figure(8)
plt.title('USGS and OGS 2010-2020 Oklahoma Event Depths vs Magnitude')
plt.ylabel('Depth(km)')
plt.xlabel('Magnitude')
plt.scatter(catdf_ogs2.magnitude,catdf_ogs2.depth_km, marker='.',c='r',label='OGS')
plt.scatter(catdf_USGS.magnitude,catdf_USGS.depth_km, marker='.',c='b', label='USGS')
plt.ylim(bottom=-2,top=40)
plt.legend(loc='best')


    #Depths vs Mags over Time
plt.figure(9, figsize=(10,8))
plt.title('USGS and OGS 2010-2020 Oklahoma Event Depths vs Magnitude')
plt.ylabel('Depth(km)')
plt.xlabel('Magnitude')
plt.scatter(catdf_ogs2.magnitude,catdf_ogs2.depth_km, marker='.',
              c=catdf_ogs2['origintime'], label='OGS')
plt.scatter(catdf_USGS.magnitude,catdf_USGS.depth_km, marker='x',
            c=catdf_USGS['origintime'],label='USGS')
plt.ylim(bottom=-2,top=40)
plt.legend(loc='best')
cbar=plt.colorbar(pad=0.05)
N_TICKS=7
indexes = [catdf_USGS['origintime'].iloc[i].strftime('%Y-%m-%d') for i in np.linspace(0,catdf_USGS.shape[0]-1,N_TICKS).astype(int)] 
cbar.ax.set_yticklabels(indexes)

#Same as fig9 but making background grey so light-colored events are easier to see
plt.figure(10, figsize=(10,8))
plt.axes(facecolor='grey')
plt.title('USGS and OGS 2010-2020 Oklahoma Event Depths vs Magnitude')
plt.ylabel('Depth(km)')
plt.xlabel('Magnitude')
plt.scatter(catdf_ogs2.magnitude,catdf_ogs2.depth_km, marker='.',
              c=catdf_ogs2['origintime'], label='OGS')
plt.scatter(catdf_USGS.magnitude,catdf_USGS.depth_km, marker='x',
            c=catdf_USGS['origintime'],label='USGS')
plt.ylim(bottom=-2,top=40)
plt.legend(loc='best')
cbar=plt.colorbar(pad=0.05)
N_TICKS=7
indexes = [catdf_USGS['origintime'].iloc[i].strftime('%Y-%m-%d') for i in np.linspace(0,catdf_USGS.shape[0]-1,N_TICKS).astype(int)] 
cbar.ax.set_yticklabels(indexes)


#Depth v Mag over Time OGS
plt.figure(11, figsize=(10,8))
plt.axes(facecolor='grey')
plt.title('OGS 2010-2020 Oklahoma Event Depths vs Magnitude')
plt.ylabel('Depth(km)')
plt.xlabel('Magnitude')
plt.scatter(catdf_ogs2.magnitude,catdf_ogs2.depth_km, marker='.',
              c=catdf_ogs2['origintime'])
plt.ylim(bottom=-2,top=40)
cbar=plt.colorbar(pad=0.05)
N_TICKS=7
indexes = [catdf_ogs2['origintime'].iloc[i].strftime('%Y-%m-%d') for i in np.linspace(0,catdf_ogs2.shape[0]-1,N_TICKS).astype(int)] 
cbar.ax.set_yticklabels(indexes)


#Depth v Mag over time USGS
plt.figure(12, figsize=(10,8))
plt.axes(facecolor='grey')
plt.title('USGS 2010-2020 Oklahoma Event Depths vs Magnitude')
plt.ylabel('Depth(km)')
plt.xlabel('Magnitude')
plt.scatter(catdf_USGS.magnitude,catdf_USGS.depth_km, marker='.',
            c=catdf_USGS['origintime'])
plt.ylim(bottom=-2,top=40)
cbar=plt.colorbar(pad=0.05)
N_TICKS=7
indexes = [catdf_USGS['origintime'].iloc[i].strftime('%Y-%m-%d') for i in np.linspace(0,catdf_USGS.shape[0]-1,N_TICKS).astype(int)] 
cbar.ax.set_yticklabels(indexes)

