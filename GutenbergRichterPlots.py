#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 14:39:41 2021

@author: kayceeschaper
"""

"""
Gutenberg–Richter Catalog Comparison
"""
#import obspy
import numpy as np
import matplotlib.pyplot as plt
#import easyQuake
import pandas as pd
#import datetime


from obspy import read_events
#from mpl_toolkits.basemap import Basemap
from easyQuake import simple_cat_df                                                            
from obspy.core import UTCDateTime

from obspy.clients.fdsn import Client


    # easyQuake Catalog_2010 
#cat_eQuake= read_events('/Users/kayceeschaper/Earthquake_Catalogs/catalog_2010_mag.xml')  
catdf_eQuake= pd.read_csv('/Users/kayceeschaper/Earthquake_Catalogs/catalog_okla_hyp_all.csv')# header=None delimiter=r"\s+")

#catdf_eQuake=simple_cat_df(cat_eQuake )

#catdf_eQuake['origintime'] = catdf_eQuake.index

catdf_eQuake['origintime'] = pd.to_datetime(catdf_eQuake['origintime'])
catdf_eQuake['magnitude'] = catdf_eQuake['magnitude'].astype(float)

eQuake2010= catdf_eQuake[(catdf_eQuake['origintime']>'2010-01-01') &
                            (catdf_eQuake['origintime']<'2011-01-01')]


    # OGS Catalog 
catdf_Wich = pd.read_csv('http://wichita.ogs.ou.edu/eq/catalog/complete/complete.csv')
catdf_Wich['origintime'] = pd.to_datetime(catdf_Wich['origintime'])
catdf_Wich['magnitude'] = catdf_Wich['magnitude'].str.replace("None", "0").astype(float)
catdf_Wich['magnitude'] = catdf_Wich['magnitude'].astype(float)

#to limit the data to only 2010
catdf_Wich2 = catdf_Wich[(catdf_Wich['origintime']>'2010-01-01') & 
                         (catdf_Wich['origintime']<'2011-01-01')]
#to go back to using the full catalog
#catdf_Wich['origintime']=pd.to_datetime(catdf_Wich['origintime'])

        # USGS 2010 Catalog
#cat_USGS2010=read_events('/Users/kayceeschaper/Earthquake_Catalogs/catalog_query_USGS_2010.xml')
#catdf_USGS2010=simple_cat_df(cat_USGS2010)


#Fetching the desired catalog from USGS via Obspy Client
starttime= UTCDateTime("2010-01-01")
endtime = UTCDateTime("2011-01-01")
client = Client("USGS")
cat_USGS2010 = client.get_events(starttime=starttime, endtime=endtime, minlatitude=33.,
                         maxlatitude=37.5, minlongitude=-104, maxlongitude=-94, 
                         minmagnitude=1)
catdf_USGS2010=simple_cat_df(cat_USGS2010)
catdf_USGS2010['origintime'] = catdf_USGS2010.index

catdf_USGS_Rev= catdf_USGS2010.reindex(index=catdf_USGS2010.index[::-1])



    #GR plot for easyQuake data #Use eQuake2010 for specifically 2010
hist, edges = np.histogram(a=eQuake2010['magnitude'].values, bins=101, range=(0,10))
chist = np.cumsum(hist[::-1])[::-1]

fig, ax = plt.subplots()

#these two lines allow subscript in axis titles
params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)
#
ax.plot(edges[:-1], hist, marker='.', color='k', linestyle='', label='Magnitude-Frequency')
ax.plot(edges[:-1], chist, marker='o', color='k', linestyle='',label='Cumulative Magnitude Frequency')
ax.set_yscale('log')
plt.xlim(-0.5,6)
plt.ylim(0,10**4)

plt.title('easyQuake 2010')
plt.xlabel('Magnitude')
plt.ylabel('Frequency')

#adding gridlines for readability
plt.grid(which='major')
plt.minorticks_on()
plt.grid(which='minor')
plt.legend()


    #GR plot for Ogs data #Use catdf_Wich2 for specific daterange
hist, edges = np.histogram(a=catdf_Wich2['magnitude'].values, bins=101, range=(0,10))
chist = np.cumsum(hist[::-1])[::-1]

fig2, ax = plt.subplots()
ax.plot(edges[:-1], hist, marker='.', color='k', linestyle='', label='Magnitude-Frequency')
ax.plot(edges[:-1], chist, marker='o', color='k', linestyle='',label='Cumulative Magnitude Frequency')
ax.set_yscale('log')
plt.xlim(-0.5,6)
plt.ylim(0,10**4)

plt.title('OGS 2010')
plt.xlabel('Magnitude')
plt.ylabel('Frequency')
plt.grid(which='major')
plt.minorticks_on()
plt.grid(which='minor')
plt.legend()


    #GR plot for USGS data #daterange selected via starttime/endtime in Client 
hist, edges = np.histogram(a=catdf_USGS_Rev['magnitude'].values, bins=101, range=(0,10))
chist = np.cumsum(hist[::-1])[::-1]

fig3, ax = plt.subplots()
ax.plot(edges[:-1], hist, marker='.', color='k', linestyle='',label='Magnitude-Frequency')
ax.plot(edges[:-1], chist, marker='o', color='k', linestyle='',label='Cumulative Magnitude Frequency')
ax.set_yscale('log')
plt.xlim(-0.5,6)
plt.ylim(0,10**4)

plt.title('USGS 2010')
plt.xlabel('Magnitude')
plt.ylabel('Frequency')
plt.grid(which='major')
plt.minorticks_on()
plt.grid(which='minor')
plt.legend()


#commented out since no Magnitude data for Keranen
#GR plot for Keranen data




#Limiting to Jones

#Injection site Coord 35.43,-97.46
#Keranen Jones zone 35.4-35.8, -97.5- -96.75


#easyQuake submodule; creates new dataframe geographic limitations
def catdf_narrowbounds(catdf=None,lat_a=None,lat_b=None,lon_a=None,lon_b=None):
    catdf = catdf[(catdf['latitude']>lat_a) & (catdf['latitude']<lat_b) & (catdf['longitude']>lon_a) & (catdf['longitude']<lon_b)]
    return catdf

catdf_eQuakeNarrow=catdf_narrowbounds(catdf=eQuake2010, lat_a=35.4, lat_b=35.8, lon_a=-97.5, lon_b=-96.75)
catdf_Wich2Narrow=catdf_narrowbounds(catdf=catdf_Wich2, lat_a=35.4, lat_b=35.8, lon_a=-97.5, lon_b=-96.75)
catdf_USGSNarrow=catdf_narrowbounds(catdf=catdf_USGS2010, lat_a=35.4, lat_b=35.8, lon_a=-97.5, lon_b=-96.75)




###After Narrowing geopgraphic bounds

    #easyQuake data #Use catdf_eQuake2 for specific daterange
hist, edges = np.histogram(a=catdf_eQuakeNarrow['magnitude'].values, bins=101, range=(0,10))
chist = np.cumsum(hist[::-1])[::-1]

fig, ax = plt.subplots()

#these two lines allow subscript in axis titles
params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)
#
ax.plot(edges[:-1], hist, marker='.', color='k', linestyle='', label='Magnitude-Frequency')
ax.plot(edges[:-1], chist, marker='o', color='k', linestyle='',label='Cumulative Magnitude Frequency')
ax.set_yscale('log')
plt.xlim(-0.5,6)
plt.ylim(0,10**4)

plt.title('easyQuake 2010 Jones')
plt.xlabel('Magnitude')
plt.ylabel('Frequency')

#adding gridlines for readability
plt.grid(which='major')
plt.minorticks_on()
plt.grid(which='minor')
plt.legend()


    #OGS data #Use catdf_Wich2 for specific daterange
hist, edges = np.histogram(a=catdf_Wich2Narrow['magnitude'].values, bins=101, range=(0,10))
chist = np.cumsum(hist[::-1])[::-1]

fig2, ax = plt.subplots()
ax.plot(edges[:-1], hist, marker='.', color='k', linestyle='', label='Magnitude-Frequency')
ax.plot(edges[:-1], chist, marker='o', color='k', linestyle='',label='Cumulative Magnitude Frequency')
ax.set_yscale('log')
plt.xlim(-0.5,6)
plt.ylim(0,10**4)

plt.title('OGS 2010 Jones')
plt.xlabel('Magnitude')
plt.ylabel('Frequency')
plt.grid(which='major')
plt.minorticks_on()
plt.grid(which='minor')
plt.legend()


    #USGS data #daterange selected via starttime/endtime in Client 
hist, edges = np.histogram(a=catdf_USGSNarrow['magnitude'].values, bins=101, range=(0,10))
chist = np.cumsum(hist[::-1])[::-1]

fig3, ax = plt.subplots()
ax.plot(edges[:-1], hist, marker='.', color='k', linestyle='',label='Magnitude-Frequency')
ax.plot(edges[:-1], chist, marker='o', color='k', linestyle='',label='Cumulative Magnitude Frequency')
ax.set_yscale('log')
plt.xlim(-0.5,6)
plt.ylim(0,10**4)

plt.title('USGS 2010 Jones')
plt.xlabel('Magnitude')
plt.ylabel('Frequency')
plt.grid(which='major')
plt.minorticks_on()
plt.grid(which='minor')
plt.legend()





