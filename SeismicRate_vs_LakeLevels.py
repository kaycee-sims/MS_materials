#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  2 09:23:24 2021

@author: kayceeschaper

Contains starter script clip created by JWalter
"""
# Goals of this code:
    #pull seismicity rates from the catalogs I have to work with
        #for this, essentially we want to look only at origintme data and ignore the rest of the columns
    #look at rates of seismic activity over time
    #pull in lake/reservoir level data
    #compare said water data to seis rate data

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime

from obspy import read_events
from easyQuake import simple_cat_df                                                            
from obspy.core import UTCDateTime 
from obspy.clients.fdsn import Client


#Apporoximate catalog magnitudes of completeness as seen using "GutenbergRichterPlots.py":
   # easyQuake:2.0
   # OGS:2.2
   # USGS:2.5
   
    #Loading relevent catalogs
 # OGS Catalog 
catdf_OGS = pd.read_csv('http://wichita.ogs.ou.edu/eq/catalog/complete/complete.csv')
catdf_OGS['origintime'] = pd.to_datetime(catdf_OGS['origintime'])
catdf_OGS['magnitude'] = catdf_OGS['magnitude'].str.replace("None", "0").astype(float)
catdf_OGS['magnitude'] = catdf_OGS['magnitude'].astype(float)
#catdf_OGS2 = catdf_OGS[(catdf_OGS['origintime']>'2009-01-01') & 
#                         (catdf_OGS['origintime']<'2015-01-01')]
   
    # easyQuake Catalog_2010 
cat_eQuake= read_events('/Users/kayceeschaper/Earthquake_Catalogs/catalog_2010_mag.xml')     
catdf_eQuake=simple_cat_df(cat_eQuake )

catdf_eQuake['origintime'] = catdf_eQuake.index

#catdf_eQuake2= catdf_eQuake[(catdf_eQuake['origintime']>'2009-01-01') & 
#                            (catdf_eQuake['origintime']<'2015-01-01')]



    # USGS downloaded catalog for 2009 to 2015 ; FDSN down as of 11am 4/8/21
catdf_USGS=pd.read_csv('/Users/kayceeschaper/Earthquake_Catalogs/USGS2009_2021Oklahoma.csv')

 #USGS Catalog
# starttime= UTCDateTime("2009-01-01")
# endtime = UTCDateTime("2015-01-01")
# client = Client("USGS")
# cat_USGS = client.get_events(starttime=starttime, endtime=endtime, minlatitude=33.,
#                          maxlatitude=37.5, minlongitude=-104, maxlongitude=-94, 
#                          minmagnitude=1)

#catdf_USGS=simple_cat_df(cat_USGS)
catdf_USGS['origintime'] = pd.to_datetime(catdf_USGS.time)

catdf_USGS['depth_km'] = catdf_USGS['depth']/10
catdf_USGS['magnitude'] = catdf_USGS['mag']

LakeLevels=pd.read_csv('/Users/kayceeschaper/Lake_Level_Data/dv.txt',sep='\t')
LakeLevels['Date']=pd.to_datetime(LakeLevels['date'])
LakeLevels['WWatDate']=LakeLevels.Date[0:1904]
LakeLevels['ShawDate']=LakeLevels.Date[1904:4096]

#Reservoir Gage Heights Over Time
plt.figure()
plt.scatter(LakeLevels.WWatDate,LakeLevels.GageMeanHeight_ft,marker='.', c='r', label='Wes Watkins Reservoir')
plt.scatter(LakeLevels.ShawDate,LakeLevels.GageMeanHeight_ft,marker='.', c='b', label= 'Shawnee Reservoir')
plt.title('Reservoir Gage Heights Over Time')
plt.xlabel('Date')
plt.ylabel('Mean Gage Height (ft)')
plt.legend() 



#Gage Heights vs Seismicity Levels

cat=catdf_USGS
min_mag=2.5
catdf = (cat)

catdf3 = catdf[catdf['magnitude']>=min_mag]
   
m3eqcount = catdf3['origintime'].groupby(catdf3.origintime.dt.to_period("M")).agg('count')
m3eqcountd = catdf3['origintime'].groupby(catdf3.origintime.dt.to_period("D")).agg('count')
alleqcount = catdf['origintime'].groupby(catdf.origintime.dt.to_period("M")).agg('count')
alleqcountd = catdf['origintime'].groupby(catdf.origintime.dt.to_period("D")).agg('count')
   
df3 = m3eqcount.to_frame()
df3d = m3eqcountd.to_frame()
dfall = alleqcount.to_frame()
dfalld = alleqcountd.to_frame()
   
df3.index = df3.index.to_timestamp()
df3d.index = df3d.index.to_timestamp()
dfall.index = dfall.index.to_timestamp()
dfalld.index = dfalld.index.to_timestamp()
   
df3 = df3.resample('MS').sum()
df3d = df3d.resample('D').sum()
dfall = dfall.resample('MS').sum()
dfalld = dfalld.resample('D').sum()



#Gage Height vs Events per Month
x=dfall.index
y1=dfall.origintime
y2=LakeLevels.GageMeanHeight_ft


fig, ax1 = plt.subplots()
plt.plot(dfall.index, dfall.origintime, color='k', label='Seismicity Rate')
ax1.set_ylabel('Earthquakes per month')
plt.legend() 
ax2 = ax1.twinx()
plt.scatter(LakeLevels.WWatDate, y2, color='r',marker='.',label='Wes Watkins Reservoir')
plt.scatter(LakeLevels.ShawDate, y2, color='b',marker='.',label='Shawnee Reservoir')

ax2.set_ylabel('Mean Gage Height (ft)')

plt.xlim([datetime.date(2009, 1, 1), datetime.date(2021,1,1)])

plt.title('Earthquakes per Month ')
plt.legend()
plt.show()


#GageHeight vs Events per day
fig, ax1 = plt.subplots()
plt.plot(dfalld.index,dfalld.origintime, color='k', label='Seismicity Rate')
ax1.set(ylabel = 'Earthquakes per day')
plt.legend()
ax2 = ax1.twinx()
plt.scatter(LakeLevels.WWatDate, y2, color='r',marker='.',label='Wes Watkins Reservoir')
plt.scatter(LakeLevels.ShawDate, y2, color='b',marker='.',label='Shawnee Reservoir')

ax2.set_ylabel('Mean Gage Height (ft)')

plt.xlim([datetime.date(2009, 1, 1), datetime.date(2015,1,1)])

plt.title('Earthquakes per Day vs Reservoir Gage Height')
plt.legend()
plt.show()
 


#Events per Month >Mc
fig, ax1=plt.subplots()
plt.plot(df3.index,df3.origintime, color='k', label='Seismicity Rate')
ax1.set(ylabel = 'Earthquakes M>'+str(min_mag)+' per month')
plt.legend()

plt.title('Earthquakes per Month >2.5')
plt.legend()
plt.show()


#Events per day >Mc
fig, ax1=plt.subplots()
plt.plot(df3d.index,df3d.origintime, color='k', label='Seismicity Rate')
ax1.set(ylabel = 'Earthquakes M>'+str(min_mag)+' per day')
plt.legend()

plt.title('Earthquakes per Day >2.5')
plt.legend()
plt.show()
  
