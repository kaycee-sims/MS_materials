#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 22:24:54 2021

@author: kayceeschaper
"""

#import obspy
import numpy as np
import matplotlib.pyplot as plt
#import easyQuake
import pandas as pd
import datetime


from obspy import read_events
from mpl_toolkits.basemap import Basemap
from easyQuake import simple_cat_df                                                            
from obspy.core import UTCDateTime

from obspy.clients.fdsn import Client

#Defining/Loading Catalogs


    # easyQuake Catalog_2010 
cat_eQuake= read_events('/Users/kayceeschaper/Earthquake_Catalogs/catalog_2010_mag.xml')     
catdf_eQuake=simple_cat_df(cat_eQuake )

#catdf_eQuake['origintime'] = catdf_eQuake.index

catdf2_eQuake= catdf_eQuake[(catdf_eQuake['origintime']>'2010-03-01') & (catdf_eQuake['origintime']<'2010-11-01')]



    # USGS 2010 Catalog
#cat_USGS2010=read_events('/Users/kayceeschaper/Earthquake_Catalogs/catalog_query_USGS_2010.xml')
#catdf_USGS2010=simple_cat_df(cat_USGS2010)


#Fetching the desired catalog from USGS via Obspy Client
starttime= UTCDateTime("2010-03-19 22:07:26.529763")
endtime = UTCDateTime("2010-11-01 04:52:31.193928")
client = Client("USGS")
cat_USGS2010 = client.get_events(starttime=starttime, endtime=endtime, minlatitude=33.,
                         maxlatitude=37.5, minlongitude=-104, maxlongitude=-94, 
                         minmagnitude=1)
catdf_USGS2010=simple_cat_df(cat_USGS2010)
catdf_USGS2010['origintime'] = catdf_USGS2010.index

catdf_USGS_Rev= catdf_USGS2010.reindex(index=catdf_USGS2010.index[::-1])



    # OGS Catalog 
catdf_OGS = pd.read_csv('http://Wichita.ogs.ou.edu/eq/catalog/complete/complete.csv')
catdf_OGS['origintime'] = pd.to_datetime(catdf_OGS['origintime'])
catdf_OGS['magnitude'] = catdf_OGS['magnitude'].str.replace("None", "0").astype(float)
catdf_OGS['magnitude'] = catdf_OGS['magnitude'].astype(float)

#to limit the data to only 2010
catdf_OGS2 = catdf_OGS[(catdf_OGS['origintime']>'2010-03-01') & 
                         (catdf_OGS['origintime']<'2010-11-01')]
#to go back to using the full catalog
#catdf_OGS['origintime']=pd.to_datetime(catdf_OGS['origintime'])

    

    # Keranen Catalog
catdf_Ker=pd.read_excel('/Users/kayceeschaper/Earthquake_Catalogs/Keranen_catalog-3.xlsx')

#catdf_Ker['origintime'] = pd.to_datetime(catdf_Ker['origintime'])
#catdf_Ker2 = catdf_Ker[(catdf_Ker['origintime']>'2010-03-01') & (catdf_Ker['origintime']<'2010-10-01')]

date1 = []
for date in catdf_Ker['Date']:
    date1.append((datetime.datetime(int(str(date)[0:2])+2000,int(str(date)[2:4]),
                                    int(str(date)[4:6]),0,0,0)))
catdf_Ker['origintime']=pd.to_datetime(date1)

catdf_Ker2 = catdf_Ker[(catdf_Ker['origintime']>'2010-03-01') & (catdf_Ker['origintime']<'2010-11-01')]



#easyQuake Big catalog

BigCat = pd.read_csv('/Users/kayceeschaper/HYPODD/bigRelocInputs/bigrelocations/eQuakeBig.csv', header= None) 
BigCat['longitude']=BigCat[2]
BigCat['latitude']=BigCat[1]
BigCat['magnitude']=BigCat[17] 
BigCat['origintime']=pd.to_datetime(BigCat[10])

BigCat2010=BigCat[(BigCat['origintime']>'2010-01-01') & 
                         (BigCat['origintime']<'2011-01-01')]


#adding locations for the stations used to collect Keranen data
KerStatLoc=pd.read_excel('/Users/kayceeschaper/KeranenEtAl2014_Supplemental_Material/KerStationLocations.xlsx')
        

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

#a double loop to make a list object for latitude and logitude data from the OGS station data
latitudes_OGS = []
longitudes_OGS = []
for net in inv:
    for sta in net:
        latitudes_OGS.append(sta.latitude)
        longitudes_OGS.append(sta.longitude)

#client = Client("IRIS")
#starttime = UTCDateTime("2020-01-01")
#endtime = UTCDateTime("2020-01-02")
#lat1=35
#lat2=38
#lon1=-96.75
#lon2=-98
#inventory = client.get_stations(minlatitude=lat1, minlongitude=lon1,maxlatitude=lat2, maxlongitude=lon2,
#                                starttime=starttime,
#                                endtime=endtime)



#Basemap plots for each catalog 


    #easyQuake Catalog_2010 
plt.figure(1, figsize=(8,5))
m = Basemap(projection='lcc', resolution='l', 
            lat_0=35.6, lon_0=-97.2,
            width=.85E5, height=.5E5)

m.scatter(catdf_eQuake.longitude.values,catdf_eQuake.latitude.values,
          #s=catdf_eQuake.magnitude.values**3*10, 
          c=catdf_eQuake.index, latlon=True)
#m.etopo()

cbar=plt.colorbar(pad=0.25)
N_TICKS=12
indexes = [catdf_eQuake['origintime'].iloc[i].strftime('%Y-%m-%d') for i in np.linspace(0,catdf_eQuake.shape[0]-1,N_TICKS).astype(int)] 
cbar.ax.set_yticklabels(indexes)

plt.title('easyQuake 2010 Jones,OK')
plt.xlabel('longitude', labelpad=(40))
plt.ylabel('latitude', labelpad=(40))
parallels = np.arange(-90.,91.,.1)
meridians = np.arange(-180.,181., .25)
m.drawmeridians(meridians,labels=[True,False,False,True], rotation=45)
m.drawparallels(parallels,labels=[False,True,True,False])

#Keranen Station Locations 
#m.scatter(KerStatLoc.Longitude.values,KerStatLoc.Latitude.values,
#         c='red', marker='^', latlon=True, label= 'Keranen Stations')
#plt.scatter(longitudes,latitudes, marker='^',c='black',label='OGS Stations')

#plt.legend()

    #USGS Catalog for 2010
plt.figure(2,figsize=(8,5)) 

m = Basemap(projection='lcc', resolution='l', 
            lat_0=35.6, lon_0=-97.2,
            width=.85E5, height=.5E5)
plt.title('USGS Catalog March-October 2010')
#m.etopo()
m.scatter(catdf_USGS_Rev.longitude.values,catdf_USGS_Rev.latitude.values,
         # s=catdf_USGS2010.magnitude.values**3*10, 
          c=catdf_USGS_Rev.index, latlon=True)
plt.xlabel('longitude', labelpad=(40))
plt.ylabel('latitude', labelpad=(10))
parallels = np.arange(-90.,91.,.1)
meridians = np.arange(-180.,181., .25)
m.drawmeridians(meridians,labels=[True,False,False,True], rotation=45)
m.drawparallels(parallels,labels=[False,True,True,False])

cbar=plt.colorbar(pad=0.25)
N_TICKS=12
indexes = [catdf_USGS_Rev['origintime'].iloc[i].strftime('%Y-%m-%d') for i in np.linspace(0,catdf_USGS_Rev.shape[0]-1,N_TICKS).astype(int)] 
cbar.ax.set_yticklabels(indexes)
#cbar.ax.invert_yaxis() #just flips the whole thing so the colors still go to the same dates which isn't what I want it to do

#Keranen Station Locations 
#m.scatter(KerStatLoc.Longitude.values,KerStatLoc.Latitude.values,
#          c='red', marker='^', latlon=True, label= 'Keranen Stations')

#plt.legend()


    #OGS Catalog 2010
plt.figure(3, figsize=(8,5))
plt.title('OGS Catalog March-October 2010')

m = Basemap(projection='lcc', resolution='l', 
            lat_0=35.6, lon_0=-97.2,
            width=.85E5, height=.5E5)
m.scatter(catdf_OGS2.longitude.values,catdf_OGS2.latitude.values,
          #s=catdf_OGS2.magnitude.values**3*10, 
          c=catdf_OGS2['origintime'], latlon=True)
#m.etopo()

cbar=plt.colorbar(pad=0.15)
N_TICKS=12
indexes = [catdf_OGS2['origintime'].iloc[i].strftime('%Y-%m-%d') for i in np.linspace(0,catdf_OGS2.shape[0]-1,N_TICKS).astype(int)] 
cbar.ax.set_yticklabels(indexes)

plt.xlabel('longitude', labelpad=(40))
plt.ylabel('latitude', labelpad=(10))
parallels = np.arange(-90.,91.,.1)
meridians = np.arange(-180.,181., .25)
m.drawmeridians(meridians,labels=[True,False,False,True], rotation=45)
m.drawparallels(parallels,labels=[False,True,True,False])

#Keranen Station Locations 
#m.scatter(KerStatLoc.Longitude.values,KerStatLoc.Latitude.values,
#          c='red', marker='^', latlon=True, label= 'Keranen Stations')
#plt.legend()


    #Keranen 
plt.figure(4, figsize=(8,5))
plt.title('Keranen Catalog March-October 2010')
m = Basemap(projection='lcc', resolution='l', 
            lat_0=35.6, lon_0=-97.2,
            width=.85E5, height=.5E5)
m.scatter(catdf_Ker2.Longitude.values,catdf_Ker2.Latitude.values,
          c=catdf_Ker2['origintime'], latlon=True)
#m.etopo()
#m.drawstates(color='gray')

cbar=plt.colorbar(pad=0.15)
N_TICKS=12
indexes = [catdf_Ker2['origintime'].iloc[i].strftime('%Y-%m-%d') for i in np.linspace(0,catdf_Ker2.shape[0]-1,N_TICKS).astype(int)] 
cbar.ax.set_yticklabels(indexes)


plt.xlabel('longitude', labelpad=(40))
plt.ylabel('latitude', labelpad=(10))
parallels = np.arange(-90.,91.,.1)
meridians = np.arange(-180.,181., .25)
m.drawmeridians(meridians,labels=[True,False,False,True], rotation=45)
m.drawparallels(parallels,labels=[False,True,True,False])

#Keranen Station Locations 
m.scatter(KerStatLoc.Longitude.values,KerStatLoc.Latitude.values,
          c='red', marker='^', latlon=True, label='Keranen Stations')
plt.legend(loc='upper left')




#Cluster of events that is visible in OGS Catalog and easyQuake Catalog_2010; Not seen in Keranen or USGS Catalogs 
plt.figure(5, figsize=(8,5))
plt.title('easyQuake, OGS, & Keranen March-October 2010')
plt.xlabel('longitude', labelpad=(40))
plt.ylabel('latitude', labelpad=(10))

m = Basemap(projection='lcc', resolution='l', 
            lat_0=35.6, lon_0=-97.2,
            width=.85E5, height=.5E5)

parallels = np.arange(-90.,91.,.1)
meridians = np.arange(-180.,181., .25)
m.drawmeridians(meridians,labels=[True,False,False,True], rotation=45)
m.drawparallels(parallels,labels=[False,True,True,False])


#easyQuake Catalog_2010
m.scatter(catdf_eQuake.longitude.values,catdf_eQuake.latitude.values,
          #s=catdf_eQuake.magnitude.values**3*10, 
          c=catdf_eQuake.index, latlon=True)

#OGS
m.scatter(catdf_OGS2.longitude.values,catdf_OGS2.latitude.values,
         # s=catdf_OGS2.magnitude.values**3*10, 
          c=catdf_OGS2['origintime'], latlon=True)

#Keranen
m.scatter(catdf_Ker2.Longitude.values,catdf_Ker2.Latitude.values,
          c=catdf_Ker2['origintime'], latlon=True)


cbar=plt.colorbar(pad=0.15)
N_TICKS=12
indexes = [catdf_Ker2['origintime'].iloc[i].strftime('%Y-%m-%d') for i in np.linspace(0,catdf_Ker2.shape[0]-1,N_TICKS).astype(int)] 
cbar.ax.set_yticklabels(indexes)

m.scatter(KerStatLoc.Longitude.values,KerStatLoc.Latitude.values,
          c='red', marker='^', latlon=True,label= 'Keranen Stations')
m.scatter(longitudes_OGS,latitudes_OGS,
         c='black', marker='^', latlon=True, label= 'OGS Stations')
plt.legend()




#closer look at March/April 

#easyQuake Catalog_2010 April
catdf_eQuake_April = catdf_eQuake[(catdf_eQuake['origintime']>'2010-03-01') & 
                                    (catdf_eQuake['origintime']<'2010-05-01')]
#OGS April
catdf_OGS2_April = catdf_OGS[(catdf_OGS['origintime']>'2010-03-01') & 
                         (catdf_OGS['origintime']<'2010-05-01')]
#Keranen April
catdf_Ker2_April = catdf_Ker[(catdf_Ker['origintime']>'2010-03-31') & 
                             (catdf_Ker['origintime']<'2010-05-01')]

plt.figure(6, figsize=(8,5))
plt.title('easyQuake, OGS, & Keranen March/April 2010')
plt.xlabel('longitude', labelpad=(40))
plt.ylabel('latitude', labelpad=(10))

m = Basemap(projection='lcc', resolution='l', 
            lat_0=35.6, lon_0=-97.2,
            width=.85E5, height=.5E5)

parallels = np.arange(-90.,91.,.1)
meridians = np.arange(-180.,181., .25)
m.drawmeridians(meridians,labels=[True,False,False,True], rotation=45)
m.drawparallels(parallels,labels=[False,True,True,False])

#easyQuake Catalog_2010 
m.scatter(catdf_eQuake_April.longitude.values,catdf_eQuake_April.latitude.values,
          #s=catdf_eQuake_April.magnitude.values**3*10, 
          c=catdf_eQuake_April.index, latlon=True)

#OGS
m.scatter(catdf_OGS2_April.longitude.values,catdf_OGS2_April.latitude.values,
          #s=catdf_OGS2_April.magnitude.values**3*10, 
          c=catdf_OGS2_April['origintime'], latlon=True)

#Keranen
m.scatter(catdf_Ker2_April.Longitude.values,catdf_Ker2_April.Latitude.values,
          c=catdf_Ker2_April['origintime'], latlon=True)

cbar=plt.colorbar(pad=0.15)
N_TICKS=12
#indexes = [catdf_OGS2_April['origintime'].iloc[i].strftime('%Y-%m-%d') for i in np.linspace(0,catdf_OGS2_April.shape[0]-1,N_TICKS).astype(int)] 
indexes = [catdf_OGS2['origintime'].iloc[i].strftime('%Y-%m-%d') for i in np.linspace(0,catdf_OGS2.shape[0]-1,N_TICKS).astype(int)] 
cbar.ax.set_yticklabels(indexes)
m.scatter(longitudes_OGS,latitudes_OGS,
         c='black', marker='^', latlon=True, label= 'OGS Stations')
#Keranen Station Locations
m.scatter(KerStatLoc.Longitude.values,KerStatLoc.Latitude.values,
          c='red', marker='^', latlon=True, label= 'Keranen Stations')

plt.legend()



#Looking at the cluster around 37.825N 97.375W; longer timeframe, ignoring magnitudes

#easyQuake Catalog_2010; just repeating label from top to get whole catalog even though it only covers 2010
cat_eQuake = read_events('/Users/kayceeschaper/Earthquake_Catalogs/catalog_2010_mag.xml')     
catdf_eQuake=simple_cat_df(cat_eQuake)

catdf_eQuake['origintime'] = catdf_eQuake.index

#OGS
catdf_OGS2_2010to2013 = catdf_OGS[(catdf_OGS['origintime']>'2010-01-01') & 
                         (catdf_OGS['origintime']<'2014-01-01')]
#Keranen
catdf_Ker2_2010to2013 = catdf_Ker[(catdf_Ker['origintime']>'2010-01-01') & 
                             (catdf_Ker['origintime']<'2013-04-27')]#last date in catalog




#easyQuake, OGS, & Keranen 2010-2013 

plt.figure(7, figsize=(8,5))
plt.title('easyQuake, OGS, & Keranen 2010-2013')
plt.xlabel('longitude', labelpad=(40))
plt.ylabel('latitude', labelpad=(10))

m = Basemap(projection='lcc', resolution='l', 
            lat_0=35.6, lon_0=-97.2,
            width=.85E5, height=.75E5)

parallels = np.arange(-90.,91.,.1)
meridians = np.arange(-180.,181., .25)
m.drawmeridians(meridians,labels=[True,False,False,True], rotation=45)
m.drawparallels(parallels,labels=[False,True,True,False])

#easyQuake Catalog_2010
m.scatter(catdf_eQuake.longitude.values,catdf_eQuake.latitude.values,
          #s=catdf_eQuake.magnitude.values**3*10, 
          c=catdf_eQuake.index, latlon=True)

#OGS
m.scatter(catdf_OGS2_2010to2013.longitude.values,catdf_OGS2_2010to2013.latitude.values,
          #s=catdf_OGS2_2010to2013.magnitude.values**3*10, 
          c=catdf_OGS2_2010to2013['origintime'], latlon=True)

#Keranen
m.scatter(catdf_Ker2_2010to2013.Longitude.values,catdf_Ker2_2010to2013.Latitude.values,
          c=catdf_Ker2_2010to2013['origintime'], latlon=True)

#Keranen Station Locations 
m.scatter(KerStatLoc.Longitude.values,KerStatLoc.Latitude.values,
          c='red', marker='^', latlon=True,label= 'Keranen Stations')
m.scatter(longitudes_OGS,latitudes_OGS,
         c='black', marker='^', latlon=True, label= 'OGS Stations')

cbar=plt.colorbar(pad=0.15)
N_TICKS=8
indexes = [catdf_Ker['origintime'].iloc[i].strftime('%Y-%m-%d') for i in np.linspace(0,catdf_Ker.shape[0]-1,N_TICKS).astype(int)] 
cbar.ax.set_yticklabels(indexes)

plt.legend(loc='lower left')



#same data as figure 7 but looking at smaller coordinate range
plt.figure(8, figsize=(8,6))
plt.title('easyQuake, OGS, & Keranen 2010-2013')
plt.xlabel('longitude', labelpad=(40))
plt.ylabel('latitude', labelpad=(10))

m = Basemap(projection='lcc', resolution='l', 
            lat_0=35.85, lon_0=-97.3,
            width=.5E5, height=.6E5)

parallels = np.arange(-90.,91.,.05)
meridians = np.arange(-180.,181., .05)
m.drawmeridians(meridians,labels=[True,False,False,True], rotation=45)
m.drawparallels(parallels,labels=[False,True,True,False])

#easyQuake Catalog_2010
m.scatter(catdf_eQuake.longitude.values,catdf_eQuake.latitude.values,
          #s=catdf_eQuake.magnitude.values**3*10, 
          c=catdf_eQuake.index, latlon=True)

#OGS
m.scatter(catdf_OGS2_2010to2013.longitude.values,catdf_OGS2_2010to2013.latitude.values,
          #s=catdf_OGS2_2010to2013.magnitude.values**3*10, 
          c=catdf_OGS2_2010to2013['origintime'], latlon=True)

#Keranen
m.scatter(catdf_Ker2_2010to2013.Longitude.values,catdf_Ker2_2010to2013.Latitude.values,
          c=catdf_Ker2_2010to2013['origintime'], latlon=True)

cbar=plt.colorbar(pad=0.15)
N_TICKS=6
indexes = [catdf_Ker['origintime'].iloc[i].strftime('%Y-%m-%d') for i in np.linspace(0,catdf_Ker.shape[0]-1,N_TICKS).astype(int)] 
cbar.ax.set_yticklabels(indexes)

#Keranen Station Locations 
m.scatter(KerStatLoc.Longitude.values,KerStatLoc.Latitude.values,
          c='red', marker='^', latlon=True,label= 'Keranen Stations')
#OGS Stations
m.scatter(longitudes_OGS,latitudes_OGS,
         c='black', marker='^', latlon=True, label= 'OGS Stations')
plt.legend(loc='best')




#Fig8 boundaries but only Keranen
plt.figure(9, figsize=(8,6))
plt.title('35.8N 97.4W Keranen 2010-2013')
plt.xlabel('longitude', labelpad=(40))
plt.ylabel('latitude', labelpad=(10))

m = Basemap(projection='lcc', resolution='l', 
            lat_0=35.85, lon_0=-97.3,
            width=.5E5, height=.6E5)

parallels = np.arange(-90.,91.,.05)
meridians = np.arange(-180.,181., .05)
m.drawmeridians(meridians,labels=[True,False,False,True], rotation=45)
m.drawparallels(parallels,labels=[False,True,True,False])

m.scatter(catdf_Ker2_2010to2013.Longitude.values,catdf_Ker2_2010to2013.Latitude.values,
          c=catdf_Ker2_2010to2013['origintime'], latlon=True)

cbar=plt.colorbar(pad=0.15)
N_TICKS=6
indexes = [catdf_Ker['origintime'].iloc[i].strftime('%Y-%m-%d') for i in np.linspace(0,catdf_Ker.shape[0]-1,N_TICKS).astype(int)] 
cbar.ax.set_yticklabels(indexes)

m.scatter(KerStatLoc.Longitude.values,KerStatLoc.Latitude.values,
          c='red', marker='^', latlon=True,label= 'Keranen Stations')
plt.legend(loc="upper right")



#Fig8 boundaries but only OGS 2010
plt.figure(10, figsize=(8,6))
plt.title('35.8N 97.4W OGS')
plt.xlabel('longitude', labelpad=(40))
plt.ylabel('latitude', labelpad=(10))

m = Basemap(projection='lcc', resolution='l', 
            lat_0=35.85, lon_0=-97.3,
            width=.5E5, height=.6E5)

parallels = np.arange(-90.,91.,.05)
meridians = np.arange(-180.,181., .05)
m.drawmeridians(meridians,labels=[True,False,False,True], rotation=45)
m.drawparallels(parallels,labels=[False,True,True,False])

m.scatter(catdf_OGS2_2010to2013.longitude.values,catdf_OGS2_2010to2013.latitude.values,
          #s=catdf_OGS2_2010to2013.magnitude.values**3*10, 
          c=catdf_OGS2_2010to2013['origintime'], latlon=True)

cbar=plt.colorbar(pad=0.15)
N_TICKS=7
indexes = [catdf_OGS2['origintime'].iloc[i].strftime('%Y-%m-%d') for i in np.linspace(0,catdf_OGS2.shape[0]-1,N_TICKS).astype(int)] 
cbar.ax.set_yticklabels(indexes)
m.scatter(longitudes_OGS,latitudes_OGS,
         c='black', marker='^', latlon=True, label= 'OGS Stations')
plt.legend()





#First half(January-July) of 2010; comparison of eQuake and OGS catalogs

#easyQuake
catdf_eQuake_JanJul= catdf2_eQuake[(catdf2_eQuake['origintime']>'2010-01-01') & 
                               (catdf2_eQuake['origintime']<'2010-07-01')]


#OGS
catdf_OGS_JanJul = catdf_OGS[(catdf_OGS['origintime']>'2010-01-01') & 
                               (catdf_OGS['origintime']<'2010-07-01')]

#full coordinate ranges/ all events (within date range) shown
plt.figure(11, figsize=(8,8))
plt.title('easyQuake and OGS Jan-July 2010')
plt.xlabel('longitude', labelpad=(45))
plt.ylabel('latitude', labelpad=(10))


m = Basemap(projection='lcc', resolution='l', 
            lat_0=35.5, lon_0=-97.2,
            width=.5E6, height=.35E6)

parallels = np.arange(-90.,91.,.1)
meridians = np.arange(-180.,181., .25)
m.drawmeridians(meridians,labels=[True,False,False,True], rotation=60)
m.drawparallels(parallels,labels=[False,True,True,False])

m.scatter(catdf_OGS_JanJul.longitude.values,catdf_OGS_JanJul.latitude.values,
          #s=catdf_OGS_JanJul.magnitude.values**3*10, 
          #c=catdf_OGS_JanJul['origintime'], 
          c='red', latlon=True, label='OGS Catalog Events')
m.scatter(catdf_eQuake_JanJul.longitude.values,catdf_eQuake_JanJul.latitude.values,
          #s=catdf_eQuake_JanJul.magnitude.values**3*10, 
          #c=catdf_eQuake_JanJul['origintime'], 
          c='b', latlon=True, label='easyQuake Catalog Events')

m.scatter(longitudes_OGS,latitudes_OGS,
         c='black', marker='^', latlon=True, label= 'OGS Stations')

plt.legend(loc="upper left")

#Jones, OK region; OGS vs eQuake (w/o magnitudes)
plt.figure(12, figsize=(8,8))
plt.title('Jones: easyQuake and OGS Jan-July 2010')
plt.xlabel('longitude', labelpad=(45))
plt.ylabel('latitude', labelpad=(10))

m = Basemap(projection='lcc', resolution='l', 
            lat_0=35.6, lon_0=-97.2,
            width=.85E5, height=.5E5)

parallels = np.arange(-90.,91.,.1)
meridians = np.arange(-180.,181., .25)
m.drawmeridians(meridians,labels=[True,False,False,True], rotation=60)
m.drawparallels(parallels,labels=[False,True,True,False])

m.scatter(catdf_OGS_JanJul.longitude.values,catdf_OGS_JanJul.latitude.values,
          #s=catdf_OGS_JanJul.magnitude.values**3*10, 
          #c=catdf_OGS_JanJul['origintime'], 
          c='green', latlon=True, label='OGS')
m.scatter(catdf_eQuake_JanJul.longitude.values,catdf_eQuake_JanJul.latitude.values,
          #s=catdf_eQuake_JanJul.magnitude.values**3*10, 
          #c=catdf_eQuake_JanJul['origintime'], 
          c='blue', marker='x',latlon=True, label='easyQuake')

m.scatter(longitudes_OGS,latitudes_OGS,
         c='black', marker='^', latlon=True, label= 'OGS Stations')

plt.legend(loc="upper right")


#eQuake and OGS cluster 35.75N 97.45W w/magnitudes
plt.figure(13, figsize=(8,8))
plt.title('easyQuake and OGS Jan-July 2010')
plt.xlabel('longitude', labelpad=(45))
plt.ylabel('latitude', labelpad=(10))

m = Basemap(projection='lcc', resolution='l', 
            lat_0=35.75, lon_0=-97.45,
            width=.85E4, height=.5E4)

parallels = np.arange(-90.,91.,.025)
meridians = np.arange(-180.,181., .025)
m.drawmeridians(meridians,labels=[True,False,False,True], rotation=60)
m.drawparallels(parallels,labels=[False,True,True,False])

m.scatter(catdf_OGS_JanJul.longitude.values,catdf_OGS_JanJul.latitude.values,
          s=catdf_OGS_JanJul.magnitude.values**3*10, 
          #c=catdf_OGS_JanJul['origintime'], 
          c='green', latlon=True, label='OGS')
m.scatter(catdf_eQuake_JanJul.longitude.values,catdf_eQuake_JanJul.latitude.values,
          s=catdf_eQuake_JanJul.magnitude.values**3*10, 
          #c=catdf_eQuake_JanJul['origintime'], 
          c='blue',latlon=True, label='easyQuake')
m.scatter(longitudes_OGS,latitudes_OGS,
         c='black', marker='^', latlon=True, label= 'OGS Stations')

plt.legend(loc="upper right")




#Jones, OK region; OGS vs eQuake (w/magnitudes)
plt.figure(14, figsize=(8,8))
plt.title('Jones: easyQuake and OGS Jan-July 2010')
plt.xlabel('longitude', labelpad=(45))
plt.ylabel('latitude', labelpad=(10))

m = Basemap(projection='lcc', resolution='l', 
            lat_0=35.6, lon_0=-97.2,
            width=.85E5, height=.5E5)

parallels = np.arange(-90.,91.,.1)
meridians = np.arange(-180.,181., .25)
m.drawmeridians(meridians,labels=[True,False,False,True], rotation=60)
m.drawparallels(parallels,labels=[False,True,True,False])

m.scatter(catdf_OGS_JanJul.longitude.values,catdf_OGS_JanJul.latitude.values,
          s=catdf_OGS_JanJul.magnitude.values**3*10, 
          #c=catdf_OGS_JanJul['origintime'], 
          c='green', latlon=True, label='OGS')
m.scatter(catdf_eQuake_JanJul.longitude.values,catdf_eQuake_JanJul.latitude.values,
          s=catdf_eQuake_JanJul.magnitude.values**3*10, 
          #c=catdf_eQuake_JanJul['origintime'], 
          c='blue',latlon=True, label='easyQuake')
m.scatter(longitudes_OGS,latitudes_OGS,
         c='black', marker='^', latlon=True, label= 'OGS Stations')
plt.legend(loc="upper right")





#Jones, OK region; USGS vs eQuakeBig (w/magnitudes)
plt.figure(15, figsize=(8,8))
plt.title('Jones: easyQuake and USGS 2010')
plt.xlabel('longitude', labelpad=(45))
plt.ylabel('latitude', labelpad=(10))

m = Basemap(projection='lcc', resolution='l', 
            lat_0=35.6, lon_0=-97.2,
            width=.85E5, height=.5E5)

parallels = np.arange(-90.,91.,.1)
meridians = np.arange(-180.,181., .25)
m.drawmeridians(meridians,labels=[True,False,False,True], rotation=60)
m.drawparallels(parallels,labels=[False,True,True,False])
m.scatter(BigCat.longitude.values,BigCat.latitude.values,
          #s=BigCat.magnitude.values**3*10, 
          #c=catdf_eQuake_JanJul['origintime'], 
          c='blue',latlon=True, label='easyQuake')
m.scatter(catdf_USGS_Rev.longitude.values,catdf_USGS_Rev.latitude.values,
          #s=catdf_USGS_Rev.magnitude.values**3*10, 
          #c=catdf_OGS_JanJul['origintime'], 
          c='green', latlon=True, label='USGS')

m.scatter(longitudes_OGS,latitudes_OGS,
         c='black', marker='^', latlon=True, label= 'OGS Stations')
plt.legend(loc="upper right")




    # USGS 2010 Catalog
#Fetching the desired catalog from USGS via Obspy Client
starttime= UTCDateTime("2010-01-01 00:00:00.00")
endtime = UTCDateTime("2011-01-01 00:00:00.00")
client = Client("USGS")
cat_USGS2010 = client.get_events(starttime=starttime, endtime=endtime, minlatitude=33.,
                         maxlatitude=37.5, minlongitude=-104, maxlongitude=-94, 
                         minmagnitude=1)
catdf_USGS2010=simple_cat_df(cat_USGS2010)


plt.figure(16, figsize=(8,8))
plt.title('Jones: easyQuake and USGS 2010')
plt.xlabel('origintime', labelpad=(45))
plt.ylabel('magnitude', labelpad=(10))

plt.scatter(catdf_USGS2010.origintime,catdf_USGS2010.magnitude.values,
          c='green', label='USGS')
plt.scatter(BigCat2010.origintime,BigCat2010.magnitude.values,
          marker='+',c='blue', label='easyQuake')


plt.legend(loc="upper right")
