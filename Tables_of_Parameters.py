#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 14:35:00 2022

@author: kayceeschaper
"""

from tabulate import tabulate
#DBSCAN trials 
table = [['eps','min_samples','Clusters', 'Noise Points'],[800,5,1412,22519],[700,5,1693,26749],[750,5,1547,24583],
         [775,5,1499,23479],[900,5,1246,19342],[900,7,845,24545],[1000,7,747,21269],[600,7,1360,41873],
         [800,10,675,36333],[1200,10,450,21052],[1300,10,425,18990],[1300,15,282,25656],
         [1300,5,872,12313],[1500,5,742,10327],[700,7,1108,34345],[700,6,1511,29376],[600,6,1776,36089],[600,5,2035,32631],
         [650,5,1826,29447]]
print(tabulate(table, headers='firstrow', tablefmt='fancy_grid'))

#vertical error dropped if greater than n
table = [['n','eps','min_samples','Clusters', 'Noise Points'],[5000,800,5,1476,22093],
         [5000,700,5,1774,26229],[5000,600,5,2050,31967],[5000,750,5,1594,24118],
         [5000,850,5,1400,20407],[5000,900,5,1300,18927]]
print(tabulate(table, headers='firstrow', tablefmt='fancy_grid'))

#Vert error trials ctd; eps850,minsamp5
table = [['n','Clusters','Noise Points'],[5000,1400,20407],[2500,1408,21093],
         [1000,1342,21405],[10000,1293,19500],[7500,1345,19944],[6000,1384,20207],
         [5500,1397,20334],[5250,1394,20384],[4500,1402,20514],[3000,1411,20878],
         [3500,1418,20710]]
print(tabulate(table, headers='firstrow', tablefmt='fancy_grid'))

#RMS filter trials Dropping values outside range
table = [['RMS','Clusters','Noise Points'],['No Filter',1848,26950],['RMS<1',1207,16855],['RMS>0',1631,24584],
         ['1>RMS>0',1196,16855],['10>RMS>0',1631,24584],['8>RMS>0', 1549,23358],['5>RMS>0',1412,20668],
         ['2>RMS>0',1299,18431],['3>RMS>0',1371,19362]]
print(tabulate(table, headers='firstrow', tablefmt='fancy_grid'))
