"""
Created on Thu Apr 7 2022

Analysis carried out from: Estimating the Open Solar Flux from In Situ Measurements

@author: Anna M Frost
         a.m.frost@pgr.reading.ac.uk
"""

import numpy as np
import xarray as xr
import datetime
from datetime import timedelta
import pandas as pd
import matplotlib.pyplot as plt
import pdb #pdb.set_trace()
from matplotlib.collections import LineCollection
plt.rcParams["font.family"] = "Times New Roman"




OSF_Pagel_txt = np.genfromtxt("TopologyTextFiles/OSF_Wind_Pagel.txt", delimiter = None)

CS_Pagel_txt = np.genfromtxt("TopologyTextFiles/N_CL_Wind_Pagel.txt", delimiter = None)

Undetermined_Pagel_txt = np.genfromtxt("TopologyTextFiles/N_U_Wind_Pagel.txt", delimiter = None)

opposite_Pagel_txt = np.genfromtxt("TopologyTextFiles/opposite_Wind_Pagel.txt", delimiter = None)

bgboundary_Pagel_txt = np.genfromtxt("TopologyTextFiles/bgboundary_Wind_Pagel.txt", delimiter = None)


x_Pagel = np.array((opposite_Pagel_txt*100))

y_Pagel = np.array((bgboundary_Pagel_txt*100))


OSF_Pagel = np.array(OSF_Pagel_txt)

CS_Pagel = np.array(CS_Pagel_txt)

Undetermined_Pagel = np.array(Undetermined_Pagel_txt)

idx_Pagel = np.lexsort((y_Pagel, x_Pagel)).reshape(41,41)

Counterstreaming_Pagel = CS_Pagel[idx_Pagel]

Opposite_Pagel = x_Pagel[idx_Pagel]

Background_Pagel = y_Pagel[idx_Pagel]

Undetermined_Pagel = Undetermined_Pagel[idx_Pagel]







###############################################################################
OSF_Best_txt = np.genfromtxt("TopologyTextFiles/OSF_Best.txt", delimiter = None)

CS_Best_txt = np.genfromtxt("TopologyTextFiles/N_CL_Best.txt", delimiter = None)

INV_Best_txt = np.genfromtxt("TopologyTextFiles/N_SS_Best.txt", delimiter = None)

UNINV_Best_txt = np.genfromtxt("TopologyTextFiles/N_AS_Best.txt", delimiter = None)


opposite_Best_txt = np.genfromtxt("TopologyTextFiles/opposite_Best.txt", delimiter = None)

bgboundary_Best_txt = np.genfromtxt("TopologyTextFiles/bgboundary_Best.txt", delimiter = None)

Undetermined_Best_txt = np.genfromtxt("TopologyTextFiles/N_U_Best.txt", delimiter = None)


x_Best = np.array((opposite_Best_txt*100))

y_Best = np.array((bgboundary_Best_txt*100))


OSF_Best = np.array(OSF_Best_txt)

CS_Best = np.array(CS_Best_txt)

Inv_Best = np.array(INV_Best_txt)

Uninv_Best = np.array(UNINV_Best_txt)

Undetermined_Best = np.array(Undetermined_Best_txt)

idx_Best = np.lexsort((y_Best, x_Best)).reshape(41,41)

Counterstreaming_Best = CS_Best[idx_Best]

Inverted_Best = Inv_Best[idx_Best]

Uninverted_Best = Uninv_Best[idx_Best]

Opposite_Best = x_Best[idx_Best]

Background_Best = y_Best[idx_Best]

Undetermined_Best = Undetermined_Best[idx_Best]




###############################################################################
#OSF_min_GoslingACE_txt = np.genfromtxt("E:/Anna Backups/ACE/Strahl Errors/OSF_min_ACE.txt", delimiter = None)
OSF_max_GoslingACE_txt = np.genfromtxt("TopologyTextFiles/OSF_ACE_Gosling.txt", delimiter = None)

#CS_min_GoslingACE_txt = np.genfromtxt("E:/Anna Backups/ACE/Strahl Errors/N_CL_min_ACE.txt", delimiter = None)
CS_max_GoslingACE_txt = np.genfromtxt("TopologyTextFiles/N_CL_ACE_Gosling.txt", delimiter = None)

opposite_GoslingACE_txt = np.genfromtxt("TopologyTextFiles/opposite_ACE_Gosling.txt", delimiter = None)

bgboundary_GoslingACE_txt = np.genfromtxt("TopologyTextFiles/bgboundary_ACE_Gosling.txt", delimiter = None)

#Undetermined_min_GoslingACE_txt = np.genfromtxt("E:/Anna Backups/ACE/Strahl Errors/N_U_min_ACE.txt", delimiter = None)
Undetermined_max_GoslingACE_txt = np.genfromtxt("TopologyTextFiles/N_U_ACE_Gosling.txt", delimiter = None)


x_GoslingACE = np.array((opposite_GoslingACE_txt*100))

y_GoslingACE = np.array((bgboundary_GoslingACE_txt*100))


#OSF_min_GoslingACE = np.array(OSF_min_GoslingACE_txt)
OSF_max_GoslingACE = np.array(OSF_max_GoslingACE_txt)


#CS_min_GoslingACE = np.array(CS_min_GoslingACE_txt)
CS_max_GoslingACE = np.array(CS_max_GoslingACE_txt)

#Undetermined_min_GoslingACE = np.array(Undetermined_min_GoslingACE_txt)
Undetermined_max_GoslingACE = np.array(Undetermined_max_GoslingACE_txt)


idx_GoslingACE = np.lexsort((y_GoslingACE, x_GoslingACE)).reshape(41,41)

#Counterstreaming_minimum_GoslingACE = CS_min_GoslingACE[idx_GoslingACE]
Counterstreaming_maximum_GoslingACE = CS_max_GoslingACE[idx_GoslingACE]

Opposite_GoslingACE = x_GoslingACE[idx_GoslingACE]

Background_GoslingACE = y_GoslingACE[idx_GoslingACE]

Undetermined_GoslingACE = Undetermined_max_GoslingACE[idx_GoslingACE]




###############################################################################

#OSF_min_GoslingWIND_txt = np.genfromtxt("E:/Anna Backups/WIND/Strahl_Errors/OSF_min_WIND.txt", delimiter = None)
OSF_max_GoslingWIND_txt = np.genfromtxt("TopologyTextFiles/OSF_Wind_Gosling.txt", delimiter = None)

#CS_min_GoslingWIND_txt = np.genfromtxt("E:/Anna Backups/WIND/Strahl_Errors/N_CL_min_WIND.txt", delimiter = None)
CS_max_GoslingWIND_txt = np.genfromtxt("TopologyTextFiles/N_CL_Wind_Gosling.txt", delimiter = None)

#Undetermined_min_GoslingWIND_txt = np.genfromtxt("E:/Anna Backups/WIND/Strahl_Errors/N_U_min_WIND.txt", delimiter = None)
Undetermined_max_GoslingWIND_txt = np.genfromtxt("TopologyTextFiles/N_U_Wind_Gosling.txt", delimiter = None)

opposite_GoslingWIND_txt = np.genfromtxt("TopologyTextFiles/opposite_Wind_Gosling.txt", delimiter = None)

bgboundary_GoslingWIND_txt = np.genfromtxt("TopologyTextFiles/bgboundary_Wind_Gosling.txt", delimiter = None)


x_GoslingWIND = np.array((opposite_GoslingWIND_txt*100))

y_GoslingWIND = np.array((bgboundary_GoslingWIND_txt*100))


#OSF_min_GoslingWIND = np.array(OSF_min_GoslingWIND_txt)
OSF_max_GoslingWIND = np.array(OSF_max_GoslingWIND_txt)

#CS_min_GoslingWIND = np.array(CS_min_GoslingWIND_txt)
CS_max_GoslingWIND = np.array(CS_max_GoslingWIND_txt)

#Undetermined_GoslingWIND_min = np.array(Undetermined_min_GoslingWIND_txt)
Undetermined_GoslingWIND_max = np.array(Undetermined_max_GoslingWIND_txt)


idx_GoslingWIND = np.lexsort((y_GoslingWIND, x_GoslingWIND)).reshape(41,41)

#Counterstreaming_minimum_GoslingWIND = CS_min_GoslingWIND[idx_GoslingWIND]
Counterstreaming_maximum_GoslingWIND = CS_max_GoslingWIND[idx_GoslingWIND]

Opposite_GoslingWIND = x_GoslingWIND[idx_GoslingWIND]

Background_GoslingWIND = y_GoslingWIND[idx_GoslingWIND]

Undetermined_GoslingWIND = Undetermined_GoslingWIND_max[idx_GoslingWIND]




###############################################################################

OSF_Anderson_txt = np.genfromtxt("TopologyTextFiles/OSF_ACE_Anderson.txt", delimiter = None)

CS_Anderson_txt = np.genfromtxt("TopologyTextFiles/N_CL_ACE_Anderson.txt", delimiter = None)

Undetermined_Anderson_txt = np.genfromtxt("TopologyTextFiles/N_U_ACE_Anderson.txt", delimiter = None)

opposite_Anderson_txt = np.genfromtxt("TopologyTextFiles/opposite_ACE_Anderson.txt", delimiter = None)

bgboundary_Anderson_txt = np.genfromtxt("TopologyTextFiles/bgboundary_ACE_Anderson.txt", delimiter = None)


x_Anderson = np.array((opposite_Anderson_txt*100))

y_Anderson = np.array((bgboundary_Anderson_txt*100))


OSF_Anderson = np.array(OSF_Anderson_txt)

CS_Anderson = np.array(CS_Anderson_txt)

Undetermined_Anderson = np.array(Undetermined_Anderson_txt)

idx_Anderson = np.lexsort((y_Anderson, x_Anderson)).reshape(41,41)

Counterstreaming_Anderson = CS_Anderson[idx_Anderson]

Opposite_Anderson = x_Anderson[idx_Anderson]

Background_Anderson = y_Anderson[idx_Anderson]

Undetermined_Anderson = Undetermined_Anderson[idx_Anderson]




###############################################################################

OSF_Skoug_txt = np.genfromtxt("TopologyTextFiles/OSF_ACE_Skoug.txt", delimiter = None)

CS_Skoug_txt = np.genfromtxt("TopologyTextFiles/N_CL_ACE_Skoug.txt", delimiter = None)

Undetermined_Skoug_txt = np.genfromtxt("TopologyTextFiles/N_U_ACE_Skoug.txt", delimiter = None)

opposite_Skoug_txt = np.genfromtxt("TopologyTextFiles/opposite_ACE_Skoug.txt", delimiter = None)

bgboundary_Skoug_txt = np.genfromtxt("TopologyTextFiles/bgboundary_ACE_Skoug.txt", delimiter = None)


x_Skoug = np.array((opposite_Skoug_txt*100))

y_Skoug = np.array((bgboundary_Skoug_txt*100))


OSF_Skoug = np.array(OSF_Skoug_txt)

CS_Skoug = np.array(CS_Skoug_txt)

Undetermined_Skoug = np.array(Undetermined_Skoug_txt)


idx_Skoug = np.lexsort((y_Skoug, x_Skoug)).reshape(41,41)

Counterstreaming_Skoug = CS_Skoug[idx_Skoug]

Opposite_Skoug = x_Skoug[idx_Skoug]

Background_Skoug = y_Skoug[idx_Skoug]

Undetermined_Skoug = Undetermined_Skoug[idx_Skoug]




###############################################################################

fig = plt.figure()
ax = fig.add_subplot(111)

ax.set_xlabel('Percentage Above Opposite', fontsize=40)
ax.set_ylabel('Percentage Above Background', fontsize=40)
ax.tick_params(axis='both', which='major', labelsize=40)    #specified size of axis numbers

ax.pcolormesh((np.array(Opposite_Best)-100), (np.array(Background_Best)-100), OSF_Best[idx_Best])  
ax.set_aspect('equal', adjustable='box')

img1 = plt.pcolormesh((np.array(Opposite_Best)-100), (np.array(Background_Best)-100), OSF_Best[idx_Best])   

###############################################################################
#           Plotting 66%
###############################################################################
"""
ax.contour(Opposite_GoslingACE, Background_GoslingACE, Counterstreaming_maximum_GoslingACE, [5,24.4], colors='r')
ax.contourf(Opposite_GoslingACE, Background_GoslingACE, Counterstreaming_maximum_GoslingACE, [0,5], colors='r', alpha=0.5)
ax.contourf(Opposite_GoslingACE, Background_GoslingACE, Counterstreaming_maximum_GoslingACE, [24.4,100], colors='r', alpha=0.5)


ax.contour(Opposite_GoslingWIND, Background_GoslingWIND, Counterstreaming_maximum_GoslingWIND, [5,24.4], colors='b')
ax.contourf(Opposite_GoslingWIND, Background_GoslingWIND, Counterstreaming_maximum_GoslingWIND, [0,5], colors='b', alpha=0.5)
ax.contourf(Opposite_GoslingWIND, Background_GoslingWIND, Counterstreaming_maximum_GoslingWIND, [24.4,100], colors='b', alpha=0.5)


ax.contour(Opposite_Anderson, Background_Anderson, Counterstreaming_Anderson, [3.4,16.6], colors='m')
ax.contourf(Opposite_Anderson, Background_Anderson, Counterstreaming_Anderson, [0,3.4], colors='m', alpha=0.5)
ax.contourf(Opposite_Anderson, Background_Anderson, Counterstreaming_Anderson, [16.6,100], colors='m', alpha=0.5)


ax.contourf(Opposite_Pagel, Background_Pagel, Undetermined_Pagel, [0,3.4], colors='g', alpha=0.5)
ax.contourf(Opposite_Pagel, Background_Pagel, Undetermined_Pagel, [16.6, 100], colors='g', alpha=0.5)
ax.contour(Opposite_Pagel, Background_Pagel, Undetermined_Pagel, [3.4,16.6], colors='g')


ax.contour(Opposite_Skoug, Background_Skoug, Counterstreaming_Skoug, [5.4,26.6], colors='k')
ax.contourf(Opposite_Skoug, Background_Skoug, Counterstreaming_Skoug, [0,5.4], colors='k', alpha=0.5)
ax.contourf(Opposite_Skoug, Background_Skoug, Counterstreaming_Skoug, [26.6,100], colors='k', alpha=0.5)
"""
###############################################################################




###############################################################################
#           Plotting 33%
###############################################################################

Opposite_Pagel = np.array(Opposite_Pagel) - 100
Background_Pagel = np.array(Background_Pagel) - 100

Opposite_GoslingACE = np.array(Opposite_GoslingACE) - 100
Background_GoslingACE = np.array(Background_GoslingACE) - 100

Opposite_GoslingWIND = np.array(Opposite_GoslingWIND) - 100
Background_GoslingWIND = np.array(Background_GoslingWIND) - 100

Opposite_Anderson = np.array(Opposite_Anderson) - 100
Background_Anderson = np.array(Background_Anderson) - 100

Opposite_Skoug = np.array(Opposite_Skoug) - 100
Background_Skoug = np.array(Background_Skoug) - 100


ax.contour(Opposite_Pagel, Background_Pagel, Undetermined_Pagel, [6.7,13.3], colors='y')
ax.contourf(Opposite_Pagel, Background_Pagel, Undetermined_Pagel, [0,6.7], colors='w', alpha=1)
ax.contourf(Opposite_Pagel, Background_Pagel, Undetermined_Pagel, [13.3,100], colors='w', alpha=1)


ax.contour(Opposite_GoslingACE, Background_GoslingACE, Counterstreaming_maximum_GoslingACE, [9.8,19.6], colors='r')
ax.contourf(Opposite_GoslingACE, Background_GoslingACE, Counterstreaming_maximum_GoslingACE, [0,9.8], colors='r', alpha=0.5)
ax.contourf(Opposite_GoslingACE, Background_GoslingACE, Counterstreaming_maximum_GoslingACE, [19.6,100], colors='r', alpha=0.5)


ax.contour(Opposite_GoslingWIND, Background_GoslingWIND, Counterstreaming_maximum_GoslingWIND, [9.8,19.6], colors='b')
ax.contourf(Opposite_GoslingWIND, Background_GoslingWIND, Counterstreaming_maximum_GoslingWIND, [0,9.8], colors='b', alpha=0.5)
ax.contourf(Opposite_GoslingWIND, Background_GoslingWIND, Counterstreaming_maximum_GoslingWIND, [19.6,100], colors='b', alpha=0.5)


ax.contour(Opposite_Anderson, Background_Anderson, Counterstreaming_Anderson, [6.7,13.3], colors='m')
ax.contourf(Opposite_Anderson, Background_Anderson, Counterstreaming_Anderson, [0,6.7], colors='m', alpha=0.5)
ax.contourf(Opposite_Anderson, Background_Anderson, Counterstreaming_Anderson, [13.3,100], colors='w', alpha=1)



ax.contour(Opposite_Skoug, Background_Skoug, Counterstreaming_Skoug, [10.7,21.3], colors='k')
ax.contourf(Opposite_Skoug, Background_Skoug, Counterstreaming_Skoug, [0,10.7], colors='w', alpha=1)
ax.contourf(Opposite_Skoug, Background_Skoug, Counterstreaming_Skoug, [21.3,100], colors='w', alpha=0.5)

###############################################################################


cbar = fig.colorbar(img1)

cbar.set_label(label='OSF (x$10^{14}$ Wb)', size=40)

cbar.ax.tick_params(labelsize=40)

plt.show()



"""
# Legend for above

plt.plot([3,2], color='y', label='Pagel')
plt.plot([2.8,2], color='r', label='Gosling (ACE data)')
plt.plot([2.8,2], color='b', label='Gosling (Wind data)')
plt.plot([2.8,2], color='m', label='Anderson')
plt.plot([2.8,2], color='k', label='Skoug')
plt.legend(fontsize=30)
"""




fig = plt.figure()
ax = fig.add_subplot(111)

#ax.set_title('Open Solar Flux [x$10^{14}$ Wb]', fontsize=35)

ax.set_xlabel('Percentage Above Opposite', fontsize=40)
ax.set_ylabel('Percentage Above Background', fontsize=40)
ax.tick_params(axis='both', which='major', labelsize=40)    #specified size of axis numbers

ax.pcolormesh((np.array(Opposite_Best)-100), (np.array(Background_Best)-100), OSF_Best[idx_Best])  
ax.set_aspect('equal', adjustable='box')

img1 = plt.pcolormesh((np.array(Opposite_Best)-100), (np.array(Background_Best)-100), OSF_Best[idx_Best])
cbar = fig.colorbar(img1)

cbar.set_label(label='OSF (x$10^{14}$ Wb)', size=40)

cbar.ax.tick_params(labelsize=40)

plt.show()

