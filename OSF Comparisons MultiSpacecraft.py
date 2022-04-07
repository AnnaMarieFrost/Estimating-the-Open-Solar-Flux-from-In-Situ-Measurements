"""
Created on Thu Apr 7 2022

Analysis carried out from: Estimating the Open Solar Flux from In Situ Measurements

@author: Anna M Frost
         a.m.frost@pgr.reading.ac.uk
"""

import matplotlib.pyplot as plt
import numpy as np

import datetime # python package for dealing with time
import pandas as pd
import OSF_Best_Estimate_with_Error_Bars_and_Error_Calc_27days as BE
plt.rcParams["font.family"] = "Times New Roman"


dates_temp = BE.dates
dates = BE.dates
OSF = BE.OSF
minimum = BE.minimum
maximum = BE.maximum
PercData = BE.PercData


best_df = pd.DataFrame(data=OSF, index=dates, columns=['OSF'])


OMNI = np.genfromtxt("Open_Flux_Values_OMNI.txt", delimiter = None)
ADAPT = np.genfromtxt("Open_Flux_Values_ADAPT_WSA.txt", delimiter = None)
EUV = np.genfromtxt("Open_Flux_Values_EUV.txt", delimiter = None)
Harvey = np.genfromtxt("Open_Flux_Values_Harvey.txt", delimiter = None)
KPVT = np.genfromtxt("Open_Flux_Values_KPVT_VSM_WSA.txt", delimiter = None)

OMNI_dates_temp = OMNI[:,1]
OMNI_OSF = (OMNI[:,2])*(1e-14)

ADAPT_dates_temp = ADAPT[:,1]
ADAPT_OSF = (ADAPT[:,2])*(1e-14)

EUV_dates_temp = EUV[:,1]
EUV_OSF = (EUV[:,2])*(1e-14)

Harvey_dates_temp = Harvey[:,1]
Harvey_OSF = (Harvey[:,2])*(1e-14)

KPVT_dates_temp = KPVT[:,1]
KPVT_OSF = (KPVT[:,2])*(1e-14)

OMNI_dates = []
ADAPT_dates = []
EUV_dates = []
Harvey_dates = []
KPVT_dates = []


for i in range(0,len(OMNI_dates_temp)):
    dates_temp_OMNI = datetime.datetime.strptime(str(OMNI_dates_temp[i]), "%Y%m%d.0")
    OMNI_dates.append(dates_temp_OMNI)
    
for i in range(0,len(ADAPT_dates_temp)):
    dates_temp_ADAPT = datetime.datetime.strptime(str(ADAPT_dates_temp[i]), "%Y%m%d.0")
    ADAPT_dates.append(dates_temp_ADAPT)

for i in range(0,len(EUV_dates_temp)):
    dates_temp_EUV = datetime.datetime.strptime(str(EUV_dates_temp[i]), "%Y%m%d.0")
    EUV_dates.append(dates_temp_EUV)

for i in range(0,len(Harvey_dates_temp)):
    dates_temp_Harvey = datetime.datetime.strptime(str(Harvey_dates_temp[i]), "%Y%m%d.0")
    Harvey_dates.append(dates_temp_Harvey)
    
for i in range(0,len(KPVT_dates_temp)):
    dates_temp_KPVT = datetime.datetime.strptime(str(KPVT_dates_temp[i]), "%Y%m%d.0")
    KPVT_dates.append(dates_temp_KPVT)


# Dataframes
OMNI_df = pd.DataFrame(data=OMNI_OSF, index=OMNI_dates, columns=['OSF'])
OMNI_OSF_yearly = []
OMNI_dates_yearly = []

for j in range(0,25):
        year = np.array(np.arange(1988,2014))
        start_date = pd.to_datetime(datetime.date(year[j], 1, 1))
        end_date = pd.to_datetime(datetime.date(year[j], 12, 31))
        OMNI_temp = OMNI_df.iloc[(OMNI_df.index > start_date) & (OMNI_df.index < end_date)]
        OMNI_OSF_mean = np.nanmean(OMNI_temp.OSF.values)
        OMNI_OSF_yearly.append(OMNI_OSF_mean)
        OMNI_dates_yearly.append(datetime.date(year[j], 1, 1))

OMNI_yearly = pd.DataFrame(data=OMNI_OSF_yearly, index=OMNI_dates_yearly, columns=['OSF'])




ADAPT_df = pd.DataFrame(data=ADAPT_OSF, index=ADAPT_dates, columns=['OSF'])
ADAPT_OSF_yearly = []
ADAPT_dates_yearly = []

for j in range(0,15):
        year = np.array(np.arange(1998, 2014))
        start_date = pd.to_datetime(datetime.date(year[j], 1, 1))
        end_date = pd.to_datetime(datetime.date(year[j], 12, 31))
        ADAPT_temp = ADAPT_df.iloc[(ADAPT_df.index > start_date) & (ADAPT_df.index < end_date)]
        ADAPT_OSF_mean = np.nanmean(ADAPT_temp.OSF.values)
        ADAPT_OSF_yearly.append(ADAPT_OSF_mean)
        ADAPT_dates_yearly.append(datetime.date(year[j], 1, 1))

ADAPT_yearly = pd.DataFrame(data=ADAPT_OSF_yearly, index=ADAPT_dates_yearly, columns=['OSF'])




EUV_df = pd.DataFrame(data=EUV_OSF, index=EUV_dates, columns=['OSF'])
EUV_OSF_yearly = []
EUV_dates_yearly = []

for j in range(0,7):
        year = np.array(np.arange(2007, 2014))
        start_date = pd.to_datetime(datetime.date(year[j], 1, 1))
        end_date = pd.to_datetime(datetime.date(year[j], 12, 31))
        EUV_temp = EUV_df.iloc[(EUV_df.index > start_date) & (EUV_df.index < end_date)]
        EUV_OSF_mean = np.nanmean(EUV_temp.OSF.values)
        EUV_OSF_yearly.append(EUV_OSF_mean)
        EUV_dates_yearly.append(datetime.date(year[j], 1, 1))

EUV_yearly = pd.DataFrame(data=EUV_OSF_yearly, index=EUV_dates_yearly, columns=['OSF'])




Harvey_df = pd.DataFrame(data=Harvey_OSF, index=Harvey_dates, columns=['OSF'])
Harvey_OSF_yearly = []
Harvey_dates_yearly = []

for j in range(0,13):
        year = np.array(np.arange(1989, 2002))
        start_date = pd.to_datetime(datetime.date(year[j], 1, 1))
        end_date = pd.to_datetime(datetime.date(year[j], 12, 31))
        Harvey_temp = Harvey_df.iloc[(Harvey_df.index > start_date) & (Harvey_df.index < end_date)]
        Harvey_OSF_mean = np.nanmean(Harvey_temp.OSF.values)
        Harvey_OSF_yearly.append(Harvey_OSF_mean)
        Harvey_dates_yearly.append(datetime.date(year[j], 1, 1))

Harvey_yearly = pd.DataFrame(data=Harvey_OSF_yearly, index=Harvey_dates_yearly, columns=['OSF'])




KPVT_df = pd.DataFrame(data=KPVT_OSF, index=KPVT_dates, columns=['OSF'])
KPVT_OSF_yearly = []
KPVT_dates_yearly = []

for j in range(0,15):
        year = np.array(np.arange(1998, 2014))
        start_date = pd.to_datetime(datetime.date(year[j], 1, 1))
        end_date = pd.to_datetime(datetime.date(year[j], 12, 31))
        KPVT_temp = KPVT_df.iloc[(KPVT_df.index > start_date) & (KPVT_df.index < end_date)]
        KPVT_OSF_mean = np.nanmean(KPVT_temp.OSF.values)
        KPVT_OSF_yearly.append(KPVT_OSF_mean)
        KPVT_dates_yearly.append(datetime.date(year[j], 1, 1))

KPVT_yearly = pd.DataFrame(data=KPVT_OSF_yearly, index=KPVT_dates_yearly, columns=['OSF'])





# Results from Owens 2017
test = np.genfromtxt("OSF_CR_Owens2017.dat", delimiter = None)
OSF_owens = test[:,16]
OSF_owens_lower = test[:,17]
OSF_owens_higher = test[:,18]
dates_owens = pd.date_range(datetime.datetime(1998,1,22,12,6,46), datetime.datetime(2011,7,30), freq='27.27D')



###############################################################################
"""
Difference between best estimate and both WSA estimates. 
(1998/02/05 - 2013/03/07)
"""
###############################################################################

ADAPT_temp = ADAPT_df.iloc[(ADAPT_df.index >= datetime.datetime(1994, 11, 5, 21, 28, 22)) & (ADAPT_df.index < datetime.datetime(2013, 3, 8))]
KPVT_temp = KPVT_df.iloc[(KPVT_df.index >= datetime.datetime(1994, 11, 5, 21, 28, 22)) & (KPVT_df.index < datetime.datetime(2013, 3, 8))]

ADAPT_temp = ADAPT_temp.reindex(KPVT_temp.index, fill_value=np.nan)

ADAPT = ADAPT_temp.resample('1D').pad()
KPVT = KPVT_temp.resample('1D').pad()



best_df = pd.DataFrame(data=OSF, index=dates, columns=['OSF'])
best_df = best_df.reindex((pd.date_range(datetime.datetime(1994, 11, 6, 21, 28, 22), datetime.datetime(2021,12,31), freq='27.27D')), fill_value=np.nan)

best_temp = best_df.iloc[(best_df.index >= datetime.datetime(1994, 11, 5, 21, 28, 22)) & (best_df.index < datetime.datetime(2013, 3, 8))]
BEST_temps = best_temp.resample('1D').pad()


Difference_bw_Best_and_WSA = np.nanmean(BEST_temps) / ((np.nanmean(KPVT) + np.nanmean(ADAPT))/2)





###############################################################################
"""
Plotting for Figure 9 in Estimating the Open Solar Flux from In-Situ Measurements Paper
"""
###############################################################################


fig = plt.figure()
ax = fig.add_subplot(111)

plt.ylabel('Open Solar Flux [x$10^{14}$ Wb]', fontsize=40)

#Axis number sizing:
plt.tick_params(axis='both', which='major', labelsize=40)

ax.fill_between(dates, OSF-minimum, OSF+maximum, color='black', alpha=0.5)

plt.plot(KPVT_df.index, KPVT_df.OSF, label='KPVT WSA', color = 'red', linewidth=2)
plt.plot(ADAPT_df.index, ADAPT_df.OSF, label='ADAPT WSA', color = 'orange', linewidth=2)
plt.plot(EUV_df.index, EUV_df.OSF, label='EUV', color = 'cyan', linewidth=2)
plt.plot(Harvey_df.index, Harvey_df.OSF, label='Harvey', color = 'blue', linewidth=2)
plt.plot(dates_owens, OSF_owens, '--', color='m', label='ACE Owens 2017', alpha=0.8, linewidth=2)
ax.plot(dates, OSF, color='k', label='OSF Strahl Method', linewidth=2)
plt.legend(loc=1, fontsize=30)
plt.show()

