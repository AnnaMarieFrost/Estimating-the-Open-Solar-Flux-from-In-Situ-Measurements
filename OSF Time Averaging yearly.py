"""
Created on Thu Apr 7 2022

Analysis carried out from: Estimating the Open Solar Flux from In Situ Measurements

@author: Anna M Frost
         a.m.frost@pgr.reading.ac.uk
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import datetime
from datetime import timedelta
import pdb #pdb.set_trace()

import OSF_Best_Estimate_with_Error_Bars_and_Error_Calc_yearly as best


def FilterWINDLow(start_date, end_date, Resolution):
        start_date = pd.to_datetime(start_date)
        end_date = pd.to_datetime(end_date)
        nc = xr.open_dataset('C:\PhD Coding\Br best times.nc')

        WIND_temp = nc.to_dataframe()
        Wind = WIND_temp.resample(str(Resolution)+'H').mean()

        WIND = Wind.iloc[(Wind.index >= start_date) & (Wind.index < end_date)]
        return(WIND)



def LoopingWIND(start_date, end_date, Resolution):
        WIND = FilterWINDLow(start_date, end_date, Resolution)
        
        Br = np.nanmean(abs(WIND.Br.values))*1e-9
        OSF = 4*np.pi*((1.496e11)**2)*Br
        
        OSF = (OSF)*1e-14
        
        return(OSF)



def LoopingFinal(Resolution):
        OSF = []
    
        start_date = datetime.datetime(1994,11,18)
        end_date = datetime.datetime(2022,12,15)
        
        
        delta = timedelta(days=366)
        #pdb.set_trace()
        while start_date <= end_date:
                start_year = start_date.year
                start_date_loop = datetime.datetime(start_year, 1, 1)
                end_date_loop = datetime.datetime(start_year, 12, 31)
                
                OSF_i = LoopingWIND(start_date_loop, end_date_loop, Resolution)
                
                OSF.append(OSF_i)
                
                #print('start date:', start_date_loop)
                print('Resolution = ', Resolution, ' hours')
                start_date += delta
        
        return(OSF)



OSF_best = best.OSF
dates_best = best.dates

df_best = pd.DataFrame(OSF_best, pd.DatetimeIndex(dates_best), columns=['OSF'])
df_best_corr = df_best.resample('1Y').mean()
df_best_corr = df_best





###############################################################################
"""
Code used to determine which time averaging produced the closest match to the 
Strahl method.
"""
###############################################################################


max_hours = 50
DATA = []
OSF_time_averages = []
for k in range(1, (max_hours+1)):
    temp_osf = LoopingFinal(k)
    DATA.append(temp_osf)
    OSF_time_averages.append(np.nanmean(temp_osf))




df_timeavs = pd.DataFrame(data=np.array(DATA).T, index=dates_best)


MIN_diff = []
for i in range(0, len(df_best)):
    Best_OSF_Compare = df_best.values[i:i+1]
    temp = abs(df_timeavs[i:i+1] - Best_OSF_Compare)
    min_diff = float(np.sort(temp.values)[0,0])
    smallest_res = temp[temp==min_diff]
    smallest_res[smallest_res>0] = 1
    MIN_diff.append(min_diff)
    if i==0:
        Best_Time_Av = smallest_res.values
        
    else:
        Best_Time_Av = np.vstack([Best_Time_Av, smallest_res.values])
    
Binned_time_av = []

for j in range(0, max_hours):
    column_temp = Best_Time_Av[:,j]
    sum_temp = np.nansum(column_temp)
    Binned_time_av.append(sum_temp)
  

MAE = []
for n in range(0, max_hours):
    temp = np.nanmean(abs(df_timeavs.values[:,n] - OSF_best))
    MAE.append(temp)


HourAv = np.arange(1,(max_hours+1))




fig = plt.figure()
ax = fig.add_subplot(111)

ax.set_xlabel('Data Time Averaging Length (Hours)', fontsize=40)
ax.set_ylabel('MAE', fontsize=40)
#plt.title('Predicted and True OSF Error', fontsize=35)

ax.tick_params(axis='both', which='major', labelsize=40)

plt.bar(HourAv, MAE)
ax.set_xlim([0, (max_hours+1)])
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')

plt.show()



fig = plt.figure()
ax = fig.add_subplot(111)

ax.set_xlabel('Data Time Averaging Length (Hours)', fontsize=40)
ax.set_ylabel('Number of CRs', fontsize=40)
#plt.title('Predicted and True OSF Error', fontsize=35)

ax.tick_params(axis='both', which='major', labelsize=40)

plt.bar(HourAv, Binned_time_av)
ax.set_xlim([0, (max_hours+1)])
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')

plt.show() 




###############################################################################
"""
Code used to produce Figure 7 in 'Estimating the Open Solar Flux from In-Situ 
Measurements'
"""
###############################################################################

"""

fig = plt.figure()
ax = fig.add_subplot(111)

ax.set_xlabel('Year', fontsize=25)
ax.set_ylabel('Open Solar Flux [x$10^{14}$ Wb]', fontsize=25)

ax.tick_params(axis='both', which='major', labelsize=20)

plt.plot(np.arange(1994,2023), df_best.OSF, label='Strahl best estimate', color='r')

plt.fill_between(np.arange(1994,2023), df_timeavs[0], color='k', alpha=0.2, label='Br 1H Averages')
plt.fill_between(np.arange(1994,2023), df_timeavs[9], color='k', alpha=0.3, label='Br 10H Averages')
plt.fill_between(np.arange(1994,2023), df_timeavs[17], color='k', alpha=0.4, label='Br 18H Averages')
plt.fill_between(np.arange(1994,2023), df_timeavs[19], color='k', alpha=0.5, label='Br 20H Averages')
plt.fill_between(np.arange(1994,2023), df_timeavs[23], color='k', alpha=0.6, label='Br 24H Averages')
plt.fill_between(np.arange(1994,2023), df_timeavs[47], color='k', alpha=0.7, label='Br 48H Averages')

ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')

plt.legend(fontsize=20, loc=3)

plt.show()

"""




###############################################################################
"""
Code used to produce Figure 8a and c in 'Estimating the Open Solar Flux from 
In-Situ Measurements'
"""
###############################################################################

"""
temp_osf = LoopingFinal(20)     # 20 hours giving the closest approximation
temp_osf = np.array(temp_osf, dtype=float)
best_osf = np.array(df_best.OSF.values, dtype=float)

Best_osf = best_osf[(np.isfinite(best_osf)) & (np.isfinite(temp_osf))]
osf_20 = temp_osf[(np.isfinite(best_osf)) & (np.isfinite(temp_osf))]


fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('1-year Averages', fontsize=40)
ax.set_ylabel('OSF, Strahl method (x$10^{14}$ Wb)', fontsize=40)
ax.set_xlabel('OSF, 20-hour Br method (x$10^{14}$ Wb)', fontsize=40)

ax.tick_params(axis='both', which='major', labelsize=40)

plt.scatter(osf_20, Best_osf, color='r', marker='x', s=400)

plt.plot([1,12], [1,12], color='k', label='Ideal fit \n (y = x)')

bestfit = np.poly1d(np.polyfit(osf_20, Best_osf, 1))
x = np.linspace(0, 13, 100)

plt.plot(x, bestfit(x), color='k', linestyle='--', label='Best fit \n (y=1.05x-0.28)')
plt.xlim(2, 11)
plt.ylim(2, 11)
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
plt.legend(fontsize=30)
plt.show()







Perc_Error = ((osf_20-Best_osf)/Best_osf)*100


fig = plt.figure()
ax = fig.add_subplot(111)

ax.set_ylabel('Frequency', fontsize=40)
ax.set_xlabel('Percentage Error (%)', fontsize=40)

ax.tick_params(axis='both', which='major', labelsize=40)

plt.hist(Perc_Error, bins=100, edgecolor='black')

ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')

plt.show()

"""