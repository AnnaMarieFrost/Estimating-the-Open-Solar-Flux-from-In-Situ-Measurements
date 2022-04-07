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
import pdb #pdb.set_trace()
import ACE_OSFCalc_27days as ACE
import WIND_OSFcalc_27days as WIND
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"



dates_ACE = ACE.dates
OSF_ACE = ACE.OSF
Percentage_ACE = ACE.Perc_Available
N_AS_ACE = ACE.N_AS
N_SS_ACE = ACE.N_SS
N_CL_ACE = ACE.N_CL
N_U_ACE = ACE.N_U



dates_WIND = WIND.dates
OSF_WIND = WIND.OSF
Percentage_WIND = WIND.Perc_Available
N_AS_WIND = WIND.N_AS
N_SS_WIND = WIND.N_SS
N_CL_WIND = WIND.N_CL
N_U_WIND = WIND.N_U


OSF_ACE = (np.array(OSF_ACE))/1e14

OSF_WIND = (np.array(OSF_WIND))/1e14



Ace_df = pd.DataFrame(data=np.array([OSF_ACE]).T, index=dates_ACE, columns=['OSF'])
Ace_perc_df = pd.DataFrame(data=np.array([Percentage_ACE]).T, index=dates_ACE, columns=['Percentage_ACE'])
Ace_df_topology = pd.DataFrame(data=np.array([N_AS_ACE, N_SS_ACE, N_CL_ACE, N_U_ACE]).T, index=dates_ACE, columns=['N_AS_ACE', 'N_SS_ACE', 'N_CL_ACE', 'N_U_ACE'])

Wind_df = pd.DataFrame(data=np.array([OSF_WIND]).T, index=dates_WIND, columns=['OSF'])
Wind_perc_df = pd.DataFrame(data=np.array([Percentage_WIND]).T, index=dates_WIND, columns=['Percentage_WIND'])
Wind_df_topology = pd.DataFrame(data=np.array([N_AS_WIND, N_SS_WIND, N_CL_WIND, N_U_WIND]).T, index=dates_WIND, columns=['N_AS_WIND', 'N_SS_WIND', 'N_CL_WIND', 'N_U_WIND'])


OSF = []
N_SS = []
N_AS = []
N_CL = []
N_U = []
dates = []

# Full time period
start_date = datetime.datetime(1994, 11, 6, 21, 28, 22)
end_date = datetime.datetime(2022,2,15)

idx = pd.date_range(datetime.datetime(1994, 11, 6, 21, 28, 22), datetime.datetime(2022,2,15), freq='27.27D')

Ace_df = Ace_df.reindex(idx, fill_value=np.nan)
Ace_perc_df = Ace_perc_df.reindex(idx, fill_value=0)
Ace_df_topology = Ace_df_topology.reindex(idx, fill_value=0)

Wind_df = Wind_df.reindex(idx, fill_value=np.nan)
Wind_perc_df = Wind_perc_df.reindex(idx, fill_value=0)
Wind_df_topology = Wind_df_topology.reindex(idx, fill_value=0)


Percentage_ACE = Ace_perc_df.Percentage_ACE.values
Percentage_WIND = Wind_perc_df.Percentage_WIND.values

PercData = []


i = 0

delta = timedelta(days=27.27)


while start_date <= end_date:
        end_date_temp = start_date + delta
        
        Start_date = start_date
        End_date = end_date_temp
        print('i = ', i)
        if Percentage_ACE[i]>50 or Percentage_WIND[i]>50:
                
                if Percentage_ACE[i] > Percentage_WIND[i]:
                        PercData.append(Percentage_ACE[i])
                        best_df = Ace_df.iloc[(Ace_df.index >= Start_date) & (Ace_df.index < End_date)]
                        best_value = best_df.OSF.values
                        topology_df = Ace_df_topology.iloc[(Ace_df_topology.index >= Start_date) & (Ace_df_topology.index < End_date)]
                        dates.append(start_date)
                        if len(best_value)==0:
                            OSF.append(np.nan)
                            
                            N_AS.append(np.nan)
                            N_SS.append(np.nan)
                            N_CL.append(np.nan)
                            N_U.append(np.nan)
                            i = i+1
                            start_date += delta
                        else:
                            OSF.append(float(best_value))
                            
                            N_AS.append(float(topology_df.N_AS_ACE.values))
                            N_SS.append(float(topology_df.N_SS_ACE.values))
                            N_CL.append(float(topology_df.N_CL_ACE.values))
                            N_U.append(float(topology_df.N_U_ACE.values))
                            
                            i = i+1
                            start_date += delta
                        
                elif Percentage_WIND[i] > Percentage_ACE[i]:
                        PercData.append(Percentage_WIND[i])
                        best_df = Wind_df.iloc[(Wind_df.index >= Start_date) & (Wind_df.index < End_date)]
                        best_value = best_df.OSF.values
                        topology_df = Wind_df_topology.iloc[(Wind_df_topology.index >= Start_date) & (Wind_df_topology.index < End_date)]
                        dates.append(start_date)
                        if len(best_value)==0:
                            OSF.append(np.nan)
                            
                            N_AS.append(np.nan)
                            N_SS.append(np.nan)
                            N_CL.append(np.nan)
                            N_U.append(np.nan)
                            i = i+1
                            start_date += delta
                        else:
                            OSF.append(float(best_value))
                            
                            N_AS.append(float(topology_df.N_AS_WIND.values))
                            N_SS.append(float(topology_df.N_SS_WIND.values))
                            N_CL.append(float(topology_df.N_CL_WIND.values))
                            N_U.append(float(topology_df.N_U_WIND.values))
                            
                            i = i+1
                            start_date += delta
                
                else:
                        dates.append(start_date)
                        PercData.append(Percentage_WIND[i])
                        best_ace_df = Ace_df.iloc[(Ace_df.index >= Start_date) & (Ace_df.index < End_date)]
                        best_ace_value = best_ace_df.OSF.values
                        topology_ace_df = Ace_df_topology.iloc[(Ace_df_topology.index >= Start_date) & (Ace_df_topology.index < End_date)]
                        
                        N_AS_ace = topology_ace_df.N_AS_ACE.values
                        N_SS_ace = topology_ace_df.N_SS_ACE.values
                        N_CL_ace = topology_ace_df.N_CL_ACE.values
                        N_U_ace = topology_ace_df.N_U_ACE.values
                        
                        
                        best_wind_df = Wind_df.iloc[(Wind_df.index >= Start_date) & (Wind_df.index < End_date)]
                        best_wind_value = best_wind_df.OSF.values
                        topology_wind_df = Wind_df_topology.iloc[(Wind_df_topology.index >= Start_date) & (Wind_df_topology.index < End_date)]
                        
                        N_AS_wind = topology_wind_df.N_AS_WIND.values
                        N_SS_wind = topology_wind_df.N_SS_WIND.values
                        N_CL_wind = topology_wind_df.N_CL_WIND.values
                        N_U_wind = topology_wind_df.N_U_WIND.values
                        
                        
                        best_value = np.nanmean([best_ace_value, best_wind_value])
                        if len(best_value)==0:
                            OSF.append(np.nan)
                            
                            N_AS.append(np.nan)
                            N_SS.append(np.nan)
                            N_CL.append(np.nan)
                            N_U.append(np.nan)
                            i = i+1
                            start_date += delta
                        else:
                            OSF.append(float(best_value))
                            
                            N_AS.append(float(np.nanmean([N_AS_ace, N_AS_wind])))
                            N_SS.append(float(np.nanmean([N_SS_ace, N_SS_wind])))
                            N_CL.append(float(np.nanmean([N_CL_ace, N_CL_wind])))
                            N_U.append(float(np.nanmean([N_U_ace, N_U_wind])))
                            
                            i = i+1
                            start_date += delta


        else:
                dates.append(start_date)
                PercData.append(np.nan)
                OSF.append(np.nan)
                
                N_AS.append(np.nan)
                N_SS.append(np.nan)
                N_CL.append(np.nan)
                N_U.append(np.nan)
                
                i = i+1
                start_date += delta


N_AS = np.array(N_AS)
N_SS = np.array(N_SS)
N_CL = np.array(N_CL)
N_U = np.array(N_U)


BEST_OSF_DF = pd.DataFrame(data=np.array([OSF]).T, index=dates, columns=['OSF'])

OSF = BEST_OSF_DF.OSF.values
dates = BEST_OSF_DF.index



###############################################################################
"""
Reading in the true errors text files, binning into 5 degree bins.
"""
###############################################################################

True_Plot_Perc_Error = np.genfromtxt("True OSF Difference (ACE greater 95).txt", delimiter = None)
True_WINDPercentageAvailable = np.genfromtxt("True OSF WIND (ACE greater 95).txt", delimiter = None)
True_ACEPercentageAvailable = np.genfromtxt("True OSF ACE Percentage (ACE greater 95).txt", delimiter = None)

Error = []
x = []
temp = np.arange(0,100.1, 0.1)

for i in range(0, len(temp)):
        error = np.nanmean(True_Plot_Perc_Error[(True_WINDPercentageAvailable<(temp[i]+0.1))&(True_WINDPercentageAvailable>=temp[i])])
        print('i=', i)
        if np.isfinite(error)==False:
            i=i+0.1
        else:
            Error.append(error)
            x.append(i)
            i=i+0.1

ERROR = []
temps = np.arange(0,105,5)
j = []
i=0
for i in range(0, len(temps)):
        error = np.nanmean(True_Plot_Perc_Error[(True_WINDPercentageAvailable<(temps[i]+5))&(True_WINDPercentageAvailable>=temps[i])])
        print('i=', i)
        if np.isfinite(error)==False:
            i+=5
        else:
            ERROR.append(error)
            j.append(temps[i])
            i+=5

            
x = np.array(x)/10

fig = plt.figure()
ax = fig.add_subplot(111)

ax.set_xlabel('Data Available (%)', fontsize=30)
ax.set_ylabel('Error in OSF (%)', fontsize=30)

ax.tick_params(axis='both', which='major', labelsize=20)

y = Error

coefficients = np.polyfit(x, y, 3)

poly = np.poly1d(coefficients)

new_x = np.linspace(x[0], x[-1])

new_y = poly(new_x)


###############################################################################
"""
Extracting best fit values for plotting error bars.
"""
###############################################################################

OSFError = []
for j in range(0, len(PercData)):
        ExtrapolatedError = np.interp(PercData[j], new_x, new_y)
        OSFError.append(float(ExtrapolatedError))


minimum = abs((np.array(OSFError)/100) * (np.array(OSF)))
maximum = abs((np.array(OSFError)/100)*(np.array(OSF)))



###############################################################################
"""
Creating text file as supplementary information for Estimating the Open Solar 
Flux from In-Situ Measurements.
"""
###############################################################################
"""
OSF_lower = OSF-minimum
OSF_higher = OSF+maximum
abs_error = np.array(OSFError)

YEARS = dates.year.values
MONTHS = dates.month.values
DAYS = dates.day.values
HOURS = dates.hour.values
MINUTES = dates.minute.values
SECONDS = dates.second.values
OSF[np.isnan(OSF)] = -9999
N_SS[np.isnan(N_SS)] = -9999
N_AS[np.isnan(N_AS)] = -9999
N_CL[np.isnan(N_CL)] = -9999
N_U[np.isnan(N_U)] = -9999
OSF_lower[np.isnan(OSF_lower)] = -9999
OSF_higher[np.isnan(OSF_higher)] = -9999
abs_error[np.isnan(abs_error)] = -9999
np.array(PercData)[np.isnan(np.array(PercData))] = -9999

data = np.array([YEARS, MONTHS, DAYS, HOURS, MINUTES, SECONDS, OSF, OSF_lower, 
                 OSF_higher, abs_error, N_AS, N_SS, N_CL, N_U, PercData])
column_names = ['Year', 'Month', 'Day', 'Hour', 'Minute', 'Second', 
                'OSF (x10^{14} Wb)', 'OSF minimum (x10^{14} Wb)', 
                'OSF maximum (x10^{14} Wb)', 'Absolute OSF error (x10^{14} Wb)', 
                'Open Flux (%)', 'Inverted Flux (%)', 'Newly-Opened Flux (%)', 
                'Disconnected Flux (%)', 'Available Data (%)']

df = pd.DataFrame(data=data.T, columns=column_names)

#np.savetxt('FrostAM_SupplementaryData.txt', df)

with pd.ExcelWriter('output.xlsx') as writer:  
    df.to_excel(writer)

"""

###############################################################################
"""
Code used to produce Figure 6a in Estimating the Open Solar Flux 
from In-Situ Measurements
"""
###############################################################################


fig = plt.figure()
ax = fig.add_subplot(111)

ax.set_ylabel('Open Solar Flux [x$10^{14}$ Wb]', fontsize=40)

ax.tick_params(axis='both', which='major', labelsize=40)
ax.tick_params(axis='x', labelrotation=45)

ax.plot(dates_ACE, OSF_ACE, color='r', label='OSF from ACE')

ax.plot(dates_WIND, OSF_WIND, color='b', label='OSF from WIND')

ax.plot(dates, OSF, color='k', label='OSF')

ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')

ax.legend(fontsize=30)

plt.show() 




###############################################################################
"""
Code used to produce Figure 6b in Estimating the Open Solar Flux 
from In-Situ Measurements
"""
###############################################################################


fig = plt.figure()
ax = fig.add_subplot(111)

ax.set_ylabel('Open Solar Flux [x$10^{14}$ Wb]', fontsize=40)

ax.tick_params(axis='both', which='major', labelsize=40)
ax.tick_params(axis='x', labelrotation=45)

ax.errorbar(dates, OSF, yerr=[minimum, maximum], fmt='o', color='k', alpha=0.5, capsize=4, label='OSF % error')

ax.plot(dates, OSF, color='k', label='OSF')

ax.scatter(dates, OSF, color='k')

ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')

ax.legend(fontsize=30)

plt.show()






###############################################################################
"""
Code used to investigate the occurrence of switchbacks with solar cycle.
"""
###############################################################################

"""
sunspot_txt = np.genfromtxt("C:\PhD Coding\CodeForPaper\SN_m_tot_V2.0.txt", delimiter = None)
year_txt = sunspot_txt[:,0]
month_txt = sunspot_txt[:,1]
#day_txt = sunspot_txt[:,2]
ssnumber_txt = sunspot_txt[:,3]

ss_datetime = []
for i in range(0, year_txt.size):
    ss_datetime.append(datetime.datetime(np.int(year_txt[i]), np.int(month_txt[i]), 1))

SUN_df = pd.DataFrame(data=ssnumber_txt, index=ss_datetime, columns=['SSN'])
sun_df = SUN_df.resample('1D').pad()

x = np.arange(0, 366)
y_temp = N_SS

y = np.array(N_SS)[np.isfinite(y_temp)]
x = np.array(x)[np.isfinite(y_temp)]

coefficients = np.polyfit(x, y, 3)

poly = np.poly1d(coefficients)

new_x = np.linspace(x[0], x[-1])

new_y = poly(new_x)

start_date = datetime.datetime(1994, 11, 6, 21, 28, 22)
end_date = datetime.datetime(2022,2,15)
fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True)
ax1.set_ylabel('% of Inversions', fontsize=30)
ax2.set_ylabel('Sunspot number', fontsize=30)

#Axis number sizing:
ax1.tick_params(axis='both', which='major', labelsize=30)   #specified size of axis numbers
ax2.tick_params(axis='both', which='major', labelsize=30)   #specified size of axis numbers
ax1.plot(dates, N_SS)
ax2.plot(sun_df.iloc[(sun_df.index >= start_date) & (sun_df.index < end_date)])
"""


###############################################################################
"""
Code used to produce Figure 2c in Estimating the Open Solar Flux 
from In-Situ Measurements
"""
###############################################################################


fig = plt.figure()
ax = fig.add_subplot(111)

ax.set_ylabel('Missing Data (%)', fontsize=30)

ax.tick_params(axis='both', which='major', labelsize=30)
ax.tick_params(axis='x', labelrotation=45)

plt.plot((100-Ace_perc_df), color='k', label='ACE')
plt.plot((100-Wind_perc_df), color='r', label='Wind')

ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')

ax.set_ylim([0, 101])

ax.legend(fontsize=30)

plt.show() 



