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


###############################################################################
"""
Magnetic field data from each spacecraft that has the highest percentage. Used
in determining the approximate time averaging.
"""
###############################################################################

def FilterACE(start_date, end_date):
        start_date = pd.to_datetime(start_date)
        end_date = pd.to_datetime(end_date)
        nc = xr.open_dataset('ACE Br Same as Electron data.nc')
        magswe_temp = nc.to_dataframe()
        magswe = magswe_temp.iloc[(magswe_temp.index > start_date) & (magswe_temp.index < end_date)]
        return(magswe)


def FilterWIND(start_date, end_date):
        start_date = pd.to_datetime(start_date)
        end_date = pd.to_datetime(end_date)
        nc = xr.open_dataset('WIND Br Same as Electron data.nc')
        WIND_temp = nc.to_dataframe()
        WIND = WIND_temp.iloc[(WIND_temp.index > start_date) & (WIND_temp.index < end_date)]
        return(WIND)


nc = xr.open_dataset('WIND Br Same as Electron data.nc')
WIND_temp = nc.to_dataframe()


start_date = datetime.datetime(1994, 11, 6, 21, 28, 22)
end_date = datetime.datetime(2022,2,15)

delta = timedelta(days=27.27)

BR = []

i = 0
        
while start_date <= end_date:
    #pdb.set_trace()
        end_date_temp = start_date + delta
        print('i = ', i)
        
        if Percentage_ACE[i]>75 or Percentage_WIND[i]>75:
                
                if Percentage_ACE[i] > Percentage_WIND[i]:
                        Br = FilterACE(start_date, end_date_temp)
                        if len(BR)==0:
                            BR = Br
                            i = i+1
                            start_date += delta
                        else:
                            BR = BR.append(Br)
                            i = i+1
                            start_date += delta
                        
                elif Percentage_WIND[i] > Percentage_ACE[i]:
                        Br = FilterWIND(start_date, end_date_temp)
                        if len(BR)==0:
                            BR = Br
                            i = i+1
                            start_date += delta
                        else:
                            BR = BR.append(Br)
                        
                            i = i+1
                            start_date += delta
                
                else:
                        Br_ace = FilterACE(start_date, end_date_temp)
                        Br_wind = FilterWIND(start_date, end_date_temp)
                        mean_Br = np.nanmean(Br_ace, Br_wind)
                        if len(BR)==0:
                            BR = Br
                            i = i+1
                            start_date += delta
                        else:
                            BR = BR.append(mean_Br)
                        
                            i = i+1
                            start_date += delta


        else:
                           
                i = i+1
                start_date += delta



BR = BR.reindex(WIND_temp.index, fill_value=np.nan)

ds = BR.to_xarray()
ds.to_netcdf('Br best times.nc')

