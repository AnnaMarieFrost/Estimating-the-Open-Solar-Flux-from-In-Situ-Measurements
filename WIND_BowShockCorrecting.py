"""
Created on Thu Apr 7 2022

Analysis carried out from: Estimating the Open Solar Flux from In Situ Measurements

@author: Anna M Frost
         a.m.frost@pgr.reading.ac.uk
"""


###############################################################################
"""
This code reads in the Bow Shock Crossings from the txt files, apply to the 
WIND corrected time file, before reindexing the entire dataset to fill the 
missing intervals with NaNs.

Before running, check that the Bow Shock Crossings txt file is the correct name
and location. Also check the target netcdf at the end is correct.
"""
###############################################################################


import numpy as np
import xarray as xr
import datetime
from datetime import timedelta
import pandas as pd
import matplotlib.pyplot as plt
import pdb #pdb.set_trace()


nc_bs = xr.open_dataset('BowShockDates.nc')
bowshock_temp = nc_bs.to_dataframe()

temps = bowshock_temp.index


nc = xr.open_dataset('WIND-L Corrected Time.nc')     # WIND Time Raw
WIND_temp = nc.to_dataframe()


WINDtime = WIND_temp.index
BowShockBool = np.full(len(WINDtime), True, bool)

for i in range(0,len(temps)):
        print('i = ', i, 'of 1123625')
        datestart = temps[i]
        dateend = temps[i] + pd.Timedelta(minutes=10)
        ind_inv = (datestart<=WINDtime) & (WINDtime<dateend)
        ind = np.invert(ind_inv)
        BowShockBool = BowShockBool & ind
  
WIND_temp[BowShockBool==False] = np.nan       #  -9999999.99
        
WIND = WIND_temp
WIND.sort_index()

idx = pd.date_range(WIND.index[0], WIND.index[-1], freq='128S')



WIND = WIND.reindex(idx, fill_value=np.nan) 

ds = WIND.to_xarray()
ds.to_netcdf('C:/PhD Coding/WIND/Data NetCDFs/WIND-L Corrected Time Bow Shock Correction.nc')
