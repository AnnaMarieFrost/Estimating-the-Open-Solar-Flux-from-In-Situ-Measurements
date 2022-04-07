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




def FilterElectron(start_date, end_date):
        start_date = pd.to_datetime(start_date)
        end_date = pd.to_datetime(end_date)
        nc = xr.open_dataset('ACE Electron Corrected Time Same Length.nc')
        swepam_temp = nc.to_dataframe()
        swepam = swepam_temp.iloc[(swepam_temp.index >= start_date) & (swepam_temp.index < end_date)]
        return(swepam)




def FilterMagnetic(start_date, end_date):
        start_date = pd.to_datetime(start_date)
        end_date = pd.to_datetime(end_date)
        nc = xr.open_dataset('ACE Magnetic Corrected Time Same Length.nc')
        magswe_temp = nc.to_dataframe()
        magswe = magswe_temp.iloc[(magswe_temp.index > start_date) & (magswe_temp.index < end_date)]
        return(magswe)




def FilterWINDLow(start_date, end_date):
        start_date = pd.to_datetime(start_date)
        end_date = pd.to_datetime(end_date)
        nc = xr.open_dataset('WIND-L Corrected Time Bow Shock Correction Same Length.nc')
        WIND_temp = nc.to_dataframe()
        WIND = WIND_temp.iloc[(WIND_temp.index > start_date) & (WIND_temp.index < end_date)]
        return(WIND)
    



swepam_min1 = FilterElectron(start_date=datetime.date(2008,1,1), end_date=datetime.date(2008,12,31))
swepam_min2 = FilterElectron(start_date=datetime.date(2019,1,1), end_date=datetime.date(2019,12,31))

swepam_max1 = FilterElectron(start_date=datetime.date(2000,1,1), end_date=datetime.date(2000,12,31))
swepam_max2 = FilterElectron(start_date=datetime.date(2014,1,1), end_date=datetime.date(2014,12,31))



magswe_min1 = FilterMagnetic(start_date=datetime.date(2008,1,1), end_date=datetime.date(2008,12,31))
magswe_min2 = FilterMagnetic(start_date=datetime.date(2019,1,1), end_date=datetime.date(2019,12,31))

magswe_max1 = FilterMagnetic(start_date=datetime.date(2000,1,1), end_date=datetime.date(2000,12,31))
magswe_max2 = FilterMagnetic(start_date=datetime.date(2014,1,1), end_date=datetime.date(2014,12,31))



WIND_min1 = FilterWINDLow(start_date=datetime.date(2008,1,1), end_date=datetime.date(2008,12,31))
WIND_min2 = FilterWINDLow(start_date=datetime.date(2019,1,1), end_date=datetime.date(2019,12,31))

WIND_max1 = FilterWINDLow(start_date=datetime.date(2000,1,1), end_date=datetime.date(2000,12,31))
WIND_max2 = FilterWINDLow(start_date=datetime.date(2014,1,1), end_date=datetime.date(2014,12,31))





swepam_min = [swepam_min1, swepam_min2]
swepam_min_df = pd.concat(swepam_min)

swepam_max = [swepam_max1, swepam_max2]
swepam_max_df = pd.concat(swepam_max)



magswe_min = [magswe_min1, magswe_min2]
magswe_min_df = pd.concat(magswe_min)

magswe_max = [magswe_max1, magswe_max2]
magswe_max_df = pd.concat(magswe_max)



WIND_min = [WIND_min1, WIND_min2]
WIND_min_df = pd.concat(WIND_min)

WIND_max = [WIND_max1, WIND_max2]
WIND_max_df = pd.concat(WIND_max)





dsmin = swepam_min_df.to_xarray()
dsmin.to_netcdf('ACE Electron Corrected Time Solar Minimum One Year.nc')
dsmin.close()

dsmax = swepam_max_df.to_xarray()
dsmax.to_netcdf('ACE Electron Corrected Time Solar Maximum One Year.nc')
dsmax.close()



dmmin = magswe_min_df.to_xarray()
dmmin.to_netcdf('ACE Magnetic Corrected Time Solar Minimum One Year.nc')
dmmin.close()

dmmax = magswe_max_df.to_xarray()
dmmax.to_netcdf('ACE Magnetic Corrected Time Solar Maximum One Year.nc')
dmmax.close()



dwmin = WIND_min_df.to_xarray()
dwmin.to_netcdf('WIND-L Corrected Time Bow Shock Correction Solar Minimum One Year.nc')
dwmin.close()

dwmax = WIND_max_df.to_xarray()
dwmax.to_netcdf('WIND-L Corrected Time Bow Shock Correction Solar Maximum One Year.nc')
dwmax.close()
