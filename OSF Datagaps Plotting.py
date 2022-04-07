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
plt.rcParams["font.family"] = "Times New Roman"


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
        magswe = magswe_temp.iloc[(magswe_temp.index >= start_date) & (magswe_temp.index < end_date)]
        return(magswe)


def FilterWINDLow(start_date, end_date):
        start_date = pd.to_datetime(start_date)
        end_date = pd.to_datetime(end_date)
        nc = xr.open_dataset('WIND-L Corrected Time Bow Shock Correction Same Length.nc')
        WIND_temp = nc.to_dataframe()
        WIND = WIND_temp.iloc[(WIND_temp.index >= start_date) & (WIND_temp.index < end_date)]
        return(WIND)
    


def Body(start_date, end_date): 

        WIND = FilterWINDLow(start_date, end_date)
                
        total_data = WIND
        
        timestamp = WIND.index
                
        Br = WIND.Bx
        BSBr = (Br == -9999999.99)
        Br[BSBr] = np.nan
        BSBrSize = np.size(Br[BSBr])
 
        sflux1 = WIND.Bin_1
        sflux2 = WIND.Bin_2
        sflux3 = WIND.Bin_3
        sflux4 = WIND.Bin_4
        sflux5 = WIND.Bin_5
        sflux6 = WIND.Bin_6
        sflux7 = WIND.Bin_7
        sflux8 = WIND.Bin_8
        
        inval1 = (sflux1 == 0.00000000e+00)
        sflux1[inval1] = np.nan
        BSsflux1 = (sflux1 == -9999999.99)
        sflux1[BSsflux1] = np.nan
        BSslux1Size = np.size(sflux1[BSsflux1])
        
        inval2 = (sflux2 == 0.00000000e+00)
        sflux2[inval2] = np.nan
        BSsflux2 = (sflux2 == -9999999.99)
        sflux2[BSsflux2] = np.nan
        BSslux2Size = np.size(sflux2[BSsflux2])

        inval3 = (sflux3 == 0.00000000e+00)
        sflux3[inval3] = np.nan
        BSsflux3 = (sflux3 == -9999999.99)
        sflux3[BSsflux3] = np.nan
        BSslux3Size = np.size(sflux3[BSsflux3])

        inval4 = (sflux4 == 0.00000000e+00)
        sflux4[inval4] = np.nan
        BSsflux4 = (sflux4 == -9999999.99)
        sflux4[BSsflux4] = np.nan
        BSslux4Size = np.size(sflux4[BSsflux4])

        inval5 = (sflux5 == 0.00000000e+00)
        sflux5[inval5] = np.nan
        BSsflux5 = (sflux5 == -9999999.99)
        sflux5[BSsflux5] = np.nan
        BSslux5Size = np.size(sflux5[BSsflux5])

        inval6 = (sflux6 == 0.00000000e+00)
        sflux6[inval6] = np.nan
        BSsflux6 = (sflux6 == -9999999.99)
        sflux6[BSsflux6] = np.nan
        BSslux6Size = np.size(sflux6[BSsflux6])
        
        inval7 = (sflux7 == 0.00000000e+00)
        sflux7[inval7] = np.nan
        BSsflux7 = (sflux7 == -9999999.99)
        sflux7[BSsflux7] = np.nan
        BSslux7Size = np.size(sflux7[BSsflux7])
 
        inval8 = (sflux8 == 0.00000000e+00)
        sflux8[inval8] = np.nan
        BSsflux8 = (sflux8 == -9999999.99)
        sflux8[BSsflux8] = np.nan
        BSslux8Size = np.size(sflux8[BSsflux8])
        
        
        fluxplot = np.array([sflux1, sflux2, sflux3, sflux4, sflux5, sflux6, sflux7, sflux8])
        
        BGflux = np.nanmean([sflux4,sflux5], axis=0)
                
        #Bins 1 and 2
        FluxStrTemp = np.nanmean([sflux1,sflux2], axis=0)  
                      
        #Bins 7 and 8
        FluxAnStrTemp = np.nanmean([sflux7,sflux8], axis=0)
        
        # Missing
        missing = (np.invert(np.isfinite(FluxStrTemp/BGflux))) | (np.invert(np.isfinite(FluxAnStrTemp/BGflux)))
        
        
        Br_missing = Br[missing]
                
        return(total_data, Br_missing)
        



###############################################################################
"""
Reading in the true errors text files.
"""
###############################################################################

True_Plot_Perc_Error = np.genfromtxt("True OSF Difference (ACE greater 95).txt", delimiter = None)
True_WINDPercentageAvailable = np.genfromtxt("True OSF WIND Percentage (ACE greater 95).txt", delimiter = None)
True_ACEPercentageAvailable = np.genfromtxt("True OSF ACE Percentage (ACE greater 95).txt", delimiter = None)




###############################################################################
"""
Binning the true errors in 5 degree bins.
"""
###############################################################################

Error = []

x = []

temp = np.arange(0,100.1,0.1)

for i in range(0, len(temp)):
        error = np.nanmean(True_Plot_Perc_Error[(True_WINDPercentageAvailable<(temp[i]+0.1))&(True_WINDPercentageAvailable>=temp[i])])
        print('i=', i)
        #pdb.set_trace()
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
            i+=1
        else:
            ERROR.append(error)
            j.append(temps[i])
            i+=1

            
x = np.array(x)/10





###############################################################################
"""
Plotting the binned errors with 3rd order best fit.
"""
###############################################################################

fig = plt.figure()
ax = fig.add_subplot(111)

ax.set_xlabel('Data Available (%)', fontsize=40)
ax.set_ylabel('Error in OSF (%)', fontsize=40)

ax.tick_params(axis='both', which='major', labelsize=30)

y = Error

coefficients = np.polyfit(x, y, 3)

poly = np.poly1d(coefficients)


new_x = np.linspace(x[0], x[-1])

new_y = poly(new_x)

ax.plot(new_x, new_y, color='k')
ax.scatter(j, ERROR, marker='x', color='k', s=300)

ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')

