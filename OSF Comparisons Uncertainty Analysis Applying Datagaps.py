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
from astropy import units as u
import requests
import os
import cdflib
import cdflib.epochs as epochs #this is the part of cdflib that is useful for converting out of 'cdf time' into something sensible
import astropy.time as as_time
import astropy
import pdb
import os.path



###############################################################################
"""
Applies datagaps from WIND to all ACE intervals with >95% data availability 
and calculates the percentage error.
"""
###############################################################################

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


def WINDPercentage(start_date, end_date): 
        WIND = FilterWINDLow(start_date, end_date)
                
        total_data = WIND
                        
        Br = WIND.Bx
        BSBr = (Br == -9999999.99)
        Br[BSBr] = np.nan
 
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
        
        inval2 = (sflux2 == 0.00000000e+00)
        sflux2[inval2] = np.nan
        BSsflux2 = (sflux2 == -9999999.99)
        sflux2[BSsflux2] = np.nan
        
        inval3 = (sflux3 == 0.00000000e+00)
        sflux3[inval3] = np.nan
        BSsflux3 = (sflux3 == -9999999.99)
        sflux3[BSsflux3] = np.nan
        
        inval4 = (sflux4 == 0.00000000e+00)
        sflux4[inval4] = np.nan
        BSsflux4 = (sflux4 == -9999999.99)
        sflux4[BSsflux4] = np.nan
        
        inval5 = (sflux5 == 0.00000000e+00)
        sflux5[inval5] = np.nan
        BSsflux5 = (sflux5 == -9999999.99)
        sflux5[BSsflux5] = np.nan
        
        inval6 = (sflux6 == 0.00000000e+00)
        sflux6[inval6] = np.nan
        BSsflux6 = (sflux6 == -9999999.99)
        sflux6[BSsflux6] = np.nan
        
        inval7 = (sflux7 == 0.00000000e+00)
        sflux7[inval7] = np.nan
        BSsflux7 = (sflux7 == -9999999.99)
        sflux7[BSsflux7] = np.nan
        
        inval8 = (sflux8 == 0.00000000e+00)
        sflux8[inval8] = np.nan
        BSsflux8 = (sflux8 == -9999999.99)
        sflux8[BSsflux8] = np.nan
        

        BGflux = np.nanmean([sflux4,sflux5], axis=0)
        
        FluxStrTemp = np.nanmean([sflux1, sflux2], axis=0)  
                              
        FluxAnStrTemp = np.nanmean([sflux7, sflux8], axis=0)
                         
        # Missing
        missing = (np.invert(np.isfinite(FluxStrTemp/BGflux))) | (np.invert(np.isfinite(FluxAnStrTemp/BGflux)))
        Br_missing = Br[missing]
        
        PercentageWIND = ((len(total_data)-len(Br_missing))/len(total_data))*100
                
        return(PercentageWIND)


def ACEPercentage(start_date, end_date): 
        swepam = FilterElectron(start_date, end_date)
        magswe = FilterMagnetic(start_date, end_date)
        
        total_data = magswe         
                
        Br = magswe.Bx

        sflux1 = swepam.Bin_1
        inval1 = (sflux1 == 0.00000000e+00)
        sflux1[inval1] = np.nan
        
        sflux2 = swepam.Bin_2
        inval2 = (sflux2 == 0.00000000e+00)
        sflux2[inval2] = np.nan
        
        sflux3 = swepam.Bin_3
        inval3 = (sflux3 == 0.00000000e+00)
        sflux3[inval3] = np.nan
        
        sflux4 = swepam.Bin_4
        inval4 = (sflux4 == 0.00000000e+00)
        sflux4[inval4] = np.nan
        
        sflux5 = swepam.Bin_5
        inval5 = (sflux5 == 0.00000000e+00)
        sflux5[inval5] = np.nan
        
        sflux6 = swepam.Bin_6
        inval6 = (sflux6 == 0.00000000e+00)
        sflux6[inval6] = np.nan
        
        sflux7 = swepam.Bin_7
        inval7 = (sflux7 == 0.00000000e+00)
        sflux7[inval7] = np.nan
        
        sflux8 = swepam.Bin_8   
        inval8 = (sflux8 == 0.00000000e+00)
        sflux8[inval8] = np.nan
        
        sflux9 = swepam.Bin_9
        inval9 = (sflux9 == 0.00000000e+00)
        sflux9[inval9] = np.nan
        
        sflux10 = swepam.Bin_10
        inval10 = (sflux10 == 0.00000000e+00)
        sflux10[inval10] = np.nan
        
        sflux11 = swepam.Bin_11
        inval11 = (sflux11 == 0.00000000e+00)
        sflux11[inval11] = np.nan
        
        sflux12 = swepam.Bin_12
        inval12 = (sflux12 == 0.00000000e+00)
        sflux12[inval12] = np.nan
        
        sflux13 = swepam.Bin_13
        inval13 = (sflux13 == 0.00000000e+00)
        sflux13[inval13] = np.nan
        
        sflux14 = swepam.Bin_14
        inval14 = (sflux14 == 0.00000000e+00)
        sflux14[inval14] = np.nan
        
        sflux15 = swepam.Bin_15
        inval15 = (sflux15 == 0.00000000e+00)
        sflux15[inval15] = np.nan
        
        sflux16 = swepam.Bin_16
        inval16 = (sflux16 == 0.00000000e+00)
        sflux16[inval16] = np.nan
        
        sflux17 = swepam.Bin_17
        inval17 = (sflux17 == 0.00000000e+00)
        sflux17[inval17] = np.nan
        
        sflux18 = swepam.Bin_18  
        inval18 = (sflux18 == 0.00000000e+00)
        sflux18[inval18] = np.nan
         
        sflux19 = swepam.Bin_19
        inval19 = (sflux19 == 0.00000000e+00)
        sflux19[inval19] = np.nan
        
        sflux20 = swepam.Bin_20
        inval20 = (sflux20 == 0.00000000e+00)
        sflux20[inval20] = np.nan
        
        BGflux = np.nanmean([sflux9, sflux10, sflux11, sflux12], axis=0)
        
        #Bins 1 and 2
        FluxStrTemp = np.nanmean([sflux1, sflux2, sflux3, sflux4], axis=0)
                      
        #Bins 7 and 8
        FluxAnStrTemp = np.nanmean([sflux17, sflux18, sflux19, sflux20], axis=0)
        
        # Missing
        missing = (np.invert(np.isfinite(FluxStrTemp/BGflux))) | (np.invert(np.isfinite(FluxAnStrTemp/BGflux)))
        
        Br_missing = Br[missing]
        
        PercentageACE = ((len(total_data)-len(Br_missing))/len(total_data))*100
        return(PercentageACE)




def BodywithError(start_date_0, end_date_0, start_date, end_date): 
        WIND = FilterWINDLow(start_date, end_date)
        Br_wind = WIND.Bx
        Opposite=2.45
        BGboundary=1.45
        
        swepam = FilterElectron(start_date_0, end_date_0)
        magswe = FilterMagnetic(start_date_0, end_date_0)
        
        total_data = magswe
        
        timestamp = magswe.index
               
        Br_temp = magswe.Bx
        
        PercentageACE = ACEPercentage(start_date_0, end_date_0)
        

        if PercentageACE > 95:
            Br_i = Br_temp[np.invert(np.isnan(Br_wind.values))]
            Br = Br_i.reindex(magswe.index, fill_value=np.nan)
        else:
            Br = Br_temp
        
        sflux1 = swepam.Bin_1
        inval1 = (sflux1 == 0.00000000e+00)
        sflux1[inval1] = np.nan
        
        sflux2 = swepam.Bin_2
        inval2 = (sflux2 == 0.00000000e+00)
        sflux2[inval2] = np.nan
        
        sflux3 = swepam.Bin_3
        inval3 = (sflux3 == 0.00000000e+00)
        sflux3[inval3] = np.nan
        
        sflux4 = swepam.Bin_4
        inval4 = (sflux4 == 0.00000000e+00)
        sflux4[inval4] = np.nan
        
        sflux5 = swepam.Bin_5
        inval5 = (sflux5 == 0.00000000e+00)
        sflux5[inval5] = np.nan
        
        sflux6 = swepam.Bin_6
        inval6 = (sflux6 == 0.00000000e+00)
        sflux6[inval6] = np.nan
        
        sflux7 = swepam.Bin_7
        inval7 = (sflux7 == 0.00000000e+00)
        sflux7[inval7] = np.nan
        
        sflux8 = swepam.Bin_8   
        inval8 = (sflux8 == 0.00000000e+00)
        sflux8[inval8] = np.nan
        
        sflux9 = swepam.Bin_9
        inval9 = (sflux9 == 0.00000000e+00)
        sflux9[inval9] = np.nan
        
        sflux10 = swepam.Bin_10
        inval10 = (sflux10 == 0.00000000e+00)
        sflux10[inval10] = np.nan
        
        sflux11 = swepam.Bin_11
        inval11 = (sflux11 == 0.00000000e+00)
        sflux11[inval11] = np.nan
        
        sflux12 = swepam.Bin_12
        inval12 = (sflux12 == 0.00000000e+00)
        sflux12[inval12] = np.nan
        
        sflux13 = swepam.Bin_13
        inval13 = (sflux13 == 0.00000000e+00)
        sflux13[inval13] = np.nan
        
        sflux14 = swepam.Bin_14
        inval14 = (sflux14 == 0.00000000e+00)
        sflux14[inval14] = np.nan
        
        sflux15 = swepam.Bin_15
        inval15 = (sflux15 == 0.00000000e+00)
        sflux15[inval15] = np.nan
        
        sflux16 = swepam.Bin_16
        inval16 = (sflux16 == 0.00000000e+00)
        sflux16[inval16] = np.nan
        
        sflux17 = swepam.Bin_17
        inval17 = (sflux17 == 0.00000000e+00)
        sflux17[inval17] = np.nan
        
        sflux18 = swepam.Bin_18  
        inval18 = (sflux18 == 0.00000000e+00)
        sflux18[inval18] = np.nan
         
        sflux19 = swepam.Bin_19
        inval19 = (sflux19 == 0.00000000e+00)
        sflux19[inval19] = np.nan
        
        sflux20 = swepam.Bin_20
        inval20 = (sflux20 == 0.00000000e+00)
        sflux20[inval20] = np.nan
        
        fluxplot = np.array([sflux1, sflux2, sflux3, sflux4, sflux5, sflux6, sflux7, sflux8, sflux9, sflux10, sflux11, sflux12, sflux13, sflux14, sflux15, sflux16, sflux17, sflux18, sflux19, sflux20])
                
        BGflux = np.nanmean([sflux9, sflux10, sflux11, sflux12], axis=0)
        
        #Bins 1 and 2
        FluxStrTemp = np.nanmean([sflux1, sflux2, sflux3, sflux4], axis=0)
                      
        #Bins 7 and 8
        FluxAnStrTemp = np.nanmean([sflux17, sflux18, sflux19, sflux20], axis=0)
        
        #Parallel Strahl
        Parallel = (FluxStrTemp/BGflux >= BGboundary) & (FluxStrTemp/FluxAnStrTemp >= Opposite) | ((FluxStrTemp/BGflux >= BGboundary) & (FluxStrTemp/FluxAnStrTemp < Opposite) & (FluxAnStrTemp/BGflux < BGboundary))
        
        #Antiparallel
        Antiparallel = (FluxAnStrTemp/BGflux >= BGboundary) & (FluxAnStrTemp/FluxStrTemp >= Opposite) | ((FluxAnStrTemp/BGflux >= BGboundary) & (FluxAnStrTemp/FluxStrTemp < Opposite) & (FluxStrTemp/BGflux < BGboundary))
        
        # No Strahl
        Unclassified = ((FluxStrTemp/BGflux < BGboundary) & (FluxAnStrTemp/BGflux < BGboundary))
        
        # Counterstreaming
        CS = (FluxStrTemp/BGflux >= BGboundary) & (FluxAnStrTemp/BGflux >= BGboundary) & (FluxStrTemp/FluxAnStrTemp < Opposite) & (FluxAnStrTemp/FluxStrTemp < Opposite) 
        
        
        # Missing
        missing = (np.invert(np.isfinite(FluxStrTemp/BGflux))) | (np.invert(np.isfinite(FluxAnStrTemp/BGflux)))
        
        uninverted = ((Parallel)&(Br>0))|((Antiparallel)&(Br<0))
        inverted = ((Parallel)&(Br<0))|((Antiparallel)&(Br>0))
        
        Br_CS = Br[CS]
        Br_unclass = Br[Unclassified]
        Br_uninv = Br[uninverted]
        Br_inv = Br[inverted]
        Br_missing = Br[missing]
        return(Br, total_data, fluxplot, timestamp, Br_uninv, Br_inv, Br_CS, Br_unclass, Br_missing)    


def Body(start_date, end_date): 
        Opposite=2.4
        BGboundary=1.45
        swepam = FilterElectron(start_date, end_date)
        magswe = FilterMagnetic(start_date, end_date)
        
        total_data = magswe
        
        timestamp = magswe.index
        
        Br = magswe.Bx

        sflux1 = swepam.Bin_1
        inval1 = (sflux1 == 0.00000000e+00)
        sflux1[inval1] = np.nan
        
        sflux2 = swepam.Bin_2
        inval2 = (sflux2 == 0.00000000e+00)
        sflux2[inval2] = np.nan
        
        sflux3 = swepam.Bin_3
        inval3 = (sflux3 == 0.00000000e+00)
        sflux3[inval3] = np.nan
        
        sflux4 = swepam.Bin_4
        inval4 = (sflux4 == 0.00000000e+00)
        sflux4[inval4] = np.nan
        
        sflux5 = swepam.Bin_5
        inval5 = (sflux5 == 0.00000000e+00)
        sflux5[inval5] = np.nan
        
        sflux6 = swepam.Bin_6
        inval6 = (sflux6 == 0.00000000e+00)
        sflux6[inval6] = np.nan
        
        sflux7 = swepam.Bin_7
        inval7 = (sflux7 == 0.00000000e+00)
        sflux7[inval7] = np.nan
        
        sflux8 = swepam.Bin_8   
        inval8 = (sflux8 == 0.00000000e+00)
        sflux8[inval8] = np.nan
        
        sflux9 = swepam.Bin_9
        inval9 = (sflux9 == 0.00000000e+00)
        sflux9[inval9] = np.nan
        
        sflux10 = swepam.Bin_10
        inval10 = (sflux10 == 0.00000000e+00)
        sflux10[inval10] = np.nan
        
        sflux11 = swepam.Bin_11
        inval11 = (sflux11 == 0.00000000e+00)
        sflux11[inval11] = np.nan
        
        sflux12 = swepam.Bin_12
        inval12 = (sflux12 == 0.00000000e+00)
        sflux12[inval12] = np.nan
        
        sflux13 = swepam.Bin_13
        inval13 = (sflux13 == 0.00000000e+00)
        sflux13[inval13] = np.nan
        
        sflux14 = swepam.Bin_14
        inval14 = (sflux14 == 0.00000000e+00)
        sflux14[inval14] = np.nan
        
        sflux15 = swepam.Bin_15
        inval15 = (sflux15 == 0.00000000e+00)
        sflux15[inval15] = np.nan
        
        sflux16 = swepam.Bin_16
        inval16 = (sflux16 == 0.00000000e+00)
        sflux16[inval16] = np.nan
        
        sflux17 = swepam.Bin_17
        inval17 = (sflux17 == 0.00000000e+00)
        sflux17[inval17] = np.nan
        
        sflux18 = swepam.Bin_18  
        inval18 = (sflux18 == 0.00000000e+00)
        sflux18[inval18] = np.nan
         
        sflux19 = swepam.Bin_19
        inval19 = (sflux19 == 0.00000000e+00)
        sflux19[inval19] = np.nan
        
        sflux20 = swepam.Bin_20
        inval20 = (sflux20 == 0.00000000e+00)
        sflux20[inval20] = np.nan
        
        fluxplot = np.array([sflux1, sflux2, sflux3, sflux4, sflux5, sflux6, sflux7, sflux8, sflux9, sflux10, sflux11, sflux12, sflux13, sflux14, sflux15, sflux16, sflux17, sflux18, sflux19, sflux20])
                
        BGflux = np.nanmean([sflux9, sflux10, sflux11, sflux12], axis=0)
        
        #Bins 1 and 2
        FluxStrTemp = np.nanmean([sflux1, sflux2, sflux3, sflux4], axis=0)
                      
        #Bins 7 and 8
        FluxAnStrTemp = np.nanmean([sflux17, sflux18, sflux19, sflux20], axis=0)
        
        #Parallel Strahl
        Parallel = (FluxStrTemp/BGflux >= BGboundary) & (FluxStrTemp/FluxAnStrTemp >= Opposite) | ((FluxStrTemp/BGflux >= BGboundary) & (FluxStrTemp/FluxAnStrTemp < Opposite) & (FluxAnStrTemp/BGflux < BGboundary))
        
        #Antiparallel
        Antiparallel = (FluxAnStrTemp/BGflux >= BGboundary) & (FluxAnStrTemp/FluxStrTemp >= Opposite) | ((FluxAnStrTemp/BGflux >= BGboundary) & (FluxAnStrTemp/FluxStrTemp < Opposite) & (FluxStrTemp/BGflux < BGboundary))
        
        # No Strahl
        Unclassified = ((FluxStrTemp/BGflux < BGboundary) & (FluxAnStrTemp/BGflux < BGboundary))
        
        # Counterstreaming
        CS = (FluxStrTemp/BGflux >= BGboundary) & (FluxAnStrTemp/BGflux >= BGboundary) & (FluxStrTemp/FluxAnStrTemp < Opposite) & (FluxAnStrTemp/FluxStrTemp < Opposite) 
        
        
        # Missing
        missing = (np.invert(np.isfinite(FluxStrTemp/BGflux))) | (np.invert(np.isfinite(FluxAnStrTemp/BGflux)))
        
        uninverted = ((Parallel)&(Br>0))|((Antiparallel)&(Br<0))
        inverted = ((Parallel)&(Br<0))|((Antiparallel)&(Br>0))
        
        Br_CS = Br[CS]
        Br_unclass = Br[Unclassified]
        Br_uninv = Br[uninverted]
        Br_inv = Br[inverted]
        Br_missing = Br[missing]
        
        return(Br, total_data, fluxplot, timestamp, Br_uninv, Br_inv, Br_CS, Br_unclass, Br_missing)




###############################################################################
"""
Initial run of data to compare with data with datagaps applied.
"""
###############################################################################

Br_z = []
total_data_z = []
fluxplot_z = []
timestamp_z = []
Br_uninv_z = []
Br_inv_z = []
Br_CS_z = []
Br_unclass_z = []
Br_missing_z = []


start_date_0 = datetime.date(1998,1,22)
end_date_0 = datetime.date(2017,12,31)
delta = timedelta(days=27)
while start_date_0 <= end_date_0:
    end_date_temp_0 = start_date_0 + delta
    Br_z_temp, total_data_z_temp, fluxplot_z_temp, timestamp_z_temp, Br_uninv_z_temp, Br_inv_z_temp, Br_CS_z_temp, Br_unclass_z_temp, Br_missing_z_temp = Body(start_date_0, end_date_temp_0)
    Br_z.append(Br_z_temp)
    total_data_z.append(total_data_z_temp)
    fluxplot_z.append(fluxplot_z_temp)
    timestamp_z.append(Br_uninv_z_temp)
    Br_uninv_z.append(Br_uninv_z_temp)
    Br_inv_z.append(Br_inv_z_temp)
    Br_CS_z.append(Br_CS_z_temp)
    Br_unclass_z.append(Br_unclass_z_temp)
    Br_missing_z.append(Br_missing_z_temp)
    start_date_0 += delta
    
    

###############################################################################
"""
Applying datagaps to all ACE intervals with >95% data availability and saving
as .nc files.
"""
###############################################################################

OSF = []
OSF_original = []
OSF_H = []
OSF_L = []
OSF_nocorr = []
Br = []
N_AS = []
N_SS = []
N_CL = []
N_U = []
Perc_Available = []
Perc_BS = []

WindPercentageAvailable = []
Plot_Perc_Error = []


start_date = datetime.date(1998,1,22)
end_date = datetime.date(2017,12,31)


delta = timedelta(days=27)


while start_date <= end_date:
    end_date_temp = start_date + delta
    Wdate = start_date.strftime("%Y%m%d")
    if os.path.isfile('ACE OSF with Wind intervals/ACE OSF Error with Wind Intervals Applied'+Wdate+'.nc'):
        print('File already created: ', start_date)
        start_date += delta
    else:
        print('WIND dates: ', start_date)
        WIND = FilterWINDLow(start_date, end_date_temp)
        Br_temp = WIND.Bx

        temp = np.zeros(np.shape(Br_temp))

        PercentageWIND = WINDPercentage(start_date, end_date_temp)
        
        start_date_0 = datetime.date(1998,1,22)
        end_date_0 = datetime.date(2017,12,31)
        
        if PercentageWIND<=100:
                WindPercentageAvailable.append(PercentageWIND)
                OSF_perc_error = []
        
                i=0
                for i in range(0, 270):
                        print('i =', i)
                        end_date_temp_0 = start_date_0 + delta
                        Br_j, total_data_j, fluxplot_j, timestamp_j, Br_uninv_j, Br_inv_j, Br_CS_j, Br_unclass_j, Br_missing_j = BodywithError(start_date_0, end_date_temp_0, start_date, end_date_temp)
                        
                        Br_z_i = Br_z[i]
                        total_data_z_i = total_data_z[i]
                        fluxplot_z_i = fluxplot_z[i]
                        timestamp_z_i = timestamp_z[i]
                        Br_uninv_z_i = Br_uninv_z[i]
                        Br_inv_z_i = Br_inv_z[i]
                        Br_CS_z_i = Br_CS_z[i]
                        Br_unclass_z_i = Br_unclass_z[i]
                        Br_missing_z_i = Br_missing_z[i]
                        
                        if np.size(timestamp_j) != 0:
                                timestamp_start = timestamp_j[0]
                                
                                date_j = timestamp_start.strftime("%Y%m%d")
                                datetime_j = datetime.date(int(date_j[0:4]), int(date_j[4:6]), int(date_j[6:8]))
                                
                                N_AS_j = len(Br_uninv_j)
                                BR_AS_j = np.nanmean(abs(Br_uninv_j)*(1e-9))
                                
                                N_SS_j = len(Br_inv_j)
                                BR_SS_j = np.nanmean(abs(Br_inv_j)*(1e-9))
                                
                                N_CL_j = len(Br_CS_j)
                                BR_CL_j = np.nanmean(abs(Br_CS_j)*(1e-9))
                                
                                N_U_j = len(Br_unclass_j)
                                BR_U_j = np.nanmean(abs(Br_unclass_j)*(1e-9))
                                
                                N_all_j = np.nansum([N_AS_j, N_SS_j, N_CL_j, N_U_j])
                                
                                
                                
                                N_original_AS_z = len(Br_uninv_z_i)
                                BR_original_AS_z = np.nanmean(abs(Br_uninv_z_i)*(1e-9))
                                
                                N_original_SS_z = len(Br_inv_z_i)
                                BR_original_SS_z = np.nanmean(abs(Br_inv_z_i)*(1e-9))
                                
                                N_original_CL_z = len(Br_CS_z_i)
                                BR_original_CL_z = np.nanmean(abs(Br_CS_z_i)*(1e-9))
                                
                                N_original_U_z = len(Br_unclass_z_i)
                                BR_original_U_z = np.nanmean(abs(Br_unclass_z_i)*(1e-9))
                                
                                N_original_all_z = np.nansum([N_original_AS_z, N_original_SS_z, N_original_CL_z, N_original_U_z])
                                
                                
                                if N_all_j==0:
                                        start_date_0 += delta
                                        i = i+1
                                else:                                   
                                    
                                        TotalMagFlux = ((4*np.pi*((1*149597871000)**2))/N_all_j)*(np.nansum([N_AS_j*BR_AS_j, N_SS_j*BR_SS_j, N_CL_j*BR_CL_j]))
                                        
                                        NotConnected = ((4*np.pi*((1*149597871000)**2))/N_all_j)*(np.nansum([N_SS_j*BR_SS_j]))
                                        
                                        OSF_j = TotalMagFlux - 2*NotConnected
                                        
                                        
                                        TotalMagFlux_original = ((4*np.pi*((1*149597871000)**2))/N_original_all_z)*(np.nansum([N_original_AS_z*BR_original_AS_z, N_original_SS_z*BR_original_SS_z, N_original_CL_z*BR_original_CL_z]))
                        
                                        NotConnected_original = ((4*np.pi*((1*149597871000)**2))/N_original_all_z)*(np.nansum([N_original_SS_z*BR_original_SS_z]))
                                        
                                        OSF_original_z = TotalMagFlux_original - 2*NotConnected_original
                                        
                                        OSF.append(OSF_j)
                                        OSF_original.append(OSF_original_z)
                
                                        
                                        OSF_perc_error_i = (abs(np.array(OSF_original_z) - np.array(OSF_j)))/(np.array(OSF_original_z))*100
                                        OSF_perc_error.append(OSF_perc_error_i)
                                        
                                        #pdb.set_trace()
                                        if start_date_0 == datetime.date(1998,1,22):
                                                    dates = datetime_j
                                                    Br_WindApplied = pd.DataFrame(Br_j)
                                                    i = i+1
                                                    start_date_0 += delta
                               
                                        else:
                                                    dates = np.hstack([dates, datetime_j])
                                                    Br_WindApplied = Br_WindApplied.append(pd.DataFrame(Br_j))
                                                    i = i+1
                                                    start_date_0 += delta
                                                                         
                        else:
                                start_date_0 += delta
                                i = i+1
                
                dates = pd.to_datetime(dates)
                Error_DataFrame = pd.DataFrame(OSF_perc_error, index=dates, columns=['OSF_error'])
                
                Wdate = start_date.strftime("%Y%m%d")
                ds = Error_DataFrame.to_xarray()
                ds.to_netcdf('ACE OSF with WIND intervals\ACE OSF Error with Wind Intervals Applied'+Wdate+'.nc')
                ds.close()
                
                print('dates:', dates)
                print('WIND %:', WindPercentageAvailable)
                print('Mean % Error:' , Plot_Perc_Error)                        

                               
                start_date += delta
        else:
                start_date += delta



