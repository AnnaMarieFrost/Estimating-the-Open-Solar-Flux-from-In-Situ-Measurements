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
        magswe = magswe_temp.iloc[(magswe_temp.index >= start_date) & (magswe_temp.index < end_date)]
        return(magswe)


def FilterWINDLow(start_date, end_date):
        start_date = pd.to_datetime(start_date)
        end_date = pd.to_datetime(end_date)
        nc = xr.open_dataset('WIND-L Corrected Time Bow Shock Correction Same Length.nc')
        WIND_temp = nc.to_dataframe()
        WIND = WIND_temp.iloc[(WIND_temp.index >= start_date) & (WIND_temp.index < end_date)]
        return(WIND)


def BodyACE(start_date, end_date): 
        swepam = FilterElectron(start_date, end_date)
        magswe = FilterMagnetic(start_date, end_date)

        total_data = magswe
        
        timestamp = magswe.index
        
        Br = magswe.Bx.values

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
        
        
        return(Br, timestamp, BGflux, FluxStrTemp, FluxAnStrTemp, total_data)


def BodyWIND(start_date, end_date):
        WIND = FilterWINDLow(start_date, end_date)
        
        total_data = WIND
        
        timestamp = WIND.index
        
        Br = WIND.Bx.values

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
                
        
        return(Br, timestamp, BGflux, FluxStrTemp, FluxAnStrTemp, total_data)


def LoopingFinal():
        BR = []
        TIMESTAMP = []
        BGFLUX = []
        FLUXSTRTEMP = []
        FLUXANSTRTEMP = []
        PercACE = []
        PercWIND = []
        dates = []
        
        start_date = datetime.datetime(1994, 12, 15)
        end_date = datetime.datetime(2022, 2, 15)
        delta = timedelta(days=27)
        
        while start_date <= end_date:
                dates.append(start_date)
                end_date_temp = start_date + delta
                
                Br_ace_i, timestamp_ace_i, BGflux_ace_i, FluxStrTemp_ace_i, FluxAnStrTemp_ace_i, total_data_ace_i = BodyACE(start_date, end_date_temp)
                Br_wind_i, timestamp_wind_i, BGflux_wind_i, FluxStrTemp_wind_i, FluxAnStrTemp_wind_i, total_data_wind_i = BodyWIND(start_date, end_date_temp)
                
                PercentageAvailable_ace = (len(Br_ace_i[np.isfinite(Br_ace_i)])/len(total_data_ace_i))*100
                PercentageAvailable_wind = (len(Br_wind_i[np.isfinite(Br_wind_i)])/len(total_data_wind_i))*100
                
                PercACE.append(PercentageAvailable_ace)
                PercWIND.append(PercentageAvailable_wind)
                
                if PercentageAvailable_ace > PercentageAvailable_wind:
                        BR.append(Br_ace_i)
                        TIMESTAMP.append(timestamp_ace_i)
                        BGFLUX.append(BGflux_ace_i)
                        FLUXSTRTEMP.append(FluxStrTemp_ace_i)
                        FLUXANSTRTEMP.append(FluxAnStrTemp_ace_i)
                        
                        print('start date:', start_date)
                        
                        start_date += delta
        
                
                elif PercentageAvailable_wind > PercentageAvailable_ace:
                        BR.append(Br_wind_i)
                        TIMESTAMP.append(timestamp_wind_i)
                        BGFLUX.append(BGflux_wind_i)
                        FLUXSTRTEMP.append(FluxStrTemp_wind_i)
                        FLUXANSTRTEMP.append(FluxAnStrTemp_wind_i)
                        
                        print('start date:', start_date)
                        
                        start_date += delta
                        
                        
                else:
                        Br_i = []
                        BGflux_i = []
                        FluxStrTemp_i = []
                        FluxAnStrTemp_i = []
                        
                        Br_i = np.nanmean([Br_ace_i, Br_wind_i], axis=0)
                        BGflux_i = np.nanmean([BGflux_ace_i, BGflux_wind_i], axis=0)
                        FluxStrTemp_i = np.nanmean([FluxStrTemp_ace_i, FluxStrTemp_wind_i], axis=0)
                        FluxAnStrTemp_i = np.nanmean([FluxAnStrTemp_ace_i, FluxAnStrTemp_wind_i], axis=0)

                        
                        BR.append(Br_i)
                        TIMESTAMP.append(timestamp_wind_i)
                        BGFLUX.append(BGflux_i)
                        FLUXSTRTEMP.append(FluxStrTemp_i)
                        FLUXANSTRTEMP.append(FluxAnStrTemp_i)
                        
                        
                        print('start date:', start_date)
                        
                        start_date += delta

        return(dates, BR,TIMESTAMP,BGFLUX,FLUXSTRTEMP,FLUXANSTRTEMP, PercACE, PercWIND)
        

dates, BR,TIMESTAMP,BGFLUX,FLUXSTRTEMP,FLUXANSTRTEMP, PercACE, PercWIND = LoopingFinal()
