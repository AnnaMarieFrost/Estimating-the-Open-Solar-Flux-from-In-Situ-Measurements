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


def FilterWINDHigh(start_date, end_date):
        start_date = pd.to_datetime(start_date)
        end_date = pd.to_datetime(end_date)
        nc = xr.open_dataset('WIND Corrected Time Valid Both.nc')
        WIND_temp = nc.to_dataframe()
        WIND = WIND_temp.iloc[(WIND_temp.index >= start_date) & (WIND_temp.index < end_date)]
        return(WIND)

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

def Body(start_date, end_date, Case, Opposite, Analyser, BGboundary):
                
        if Analyser==1:
            WIND = FilterWINDHigh(start_date, end_date)
        if Analyser==2:
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
                                
        #Parallel
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
        
        ts = np.zeros(inverted.shape)
        ts[inverted] = -1
        ts[uninverted] = 1
        ts[CS] = 2 
        ts[Unclassified] = -2
        ts[missing] = -3
        
        return(Br, total_data, fluxplot, timestamp, ts, Br_uninv, Br_inv, Br_CS, Br_unclass, Br_missing, BSslux1Size)

     


OSF = []
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
year=[]
month=[]
day=[]
dates = []

# Full time period
start_date = datetime.datetime(1994, 11, 6, 21, 28, 22)
end_date = datetime.datetime(2022,2,15)


Case = 3
Analyser = 2   #1=High, 2=Low
Opposite=2.4
BGboundary=1.45

start_date_0 = start_date


delta = timedelta(days=366)

      
while start_date <= end_date:
        start_year = start_date.year
        start_date = datetime.datetime(start_year, 1, 1)
        end_date_temp = datetime.datetime(start_year, 12, 31)
        print(start_date)
        
        
        Br_i, total_data, fluxplot_i, timestamp_temp, ts_i, Br_uninv_i, Br_inv_i, Br_CS_i, Br_unclass_i, Br_missing_i, BowShock_i = Body(start_date, end_date_temp, Case, Opposite, Analyser, BGboundary)
        
        if len(total_data)>0:
        
                PercentageWINDAvailable = WINDPercentage(start_date, end_date_temp)
                
                if PercentageWINDAvailable>0:
                
                        if np.size(timestamp_temp) != 0:
                                timestamp_start = timestamp_temp[0]
                        
                                date_temp = timestamp_start.strftime("%Y%m%d")
                                datetime_i = datetime.datetime(int(date_temp[0:4]), int(date_temp[4:6]), int(date_temp[6:8]))
                                
                                N_AS_i = len(Br_uninv_i)
                                BR_AS_i = np.mean(abs(Br_uninv_i)*(1e-9))
                                
                                N_SS_i = len(Br_inv_i)
                                BR_SS_i = np.mean(abs(Br_inv_i)*(1e-9))
                                
                                N_CL_i = len(Br_CS_i)
                                BR_CL_i = np.mean(abs(Br_CS_i)*(1e-9))
                                
                                N_U_i = len(Br_unclass_i)
                                BR_U_i = np.mean(abs(Br_unclass_i)*(1e-9))
                                
                                N_all_i = np.nansum([N_AS_i,N_SS_i,N_CL_i,N_U_i])
                                
                                
                                if N_all_i==0:
                                        start_date += delta
                                else:
                                        TotalMagFlux = ((4*np.pi*((1*149597871000)**2))/N_all_i)*(np.nansum([N_AS_i*BR_AS_i, N_SS_i*BR_SS_i, N_CL_i*BR_CL_i]))
                                        
                                        NotConnected = ((4*np.pi*((1*149597871000)**2))/N_all_i)*((N_SS_i*BR_SS_i))
                                        
                                        
                                        OSF_i = TotalMagFlux - 2*NotConnected
                                        
                                        
                                        N_AS.append((N_AS_i/N_all_i)*100)
                                        N_SS.append((N_SS_i/N_all_i)*100)
                                        N_CL.append((N_CL_i/N_all_i)*100)
                                        N_U.append((N_U_i/N_all_i)*100)
                                        OSF.append(OSF_i)
                                        
                                        OSF_nocorr.append(TotalMagFlux)
                                        
                                        Br_temp = np.mean(abs(Br_i))
                                        Br.append(Br_temp)
                                        
                                    
                                        Percentage_BS = (BowShock_i/len(timestamp_temp))*100 
                                        Perc_BS.append(Percentage_BS)
                                        
                                        
                                        Perc_Available.append(PercentageWINDAvailable)
                                        
                                        year.append(start_date.year)
                                        month.append(start_date.month)
                                        day.append(start_date.day)
                                        
                                        if np.size(dates)==0:
                                                    dates = datetime_i
                                        else:
                                                    dates = np.hstack([dates, datetime_i])
                                        
                                        start_date += delta
                        else:
                                start_date += delta                       
                else:
                        start_date += delta
        
        else:
                start_date += delta
        
        
    