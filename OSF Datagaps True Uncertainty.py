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



###############################################################################
"""
Computes the percentage error between each CR of overlap between ACE and WIND.
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


def BodyACE(start_date, end_date): 
        BGboundary = 1.45
        Opposite = 2.4
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


def BodyWIND(start_date, end_date):
        BGboundary = 1.45
        Opposite = 2.4
        
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





ACEPercentageAvailable = []
WINDPercentageAvailable = []
Plot_Perc_Error = []
dates = []

start_date = datetime.date(1998,1,22)
end_date = datetime.date(2017,12,31)


delta = timedelta(days=27)

while start_date <= end_date:
        end_date_temp = start_date + delta
        #print('WIND dates: ', start_date)
        
        Br_ace, total_data_ace, fluxplot_ace, timestamp_ace, Br_uninv_ace, Br_inv_ace, Br_CS_ace, Br_unclass_ace, Br_missing_ace = BodyACE(start_date, end_date_temp)
        Br_wind, total_data_wind, fluxplot_wind, timestamp_wind, Br_uninv_wind, Br_inv_wind, Br_CS_wind, Br_unclass_wind, Br_missing_wind = BodyWIND(start_date, end_date_temp)
        
        PercentageACEAvailable = ((len(total_data_ace)-len(Br_missing_ace))/len(total_data_ace))*100
        PercentageWINDAvailable = ((len(total_data_wind)-len(Br_missing_wind))/len(total_data_wind))*100
        
        if PercentageACEAvailable>95:
            if PercentageWINDAvailable<=100:
                               
                
                N_AS_ace = len(Br_uninv_ace)
                BR_AS_ace = np.nanmean(abs(Br_uninv_ace)*(1e-9))
                
                N_SS_ace = len(Br_inv_ace)
                BR_SS_ace = np.nanmean(abs(Br_inv_ace)*(1e-9))
                
                N_CL_ace = len(Br_CS_ace)
                BR_CL_ace = np.nanmean(abs(Br_CS_ace)*(1e-9))
                
                N_U_ace = len(Br_unclass_ace)
                BR_U_ace = np.nanmean(abs(Br_unclass_ace)*(1e-9))
                
                N_all_ace = N_AS_ace + N_SS_ace + N_CL_ace + N_U_ace
                
                
                N_AS_wind = len(Br_uninv_wind)
                BR_AS_wind = np.nanmean(abs(Br_uninv_wind)*(1e-9))
                
                N_SS_wind = len(Br_inv_wind)
                BR_SS_wind = np.nanmean(abs(Br_inv_wind)*(1e-9))
                
                N_CL_wind = len(Br_CS_wind)
                BR_CL_wind = np.nanmean(abs(Br_CS_wind)*(1e-9))
                
                N_U_wind = len(Br_unclass_wind)
                BR_U_wind = np.nanmean(abs(Br_unclass_wind)*(1e-9))
                
                N_all_wind = np.nansum([N_AS_wind, N_SS_wind, N_CL_wind, N_U_wind])
                
                
                if N_all_ace==0:
                        start_date += delta
                        
                elif N_all_wind==0:
                        start_date += delta
                
                else:
                                            
                        TotalMagFlux_ace = ((4*np.pi*((1*149597871000)**2))/N_all_ace)*(np.nansum([N_AS_ace*BR_AS_ace, N_SS_ace*BR_SS_ace, N_CL_ace*BR_CL_ace]))
                        NotConnected_ace = ((4*np.pi*((1*149597871000)**2))/N_all_ace)*(np.nansum([N_SS_ace*BR_SS_ace]))              
                        OSF_ace = TotalMagFlux_ace - 2*NotConnected_ace
                        
                        TotalMagFlux_wind = ((4*np.pi*((1*149597871000)**2))/N_all_wind)*(np.nansum([N_AS_wind*BR_AS_wind, N_SS_wind*BR_SS_wind, N_CL_wind*BR_CL_wind]))
                        NotConnected_wind = ((4*np.pi*((1*149597871000)**2))/N_all_wind)*(np.nansum([N_SS_wind*BR_SS_wind]))              
                        OSF_wind = TotalMagFlux_wind - 2*NotConnected_wind
                        
                        OSF_difference = ((OSF_ace - OSF_wind)/OSF_ace)*100
                        print('% error: ', OSF_difference)
                        print('WIND dates: ', start_date)
                        Plot_Perc_Error.append(abs(OSF_difference))
                           
                        WINDPercentageAvailable.append(PercentageWINDAvailable)
                        ACEPercentageAvailable.append(PercentageACEAvailable)
                        dates_temp = start_date.strftime("%Y%m%d")
                        dates.append(dates_temp)
                        start_date += delta

            else:
                start_date+=delta

        else:
                start_date+=delta
                



###############################################################################
"""
Saves text files of the true uncertainties.
"""
###############################################################################
"""
np.savetxt('True OSF Difference (ACE greater 95).txt', np.array(Plot_Perc_Error))
np.savetxt('True OSF WIND Percentage (ACE greater 95).txt', np.array(WINDPercentageAvailable))
np.savetxt('True OSF ACE Percentage (ACE greater 95).txt', np.array(ACEPercentageAvailable))
"""