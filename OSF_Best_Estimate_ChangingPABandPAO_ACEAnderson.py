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
    
    

def BodyACE(start_date, end_date): 
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
                
        
        return(Br, timestamp, BGflux, FluxStrTemp, FluxAnStrTemp, total_data)
    

def Body(Br, BGflux, FluxStrTemp, FluxAnStrTemp, BGboundary, Opposite):
        Br = np.array(Br)
        BGflux = np.array(BGflux)
        FluxStrTemp = np.array(FluxStrTemp)
        FluxAnStrTemp = np.array(FluxAnStrTemp)
    
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
        
        return(Br, Br_uninv, Br_inv, Br_CS, Br_unclass, Br_missing)




def Looping(Br, BGflux, FluxStrTemp, FluxAnStrTemp, BGboundary, Opposite):
        OSF_j = []
        N_AS_j = []
        N_SS_j = []
        N_CL_j = []
        N_U_j = []

        Br_i, Br_uninv_i, Br_inv_i, Br_CS_i, Br_unclass_i, Br_missing_i = Body(Br, BGflux, FluxStrTemp, FluxAnStrTemp, BGboundary, Opposite)
        
        #N_MISSING_i = len(Br_missing_i)
                
        N_AS_i = len(Br_uninv_i)
        BR_AS_i = np.nanmean(abs(Br_uninv_i)*(1e-9))
        
        N_SS_i = len(Br_inv_i)
        BR_SS_i = np.nanmean(abs(Br_inv_i)*(1e-9))
        
        N_CL_i = len(Br_CS_i)
        BR_CL_i = np.nanmean(abs(Br_CS_i)*(1e-9))
        
        N_U_i = len(Br_unclass_i)
        BR_U_i = np.nanmean(abs(Br_unclass_i)*(1e-9))
        
        N_all_i = np.nansum([N_AS_i, N_SS_i, N_CL_i, N_U_i])
        
        #N_all_inclMissing_i = np.nansum([N_AS_i, N_SS_i, N_CL_i, N_U_i, N_MISSING_i])
                
        TotalMagFlux = ((4*np.pi*((1*149597871000)**2))/N_all_i)*(np.nansum([N_AS_i*BR_AS_i, N_SS_i*BR_SS_i, N_CL_i*BR_CL_i]))
                                        
        NotConnected = ((4*np.pi*((1*149597871000)**2))/N_all_i)*((N_SS_i*BR_SS_i))
        
        OSF_i = TotalMagFlux - 2*NotConnected
        
        OSF_j.append(OSF_i)
        N_AS_j.append((N_AS_i/N_all_i)*100)
        N_SS_j.append((N_SS_i/N_all_i)*100)
        N_CL_j.append((N_CL_i/N_all_i)*100)
        N_U_j.append((N_U_i/N_all_i)*100)
        #N_missing_j.append((N_MISSING_i/N_all_inclMissing_i)*100)
                                        
        OSF = (np.nanmean(OSF_j))*10**(-14) 
        N_AS = np.nanmean(N_AS_j)
        N_SS = np.nanmean(N_SS_j)
        N_CL = np.nanmean(N_CL_j)
        N_U = np.nanmean(N_U_j) 
        #N_missing = np.nanmean(N_missing_j)       
        
        PercentageAvailable = (N_all_i/18408)*100
        Perc_missing = 100 - PercentageAvailable
        
        return(OSF, N_AS, N_SS, N_CL, N_U, PercentageAvailable, Perc_missing)




def LoopingFinal(Br_best, BGFlux_best, FluxStrTemp_best, FluxAnStrTemp_best, BGboundary, Opposite):
        OSF = []
        N_AS = []
        N_SS = []
        N_CL = []
        N_U = []
        Perc_missing = []
        
        OSF_i, N_AS_i, N_SS_i, N_CL_i, N_U_i, PercentageAvailable_i, Perc_missing_i = Looping(Br_best, BGFlux_best, FluxStrTemp_best, FluxAnStrTemp_best, BGboundary, Opposite)
        
        OSF.append(OSF_i)
        N_AS.append(N_AS_i)
        N_SS.append(N_SS_i)
        N_CL.append(N_CL_i)
        N_U.append(N_U_i)
        Perc_missing.append(Perc_missing_i)
        
        osf = np.nanmean(np.array(OSF))
        n_as = np.nanmean(np.array(N_AS))
        n_ss = np.nanmean(np.array(N_SS))
        n_cl = np.nanmean(np.array(N_CL))
        n_u = np.nanmean(np.array(N_U))
        #perc_missing = np.nanmean(np.array(Perc_missing))
        
        #return(osf, n_as, n_ss, n_cl, n_u, perc_missing, OSF, N_AS, N_SS, N_CL, N_U, Perc_missing, dates, Perc_available_ace, Perc_available_wind)
        return(osf, n_as, n_ss, n_cl, n_u)
        

###############################################################################
"""
Using the same time period as Anderson et al. (2012).
"""
###############################################################################

start_date = datetime.date(1998,1,1)
end_date = datetime.date(2002,12,31)
Br, timestamp, BGflux, FluxStrTemp, FluxAnStrTemp, total_data = BodyACE(start_date, end_date)


opposite = []
bgboundary = []
OSF = []
N_AS = []
N_SS = [] 
N_CL = [] 
N_U = [] 


for i in range(0,41):
    for j in range(0,41):
        print('ij = ', i,j)
        Opposite_temp = np.arange(1.0,3.05,0.05)
        BGboundary_temp = np.arange(1.0,3.05,0.05)
        Opposite = Opposite_temp[i]
        BGboundary = BGboundary_temp[j]
        opposite.append(Opposite)
        bgboundary.append(BGboundary)
        
        OSF_i, N_AS_i, N_SS_i, N_CL_i, N_U_i = LoopingFinal(Br, BGflux, FluxStrTemp, FluxAnStrTemp, BGboundary, Opposite)
        
        OSF.append(OSF_i)
        N_AS.append(N_AS_i)
        N_SS.append(N_SS_i)
        N_CL.append(N_CL_i)
        N_U.append(N_U_i)



###############################################################################
"""
Code used to save variables as text files.
"""
###############################################################################

"""
np.savetxt('TopologyTextFiles/N_AS_ACE_Anderson.txt', N_AS)   
np.savetxt('TopologyTextFiles/N_SS_ACE_Anderson.txt', N_SS)
np.savetxt('TopologyTextFiles/N_CL_ACE_Anderson.txt', N_CL)   
np.savetxt('TopologyTextFiles/N_U_ACE_Anderson.txt', N_U)   
np.savetxt('TopologyTextFiles/OSF_ACE_Anderson.txt', OSF) 
  
np.savetxt('TopologyTextFiles/opposite_ACE_Anderson.txt', opposite)   
np.savetxt('TopologyTextFiles/bgboundary_ACE_Anderson.txt', bgboundary)

"""