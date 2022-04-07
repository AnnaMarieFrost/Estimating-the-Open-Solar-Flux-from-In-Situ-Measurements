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



def FilterOSFErrors(start_date, end_date, WindDateString):
        start_date = pd.to_datetime(start_date)
        end_date = pd.to_datetime(end_date)
        nc = xr.open_dataset('ACE OSF with Wind intervals\ACE OSF Error with Wind Intervals Applied'+WindDateString+'.nc')
        OSF_Error_temp = nc.to_dataframe()
        OSF_Error = OSF_Error_temp.iloc[(OSF_Error_temp.index >= start_date) & (OSF_Error_temp.index < end_date)]
        return(OSF_Error)


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
#OSF_perc_error = []
Plot_Perc_Error_temp = []
Plot_5Percentile = []
Plot_95Percentile = []
Plot_5Percentile_temp = []
Plot_95Percentile_temp = []

Plot_25Percentile = []
Plot_75Percentile = []
Plot_25Percentile_temp = []
Plot_75Percentile_temp = []

Plot_33Percentile = []
Plot_66Percentile = []
Plot_33Percentile_temp = []
Plot_66Percentile_temp = []

min_temp = []
max_temp = []

dates = []



start_date = datetime.datetime(1998,1,22)
end_date = datetime.datetime(2017,12,31)

#start_date = datetime.datetime(1998,1,22)
#end_date = datetime.datetime(2017,12,31)

delta = timedelta(days=27)


while start_date <= end_date:
        end_date_temp = start_date + delta
        print('WIND dates: ', start_date)
        
        PercentageWINDAvailable = WINDPercentage(start_date, end_date_temp)
        
        
        if PercentageWINDAvailable<=100:
        #if 88.295<=PercentageAvailable<=88.3:
                OSF_perc_error = []
        
                WindDateString = start_date.strftime("%Y%m%d")
                
                WindPercentageAvailable.append(float(PercentageWINDAvailable))
                
                nc = xr.open_dataset('ACE OSF with WIND intervals\ACE OSF Error with Wind Intervals Applied'+WindDateString+'.nc')
                OSF_perc_error_temp = nc.to_dataframe()
                
                #pdb.set_trace()
                OSF_perc_error = abs(np.array(OSF_perc_error_temp))
                #OSF_perc_error = np.array(OSF_perc_error_temp)
                temp = np.nanpercentile(OSF_perc_error, 50)
                error_5percentile = np.nanpercentile(OSF_perc_error, 5)
                error_95percentile = np.nanpercentile(OSF_perc_error, 95)
                
                error_25percentile = np.nanpercentile(OSF_perc_error, 25)
                error_75percentile = np.nanpercentile(OSF_perc_error, 75)
                
                error_33percentile = np.nanpercentile(OSF_perc_error, 33)
                error_66percentile = np.nanpercentile(OSF_perc_error, 66)
                
                minimum = np.nanpercentile(OSF_perc_error, 0)
                maximum = np.nanpercentile(OSF_perc_error, 100)
                
                Plot_5Percentile_temp.append(error_5percentile)
                Plot_95Percentile_temp.append(error_95percentile)
                Plot_25Percentile_temp.append(error_25percentile)
                Plot_75Percentile_temp.append(error_75percentile)
                Plot_33Percentile_temp.append(error_33percentile)
                Plot_66Percentile_temp.append(error_66percentile)
                min_temp.append(minimum)
                max_temp.append(maximum)
                
                Plot_Perc_Error_temp.append(temp)                
                Plot_5Percentile = abs(np.array(Plot_5Percentile_temp) - np.array(Plot_Perc_Error_temp))
                Plot_95Percentile = abs(np.array(Plot_Perc_Error_temp) - np.array(Plot_95Percentile_temp))
                Plot_25Percentile = abs(np.array(Plot_25Percentile_temp) - np.array(Plot_Perc_Error_temp))
                Plot_75Percentile = abs(np.array(Plot_Perc_Error_temp) - np.array(Plot_75Percentile_temp))
                Plot_33Percentile = abs(np.array(Plot_33Percentile_temp) - np.array(Plot_Perc_Error_temp))
                Plot_66Percentile = abs(np.array(Plot_Perc_Error_temp) - np.array(Plot_66Percentile_temp))
                Plot_Minimum = abs(np.array(min_temp) - np.array(Plot_Perc_Error_temp))
                Plot_Maximum = abs(np.array(Plot_Perc_Error_temp) - np.array(max_temp))
                                
                dates.append(start_date)
                
                print('Median % Error:', temp)                                    
                start_date += delta
        else:
                dates.append(start_date)
                start_date += delta



plot_5percentile_temp = np.array(Plot_5Percentile_temp)[np.invert(np.array(Plot_Perc_Error_temp)==0) & np.invert(np.array(WindPercentageAvailable)==0)]
plot_95percentile_temp = np.array(Plot_95Percentile_temp)[np.invert(np.array(Plot_Perc_Error_temp)==0) & np.invert(np.array(WindPercentageAvailable)==0)]
plot_25percentile_temp = np.array(Plot_25Percentile_temp)[np.invert(np.array(Plot_Perc_Error_temp)==0) & np.invert(np.array(WindPercentageAvailable)==0)]
plot_75percentile_temp = np.array(Plot_75Percentile_temp)[np.invert(np.array(Plot_Perc_Error_temp)==0) & np.invert(np.array(WindPercentageAvailable)==0)]
plot_33percentile_temp = np.array(Plot_33Percentile_temp)[np.invert(np.array(Plot_Perc_Error_temp)==0) & np.invert(np.array(WindPercentageAvailable)==0)]
plot_66percentile_temp = np.array(Plot_66Percentile_temp)[np.invert(np.array(Plot_Perc_Error_temp)==0) & np.invert(np.array(WindPercentageAvailable)==0)]

Plot_Perc_Error = np.array(Plot_Perc_Error_temp)[np.invert(np.array(Plot_Perc_Error_temp)==0) & np.invert(np.array(WindPercentageAvailable)==0)]
Plot_5Percentile = np.array(Plot_5Percentile)[np.invert(np.array(Plot_Perc_Error_temp)==0) & np.invert(np.array(WindPercentageAvailable)==0)]
Plot_95Percentile = np.array(Plot_95Percentile)[np.invert(np.array(Plot_Perc_Error_temp)==0) & np.invert(np.array(WindPercentageAvailable)==0)]
Plot_25Percentile = np.array(Plot_25Percentile)[np.invert(np.array(Plot_Perc_Error_temp)==0) & np.invert(np.array(WindPercentageAvailable)==0)]
Plot_75Percentile = np.array(Plot_75Percentile)[np.invert(np.array(Plot_Perc_Error_temp)==0) & np.invert(np.array(WindPercentageAvailable)==0)]
Plot_33Percentile = np.array(Plot_33Percentile)[np.invert(np.array(Plot_Perc_Error_temp)==0) & np.invert(np.array(WindPercentageAvailable)==0)]
Plot_66Percentile = np.array(Plot_66Percentile)[np.invert(np.array(Plot_Perc_Error_temp)==0) & np.invert(np.array(WindPercentageAvailable)==0)]
Plot_Minimum = np.array(Plot_Minimum)[np.invert(np.array(Plot_Perc_Error_temp)==0) & np.invert(np.array(WindPercentageAvailable)==0)]
Plot_Maximum = np.array(Plot_Maximum)[np.invert(np.array(Plot_Perc_Error_temp)==0) & np.invert(np.array(WindPercentageAvailable)==0)]
WindPercentageAvailable = np.array(WindPercentageAvailable)[np.invert(np.array(Plot_Perc_Error_temp)==0) & np.invert(np.array(WindPercentageAvailable)==0)]





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
Plotting the true and predicted errors.
"""
###############################################################################

fig = plt.figure()
ax = fig.add_subplot(111)

ax.set_xlabel('WIND Data Available (%)', fontsize=40)
ax.set_ylabel('Error in OSF (%)', fontsize=40)
#plt.title('Predicted and True OSF Error', fontsize=35)

ax.tick_params(axis='both', which='major', labelsize=30)

ax.errorbar(WindPercentageAvailable, Plot_Perc_Error, yerr=[Plot_5Percentile,Plot_95Percentile], fmt='o', color='b', capsize=4, label='5th and 95th percentile')
ax.errorbar(WindPercentageAvailable, Plot_Perc_Error, yerr=[Plot_25Percentile,Plot_75Percentile], fmt='o', color='r', capsize=4, label='25th and 75th percentile')
ax.errorbar(WindPercentageAvailable, Plot_Perc_Error, yerr=[Plot_33Percentile,Plot_66Percentile], fmt='o', color='k', capsize=4, label='33rd and 66th percentile')


ax.scatter(True_WINDPercentageAvailable, True_Plot_Perc_Error, marker='x', color='k', label='Actual Error', s=300)


ax.set_xlim([0, 103])
ax.set_ylim([-1, 95])

ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')

ax.legend(fontsize=30)

plt.show() 


