"""
Created on Thu Apr 7 2022

Analysis carried out from: Estimating the Open Solar Flux from In Situ Measurements

@author: Anna M Frost
         a.m.frost@pgr.reading.ac.uk
"""


import numpy as np
import xarray as xr
import datetime # python package for dealing with time
from datetime import timedelta
import pandas as pd


def FilterMAGSWE(Year, DOY):
        YEAR = str(Year)
        nc_partial = xr.open_dataset('.\MAGSWE\ACE_MERGED_MAGSWE_'+YEAR+'.nc')
        return(nc_partial)

def MAGSWE(Year, DOY):
        nc_partial = FilterMAGSWE(Year, DOY)
        Bx_GSE = nc_partial.B_gse_x.values
        By_GSE = nc_partial.B_gse_y.values
        Bz_GSE = nc_partial.B_gse_z.values
        year = nc_partial.year.values
        frac_day = nc_partial.fp_day.values
        Bx = - Bx_GSE
        By = - By_GSE
        Bz = - Bz_GSE
        t_datetime_MAGSWE = []
        for i in range(0, year.size):
                t_datetime_MAGSWE.append(datetime.datetime(np.int(year[i]), 1, 1) +datetime.timedelta(frac_day[i] - 1)) 
        invalx = (Bx == 9.99990000e+03)
        Bx[invalx] = np.nan
        invaly = (By == 9.99990000e+03)
        By[invaly] = np.nan
        invalz = (Bz == 9.99990000e+03)
        Bz[invalz] = np.nan
        return(Bx, By, Bz, t_datetime_MAGSWE)


###############################################################################
"""
The following part of code reads in the PA distributions.
"""
###############################################################################

def FilterSWEPAM(Year, DOY):
        YEAR = str(Year)
        temps_partial = np.genfromtxt(".\SWEPAM\SWEPAM_"+YEAR+".txt", delimiter = None)
        return(temps_partial)

def SWEPAM(Year, DOY):
        temps_partial = FilterSWEPAM(Year, DOY)

        YEAR = temps_partial[:,0]
        DOY = temps_partial[:,1]
        HOUR = temps_partial[:,2]
        MIN = temps_partial[:,3]
        SEC = temps_partial[:,4]
        
        frac_doy = DOY + HOUR/24 + MIN/1440 + SEC/86400
        
        t_datetime_SWEPAM = []
        for i in range(0, YEAR.size):
                t_datetime_SWEPAM.append(datetime.datetime(np.int(YEAR[i]), 1, 1) +datetime.timedelta(frac_doy[i] - 1)) 
        PA_bin_1 = temps_partial[:,6]
        PA_bin_2 = temps_partial[:,7]
        PA_bin_3 = temps_partial[:,8]
        PA_bin_4 = temps_partial[:,9]
        PA_bin_5 = temps_partial[:,10]
        PA_bin_6 = temps_partial[:,11]
        PA_bin_7 = temps_partial[:,12]
        PA_bin_8 = temps_partial[:,13]
        PA_bin_9 = temps_partial[:,14]
        PA_bin_10 = temps_partial[:,15]
        PA_bin_11 = temps_partial[:,16]
        PA_bin_12 = temps_partial[:,17]
        PA_bin_13 = temps_partial[:,18]
        PA_bin_14 = temps_partial[:,19]
        PA_bin_15 = temps_partial[:,20]
        PA_bin_16 = temps_partial[:,21]
        PA_bin_17 = temps_partial[:,22]
        PA_bin_18 = temps_partial[:,23]
        PA_bin_19 = temps_partial[:,24]
        PA_bin_20 = temps_partial[:,25]
        return(t_datetime_SWEPAM, PA_bin_1, PA_bin_2, PA_bin_3, PA_bin_4, PA_bin_5, PA_bin_6, PA_bin_7, PA_bin_8, PA_bin_9, PA_bin_10, PA_bin_11, PA_bin_12, PA_bin_13, PA_bin_14, PA_bin_15, PA_bin_16, PA_bin_17, PA_bin_18, PA_bin_19, PA_bin_20)




def CorrectedDataByDay(Year, Month, Day):
        year = int(Year)
        month = int(Month)
        day = int(Day)
        
        DOY = datetime.date(year, month, day).timetuple().tm_yday
        
        Bx, By, Bz, t_datetime_MAGSWE = MAGSWE(year, DOY)
        t_datetime_SWEPAM, PA_bin_1, PA_bin_2, PA_bin_3, PA_bin_4, PA_bin_5,\
        PA_bin_6, PA_bin_7, PA_bin_8, PA_bin_9, PA_bin_10, PA_bin_11,\
        PA_bin_12, PA_bin_13, PA_bin_14, PA_bin_15, PA_bin_16, PA_bin_17,\
        PA_bin_18, PA_bin_19, PA_bin_20 = SWEPAM(year, DOY)
        
        # MAGSWE time correction
        dat_B = np.array([Bx, By, Bz])
        columns_B = ['Bx', 'By', 'Bz']
        t_datetime_MAGSWE = pd.DatetimeIndex(t_datetime_MAGSWE)
        dataframe_B = pd.DataFrame(dat_B.T, index=t_datetime_MAGSWE, columns=columns_B)
        OneSecRes_B = dataframe_B.resample('1S', axis=0).pad()
        Magswe = OneSecRes_B.resample('128S', axis=0, loffset='-64S').mean()
        
        # SWEPAM time correction
        dat_Electron = np.array([PA_bin_1, PA_bin_2, PA_bin_3, PA_bin_4, PA_bin_5,\
                                 PA_bin_6, PA_bin_7, PA_bin_8, PA_bin_9, PA_bin_10, PA_bin_11,\
                                 PA_bin_12, PA_bin_13, PA_bin_14, PA_bin_15, PA_bin_16, PA_bin_17,\
                                 PA_bin_18, PA_bin_19, PA_bin_20])
        columns_Electron = ['Bin_1', 'Bin_2', 'Bin_3', 'Bin_4', 'Bin_5', 'Bin_6', 'Bin_7',\
                            'Bin_8', 'Bin_9', 'Bin_10', 'Bin_11', 'Bin_12', 'Bin_13', 'Bin_14',\
                            'Bin_15', 'Bin_16','Bin_17', 'Bin_18', 'Bin_19', 'Bin_20']
        t_datetime_SWEPAM = pd.DatetimeIndex(t_datetime_SWEPAM)
        dataframe_Electron = pd.DataFrame(dat_Electron.T, index=t_datetime_SWEPAM, columns=columns_Electron)
        OneSecRes_Electron = dataframe_Electron.resample('1S', axis=0).pad()
        Swepam = OneSecRes_Electron.resample('128S', axis=0, loffset='-64S').mean()
        
        idx = pd.date_range(Magswe.index[0], Magswe.index[-1], freq='128S')
        Magswe = Magswe.reindex(idx, fill_value=np.nan)
        Swepam = Swepam.reindex(idx, fill_value=np.nan)
        return(dataframe_B, dataframe_Electron, Magswe, Swepam)



start_date = datetime.date(1998,1,1)
#end_date = datetime.date(1999,12,31)

end_date = datetime.date(2017,12,31)


start_date_0 = start_date

delta = timedelta(days=366)

while start_date <= end_date:
    a = start_date.strftime("%Y%m%d")
    Year = int(a[0:4])
    Month = int(a[4:6])
    Day = int(a[6:8])
    date_target = [Year, Month, Day]
    Raw_Magswe, Raw_Swepam, Magswe, Swepam = CorrectedDataByDay(Year, Month, Day)
    if start_date==start_date_0:
        ds = Magswe.to_xarray()
        ds.to_netcdf('ACE Magnetic Corrected Time.nc')
        ds.close()
        dsE = Swepam.to_xarray()
        dsE.to_netcdf('ACE Electron Corrected Time.nc')
        dsE.close()
        
        ds1 = Raw_Magswe.to_xarray()
        ds1.to_netcdf('ACE Magnetic Raw Time.nc')
        ds1.close()
        dsE1 = Raw_Swepam.to_xarray()
        dsE1.to_netcdf('ACE Electron Raw Time.nc')
        dsE1.close()
    else:
        ds = xr.open_dataset('ACE Magnetic Corrected Time.nc')
        df = ds.to_dataframe()
        ds.close()
        df = df.append(Magswe)
        dset = df.to_xarray() 
        dset.to_netcdf('ACE Magnetic Corrected Time.nc')
        dset.close()
        
        dsE = xr.open_dataset('ACE Electron Corrected Time.nc')
        dfE = dsE.to_dataframe()
        dsE.close()
        dfE = dfE.append(Swepam)
        dsetE = dfE.to_xarray() 
        dsetE.to_netcdf('ACE Electron Corrected Time.nc')
        dsetE.close()
        
        
        ds1 = xr.open_dataset('ACE Magnetic Raw Time.nc')
        df1 = ds1.to_dataframe()
        ds1.close()
        df1 = df1.append(Raw_Magswe)
        dset1 = df1.to_xarray() 
        dset1.to_netcdf('ACE Magnetic Raw Time.nc')
        dset1.close()
        
        dsE1 = xr.open_dataset('ACE Electron Raw Time.nc')
        dfE1 = dsE1.to_dataframe()
        dsE1.close()
        dfE1 = dfE1.append(Raw_Swepam)
        dsetE1 = dfE1.to_xarray() 
        dsetE1.to_netcdf('ACE Electron Raw Time.nc')
        dsetE1.close()
                    
    start_date += delta