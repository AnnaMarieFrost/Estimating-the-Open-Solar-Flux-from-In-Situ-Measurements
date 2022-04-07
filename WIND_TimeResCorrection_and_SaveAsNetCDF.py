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
from astropy import units as u
import pandas as pd
import os
import cdflib
import cdflib.epochs as epochs #this is the part of cdflib that is useful for converting out of 'cdf time' into something sensible
import astropy.time as as_time
import os.path
from pandas.tseries.frequencies import to_offset



    
###############################################################################
"""
This function takes an input date ([YYYY, MM, DD]), the analyser selection 
(high ('h') or low ('l') energy) from the EESA detectors and the filename for
the cdf files as given from DL_3DP. The output is the target subfolder to save 
the data for the different energies in (cdf_file) and a dictionary giving the 
variables contained in the cdf file (inf). 
"""
###############################################################################

def get_cdf_Inf(date, eesa, fname):
        if eesa=='h':
                cdf_file = cdflib.CDF(fname)    #targets the subfolder "EESA_H_CDFs" to save the eesa='h' data in
                #inf = cdf_file.cdf_info()    #makes a dictionary which tells you all of the variables contained in the CDF file
        elif eesa=='l':
                cdf_file = cdflib.CDF(fname)    #targets the subfolder "EESA_L_CDFs" to save the eesa='l' data in
                #inf = cdf_file.cdf_info()    #makes a dictionary which tells you all of the variables contained in the CDF file
        return(cdf_file)



###############################################################################
"""
This function takes the input of cdf_file (which was the output from get_cdf_Inf)
and outputs a split cdf file, separated into the variables used in the 
remainder of the program.
"""
###############################################################################

def split_CDF(cdf_file):
        E = np.flip(cdf_file.varget('ENERGY'),1)     # gets the bin energies for each time step. 15 energy bins, E gives the 
                                                        # central energy for each (spaced logarithmically).
                                                     # np.flip(something, 1) flips array horizontally (axis = 1)
        PA = cdf_file.varget('PANGLE')    # there are 8 pitch angle bins, so there are 15 x 8 bins for each time step
        vsw_3dp = cdf_file.varget('VSW')    # this is the background solar wind bulk velocity (protons + ions)
        B = cdf_file.varget('MAGF')    # this is 3 components of the magnetic field
        flux = np.flip(cdf_file.varget('FLUX'),2)    # this gets the differential number flux
        flux[flux<0] = np.nan # fill values of the flux are set as <0, but it works better in python if we make them a nan
        return(E, PA, vsw_3dp, B, flux)




###############################################################################
"""
This function takes the input of cdf_file (which was the output from get_cdf_Inf)
and outputs the time variables contained in the cdf file for use later.
"""
###############################################################################

def time_Julian(cdf_file):
        t_cdf = cdf_file.varget('Epoch')    #cdf time, Epoch, defined as the number of milliseconds since 01-Jan-0000 00:00:00.000
        t_unix = cdf_file.varget('TIME')    # unix time is seconds since jan 1st 1970
        t_dates = epochs.CDFepoch.breakdown(t_cdf,to_np=True)    # cdflib gives a way to convert from cdf time
        if np.size(t_dates)==7:
                t_datetime = np.array([datetime.datetime(t_dates[0], t_dates[1], t_dates[2], t_dates[3], t_dates[4], t_dates[5], t_dates[6])])   # gives a list of datetime objects to give the time
                t_jd = as_time.Time(t_datetime).jd
        else:
                t_datetime = np.array([datetime.datetime(*x) for x in t_dates])   # gives a list of datetime objects to give the time
                t_jd = as_time.Time(t_datetime).jd    # can use the astropy time to convert to julian date
        return(t_unix, t_datetime, t_jd, t_dates)



def CorrectedDataByDay(Year, Month, Day):
        eesa = 'l'
        
        # Year/monthday of the target date
        year = int(Year)
        month = int(Month)
        day = int(Day)
        
        # Datetime for target date and one day before and after
        target_day = datetime.date(year, month, day)
        Day_before = target_day - timedelta(days=1)
        Day_after = target_day + timedelta(days=1)
        
        # Target day
        a_target = target_day.strftime("%Y%m%d")
        year_target = a_target[0:4]
        month_target = a_target[4:6]
        day_target = a_target[6:8]
        date_target = [year_target, month_target, day_target]
        
        # Day before target
        a_before = Day_before.strftime("%Y%m%d")
        year_before = a_before[0:4]
        month_before = a_before[4:6]
        day_before = a_before[6:8]
        date_before = [year_before, month_before, day_before]
        
        # Day after target
        a_after = Day_after.strftime("%Y%m%d")
        year_after = a_after[0:4]
        month_after = a_after[4:6]
        day_after = a_after[6:8]
        date_after = [year_after, month_after, day_after]
        
        
        sflux1_temp_before = np.array([])
        sflux1_temp_target = np.array([])
        sflux1_temp_after = np.array([])
        sflux2_temp_before = np.array([])
        sflux2_temp_target = np.array([])
        sflux2_temp_after = np.array([])
        sflux3_temp_before = np.array([])
        sflux3_temp_target = np.array([])
        sflux3_temp_after = np.array([])
        sflux4_temp_before = np.array([])
        sflux4_temp_target = np.array([])
        sflux4_temp_after = np.array([])
        sflux5_temp_before = np.array([])
        sflux5_temp_target = np.array([])
        sflux5_temp_after = np.array([])
        sflux6_temp_before = np.array([])
        sflux6_temp_target = np.array([])
        sflux6_temp_after = np.array([])
        sflux7_temp_before = np.array([])
        sflux7_temp_target = np.array([])
        sflux7_temp_after = np.array([])
        sflux8_temp_before = np.array([])
        sflux8_temp_target = np.array([])
        sflux8_temp_after = np.array([])
        
        Bx_temp_before = np.array([])
        Bx_temp_target = np.array([])
        Bx_temp_after = np.array([])
        
        By_temp_before = np.array([])
        By_temp_target = np.array([])
        By_temp_after = np.array([])
        
        Bz_temp_before = np.array([])
        Bz_temp_target = np.array([])
        Bz_temp_after = np.array([])
        
        t_datetime_before = np.array([])
        t_datetime_target = np.array([])
        t_datetime_after = np.array([])
        Val_before = np.array([], dtype=bool)
        Val_target = np.array([], dtype=bool) 
        Val_after = np.array([], dtype=bool)
        
        
        #pdb.set_trace()
        
        fname_target = 'EESA_L_CDFs\wi_elpd_3dp_'+str(date_target[0])+str(date_target[1])+str(date_target[2])+'.cdf' 
        if os.path.isfile(fname_target):
                cdf_file_target = get_cdf_Inf(date_target, eesa, fname_target)
                E_target, PA_target, vsw_3dp_target, B_target, flux_target = split_CDF(cdf_file_target)
                
                t_unix_target, t_datetime_target, t_jd_target, t_dates_target = time_Julian(cdf_file_target)
                
                ds_target = xr.Dataset({'flux_h':(['t','y','x'],flux_target), 'vion_3dp':(['t','z'],vsw_3dp_target),\
                                 'B_3dp':(['t','z'],B_target),'t_unix':('t',t_unix_target),\
                                 't_datetime':('t',t_datetime_target)}, coords =  {'t_jd':('t',t_jd_target),\
                                 'E_h':(['t','x'],E_target), 'PA_h':(['t','y'],PA_target)} ) 
                # This uses x_array to make a dataset out of all our variables. For each variable I'm giving it
                # a name in my dataset, ds, and saying what dimensions it varies along. I also define coordiantes
                # as time, energy, and pitch angle, since I have a unique flux measurement for each of those.
                                                                                   
                #This is kind of optional but it provides a neat way to keep the data together 
                ds_target.flux_h.attrs['units'] = '{0}'.format(u.cm**-2*u.sr**-1*u.eV**-1*u.s**-1) # xarray lets us assign arbitrary attributes that get stored for each parameter. here I'm using astropy units but this can be a string
                ds_target.E_h.attrs['units']='{0}'.format(u.eV)
                ds_target.PA_h.attrs['units'] = '{0}'.format(u.deg)
                ds_target.vion_3dp.attrs['units'] ='{0}'.format(u.km/u.s)
                ds_target.B_3dp.attrs['units']='{0}'.format(u.nT)
                ds_target.t_unix.attrs['units']='{0}'.format(u.s)
                ds_target.t_jd.attrs['units']='{0}'.format(u.day)
                ds_target.vion_3dp.attrs['coords'] = 'GSE' # magnetic field is in GSE coordinates
                ds_target.B_3dp.attrs['coords'] = 'GSE'
                
                sort_dat_target = np.sort(t_datetime_target)
                time_diff_target = np.diff(sort_dat_target)
                val_target = (time_diff_target > timedelta(seconds=0))
                valid_target = np.insert(val_target, 0, True, axis=0)
                #BowShockBool_target = BowShockFiltering(date_target, t_datetime_target)
                Val_target = valid_target #& BowShockBool_target
                
                Bx_temp_target = -B_target[:, 0]
                By_temp_target = -B_target[:, 1]
                Bz_temp_target = -B_target[:, 2]
                sflux1_temp_target = flux_target[:,0,11]  # 2 = Energy of 292eV
                sflux2_temp_target = flux_target[:,1,11]
                sflux3_temp_target = flux_target[:,2,11]
                sflux4_temp_target = flux_target[:,3,11]
                sflux5_temp_target = flux_target[:,4,11]
                sflux6_temp_target = flux_target[:,5,11]
                sflux7_temp_target = flux_target[:,6,11]
                sflux8_temp_target = flux_target[:,7,11]
                
        
        # Calling for data from day before
        fname_before = 'EESA_L_CDFs\wi_elpd_3dp_'+str(date_before[0])+str(date_before[1])+str(date_before[2])+'.cdf' 
        if os.path.isfile(fname_before):
                cdf_file_before = get_cdf_Inf(date_before, eesa, fname_before)
                E_before, PA_before, vsw_3dp_before, B_before, flux_before = split_CDF(cdf_file_before)
                t_unix_before, t_datetime_before, t_jd_before, t_dates_before = time_Julian(cdf_file_before)
                
                ds_before = xr.Dataset({'flux_h':(['t','y','x'],flux_before), 'vion_3dp':(['t','z'],vsw_3dp_before),\
                                 'B_3dp':(['t','z'],B_before),'t_unix':('t',t_unix_before),\
                                 't_datetime':('t',t_datetime_before)}, coords =  {'t_jd':('t',t_jd_before),\
                                 'E_h':(['t','x'],E_before), 'PA_h':(['t','y'],PA_before)} ) 
                # This uses x_array to make a dataset out of all our variables. For each variable I'm giving it
                # a name in my dataset, ds, and saying what dimensions it varies along. I also define coordiantes
                # as time, energy, and pitch angle, since I have a unique flux measurement for each of those.
                                                                                   
                #This is kind of optional but it provides a neat way to keep the data together 
                ds_before.flux_h.attrs['units'] = '{0}'.format(u.cm**-2*u.sr**-1*u.eV**-1*u.s**-1) # xarray lets us assign arbitrary attributes that get stored for each parameter. here I'm using astropy units but this can be a string
                ds_before.E_h.attrs['units']='{0}'.format(u.eV)
                ds_before.PA_h.attrs['units'] = '{0}'.format(u.deg)
                ds_before.vion_3dp.attrs['units'] ='{0}'.format(u.km/u.s)
                ds_before.B_3dp.attrs['units']='{0}'.format(u.nT)
                ds_before.t_unix.attrs['units']='{0}'.format(u.s)
                ds_before.t_jd.attrs['units']='{0}'.format(u.day)
                ds_before.vion_3dp.attrs['coords'] = 'GSE' # magnetic field is in GSE coordinates
                ds_before.B_3dp.attrs['coords'] = 'GSE'
                
                sort_dat_before = np.sort(t_datetime_before)
                time_diff_before = np.diff(sort_dat_before)
                val_before = (time_diff_before > timedelta(seconds=0))
                valid_before = np.insert(val_before, 0, True, axis=0)
                #BowShockBool_before = BowShockFiltering(date_before, t_datetime_before)
                Val_before = valid_before #& BowShockBool_before
        
                Bx_temp_before = -B_before[:, 0]
                By_temp_before = -B_before[:, 1]
                Bz_temp_before = -B_before[:, 2]
                sflux1_temp_before = flux_before[:,0,11]  # 2 = Energy of 292eV
                sflux2_temp_before = flux_before[:,1,11]
                sflux3_temp_before = flux_before[:,2,11]
                sflux4_temp_before = flux_before[:,3,11]
                sflux5_temp_before = flux_before[:,4,11]
                sflux6_temp_before = flux_before[:,5,11]
                sflux7_temp_before = flux_before[:,6,11]
                sflux8_temp_before = flux_before[:,7,11]
        
        
        # Calling for data from day after
        fname_after = 'EESA_L_CDFs\wi_elpd_3dp_'+str(date_after[0])+str(date_after[1])+str(date_after[2])+'.cdf' 
        if os.path.isfile(fname_after):
                cdf_file_after = get_cdf_Inf(date_after, eesa, fname_after)
                E_after, PA_after, vsw_3dp_after, B_after, flux_after = split_CDF(cdf_file_after)
                t_unix_after, t_datetime_after, t_jd_after, t_dates_after = time_Julian(cdf_file_after)
                
                ds_after = xr.Dataset({'flux_h':(['t','y','x'],flux_after), 'vion_3dp':(['t','z'],vsw_3dp_after),\
                                 'B_3dp':(['t','z'],B_after),'t_unix':('t',t_unix_after),\
                                 't_datetime':('t',t_datetime_after)}, coords =  {'t_jd':('t',t_jd_after),\
                                 'E_h':(['t','x'],E_after), 'PA_h':(['t','y'],PA_after)} ) 
                # This uses x_array to make a dataset out of all our variables. For each variable I'm giving it
                # a name in my dataset, ds, and saying what dimensions it varies along. I also define coordiantes
                # as time, energy, and pitch angle, since I have a unique flux measurement for each of those.
                                                                                   
                #This is kind of optional but it provides a neat way to keep the data together 
                ds_after.flux_h.attrs['units'] = '{0}'.format(u.cm**-2*u.sr**-1*u.eV**-1*u.s**-1) # xarray lets us assign arbitrary attributes that get stored for each parameter. here I'm using astropy units but this can be a string
                ds_after.E_h.attrs['units']='{0}'.format(u.eV)
                ds_after.PA_h.attrs['units'] = '{0}'.format(u.deg)
                ds_after.vion_3dp.attrs['units'] ='{0}'.format(u.km/u.s)
                ds_after.B_3dp.attrs['units']='{0}'.format(u.nT)
                ds_after.t_unix.attrs['units']='{0}'.format(u.s)
                ds_after.t_jd.attrs['units']='{0}'.format(u.day)
                ds_after.vion_3dp.attrs['coords'] = 'GSE' # magnetic field is in GSE coordinates
                ds_after.B_3dp.attrs['coords'] = 'GSE'
                
                sort_dat_after = np.sort(t_datetime_after)
                time_diff_after = np.diff(sort_dat_after)
                val_after = (time_diff_after > timedelta(seconds=0))
                valid_after = np.insert(val_after, 0, True, axis=0)
                
                Val_after = valid_after
                
                Bx_temp_after =  -B_after[:, 0]
                By_temp_after =  -B_after[:, 1]
                Bz_temp_after =  -B_after[:, 2]
                sflux1_temp_after = flux_after[:,0,11]  # 2 = Energy of 292eV
                sflux2_temp_after = flux_after[:,1,11]
                sflux3_temp_after = flux_after[:,2,11]
                sflux4_temp_after = flux_after[:,3,11]
                sflux5_temp_after = flux_after[:,4,11]
                sflux6_temp_after = flux_after[:,5,11]
                sflux7_temp_after = flux_after[:,6,11]
                sflux8_temp_after = flux_after[:,7,11]
        
        
        flux_sf1 = np.concatenate([sflux1_temp_before, sflux1_temp_target, sflux1_temp_after])
        flux_sf2 = np.concatenate([sflux2_temp_before, sflux2_temp_target, sflux2_temp_after])
        flux_sf3 = np.concatenate([sflux3_temp_before, sflux3_temp_target, sflux3_temp_after])
        flux_sf4 = np.concatenate([sflux4_temp_before, sflux4_temp_target, sflux4_temp_after])
        flux_sf5 = np.concatenate([sflux5_temp_before, sflux5_temp_target, sflux5_temp_after])
        flux_sf6 = np.concatenate([sflux6_temp_before, sflux6_temp_target, sflux6_temp_after])
        flux_sf7 = np.concatenate([sflux7_temp_before, sflux7_temp_target, sflux7_temp_after])
        flux_sf8 = np.concatenate([sflux8_temp_before, sflux8_temp_target, sflux8_temp_after])
        Bx_all = np.concatenate([Bx_temp_before, Bx_temp_target, Bx_temp_after])
        By_all = np.concatenate([By_temp_before, By_temp_target, By_temp_after])
        Bz_all = np.concatenate([Bz_temp_before, Bz_temp_target, Bz_temp_after])
        t_datetime_all = np.concatenate([t_datetime_before, t_datetime_target, t_datetime_after])
        
        sort_dates = np.sort(t_datetime_all)
        time_difference = np.diff(sort_dates)
        valid_time_diff = (time_difference > timedelta(seconds=0))
        valid_bool = np.insert(valid_time_diff, 0, True, axis=0)
        
        data = np.array([Bx_all, By_all, Bz_all, flux_sf1, flux_sf2, flux_sf3, flux_sf4, flux_sf5, flux_sf6, flux_sf7, flux_sf8])
        columns = ['Bx', 'By', 'Bz', 'Bin_1', 'Bin_2', 'Bin_3', 'Bin_4', 'Bin_5', 'Bin_6', 'Bin_7', 'Bin_8']
        dataframe = pd.DataFrame(data.T, index=t_datetime_all, columns=columns)
        data_sort = dataframe.sort_index()
        
        Valid = data_sort[valid_bool]
        OneSecRes = Valid.resample('0.5S', axis=0).pad()      
        
        CorrRes = OneSecRes.resample('128S', axis=0).mean()
        CorrRes.index = CorrRes.index + to_offset('-64S')
        
        idx = pd.date_range(CorrRes.index[0], CorrRes.index[-1], freq='128S')
        CorrRes = CorrRes.reindex(idx, fill_value=np.nan)
                

        # Returning data from target date only
        Val_dates = (CorrRes.index.date >= pd.to_datetime(target_day)) & (CorrRes.index.date <= pd.to_datetime(target_day))
        WIND = CorrRes.loc[Val_dates]
        
        Val_datesraw = (dataframe.index.date >= pd.to_datetime(target_day)) & (dataframe.index.date <= pd.to_datetime(target_day))
        Rawdata = dataframe.loc[Val_datesraw]
        
        return(WIND, Rawdata)
        



start_date = datetime.date(1994,12,20)
start_date_0 = datetime.date(1994,12,20)

delta = timedelta(days=1)

end_date = datetime.date(2022,2,16)

while start_date <= end_date:
    a = start_date.strftime("%Y%m%d")
    Year = (a[0:4])
    Month = (a[4:6])
    Day = (a[6:8])
    date_target = [Year, Month, Day]
    fname_target = 'C:\PhD Coding\WIND\EESA_L_CDFs\wi_elpd_3dp_'+str(date_target[0])+str(date_target[1])+str(date_target[2])+'_v02.cdf' 
    print(start_date)
    if os.path.isfile(fname_target):
            WIND, Rawdata = CorrectedDataByDay(Year, Month, Day)
            
            if start_date==start_date_0:
                    ds = Rawdata.to_xarray()
                    ds.to_netcdf('C:\PhD Coding\WIND\Data NetCDFs\WIND-L Raw Time.nc')
                    ds.close()
                    ds1 = WIND.to_xarray()
                    ds1.to_netcdf('C:\PhD Coding\WIND\Data NetCDFs\WIND-L Corrected Time.nc')
                    ds1.close()
            else:
                    
                    ds = xr.open_dataset('C:\PhD Coding\WIND\Data NetCDFs\WIND-L Raw Time.nc')
                    df = ds.to_dataframe()
                    ds.close()
                    df = df.append(Rawdata)
                    dset = df.to_xarray()
                    dset.to_netcdf('C:\PhD Coding\WIND\Data NetCDFs\WIND-L Raw Time.nc')
                    dset.close()
                    
                    ds1 = xr.open_dataset('C:\PhD Coding\WIND\Data NetCDFs\WIND-L Corrected Time.nc')
                    df1 = ds1.to_dataframe()
                    ds1.close()
                    df1 = df1.append(WIND)
                    dset1 = df1.to_xarray()
                    dset1.to_netcdf('C:\PhD Coding\WIND\Data NetCDFs\WIND-L Corrected Time.nc')
                    dset1.close()
                   
            start_date += delta
    else:
            start_date += delta
    