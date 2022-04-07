###############################################################################
"""
This code reads in the cdfs of the orbit data and returns the array 'Orbittime'
which can be saved as a txt file by running 
"np.savetxt('BowShockv1.txt', Orbittime)" in the console.
""" 
###############################################################################


import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import datetime # python package for dealing with time
from datetime import timedelta
from astropy import units as u
import pandas as pd
import requests
import os
import cdflib
import cdflib.epochs as epochs #this is the part of cdflib that is useful for converting out of 'cdf time' into something sensible
import astropy.time as as_time
import astropy
import pdb #pdb.set_trace()
import sys
from importlib import reload
reload(sys)


###############################################################################
"""
This function takes an input date ([YYYY, MM, DD]) and downloads the data files 
containing the oribt data for the WIND instrument. It outputs the status 
(whether there is a file at the URL (200) or not (404)) as well as the filename 
of the file so it can be found.
"""
###############################################################################

def DL_Orbit(date):
    s = requests.Session()
    windurl = 'https://cdaweb.gsfc.nasa.gov/pub/data/wind/orbit/pre_or/'
    for i in range(1,12):
        fname = 'wi_or_pre_'+str(date[0])+str(date[1])+str(date[2])+'_v0'+str(i)+'.cdf'
        direct = 'C:/PhD Coding/WIND/Orbit/ORBIT_CDFs/'
        url = windurl+date[0]+'/'+fname #put the whole thing together to find the individual data files for each specified date
        fnamesav = direct+fname
        
        if os.path.isfile(''+fnamesav):
            print('File already downloaded!')
            status = 200
        else:
            r = s.get(url) 
            if r.status_code == 200:
                with open(fnamesav, 'wb') as out:
                    for bits in r.iter_content(): 
                        out.write(bits)
                print('Saved File: ', date)
                out.close()
            if r.status_code == 404:
                if os.path.isfile(fnamesav):	os.remove(fnamesav)
            status = r.status_code
    return(status, fname)


###############################################################################
"""
Below is the code used to read in the orbit data and filter according to the 
radial distance from Earth.
"""
###############################################################################

def OrbitFile(date):
    for i in range(1,12):
        fname = 'wi_or_pre_'+str(date[0])+str(date[1])+str(date[2])+'_v0'+str(i)+'.cdf'
        fnamesav = 'C:/PhD Coding/WIND/Orbit/ORBIT_CDFs/' + fname
        if os.path.isfile(''+fnamesav):
            cdf_file = cdflib.CDF(fnamesav)
    return(cdf_file)


def OrbitFilter(date):
        cdf_file_orbit = OrbitFile(date)
        Pos_gse = cdf_file_orbit.varget('GSE_POS')
        time_ms = cdf_file_orbit.varget('Time_PB5')    # [Year, DOY, Milliseconds]
        Orbittime = []
        for i in range(0, time_ms[:,0].size):
                Orbit_temp = (datetime.datetime(np.int(time_ms[i,0]), 1, 1)) + datetime.timedelta(days = np.int(time_ms[i,1] - 1), milliseconds = np.int(time_ms[i,2]))
                Orbittime.append(Orbit_temp)
        Pos_X = Pos_gse[:,0]
        Pos_Y = Pos_gse[:,1]
        Pos_Z = Pos_gse[:,2]
        Orbitradius = np.sqrt((Pos_X**2)+(Pos_Y**2)+(Pos_Z**2))
        Radius = np.sqrt((Pos_Y**2)+(Pos_Z**2))
        Bow_shock = ((Orbitradius > 191130) & (Pos_X > 0)) | ((Pos_X < 0) & (Radius > 318550))
        BS = np.invert(Bow_shock)           
        """                    
        Orbittime = np.array(Orbittime)[BS]
        Pos_X = Pos_X[BS]
        Pos_Y = Pos_Y[BS]
        Pos_Z = Pos_Z[BS]
        """
        return(Orbittime, Pos_X, Pos_Y, Pos_Z)



###############################################################################
"""
Below is the main code being called and looping over multiple days.
"""
###############################################################################

start_date = datetime.date(2000,1,1)
end_date = datetime.date(2021,11,1)


start_date_0 = start_date
delta = timedelta(days=1)

Orbittime = []
Pos_X = [] 
Pos_Y = []
Pos_Z = []

while start_date <= end_date:
    a = start_date.strftime("%Y%m%d")
    year = a[0:4]
    month = a[4:6]
    day = a[6:8]
    date = [year, month, day]
    i = 0
    for i in range(1,12):
        fname_target = 'C:/PhD Coding/WIND/Orbit/ORBIT_CDFs/wi_or_pre_'+str(date[0])+str(date[1])+str(date[2])+'_v0'+str(i)+'.cdf'
        i += 1
        if os.path.isfile(fname_target):
            Orbittime_i, Pos_X_i, Pos_Y_i, Pos_Z_i = OrbitFilter(date) 
            Orbittime = np.hstack([Orbittime, Orbittime_i])
            Pos_X = np.hstack([Pos_X, Pos_X_i])
            Pos_Y = np.hstack([Pos_Y, Pos_Y_i])
            Pos_Z = np.hstack([Pos_Z, Pos_Z_i])
            start_date += delta
        else:
            start_date += delta
        

data = np.array([Pos_X, Pos_Y, Pos_Z])
columns = ['Pos_X', 'Pos_Y', 'Pos_Z']
dataframe = pd.DataFrame(data.T, index=Orbittime, columns=columns)


ds = dataframe.to_xarray()
ds.to_netcdf('RadialDist.nc')
ds.close()
"""
ds = xr.open_dataset('RadialDist.nc')
df = ds.to_dataframe()
ds.close()
df = df.append(dataframe)
dset = df.to_xarray()
dset.to_netcdf('RadialDist.nc')
dset.close()
"""