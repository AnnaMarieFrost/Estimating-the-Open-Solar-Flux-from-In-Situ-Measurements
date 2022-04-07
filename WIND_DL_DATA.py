"""
Created on Thu Apr 7 2022

Analysis carried out from: Estimating the Open Solar Flux from In Situ Measurements

@author: Anna M Frost
         a.m.frost@pgr.reading.ac.uk
"""


###############################################################################
"""
This script downloads the Wind 3dp data files and saves them as a common name. 
Lines 21 and 25 are specific subfolders that keep everything tidy.
"""
###############################################################################

import datetime
from datetime import timedelta
import requests
import os
import os.path


def DL_3DP(date, eesa):#eesa as an input so you can specify the higher or lower energy analyser
	s = requests.Session()  #requests library is the standard for making HTTP requests in Python
	
	if eesa == 'h':  #for the high energy analyser
		windurl = 'https://cdaweb.gsfc.nasa.gov/pub/data/wind/3dp/3dp_ehpd/' #url for the wind instrument
		for i in range(1,5):
			fname = 'wi_ehpd_3dp_'+str(date[0])+date[1]+date[2]+'_v0'+str(i)+'.cdf' #filenames as found here: https://cdaweb.gsfc.nasa.gov/pub/data/wind/3dp/3dp_ehpd/1994/ with variable date inputs
		direct = './EESA_H_CDFs/' #this is the subdirectory for saving data, make sure folders are set up
	elif eesa == 'l':  #for the low energy analyser
		windurl = 'https://cdaweb.gsfc.nasa.gov/pub/data/wind/3dp/3dp_elpd/' #url for the wind instrument
		fname = 'wi_elpd_3dp_'+str(date[0])+str(date[1])+str(date[2])+'_v02.cdf' #filenames as found here: https://cdaweb.gsfc.nasa.gov/pub/data/wind/3dp/3dp_ehpd/1994/ with variable date inputs
		direct = './EESA_L_CDFs/'	
	url = windurl+str(date[0])+'/'+fname #put the whole thing together to find the individual data files for each specified date
	fnamesav=direct+'wi_elpd_3dp_'+str(date[0])+str(date[1])+str(date[2])+'.cdf' #specifies the sub folder and filename for saving
    
	if os.path.isfile(''+fnamesav): # this checks if we've already got the file so we don't download it
		print('File already downloaded!')
		status = 200  # Means the server successfully answered the http request
	else:
		r = s.get(url) #this checks if something is at the url
		if r.status_code == 200:  #if there is then we download it and save it
			with open(fnamesav, 'wb') as out:  #Which opens a file in the chosen directory and name for the data
													#'wb' indicates that the file is opened for writing in binary mode
													# out is used to reference that file
				for bits in r.iter_content(): #takes the info from web address in r and loops over different parts of it
												# "bits" is the bit of data for each loop
					out.write(bits)   #takes file referenced by "out" and writes in the data contained in "bits"
			print('Saved File: ', date)
			out.close()   #closes the file referenced by out
		if r.status_code == 404:    # if there's not then we're sad
			if os.path.isfile(fnamesav):	os.remove(fnamesav)    #and we remove the empty file we created
		status = r.status_code

	return(status,fname)    #return info on if we got the file, and what it's name is so we can find it



start_date = datetime.date(1994,11,5)
end_date = datetime.date(2022,3,8)
 

delta = timedelta(days=1)

while start_date <= end_date:
    a = start_date.strftime("%Y%m%d")
    Year = (a[0:4])
    Month = (a[4:6])
    Day = (a[6:8])
    date = [Year, Month, Day]
    status, fname = DL_3DP(date, 'l')
    start_date+=delta

