"""
Created on Thu Apr 7 2022

Analysis carried out from: Estimating the Open Solar Flux from In Situ Measurements

@author: Anna M Frost
         a.m.frost@pgr.reading.ac.uk
"""


import numpy as np
import datetime
from datetime import timedelta
import OSF_Best_Estimate_ChangingPABandPAO_Best_Fluxes as OBEChange

Br_best = OBEChange.BR
BGFlux_best = OBEChange.BGFLUX
FluxStrTemp_best = OBEChange.FLUXSTRTEMP
FluxAnStrTemp_best = OBEChange.FLUXANSTRTEMP


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
        #N_missing_j = []

        #start_date = datetime.date(1998,1,22)
        #end_date = datetime.date(2017,12,31)

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
        dates = []
        Perc_available_ace = []
        Perc_available_wind = []
        start_date = datetime.datetime(1994, 12, 4, 3, 57, 10)
        end_date = datetime.datetime(2022,2,15)
        delta = timedelta(days=27.27)
        
        i = 0
        while start_date <= end_date:
                #end_date_temp = start_date + delta
                if len(Br_best[i])>0:
                        OSF_i, N_AS_i, N_SS_i, N_CL_i, N_U_i, PercentageAvailable_i, Perc_missing_i = Looping(Br_best[i], BGFlux_best[i], FluxStrTemp_best[i], FluxAnStrTemp_best[i], BGboundary, Opposite)
                        
                        OSF.append(OSF_i)
                        N_AS.append(N_AS_i)
                        N_SS.append(N_SS_i)
                        N_CL.append(N_CL_i)
                        N_U.append(N_U_i)
                        Perc_missing.append(Perc_missing_i)
                        
                        dates.append(start_date)
                                
                        print('start date:', start_date)
                        
                        i = i+1
                        start_date += delta
                else:
                        OSF.append(np.nan)
                        N_AS.append(np.nan)
                        N_SS.append(np.nan)
                        N_CL.append(np.nan)
                        N_U.append(np.nan)
                        Perc_missing.append(np.nan)
                        
                        dates.append(start_date)
                                
                        print('start date:', start_date)
                        
                        i = i+1
                        start_date += delta


        osf = np.nanmean(np.array(OSF))
        n_as = np.nanmean(np.array(N_AS))
        n_ss = np.nanmean(np.array(N_SS))
        n_cl = np.nanmean(np.array(N_CL))
        n_u = np.nanmean(np.array(N_U))
        #perc_missing = np.nanmean(np.array(Perc_missing))
        
        #return(osf, n_as, n_ss, n_cl, n_u, perc_missing, OSF, N_AS, N_SS, N_CL, N_U, Perc_missing, dates, Perc_available_ace, Perc_available_wind)
        return(osf, n_as, n_ss, n_cl, n_u)
        


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
        
        OSF_i, N_AS_i, N_SS_i, N_CL_i, N_U_i = LoopingFinal(Br_best, BGFlux_best, FluxStrTemp_best, FluxAnStrTemp_best, BGboundary, Opposite)
        
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
np.savetxt('TopologyTextFiles/N_AS_Best.txt', N_AS)   
np.savetxt('TopologyTextFiles/N_SS_Best.txt', N_SS)
np.savetxt('TopologyTextFiles/N_CL_Best.txt', N_CL)   
np.savetxt('TopologyTextFiles/N_U_Best.txt', N_U)   
np.savetxt('TopologyTextFiles/OSF_Best.txt', OSF) 
  
np.savetxt('TopologyTextFiles/opposite_Best.txt', opposite)   
np.savetxt('TopologyTextFiles/bgboundary_Best.txt', bgboundary)
"""
