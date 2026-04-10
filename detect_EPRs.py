#makes a csv file containing EPRs
import pandas as pd
import numpy as np
import xarray as xr
import datetime

pctile = 70
min_len = 5 #minimum length of regime
#percentiles = pd.read_csv('/storage/yeechian/data/data_csv/namer_pcp_vol_norm80_ERA5_pctiles_1991-2020.csv',delimiter=':',names=("pctile","value"),header=None)
#threshold = float(percentiles[percentiles["pctile"] == pctile]["value"])
pctile_ds = xr.open_dataset("/storage/yeechian/data/precip_pctiles_ERA5.nc")
threshold = pctile_ds["pctile_djf_pcpv1d_ERA5"].sel(percentile = pctile)

EPR_csv = open('/storage/yeechian/data/data_csv/EPRs_ERA5.csv','w')
first_year = 1941
last_year = 2025
for year in range(first_year,last_year+1):
    num_athres = 0 #number of days above threshold
    num_bthres = 0 #number of days below threshold
    num_bthres_cs = 0 #number of days below threshold consecutively
    num_days = 0 #duration of potential EPR
    new_regime_trigger = True
    is_regime = False

    pcp_vol1 = xr.open_dataset('/storage/yeechian/data/pcp_vol_ERA5_'+str(year-1)+'.nc')['pcp_vol_norm80']
    pcp_vol2 = xr.open_dataset('/storage/yeechian/data/pcp_vol_ERA5_'+str(year)+'.nc')['pcp_vol_norm80']

    if year%4 == 0: #current year is leap year
        pcp_vol = np.concatenate((pcp_vol1[334:365],pcp_vol2[:60]))
    elif (year-1)%4 == 0: #previous year is leap year
        pcp_vol = np.concatenate((pcp_vol1[335:366],pcp_vol2[:59]))
    else:
        pcp_vol = np.concatenate((pcp_vol1[334:365],pcp_vol2[:59]))

    length = len(pcp_vol)
    index = 0
    while index <= length-1: #until end of season
        '''

        if index <= length-min_len: #can use full moving window
            moving_window = pcp_vol[index:index+min_len] #moving window of min_len number of days
        else: #moving window cut off by end of data
            moving_window = pcp_vol[index:] #shortened moving window

        #calculate number of days in moving window above threshold
        num_window_athres = 0 
        for elem in moving_window:
            if elem > threshold:
                num_window_athres += 1

        if not new_regime_trigger and num_window_athres < min_len-2 and index >= min_len-1 and index <= length-min_len: #end of regime
            new_regime_trigger = True #whatever existing regime has ended, turn trigger on
            endtime = datetime.datetime(year-1,12,1)+datetime.timedelta(days=index+5)
            endtime = '%04d%02d%02d'%(endtime.year,endtime.month,endtime.day)

        if index <= length-min_len and new_regime_trigger and num_window_athres >= min_len-2: #new regime: period is not in an existing regime and entire period is above 90th percentile       
            starttime = datetime.datetime(year-1,12,1)+datetime.timedelta(days=index)
            starttime = '%04d%02d%02d'%(starttime.year,starttime.month,starttime.day)
            new_regime_trigger = False #turn off trigger for regime detection
            index = index + 1 #skip until after period
        else:
            index = index + 1 #increment by one timestamp
        '''

        if pcp_vol[index] >= threshold+0.01:
            if new_regime_trigger and index < length - min_len:
                starttime = datetime.datetime(year-1,12,1)+datetime.timedelta(days=index)
                starttime = '%04d%02d%02d'%(starttime.year,starttime.month,starttime.day)
            if num_days >= min_len and (num_athres/num_days >= 2.0/3): #already qualifies as a regime
                is_regime = True
            new_regime_trigger = False #turn off trigger for regime detection
            num_athres = num_athres + 1
            num_bthres_cs = 0
            num_days = num_days + 1
        elif not new_regime_trigger:
            num_bthres = num_bthres + 1
            num_bthres_cs = num_bthres_cs + 1
            if num_bthres_cs >= 2 or (is_regime and num_athres/num_days < 2.0/3):
                new_regime_trigger = True
                if (num_days-1 >= min_len and num_athres/(num_days-1) >= 2.0/3) or (num_days >= min_len and is_regime): #lasts long enough and at least 2/3 of days meet threshold
                    if num_bthres_cs >= 2:
                        endtime = datetime.datetime(year-1,12,1)+datetime.timedelta(days=index-2)
                        endtime = '%04d%02d%02d'%(endtime.year,endtime.month,endtime.day)
                    else:
                        endtime = datetime.datetime(year-1,12,1)+datetime.timedelta(days=index-1)
                        endtime = '%04d%02d%02d'%(endtime.year,endtime.month,endtime.day)

                    EPR_csv.write(starttime+'-'+endtime+'\n') #output time period of regime

                num_athres = 0
                num_bthres = 0
                num_bthres_cs = 0
                num_days = 0
                is_regime = False
            else:
                num_days = num_days + 1
        index = index + 1

    if not new_regime_trigger and ((num_days >= min_len and num_athres/num_days >= 2.0/3 and num_bthres_cs == 0) or (num_days-1 >= min_len and num_athres/(num_days-1) >= 2.0/3 and num_bthres_cs == 1)): #end of regime not set because of end of season
        if num_bthres_cs == 0:
            endtime = datetime.datetime(year-1,12,1)+datetime.timedelta(days=length-1) #make end of regime equal to end of season
            endtime = '%04d%02d%02d'%(endtime.year,endtime.month,endtime.day)
        else: #last day is not part of regime since there was no heavy precip
            endtime = datetime.datetime(year-1,12,1)+datetime.timedelta(days=length-2)
            endtime = '%04d%02d%02d'%(endtime.year,endtime.month,endtime.day)

        new_regime_trigger = True
        EPR_csv.write(starttime+'-'+endtime+'\n') #output time period of regime
        num_athres = 0
        num_bthres = 0
        num_bthres_cs = 0
        num_days = 0
        is_regime = False

EPR_csv.close()

