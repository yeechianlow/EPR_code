#creates file with daily precipitation from ERA5 dataset (assumed to have hourly temporal resolution for entire year)
import xarray as xr
import numpy as np

first_year = 1940
last_year = 2025

for year in range(first_year,last_year+1):
    pcp_ds = xr.open_dataset('/storage/yeechian/data/pcp_NA_ERA5_'+str(year)+'.nc').load()
    precip_year = pcp_ds['tp']

    #calculate number of days in current year
    if year%4 == 0: #leap year
        num_days = 366
    else: #non-leap year
        num_days = 365

    pcp_day = np.zeros((num_days,len(pcp_ds['latitude']),len(pcp_ds['longitude'])))
    for day in range(num_days-1):
        if year == 1940: #first few hours of 1940 are missing
            if day != 0: #leave first day of 1940 blank due to missing data
                precip = precip_year[day*24-6:day*24+18,:,:]
            else:
                continue
        else:
            precip = precip_year[day*24+1:day*24+25,:,:]
        pcp_day[day,:,:] = 1000*np.sum(precip,axis=0)

    if year != last_year: #take care of last day of year (if not last year)

        precip_next = xr.open_dataset('/storage/yeechian/data/pcp_NA_ERA5_'+str(year+1)+'.nc')['tp'][:1,:,:]
        if year == 1940: #first few hours of 1940 are missing
            precip = precip_year[num_days*24-30:,:,:]
        else:
            precip = precip_year[num_days*24-23:,:,:]
            
        precip = np.concatenate((precip,precip_next),axis=0)
        pcp_day[num_days-1,:,:] = 1000*(np.sum(precip,axis=0))

    #output daily total precipitation to file
    new_ds = xr.Dataset(data_vars={'pcp_1day':(('day','latitude','longitude'),pcp_day)},coords={'day':range(1,num_days+1),'latitude':pcp_ds['latitude'],'longitude':pcp_ds['longitude']})
    new_ds.pcp_1day.attrs["units"]="mm"
    new_ds.pcp_1day.attrs["long_name"]="Daily precipitation"
    new_ds.to_netcdf(path="/storage/yeechian/data/precip_daily_"+str(year)+"_ERA5.nc",mode="w",format="NETCDF4")    

    
