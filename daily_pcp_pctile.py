#calculates and stores winter daily precipitation percentiles

import xarray as xr
import numpy as np

start_year = 1941
end_year = 2022

for year in range(start_year,end_year+1):
    precip_ds = xr.open_dataset('/storage/yeechian/data/precip_daily_'+str(year)+'_ERA5.nc')
    
    if year%4 == 0: #leap year, exclude leap day from data
        precip = np.concatenate((precip_ds["pcp_1day"][0:59,:,:],precip_ds["pcp_1day"][336:366,:,:]),axis=0)
    else:
        precip = np.concatenate((precip_ds["pcp_1day"][0:59,:,:],precip_ds["pcp_1day"][335:365,:,:]),axis=0)

    if year == start_year:
        cum_precip = np.copy(precip)
    else:
        cum_precip = np.concatenate((cum_precip,precip),axis=0) #update cumulative precipitation

#calculate daily precipitation percentiles at each grid point
pctile = np.zeros((100,precip.shape[1],precip.shape[2]),dtype=np.float32)
for i in range(precip.shape[1]):
    for j in range(precip.shape[2]):
        #filter out zero or non-measurable precipitation amounts
        filtered_cum_precip = []
        for k in range(cum_precip.shape[0]):
            if cum_precip[k,i,j] >= 0.2: 
                filtered_cum_precip.append(cum_precip[k,i,j])
 
        if len(filtered_cum_precip) > 0:#number of days that actually had measurable precipitation
            for percentile in range(1,101):
                pctile[percentile-1,i,j] = np.percentile(filtered_cum_precip,percentile)
        else: #percentiles are NaN because no data
            pctile[:,i,j] = float("NaN")

#output percentiles into file
new_ds = xr.Dataset(data_vars={'pctile_djf_pcp1d_ERA5':(('percentile','latitude','longitude'),pctile)},coords={'percentile':range(1,101),'latitude':precip_ds['latitude'],'longitude':precip_ds['longitude']})           
new_ds.pctile_djf_pcp1d_ERA5.attrs["units"]="mm"
new_ds.pctile_djf_pcp1d_ERA5.attrs["long_name"]="Percentiles of winter daily precipitation"
new_ds.to_netcdf(path="/storage/yeechian/data/precip_pctiles_ERA5.nc",mode="w",format="NETCDF4")    


