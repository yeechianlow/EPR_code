#calculates and stores daily precipitation volume divided by 80th percentile of daily precipitation
import xarray as xr
import numpy as np
import math
import metpy

land_sea_ds = xr.open_dataset("/storage/yeechian/data/mask_file_NA.nc") #get file containing land/sea data at 0.25 x 0.25 degree resolution
landsea_mask = land_sea_ds["lsmask"].sel(lat = slice(50,25),lon = slice(-95,-50))

pi = 3.141592654
pctile_ds = xr.open_dataset("/storage/yeechian/data/precip_pctiles_ERA5.nc")
pct80_pcp = pctile_ds["pctile_djf_pcp1d_ERA5"].sel(percentile = 80, latitude = slice(50,25),longitude = slice(-95,-50))

lat_res = 0.25
lon_res = 0.25

first_year = 1940
last_year = 2025

for year in range(first_year,last_year+1):
    pcp_ds = xr.open_dataset("/storage/yeechian/data/precip_daily_"+str(year)+"_ERA5.nc") 
    pcp = pcp_ds["pcp_1day"].sel(latitude = slice(50,25),longitude = slice(-95,-50))

    if (year%4 == 0): #leap year
        num_days = 366
    else:
        num_days = 365
            
    volume = np.zeros((num_days),dtype=np.float32)

    #calculate precipitation volume by summing precipitation at each grid point multiplied by each grid point's area and by 80th percentile of daily precipitation at each grid point
    for y in range(pcp.shape[1]):
        print("y = "+str(y))
        weight = math.cos((-y*lat_res + 50)*pi/180) #weight grid point based on its latitude, since grid points are denser in higher latitudes
        for x in range(pcp.shape[2]):
            if landsea_mask[y,x] >= 1: #not in ocean
                volume += pcp[:,y,x]*weight*(111.0/4)**2/pct80_pcp[y,x]
           

    new_ds = xr.Dataset(data_vars={'pcp_vol_norm80':(('day'),np.asarray(volume))},coords={'day':range(num_days)})
    new_ds.pcp_vol_norm80.attrs["units"]="km^2"
    new_ds.pcp_vol_norm80.attrs["long_name"]="Daily precipitation volume over eastern North America normalized by 80th percentile of daily precipitation"
    new_ds.to_netcdf(path="/storage/yeechian/data/pcp_vol_ERA5_"+str(year)+".nc",mode="w",format="NETCDF4")
