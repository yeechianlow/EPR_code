#Compute maximum equivalent potential temperature (theta-e) between 700 hPa and 900 hPa inclusively

import metpy.calc as mpcalc
import metpy.interpolate
from metpy.units import units
import numpy as np
import xarray as xr
levels = [700,750,775,800,825,850,875,900]

for year in range(1977,1978):
    ds = xr.open_dataset('/storage/yeechian/data/tr_ERA5_%d.nc'%year) #Open ERA5 data file containing temperature and relative humidity with coordinates of (time, pressure level, latitude, longitude)
    thetae = np.zeros(ds['r'].shape,dtype=np.float32)
    max_thetae = np.zeros(ds['r'].metpy.sel(level=700).shape,dtype=np.float32) #Maximum theta-e values
    max_thetae_pres = np.zeros(ds['r'].metpy.sel(level=700).shape,dtype=np.int64) #Pressure at which maximum theta-e is achieved

    #First calculate theta-e at all levels between 700 hPa and 900 hPa inclusively
    for i in range(len(levels)):
        t = ds['t'].metpy.sel(level=levels[i])
        rh = np.asarray(ds['r'].metpy.sel(level=levels[i]))
        rh[rh<=0.0] = 0.001 #Remove any negative RH values
        td = mpcalc.dewpoint_from_relative_humidity(units('K')*t,rh*0.01) 
        thetae[:,i,:,:] = mpcalc.equivalent_potential_temperature(units('hPa')*levels[i],units('K')*t,td)

    #Now calculate maximum theta-e and pressure at which maximum theta-e occurs
    max_ind = np.argmax(thetae,axis=1)
    for t in range(len(thetae)):
        for y in range(thetae.shape[2]):
            for x in range(thetae.shape[3]):
                max_thetae[t,y,x] = thetae[t,max_ind[t,y,x],y,x]
                max_thetae_pres[t,y,x] = levels[max_ind[t,y,x]]
           
    new_ds = xr.Dataset(data_vars={'max_700_900_thetae':(('time','latitude','longitude'),np.asarray(max_thetae).astype('f4')),'max_thetae_pres':(('time','latitude','longitude'),np.asarray(max_thetae_pres))},coords={'time':ds['time'],'latitude':ds['latitude'],'longitude':ds['longitude']})           
    new_ds.max_700_900_thetae.attrs["units"]="K"
    new_ds.max_700_900_thetae.attrs["long_name"]="Equivalent potential temperature"
    new_ds.max_thetae_pres.attrs["units"]="hPa"
    new_ds.max_thetae_pres.attrs["long_name"]="Pressure at which maximum theta-e occurs"
    new_ds.to_netcdf(path="/storage/yeechian/data/data_csv/max_thetae_ERA5_%d.nc"%(year),mode="w",format="NETCDF4")
