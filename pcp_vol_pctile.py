#calculates percentiles of normalized winter daily precipitation volume in eastern North America
import xarray as xr
import numpy as np

all_pcp_vol = []
for year in range(1991,2021): #use 1991-2020 climatology
        
    ds = xr.open_dataset('/storage/yeechian/data/pcp_vol_ERA5_'+str(year)+'.nc').load()
    pcp_vol = ds['pcp_vol_norm80']

    #obtain winter daily precipitation volume
    if (year%4 == 0): #leap year
        pcp_vol1 = pcp_vol[1:61]
        pcp_vol2 = pcp_vol[336:367]
    else:
        pcp_vol1 = pcp_vol[1:60]
        pcp_vol2 = pcp_vol[335:366]

    all_pcp_vol += list(np.concatenate((pcp_vol1,pcp_vol2),axis=0))

#output percentiles into file
new_ds = xr.Dataset(data_vars={'pctile_djf_pcpv1d_ERA5':(('percentile'),np.percentile(all_pcp_vol,range(1,101)))},coords={'percentile':range(1,101)})           
new_ds.pctile_djf_pcpv1d_ERA5.attrs["units"]="km^2"
new_ds.pctile_djf_pcpv1d_ERA5.attrs["long_name"]="Percentiles of winter daily precipitation volume over eastern North America normalized by 80th percentile of daily precipitation"
new_ds.to_netcdf(path="/storage/yeechian/data/precip_pctiles_ERA5.nc",mode="a",format="NETCDF4")    


