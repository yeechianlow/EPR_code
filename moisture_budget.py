#Calculate moisture budget for eastern North America as defined in Low et al. 2022
#Assumes using land/sea data at 0.5 x 0.5 degree resolution and ERA5 data at 0.25 x 0.25 degree resolution

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib as mpl
mpl.use('agg')
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import metpy
import numpy as np
import xarray as xr
import datetime

def get_distance(point1, point2):
    '''
    Calculates distance between two points on Earth.

    Parameters:
    ----------
    point1 : tuple or list
        (latitude, longitude) of first point in degrees
    point2 : tuple or list
        (latitude, longitude) of second point in degrees

    Returns:
    ----------
    float
        distance between point1 and point2 in kilometers
    '''
    R = 6370 #radius of Earth in km
    lat1 = point1[0]*np.pi/180.0
    lon1 = point1[1]*np.pi/180.0 
    lat2 = point2[0]*np.pi/180.0 
    lon2 = point2[1]*np.pi/180.0 

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    distance = R * c
    return distance

min_lat = 25.0
max_lat = 50.0
min_lon = -95.0
max_lon = -50.0

lat_res = 0.25
lon_res = 0.25

first_year = 1940
last_year = 2025

for year in range(first_year,last_year+1):

        syear = year-1
        start_time = datetime.datetime(year-1,12,1)
        end_time = datetime.datetime(year,3,1)-datetime.timedelta(hours=1)
        
        #load evaporation, IVT, precipitation, and precipitable water data from ERA5
        #note that this code calculates moisture budget every 6 hours, hence IVT and precipitable water are assumed to be given every 6 hours
        #however, evaporation and precipitation from ERA5 are given hourly, and moisture budget requires those values to be added hourly, so those datasets are assumed to be given hourly
        evap_ds=xr.open_dataset('/storage/yeechian/data/evap_ERA5_'+str(year)+'.nc').load()  
        ivt_ds=xr.open_dataset('/storage/yeechian/data/IVT_ERA5_'+str(year)+'.nc').load()  
        pcp_ds=xr.open_dataset('/storage/yeechian/data/pcp_NA_ERA5_'+str(year)+'.nc').load()
        pw_ds=xr.open_dataset('/storage/yeechian/data/pw_ERA5_'+str(year)+'.nc').load()
        if year != last_year:
            evap_ds2=xr.open_dataset('/storage/yeechian/data/evap_ERA5_'+str(year+1)+'.nc').load()  
            ivt_ds2=xr.open_dataset('/storage/yeechian/data/IVT_ERA5_'+str(year+1)+'.nc').load()  
            pcp_ds2=xr.open_dataset('/storage/yeechian/data/pcp_NA_ERA5_'+str(year+1)+'.nc').load()
            pw_ds2=xr.open_dataset('/storage/yeechian/data/pw_ERA5_'+str(year+1)+'.nc').load()
            
        if year == first_year:
            lat = ivt_ds.latitude.data[:]
            lon = ivt_ds.longitude.data[:]
            lat2d = np.zeros((len(lat),len(lon)),dtype=np.float32)
            for i in range(len(lon)):
                lat2d[:,i] = lat
            lon_360 = np.copy(lon)
            lon_360[lon<0]=lon_360[lon<0]+360
            
            land_mask = np.zeros((len(lat),len(lon)),dtype=np.float32) #create land mask for eastern North America
        
        time = ivt_ds.time.data
        
        if year == 1940:
            t_len = len(time)
        else:
            t_len = len(time)+1

        mflux_northern = np.zeros((len(time)),dtype=np.float32) 
        mflux_western = np.zeros((len(time)),dtype=np.float32) 
        mflux_northern = np.zeros((len(time)),dtype=np.float32) 
        mflux_western = np.zeros((len(time)),dtype=np.float32) 
        mflux_gofm = np.zeros((len(time)),dtype=np.float32) 
        mflux_ec = np.zeros((len(time)),dtype=np.float32) 
        mflux_total = np.zeros((len(time)),dtype=np.float32) 
        net_mbudget = np.zeros((len(time)),dtype=np.float32)
        residual = np.zeros((len(time)),dtype=np.float32)
        evap_land = np.zeros((len(time)),dtype=np.float32) 
        pcp_land = np.zeros((len(time)),dtype=np.float32) 
        storage = np.zeros((len(time)),dtype=np.float32) 

        #To calculate moisture budget for full year, need to also have data for last time step of previous year and first time step of next year
        if year == 1940: #data is missing for first few hours of 1940
            sindex = 12
            sindex2 = 6
            sindex3 = 6
            net_mbudget[0:2] = np.nan
            residual[0:2] = np.nan
            evap_land[0:2] = np.nan
            pcp_land[0:2] = np.nan
            storage[0:1] = np.nan
            evap = np.concatenate((np.asarray(evap_ds['e'].sel(latitude = slice(max_lat,min_lat), longitude = slice(min_lon,max_lon))),np.asarray(evap_ds2['e'].sel(latitude = slice(max_lat,min_lat), longitude = slice(min_lon,max_lon))[:1,:,:])),axis=0)
            pcp = np.concatenate((np.asarray(pcp_ds['tp'].sel(latitude = slice(max_lat,min_lat), longitude = slice(min_lon,max_lon))),np.asarray(pcp_ds2['tp'].sel(latitude = slice(max_lat,min_lat), longitude = slice(min_lon,max_lon))[:1,:,:])),axis=0)
            pw = np.concatenate((np.asarray(pw_ds['tcwv'].sel(latitude = slice(max_lat,min_lat), longitude = slice(min_lon,max_lon))),np.asarray(pw_ds2['tcwv'].sel(latitude = slice(max_lat,min_lat), longitude = slice(min_lon,max_lon))[:1,:,:])),axis=0)
        else:
            sindex = 0
            sindex2 = 0
            sindex3 = 1
            evap_ds3=xr.open_dataset('/storage/yeechian/data/evap_ERA5_'+str(year-1)+'.nc').load()
            pcp_ds3=xr.open_dataset('/storage/yeechian/data/pcp_NA_ERA5_'+str(year-1)+'.nc').load()
            pw_ds3=xr.open_dataset('/storage/yeechian/data/pw_ERA5_'+str(year-1)+'.nc').load()
            evap = np.concatenate((np.asarray(evap_ds3['e'].sel(latitude = slice(max_lat,min_lat), longitude = slice(min_lon,max_lon))[-5:,:,:]),np.asarray(evap_ds['e'].sel(latitude = slice(max_lat,min_lat), longitude = slice(min_lon,max_lon))),np.asarray(evap_ds2['e'].sel(latitude = slice(max_lat,min_lat), longitude = slice(min_lon,max_lon))[:1,:,:])),axis=0)
            pcp = np.concatenate((np.asarray(pcp_ds3['tp'].sel(latitude = slice(max_lat,min_lat), longitude = slice(min_lon,max_lon))[-5:,:,:]),np.asarray(pcp_ds['tp'].sel(latitude = slice(max_lat,min_lat), longitude = slice(min_lon,max_lon))),np.asarray(pcp_ds2['tp'].sel(latitude = slice(max_lat,min_lat), longitude = slice(min_lon,max_lon))[:1,:,:])),axis=0)
            pw = np.concatenate((np.asarray(pw_ds3['tcwv'].sel(latitude = slice(max_lat,min_lat), longitude = slice(min_lon,max_lon))[-6::6,:,:]),np.asarray(pw_ds['tcwv'].sel(latitude = slice(max_lat,min_lat), longitude = slice(min_lon,max_lon))),np.asarray(pw_ds2['tcwv'].sel(latitude = slice(max_lat,min_lat), longitude = slice(min_lon,max_lon))[:1,:,:])),axis=0)
        evap = np.asarray([np.sum(evap[i:i+12,:,:],axis=0) for i in range(0,len(evap)-6,6)])
        ivt_u = np.concatenate((np.asarray(ivt_ds['IVT_u'].sel(latitude = slice(max_lat,min_lat), longitude = slice(min_lon,max_lon))),np.asarray(ivt_ds2['IVT_u'].sel(latitude = slice(max_lat,min_lat), longitude = slice(min_lon,max_lon))[0:1,:,:])),axis=0)
        ivt_v = np.concatenate((np.asarray(ivt_ds['IVT_v'].sel(latitude = slice(max_lat,min_lat), longitude = slice(min_lon,max_lon))),np.asarray(ivt_ds2['IVT_v'].sel(latitude = slice(max_lat,min_lat), longitude = slice(min_lon,max_lon))[0:1,:,:])),axis=0)

        pcp = np.asarray([np.sum(pcp[i:i+12,:,:],axis=0) for i in range(0,len(pcp)-6,6)])
        
        #Points of coastline
        points = [(50.0,-53.0),(50.0,-95.01),(30.0,-95.0),(30.0,-84.0),(35.0,-76.5),(40.5,-74.0),(41.5,-70.5),(43.5,-70.49),(43.5,-66.0),(47.5,-53.0),(50.0,-52.99)]
        
        #Calculate moisture flux through each coastline segment separately
        for i in range(len(points)-1,0,-1):
       
            #Make sure we are always moving eastward on piecewise trajectory
            if points[i][1] >= points[i-1][1]:
                prev_point = points[i-1]
                curr_point = points[i]
            else:
                prev_point = points[i]
                curr_point = points[i-1]

            #Calculate moisture flux through each grid point
            slope = (curr_point[0]-prev_point[0])/(curr_point[1]-prev_point[1])
            theta = np.arctan(slope)-np.pi/2
            if i == 1 or i == 6 or i == 15: #reverse orientation of theta by 180 degrees due to coastline orientation
                theta += np.pi
            if i == 6 or i == 14 or i == 15: #reverse traveling direction due to coastline orientation
                x_range = range(int((curr_point[1]-min_lon)*int(1.0/lon_res)),int((prev_point[1]-min_lon)*int(1.0/lon_res))-1,-1)
            else:
                x_range = range(int((prev_point[1]-min_lon)*int(1.0/lon_res)),int((curr_point[1]-min_lon)*int(1.0/lon_res))+1)
            
            #Create land-ocean mask using line segments
            if i >= 4 and slope != 0 and year == first_year:
                if i == 14:
                    y_range = range(int(round((max_lat-points[i-1][0])*int(1.0/lat_res),0)),int(round((max_lat-points[i][0])*int(1.0/lat_res),0))+1)
                else:
                    y_range = range(int(round((max_lat-points[i-1][0])*int(1.0/lat_res),0)),int(round((max_lat-points[i][0])*int(1.0/lat_res),0))-1,-1)
                for y in y_range:
                    for x in range(int(round(((max_lat-points[i-1][0])*int(1.0/lat_res)-y)/slope+(points[i-1][1]-min_lon)*int(1.0/lon_res),0))+1):
                        if i == 14: #special case due to coastline curvature
                            land_mask[y+int((lat[0]-max_lat)*int(1.0/lat_res)),x+int((min_lon-lon[0])*int(1.0/lon_res))] = 0
                        else:
                            land_mask[y+int((lat[0]-max_lat)*int(1.0/lat_res)),x+int((min_lon-lon[0])*int(1.0/lon_res))] = 1
                       
            prev_x = -999
            prev_y = -999
            for x in x_range:
                if slope > 1:
                    max_y = min(int(round(-slope*(x-0.5-(prev_point[1]-min_lon)*int(1.0/lon_res))+(max_lat-prev_point[0])*int(1.0/lat_res),0)),int(round((max_lat-prev_point[0])*int(1.0/lat_res),0)))
                    min_y = max(int(round(-slope*(x+0.5-(prev_point[1]-min_lon)*int(1.0/lon_res))+(max_lat-prev_point[0])*int(1.0/lat_res),0)),int(round((max_lat-curr_point[0])*int(1.0/lat_res),0)))-1              
                elif slope < -1:
                    min_y = max(int(round(-slope*(x-0.5-(prev_point[1]-min_lon)*int(1.0/lon_res))+(max_lat-prev_point[0])*int(1.0/lat_res),0)),int(round((max_lat-prev_point[0])*int(1.0/lat_res),0)))
                    max_y = min(int(round(-slope*(x+0.5-(prev_point[1]-min_lon)*int(1.0/lon_res))+(max_lat-prev_point[0])*int(1.0/lat_res),0)),int(round((max_lat-curr_point[0])*int(1.0/lat_res),0)))+1
                else:
                    max_y = int(round(-slope*(x-(prev_point[1]-min_lon)*int(1.0/lon_res))+(max_lat-prev_point[0])*int(1.0/lat_res),0))
                    min_y = max_y-1
                for y in range(max_y,min_y,-1):
                    if (prev_x != -999) and (prev_y != -999) and (np.abs(slope) <= 1 or int(round(((max_lat-prev_point[0])*int(1.0/lat_res)-y)/slope+(prev_point[1]-min_lon)*int(1.0/lon_res),0)) == x): #use IVT at end of line segment to calculate flux through line segment, and include as small of line segments as possible
                        mflux = -1000*(ivt_u[:,y,x]*np.cos(theta)+ivt_v[:,y,x]*np.sin(theta))*get_distance((max_lat-prev_y*lat_res,prev_x*lon_res+min_lon),(max_lat-y*lat_res,x*lon_res+min_lon))
                        
                        #Split flux into flux through northern boundary, western boundary, Gulf of Mexico, and Atlantic
                        if i == 1 and x*0.25+min_lon != -50.0:
                            mflux_northern += np.asarray([(mflux[j]+mflux[j+1])/2.0 for j in range(len(mflux)-1)])
                        elif i == 2:
                            mflux_western += np.asarray([(mflux[j]+mflux[j+1])/2.0 for j in range(len(mflux)-1)])
                        elif i == 3:
                            mflux_gofm += np.asarray([(mflux[j]+mflux[j+1])/2.0 for j in range(len(mflux)-1)])
                        else:
                            mflux_ec += np.asarray([(mflux[j]+mflux[j+1])/2.0 for j in range(len(mflux)-1)])
                        mflux_total += np.asarray([(mflux[j]+mflux[j+1])/2.0 for j in range(len(mflux)-1)])

                        prev_x = x
                        prev_y = y

                    if prev_x == -999 and (np.abs(slope) <= 1 or int(round(((max_lat-prev_point[0])*int(1.0/lat_res)-y)/slope+(prev_point[1]-min_lon)*int(1.0/lon_res),0)) == x): #start of line segment                       
                        prev_x = x
                        prev_y = y


        weighted_land_mask = land_mask[int((lat[0]-max_lat)*int(1.0/lat_res)):int((lat[0]-min_lat)*int(1.0/lat_res))+1,int((min_lon-lon[0])*int(1.0/lon_res)):int((max_lon-lon[0])*int(1.0/lon_res))+1]*np.cos(lat2d[int((lat[0]-max_lat)*int(1.0/lat_res)):int((lat[0]-min_lat)*int(1.0/lat_res))+1,int((min_lon-lon[0])*int(1.0/lon_res)):int((max_lon-lon[0])*int(1.0/lon_res))+1]*np.pi/180) #land mask weighted by latitude-dependent area of each grid cell
        area = np.nansum((111000/4)**2*weighted_land_mask)
        evap_land[sindex//6:] = -1.0/43.2*np.nansum((111000/4)**2*evap*weighted_land_mask,axis=(1,2))
        pcp_land[sindex//6:] = 1.0/43.2*np.nansum((111000/4)**2*pcp*weighted_land_mask,axis=(1,2))
        pw_land = np.nansum((111000/4)**2*pw*weighted_land_mask,axis=(1,2))
        net_mbudget[sindex//6:] = mflux_total[sindex//6:]+evap_land[sindex//6:]-pcp_land[sindex//6:]
        if year == 1940: #First few hours of 1940 are missing
            storage[1:] = np.asarray([1.0/43200*(pw_land[i+1] - pw_land[i-1]) for i in range(1,len(pw_land)-1)])
        else:
            storage = np.asarray([1.0/43200*(pw_land[i+1] - pw_land[i-1]) for i in range(1,len(pw_land)-1)])
        residual[sindex//6:] = net_mbudget[sindex//6:] - storage[sindex//6:]

        evap_ds.close()
        ivt_ds.close()
        pcp_ds.close()
        
        new_ds = xr.Dataset(data_vars={'area':(('t'),area*np.ones(storage.shape)),'storage':(('t'),storage),'evaporation':(('t'),evap_land),'precipitation':(('t'),pcp_land),'residual':(('t'),residual),'mflux_N':(('t'),mflux_northern),'mflux_W':(('t'),mflux_western),'mflux_GofM':(('t'),mflux_gofm),'mflux_EC':(('t'),mflux_ec)},coords={'time':time,'latitude':lat,'longitude':lon})
        new_ds.area.attrs["units"]="m^2"
        new_ds.area.attrs["long_name"]="Total area of eastern North America over which moisture budget is calculated"
        new_ds.storage.attrs["units"]="m^2 mm day^-1"
        new_ds.storage.attrs["long_name"]="12-hr averaged rate of moisture storage or increase in moisture over eastern North America"
        new_ds.evaporation.attrs["units"]="m^2 mm day^-1"
        new_ds.evaporation.attrs["long_name"]="12-hr averaged evaporation rate over eastern North America"
        new_ds.precipitation.attrs["units"]="m^2 mm day^-1"
        new_ds.precipitation.attrs["long_name"]="12-hr averaged precipitation rate over eastern North America"
        new_ds.residual.attrs["units"]="m^2 mm day^-1"
        new_ds.residual.attrs["long_name"]="Residual term of moisture budget (difference between calculated moisture budget and actual storage term)"
        new_ds.mflux_N.attrs["units"]="m^2 mm day^-1"
        new_ds.mflux_N.attrs["long_name"]="12-hr averaged moisture flux through northern boundary"
        new_ds.mflux_W.attrs["units"]="m^2 mm day^-1"
        new_ds.mflux_W.attrs["long_name"]="12-hr averaged moisture flux through western boundary"
        new_ds.mflux_GofM.attrs["units"]="m^2 mm day^-1"
        new_ds.mflux_GofM.attrs["long_name"]="12-hr averaged moisture flux through Gulf of Mexico boundary"
        new_ds.mflux_EC.attrs["units"]="m^2 mm day^-1"
        new_ds.mflux_EC.attrs["long_name"]="12-hr averaged moisture flux through East Coast boundary"
        new_ds.to_netcdf(path="/storage/yeechian/data/moisture_budget_ERA5_"+str(year)+".nc",mode="w",format="NETCDF4")

        #plot land-ocean mask and line segments through which moisture flux is calculated
        if year == first_year:
            land = cfeature.NaturalEarthFeature(category='cultural', name='admin_0_countries',scale='50m',facecolor='none')
            states_provinces = cfeature.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lakes',scale='50m',facecolor='none') # Grab data for plotting state boundaries
            gs = gridspec.GridSpec(2, 1, height_ratios=[1, .02], bottom=.07, top=.99,hspace=0.01, wspace=0.01)

            fig = plt.figure(figsize=(17., 11.))
            ax = plt.subplot(111, projection=ccrs.PlateCarree())
            gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1.5, color='black', alpha=0.5, linestyle='--')
            gl.xlocator = mticker.FixedLocator(np.arange(-105,45,5))
            gl.ylocator = mticker.FixedLocator(np.arange(20,60,5))
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER
            gl.xlabel_style = {'size': 24,'weight': 'bold'}
            gl.ylabel_style = {'size': 24,'weight': 'bold'}
            gl.xlabels_top = False
            plt.title('Land-ocean mask \n \n' ,{'fontsize':24})

            # Set extent and plot map lines
            ax.set_extent([-105.0,-50.0,20.0,55.0], ccrs.PlateCarree())
            ax.add_feature(land, edgecolor=(0.2,0.2,0.2), linewidth=1.5)
            ax.add_feature(states_provinces, edgecolor=(0.2,0.2,0.2), linewidth=1.5)
            cf1 = ax.contourf(lon_360, lat, land_mask, [0.5,1.0],extend='both', colors=['white','brown','brown'],alpha=0.45,transform=ccrs.PlateCarree())
            lcolors = ['b-','m-','c-','y-','y-','y-','y-','y-','y-','y-','y-']
            for pindex in range(-1,len(points)-1):
                plt.plot([points[pindex][1],points[pindex+1][1]],[points[pindex][0],points[pindex+1][0]],lcolors[pindex],linewidth=5.0)

            plt.savefig('/storage/yeechian/figures/land-ocean_mask.png')
            plt.close()
        
