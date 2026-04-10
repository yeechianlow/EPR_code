#Run a self-organizing map (SOM) on 1000-500 hPa thickness and SLP maps during EPRs using ERA5 data

from minisom import MiniSom
import pandas as pd
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib as mpl
import metpy
import math
mpl.use('agg')
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.colors as colors
import datetime

#All input and climatology files are assumed to have 0.5 x 0.5 degree resolution over entire Northern Hemisphere with data every 6 hours in January to March and November to December
#slp_NH_ERA5_yyyy.nc must have sea-level pressure
#z_NH_ERA5_yyyy.nc must have geopotential height at 500-hPa and 1000-hPa
#NH_climo_ERA5_1991-2020.nc must have SLP and 1000-500 hPa thickness climatology for 1991-2020

#Set coordinates to area to calculate SOM based on
ds2 = xr.open_dataset('/storage/yeechian/data/slp_NH_ERA5_1979.nc')
lat = ds2.latitude.data[60:131]
lon = ds2.longitude.data[2:261]
lon_360 = np.copy(lon)
lon_360[lon<0]=lon_360[lon<0]+360
min_lon = -179
max_lat = 60
max_lon = -50 
min_lat = 25 

#Define contour levels and grid initializations
clevs_thck = np.asarray([-42,-36,-30,-24,-18,-12,-6,6,12,18,24,30,36,42])*1
clevs_slp = np.asarray([-22,-20,-18,-16,-14,-12,-10,-8,-6,-4,-2,2,4,6,8,10,12,14,16,18,20,22])*2
clevs_pcp = np.asarray([2,5,10,15,20,30,40,50,75,100])
gs = gridspec.GridSpec(2, 1, height_ratios=[1, .02], bottom=.07, top=.99,hspace=0.01, wspace=0.01)
land = cfeature.NaturalEarthFeature(category='cultural', name='admin_0_countries_lakes',scale='50m',facecolor='none')
states_provinces = cfeature.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lakes',scale='50m',facecolor='none')

climo_ds = xr.open_dataset('/storage/yeechian/data/NH_climo_ERA5_1991-2020.nc') #File containing climatology to calculate anomaly from

#Loop over various timestamps of EPR
t_str = ["Start-3","Start-2","Start-1","Start","Start+1","Start+2","Mid-1","Mid","Mid+1","End-2","End-1","End","End+1","End+2"]
t_str2 = ["s-3","s-2","s-1","s","s+1","s+2","m-1","m","m+1","e-2","e-1","e","e+1","e+2"]
for met in range(14):
    fi = open("/storage/yeechian/data/data_csv/"+(t_str[met]).lower()+"_pvol_SE5_70_ERA5_1941.csv","r") #csv file contains dates of timestamps in EPRs, with each line being in yyyymmddhh format
    lines = fi.readlines()[1:]

    #Initialize arrays of maps of each EPR
    size = (len(lines),71,259)
    thck_ano = np.zeros(size)
    slp_ano = np.zeros(size)

    #Flatten maps to 1D arrays for SOM
    thck_ano_vec = np.zeros((size[0],(int(size[1]/2)+1)*(int(size[2]/2)+1)))
    slp_ano_vec = np.zeros((size[0],(int(size[1]/2)+1)*(int(size[2]/2)+1)))

    tpcp = np.zeros(size)
    tpcp_vec = np.zeros((len(lines),size[1]*size[2]))

    for i in range(size[0]): #Loop over each EPR

        line = lines[i]
        syear = int(line[0:4])
        smonth = int(line[4:6])
        sday = int(line[6:8])
        shour = int(line[8:10])
        end_time = datetime.datetime(syear,smonth,sday,shour,0,0)+datetime.timedelta(hours=6)

        #Calculate time indices in climatological and actual data files
        if smonth >= 11:
            cindex = int(604-(datetime.datetime(1980,1,1)-datetime.datetime(1979,smonth,sday,shour,0,0)).total_seconds()/21600)
            index = int((datetime.datetime(syear,3,1)-datetime.datetime(syear,2,1)).total_seconds()/21600+492-(datetime.datetime(1980,1,1)-datetime.datetime(1979,smonth,sday,shour,0,0)).total_seconds()/21600)
        else:
            cindex = int((datetime.datetime(1980,smonth,sday,shour,0,0)-datetime.datetime(1980,1,1)).total_seconds()/21600)
            index = cindex

        if end_time.month >= 11:
            cindex2 = int(604-(datetime.datetime(1980,1,1)-datetime.datetime(1979,end_time.month,end_time.day,end_time.hour,0,0)).total_seconds()/21600)
            index2 = int((datetime.datetime(end_time.year,3,1)-datetime.datetime(end_time.year,2,1)).total_seconds()/21600+492-(datetime.datetime(1980,1,1)-datetime.datetime(1979,end_time.month,end_time.day,end_time.hour,0,0)).total_seconds()/21600)
        else:
            cindex2 = int((datetime.datetime(1980,end_time.month,end_time.day,end_time.hour,0,0)-datetime.datetime(1980,1,1)).total_seconds()/21600)
            index2 = cindex2

        slp_ds = xr.open_dataset('/storage/yeechian/data/slp_NH_ERA5_%d.nc'%syear)
        thck_ds = xr.open_dataset('/storage/yeechian/data/z_NH_ERA5_%d.nc'%syear)
        
        if cindex2 < cindex: #EPR spans two calendar years
            thck_climo = np.mean(np.concatenate((np.asarray(climo_ds['thck_climo'].metpy.sel(latitude=slice(max_lat,min_lat+0.1),longitude=slice(min_lon,max_lon+0.1))[cindex:,:,:]),np.asarray(climo_ds['thck_climo'].metpy.sel(latitude=slice(max_lat,min_lat+0.1),longitude=slice(min_lon,max_lon+0.1))[:cindex2,:,:])),axis=0),axis=0)
            slp_climo = np.mean(np.concatenate((np.asarray(climo_ds['slp_climo'].metpy.sel(latitude=slice(max_lat,min_lat+0.1),longitude=slice(min_lon,max_lon+0.1))[cindex:,:,:]),np.asarray(climo_ds['slp_climo'].metpy.sel(latitude=slice(max_lat,min_lat+0.1),longitude=slice(min_lon,max_lon+0.1))[:cindex2,:,:])),axis=0),axis=0)
            slp_ds2 = xr.open_dataset('/storage/yeechian/data/slp_NH_ERA5_%d.nc'%eyear)
            thck_ds2 = xr.open_dataset('/storage/yeechian/data/z_NH_ERA5_%d.nc'%eyear)

            thck_ano[i,:,:] = np.mean(np.concatenate(((np.asarray(thck_ds['z'].metpy.sel(level=500,latitude=slice(max_lat,min_lat+0.1),longitude=slice(min_lon,max_lon+0.1))[index:,:,:]-thck_ds['z'].metpy.sel(level=1000,latitude=slice(max_lat,min_lat+0.1),longitude=slice(min_lon,max_lon+0.1))[index:,:,:])/9.81-thck_climo),(np.asarray(thck_ds2['z'].metpy.sel(level=500,latitude=slice(max_lat,min_lat+0.1),longitude=slice(min_lon,max_lon+0.1))[:index2,:,:]-thck_ds2['z'].metpy.sel(level=1000,latitude=slice(max_lat,min_lat+0.1),longitude=slice(min_lon,max_lon+0.1))[:index2,:,:])/9.81-thck_climo)),axis=0)*0.1,axis=0)
            slp_ano[i,:,:] = np.mean(np.concatenate(((np.asarray(slp_ds['msl'].metpy.sel(latitude=slice(max_lat,min_lat+0.1),longitude=slice(min_lon,max_lon+0.1))[index:,:,:])-slp_climo),(np.asarray(slp_ds2['msl'].metpy.sel(latitude=slice(max_lat,min_lat+0.1),longitude=slice(min_lon,max_lon+0.1))[:index2,:,:])-slp_climo)),axis=0)*0.01,axis=0)
            thck_ano_vec[i,:] = np.ravel(thck_ano[i,::2,::2]) #To reduce computation time, SOM is run with 1 x 1 degree instead of 0.5 x 0.5 degree fields, which suffices for synoptic-scale features
            slp_ano_vec[i,:] = np.ravel(slp_ano[i,::2,::2])
        else:   
            thck_climo = np.mean(np.asarray(climo_ds['thck_climo'][cindex:cindex2,60:131,2:261]),axis=0)
            slp_climo = np.mean(np.asarray(climo_ds['slp_climo'][cindex:cindex2,60:131,2:261]) ,axis=0)    
            thck_ano[i,:,:] = np.mean((np.asarray(thck_ds['z'].metpy.sel(level=500,latitude=slice(max_lat,min_lat+0.1),longitude=slice(min_lon,max_lon+0.1))[index:index2,:,:]-thck_ds['z'].metpy.sel(level=1000,latitude=slice(max_lat,min_lat+0.1),longitude=slice(min_lon,max_lon+0.1))[index:index2,:,:])/9.81-thck_climo)*0.1,axis=0)
            slp_ano[i,:,:] = np.mean((np.asarray(slp_ds['msl'].metpy.sel(latitude=slice(max_lat,min_lat+0.1),longitude=slice(min_lon,max_lon+0.1))[index:index2,:,:])-slp_climo)*0.01,axis=0)
            thck_ano_vec[i,:] = np.ravel(thck_ano[i,::2,::2])
            slp_ano_vec[i,:] = np.ravel(slp_ano[i,::2,::2])

    
    tfiles = [[],[],[]]
    sfiles = [[],[],[]]
    
    #3x3 SOM
    for n_neurons in range(3,4): #n_neurons = 3
        for m_neurons in range(3,4): #m_neurons = 3
            
            #Training
            thck_som = MiniSom(n_neurons, m_neurons, thck_ano_vec.shape[1], sigma=1.5, learning_rate=.5, neighborhood_function='gaussian', random_seed=0)
            thck_som.pca_weights_init(thck_ano_vec)
            thck_som.train(thck_ano_vec, 1000, verbose=True) 
            thck_win_coords = []
            
            slp_som = MiniSom(n_neurons, m_neurons, slp_ano_vec.shape[1], sigma=1.5, learning_rate=.5, neighborhood_function='gaussian', random_seed=0)
            slp_som.pca_weights_init(slp_ano_vec)
            slp_som.train(slp_ano_vec, 1000, verbose=True)
            slp_win_coords = []

            #Initialize arrays to store and plot composites of different SOM nodes
            thck_som_comp = np.zeros((n_neurons, m_neurons, size[1], size[2]))
            thck_som_slp_comp = np.zeros((n_neurons, m_neurons, size[1], size[2]))
            slp_som_comp = np.zeros((n_neurons, m_neurons, size[1], size[2]))
            slp_som_thck_comp = np.zeros((n_neurons, m_neurons, size[1], size[2]))
            thck_num_cases = np.zeros((n_neurons, m_neurons))
            slp_num_cases = np.zeros((n_neurons, m_neurons))

            #Put cases "belonging" to each node in csv files
            for i in range(n_neurons):
                for j in range(m_neurons):
                    tfiles[i].append(open("/storage/yeechian/data/data_csv/EPR_%s_NH_thck_%d%d.csv"%(t_str2[met],i,j),"w"))
                    sfiles[i].append(open("/storage/yeechian/data/data_csv/EPR_%s_NH_slp_%d%d.csv"%(t_str2[met],i,j),"w"))
                    
            #For each case        
            for i in range(size[0]):
                line = lines[i]

                #Add case to cases's "winner node" 1000-500 hPa thickness composite
                win = thck_som.winner(thck_ano_vec[i])
                thck_som_comp[win[0],win[1],:,:] += thck_ano[i,:,:]
                thck_som_slp_comp[win[0],win[1],:,:] += slp_ano[i,:,:]
                thck_num_cases[win[0],win[1]] += 1.0
                thck_win_coords.append(win)
                tfiles[win[0]][win[1]].write(line)
                
                #Add case to cases's "winner node" SLP composite
                win = slp_som.winner(slp_ano_vec[i]) #winner node
                slp_som_comp[win[0],win[1],:,:] += slp_ano[i,:,:]
                slp_som_thck_comp[win[0],win[1],:,:] += thck_ano[i,:,:]
                slp_num_cases[win[0],win[1]] += 1.0
                slp_win_coords.append(win)
                sfiles[win[0]][win[1]].write(line)

            for i in range(n_neurons):
                for j in range(m_neurons):
                    thck_som_comp[i,j,:,:] = thck_som_comp[i,j,:,:]/thck_num_cases[i,j]
                    thck_som_slp_comp[i,j,:,:] = thck_som_slp_comp[i,j,:,:]/thck_num_cases[i,j]
                    slp_som_comp[i,j,:,:] = slp_som_comp[i,j,:,:]/slp_num_cases[i,j]
                    slp_som_thck_comp[i,j,:,:] = slp_som_thck_comp[i,j,:,:]/slp_num_cases[i,j]
                    
            thck_win_coords = np.asarray(thck_win_coords)
            slp_win_coords = np.asarray(slp_win_coords)

            for i in range(n_neurons):
                for j in range(m_neurons):
                        #Plot composite using 1000-500 hPa SOM
                        fig = plt.figure(figsize=(17., 6.1))
                        fig.tight_layout(rect=[0,0.01,1,0.99])
                        ax = plt.subplot(111, projection=ccrs.PlateCarree())
                        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1.5, color='black', alpha=0.5, linestyle='--')
                        gl.xlocator = mticker.FixedLocator(np.arange(-180,-40,10))
                        gl.ylocator = mticker.FixedLocator(np.arange(20,80,10))
                        gl.xformatter = LONGITUDE_FORMATTER
                        gl.yformatter = LATITUDE_FORMATTER
                        gl.xlabel_style = {'size': 18,'weight': 'bold'}
                        gl.ylabel_style = {'size': 18,'weight': 'bold'}
                        gl.xlabels_top = False
                        gl.xlabels_bottom = True
                        gl.ylabels_left = True
                        gl.ylabels_right = False
                        plt.title("%s EPR 1000-500 hPa thickness anomaly %dx%d SOM, (%d, %d) node, n = %d \n 1000-500 hPa thickness anomaly (dam) and SLP anomaly (hPa)" % (t_str[met],n_neurons,m_neurons,i,j,int(thck_num_cases[i][j])),{'fontsize':18, 'fontweight':"bold"})

                        # Set extent and plot map lines
                        ax.set_extent([min_lon,max_lon,min_lat,max_lat], ccrs.PlateCarree())
                        ax.add_feature(land, edgecolor=(0.2,0.2,0.2), linewidth=2.0)
                        ax.add_feature(states_provinces, edgecolor=(0.2,0.2,0.2), linewidth=1.0)

                        cs2 = ax.contour(lon_360, lat, thck_som_slp_comp[i,j,:,:], clevs_slp, linewidths=2.5, linestyles=['dashed','dashed','dashed','dashed','dashed','dashed','dashed','dashed','dashed','dashed','dashed','solid','solid','solid','solid','solid','solid','solid','solid','solid','solid','solid'], colors='black', transform=ccrs.PlateCarree())
                        plt.clabel(cs2, cs2.levels,fontsize=18, inline=1, inline_spacing=10, fmt='%d',rightside_up=True, use_clabeltext=False)
                        cf1 = ax.contourf(lon_360, lat,thck_som_comp[i,j,:,:], clevs_thck, extend='both',cmap='bwr', transform=ccrs.PlateCarree())
                        cax1 = plt.subplot(gs[1])
                        cb = plt.colorbar(cf1, cax=cax1, orientation='horizontal', pad=0.2, aspect=100,fraction=0.1, extendrect='True', ticks=clevs_thck)
                        cb.set_label('dam', size=18)
                        
                        plt.savefig('/storage/yeechian/figures/SOM_EPR_%s_NH_%dx%d_%d%d_thck_ano.png' % (t_str2[met],n_neurons,m_neurons,i,j))
                        plt.close()
            
                        #Plot composite using SLP SOM
                        fig = plt.figure(figsize=(17., 6.1))
                        fig.tight_layout(rect=[0,0.01,1,0.99])
                        ax = plt.subplot(111, projection=ccrs.PlateCarree())
                        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1.5, color='black', alpha=0.5, linestyle='--')
                        gl.xlocator = mticker.FixedLocator(np.arange(-180,-40,10))
                        gl.ylocator = mticker.FixedLocator(np.arange(20,80,10))
                        gl.xformatter = LONGITUDE_FORMATTER
                        gl.yformatter = LATITUDE_FORMATTER
                        gl.xlabel_style = {'size': 18,'weight': 'bold'}
                        gl.ylabel_style = {'size': 18,'weight': 'bold'}
                        gl.xlabels_top = False
                        gl.xlabels_bottom = True
                        gl.ylabels_left = True
                        gl.ylabels_right = False
                        plt.title("%s EPR SLP anomaly %dx%d SOM, (%d, %d) node, n = %d \n 1000-500 hPa thickness anomaly (dam) and SLP anomaly (hPa)" % (t_str[met],n_neurons,m_neurons,i,j,int(slp_num_cases[i][j])),{'fontsize':18, 'fontweight':"bold"})

                        # Set extent and plot map lines
                        ax.set_extent([min_lon,max_lon,min_lat,max_lat], ccrs.PlateCarree())
                        ax.add_feature(land, edgecolor=(0.2,0.2,0.2), linewidth=2.0)
                        ax.add_feature(states_provinces, edgecolor=(0.2,0.2,0.2), linewidth=1.0)

                        cs2 = ax.contour(lon_360, lat, slp_som_comp[i,j,:,:], clevs_slp, linewidths=2.5, linestyles=['dashed','dashed','dashed','dashed','dashed','dashed','dashed','dashed','dashed','dashed','dashed','solid','solid','solid','solid','solid','solid','solid','solid','solid','solid','solid'], colors='black', transform=ccrs.PlateCarree())
                        plt.clabel(cs2, cs2.levels,fontsize=18, inline=1, inline_spacing=10, fmt='%d',rightside_up=True, use_clabeltext=False)
                        cf1 = ax.contourf(lon_360, lat,slp_som_thck_comp[i,j,:,:], clevs_thck, extend='both',cmap='bwr', transform=ccrs.PlateCarree())
                        cax1 = plt.subplot(gs[1])
                        cb = plt.colorbar(cf1, cax=cax1, orientation='horizontal', pad=0.2, aspect=100,fraction=0.1, extendrect='True', ticks=clevs_thck)
                        cb.set_label('dam', size=18)
                        
                        plt.savefig('/storage/yeechian/figures/SOM_EPR_%s_NH_%dx%d_%d%d_slp_ano.png' % (t_str2[met],n_neurons,m_neurons,i,j))
                        plt.close()
           
