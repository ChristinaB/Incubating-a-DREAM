"""
Script to plot and compare modeled historical climate data and observed and 
modeled streamflow

Created starting Jan 13 2017 by cbev
"""
#%%
# Import required modules
import os
import csv
import pandas as pd
import datetime
import matplotlib.pyplot as plt
import numpy as np
from collections import OrderedDict
import ulmo
import geopandas as gpd
from shapely.geometry import Point
import xarray as xr

# IMPORT DATA
# INPUT: Location Name and drainage area (m2)
loc_name='Elwha Watershed'
drainage_area_all=832944265 #m2
drainage_area_mcd=694527755 #m2
drainage_area_mills=512954011 #m2

# General Information
homedir='D:/GitHub/Elwha_Landlab/dhsvm' 
forcsdir='D:/GitHub/Elwha_Landlab/dhsvm/forcings' 
obsdir='D:/GitHub/Elwha_Landlab/observations' 
liv2013dir='D:/GitHub/Elwha_Landlab/dhsvm/forcings/raw/livneh2013'
liv2015dir='D:/GitHub/Elwha_Landlab/dhsvm/forcings/raw/livneh2015'
wrf2014dir='D:/GitHub/Elwha_Landlab/dhsvm/forcings/raw/wrf2014'
mappingfile='ElwhaClimatePoints_UTM_Elev.csv' # csv file with lat, long and elevation. Header should include 'LAT','LONG_', and 'ELEV'
forcingfile='met_forcing_gis.xlsx'
streamflow_mills_file='12044900_ElwhaMills.xls'
streamflow_mcd_file='12045500_ElwhaMcdonald.xls'

# Observations information
obs_names=['buckinghorse','waterhole','elwha_rs','sutherland','port_angeles'] # station names
obs_title=['Buckinghorse','Waterhole','Ranger Station','Sutherland','Port Angeles'] # names for plots
obs_elevation=[1527, 1484, 110, 174, 27] # station elevations
obs_met_station=[7, 26, 43, 51, 47] # met stations corresponding to the observation statations
obs_start_date=[datetime.datetime(2008, 9,10),datetime.datetime(1999,10,1),
                datetime.datetime(1942,7,1),datetime.datetime(1929,1,1),
                datetime.datetime(1917,1,1)]
obs_end_date=[datetime.datetime(2010,12,31), datetime.datetime(2010,12,31),
              datetime.datetime(2010,12,31),datetime.datetime(1976,6,1),
              datetime.datetime(2008,3,1)]

# General reference variables
# For computations
forcing_list=['liv2013','liv2015','liv2015bc','wrf2014','wrf2014bc']
variables=['precip','tmax', 'tmin', 'wind']
# For plots
wy_index=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
wy_numbers=[9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8]
month_strings=[ 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sept']
days_per_month=np.array([31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])

# Comparison Decisions
# "_c" is for comparison
compare_criteria='mills' # 'elevation' or 'mills' or 'mcdonald'
plot_title='Lake Mills Tributary'

# Elevation range- okay to be beyond the actual elevations
min_elev_c=1500 # default is min elevation: met_df.elevation.min()
max_elev_c=2500 # default is max elevation: met_df.elevation.max()

# Date range
start_date_c=datetime.date(1950,1,1)
end_date_c=datetime.date(2010,12,31)
# 
import_obs='no' # if yes, import observations


#%% Import observation stations data (COOP and SNOTEL)
#st = ulmo.ncdc.ghcn_daily.get_stations(country='US', as_dataframe=True)
#st[st.name.str.contains('ELWHA')]

if import_obs=='yes':
    # SNOTEL sites
    buckinghorse = ulmo.ncdc.ghcn_daily.get_data('USS0023B18S', as_dataframe=True)
    waterhole = ulmo.ncdc.ghcn_daily.get_data('USS0023B17S', as_dataframe=True)
    
    # COOP stations used by Livneh
    elwha_rs = ulmo.ncdc.ghcn_daily.get_data('USC00452548', as_dataframe=True)
    sutherland = ulmo.ncdc.ghcn_daily.get_data('USC00454438', as_dataframe=True)
    port_angeles = ulmo.ncdc.ghcn_daily.get_data('USC00456624', as_dataframe=True)
    
    # COOP Stations that have incomplete records and are in Livneh data file but probably not used by livneh
    #port_crescent = ulmo.ncdc.ghcn_daily.get_data('USC00456642', as_dataframe=True)
    #mount_pleasant = ulmo.ncdc.ghcn_daily.get_data('USC00455669', as_dataframe=True)
    #port_angeles_wbap = ulmo.ncdc.ghcn_daily.get_data('USW00024228', as_dataframe=True)
    
    obs_data={'elevation':obs_elevation, 'met_station':obs_met_station,
              'title':obs_title, 'start_date':obs_start_date, 'end_date':obs_end_date}
    obs_stations_df=pd.DataFrame(obs_data, index=obs_names)

#%% Import gridded climate products
os.chdir(forcsdir)

# Read in reference locations and elevations:
def read_in_longlats(mappingfile):
    maptable=[]
    with open(mappingfile, 'r') as csvfile: 
        longlat = csv.reader(csvfile, delimiter=',')
        for row in longlat:
            maptable.append(row)
    csvfile.close()
    return(maptable)
    
maptable = read_in_longlats(mappingfile) # maptable is a list of the rows/columns from the mapping file. The top row is the header of each column. 

station=maptable[0].index('FID')
latitude=maptable[0].index('LAT')
longitude=maptable[0].index('LONG')
elevation=maptable[0].index('ELEV')

station_list=[]; lat_list=[]; long_list=[]; elev_list=[]

for row in maptable:
    if maptable.index(row)!=0:
        station_list.append(row[station])
        lat_list.append(row[latitude])
        long_list.append(row[longitude])
        elev_list.append(int(row[elevation]))

# Create a data frame and geodataframe with station number, lat, long, and elevation values
met_df=pd.DataFrame({"station": station_list,"latitude": lat_list,
                     "longitude": long_list,"elevation": elev_list})
met_geometry=[Point(xy) for xy in zip(pd.to_numeric(met_df['longitude']),
              pd.to_numeric(met_df['latitude']))]
met_gdf=gpd.GeoDataFrame(met_df,geometry=met_geometry)
n_stations=len(met_df) # number of stations

ws_bnd=gpd.read_file('D:/GitHub/Elwha_Landlab/dhsvm/forcings/GIS/elwha_ws_bnd_wgs84.shp')
met_grid_poly=gpd.read_file('D:/GitHub/Elwha_Landlab/dhsvm/forcings/GIS/climate_station_poly_wgs84.shp')

# Plot watershed boundary and climate points
f, ax = plt.subplots(1, figsize=(8, 6))
ax.set_title('Gridded Meteorology Points: Elwha Watershed')
ws_bnd.plot(ax=ax, cmap='Set3')
met_gdf.plot(ax=ax, marker='^', color='black', edgecolor='white')
met_grid_poly.plot(ax=ax, color='none', edgecolor='black')
ax.set_xlim([-123.7, -123.1])
ax.set_ylim([47.6, 48.2])
plt.axis('equal')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
#%%                                                                           
# Upload Livneh 2013 meteorology data                                      
# import the data
# INPUT directory to 2013 met data
os.chdir(forcsdir+'/raw/livneh2013')                                   

# compile name of Livneh 2013 met data files
file_names=[]
for i in range(0, n_stations):
    file_names.append('_'.join(['Meteorology_Livneh_CONUSExt_v.1.2_2013',
                                met_df.latitude[i], 
                                met_df.longitude[i]]))
    # This is a list of file names of climate data, 1 file per location

# create date range for livneh 2013 daily met data since not included in the data file
liv2013_start_date=datetime.date(1915,1,1)
liv2013_end_date=datetime.date(2011,12,31)
liv2013_dates=pd.date_range(liv2013_start_date,liv2013_end_date)

# compile data from text files into a list of data frames where each row 
# corresponds to a location. The index of the row is the index in the 
# 'met_df' (lat, long, elev)

liv2013=[]
for i in range(0, n_stations):
    liv2013.append(pd.read_table(file_names[i], header=None))
    liv2013[i].columns=['precip','tmax', 'tmin', 'wind']
    liv2013[i].set_index(liv2013_dates, inplace=True) # set the index as the date

# Upload Livneh 2015 BC meteorology data                                      
os.chdir(forcsdir+'/raw/livneh2015_BC')                                      
file_names=[]
for i in range(0, n_stations):
    file_names.append('_'.join(['data', met_df.latitude[i], 
                                met_df.longitude[i]]))
liv2015bc_start_date=datetime.date(1950,1,1)
liv2015bc_end_date=datetime.date(2013,12,31)
liv2015bc_dates=pd.date_range(liv2015bc_start_date,liv2015bc_end_date)
liv2015bc=[]
for i in range(0, n_stations):
    liv2015bc.append(pd.read_table(file_names[i], delimiter=' ',header=None))
    liv2015bc[i].columns=['precip','tmax', 'tmin', 'wind']
    liv2015bc[i].set_index(liv2015bc_dates, inplace=True) # set the index as the date

# Upload Livneh 2015 meteorology data   
os.chdir(forcsdir+'/raw/livneh2015')    
file_names=[]
for i in range(0, n_stations):
    file_names.append('_'.join(['Meteorology_Livneh_NAmerExt_15Oct2014',
                                met_df.latitude[i], 
                                met_df.longitude[i]]))
liv2015_start_date=datetime.date(1950,1,1)
liv2015_end_date=datetime.date(2013,12,31)
liv2015_dates=pd.date_range(liv2015_start_date,liv2015_end_date)
liv2015=[]
for i in range(0, n_stations):
    liv2015.append(pd.read_table(file_names[i], delimiter=' ', header=None))
    liv2015[i].columns=['precip','tmax', 'tmin', 'wind']
    liv2015[i].set_index(liv2015_dates, inplace=True) # set the index as the date
    
# Upload Salathe 2014 meteorology data   
os.chdir(forcsdir+'/raw/wrf2014')   
file_names=[]
for i in range(0, n_stations):
    file_names.append('_'.join(['data',
                                met_df.latitude[i], 
                                met_df.longitude[i]]))
wrf2014_start_date=datetime.date(1950,1,1)
wrf2014_end_date=datetime.date(2010,12,31)
wrf2014_dates=pd.date_range(wrf2014_start_date,wrf2014_end_date)
wrf2014=[]
for i in range(0, n_stations):
    wrf2014.append(pd.read_table(file_names[i], delimiter=' ', header=None))
    wrf2014[i].columns=['precip','tmax', 'tmin', 'wind']
    wrf2014[i].set_index(wrf2014_dates, inplace=True) # set the index as the date

# Upload Salathe 2014 BC meteorology data
os.chdir(forcsdir+'/raw/wrf2014_BC') 
file_names=[]
for i in range(0, n_stations):
    file_names.append('_'.join(['data',
                                met_df.latitude[i], 
                                met_df.longitude[i]]))
    
wrf2014bc_dates=pd.date_range(wrf2014_start_date,wrf2014_end_date)
wrf2014bc=[]
for i in range(0, len(file_names)):
    wrf2014bc.append(pd.read_table(file_names[i], delimiter=' ', header=None))
    wrf2014bc[i].columns=['precip','tmax', 'tmin', 'wind']
    wrf2014bc[i].set_index(wrf2014_dates, inplace=True) # set the index as the date

#%%
# Upload prism
os.chdir(forcsdir+'/raw/prism/') 
file_names=[]
for i in range(0, n_stations):
    file_names.append('.'.join([('_'.join(['data',met_df.latitude[i],met_df.longitude[i]])),'prism']))
   
prism=[]
for i in range(0, len(file_names)):
    prism.append(pd.read_table(file_names[i], delimiter=' ', header=None))
    prism[i].columns=[1,2,3,4,5,6,7,8,9,10,11,12]
    
prism_ds=xr.concat([df.to_xarray() for df in prism], dim='dataset')   
prism_ds.rename({'dataset':'station'}, inplace=True)

#%%
# Upload streamflow data and create data frames
def create_q_obs_df(file_name, drainage_area):
    q=pd.read_excel(file_name, sheetname='data')
    q.columns=['year','month','day','flow_cfs']
    q_dates=pd.to_datetime(q[[0, 1, 2]])
    q.set_index(q_dates, inplace=True)
    q.drop(['year','month','day'],axis=1, inplace=True)
    q_cms=q.flow_cfs/(3.28084**3)
    q_mmday=q_cms*1000*3600*24/drainage_area
    q=pd.concat([q, q_cms, q_mmday],axis=1)
    q.columns=['flow_cms','flow_cfs', 'flow_mmday']
    return q

os.chdir(obsdir)
qobs_mills=create_q_obs_df(streamflow_mills_file, drainage_area_mills)
qobs_mcd=create_q_obs_df(streamflow_mcd_file, drainage_area_mcd)

#%%
# Upload forcing cell geometry and create station weight datasets
os.chdir(forcsdir)
station_geom_all= pd.read_excel(forcingfile, sheetname='all')
station_geom_all.set_index('Cell_ID', inplace=True)
station_geom_mills=pd.read_excel(forcingfile, sheetname='lakemills')
station_geom_mills.set_index('Cell_ID', inplace=True)
station_geom_mcd=pd.read_excel(forcingfile, sheetname='mcdonald')
station_geom_mcd.set_index('Cell_ID', inplace=True)

station_area_weights_mills_ds=xr.DataArray(station_geom_mills.Area_m2/drainage_area_mills, 
                                           coords=[station_geom_mills.index], dims=['station'])

station_area_weights_mcd_ds=xr.DataArray(station_geom_mcd.Area_m2/drainage_area_mcd, 
                                           coords=[station_geom_mcd.index], dims=['station'])

#%% Develop datasets
# Create xarray of station data
stations_ds=xr.Dataset(coords={'longitude':pd.to_numeric(met_df['longitude']).values,
                               'latitude':pd.to_numeric(met_df['latitude']).values,
                               'elevation':pd.to_numeric(met_df['elevation']).values})
    
# Create data sets and data frames of forcing data
i=0   
for name in forcing_list:
    globals()[name+'_ds']=xr.concat([df.to_xarray() for df in eval(str(forcing_list[i]))], dim="dataset")
    globals()[name+'_ds'].rename({'dataset':'station', 'index':'time'}, inplace=True)
    globals()[name+'_ds']=eval(name+'_ds').merge(stations_ds.coords)
    print(name+'_ds')
    i=i+1

#%%
if compare_criteria=='elevation':    
    # Determine station selection per desired elevation range
    stations_info_c=met_df[(met_df.elevation >= min_elev_c) & (met_df.elevation <= max_elev_c)] # Find points in elevation range of interest
    elev_min_station_c=stations_info_c.elevation.idxmin() # station of minimum elevation in analysis
    elev_min_c=stations_info_c.elevation.min() # minimum elevaiton of stations in analysis
    elev_max_station_c=stations_info_c.elevation.idxmax() # station of maximum elevation in analysis
    elev_max_c=stations_info_c.elevation.max() # maximum elevaiton of stations in analysis
    # Extract list of station numbers for indexing. Alternative, you can set the list of stations manually!
    stations_list_c=stations_info_c.station.values.astype('int64') # default is stations in desired elevation band: stations_info_c.station.values.astype('int64')
    
    # Compute monthly values in time frame and elevation band of interest
    for name in forcing_list:
        for var in variables:
            globals()[var+'_'+name+'_ds_c']=eval(name+'_ds')[var].sel(time=slice(start_date_c, end_date_c), station=stations_list_c)
            print(var+'_'+name+'_ds_c')


if compare_criteria=='mills':    
    # Compute monthly values in time frame and elevation band of interest
    elev_min_station_c=station_geom_mills.Elev_.idxmin() # station of minimum elevation in analysis
    elev_min_c=station_geom_mills.Elev_.min() # minimum elevaiton of stations in analysis
    elev_max_station_c=station_geom_mills.Elev_.idxmax() # station of maximum elevation in analysis
    elev_max_c=station_geom_mills.Elev_.max() # maximum elevaiton of stations in analysis

    # Extract list of station numbers for indexing. Alternative, you can set the list of stations manually!
    stations_list_c=station_geom_mills.index.values.astype('int64') # default is stations in desired elevation band: stations_info_c.station.values.astype('int64')
    for name in forcing_list:
        for var in variables:
            globals()[var+'_'+name+'_ds_c']=eval(name+'_ds')[var].sel(time=slice(start_date_c, end_date_c), station=station_geom_mills.index)
            print(var+'_'+name+'_ds_c')       

if compare_criteria=='mcdonald':    
    # Compute monthly values in time frame and elevation band of interest
    elev_min_station_c=station_geom_mcd.Elev_.idxmin() # station of minimum elevation in analysis
    elev_min_c=station_geom_mcd.Elev_.min() # minimum elevaiton of stations in analysis
    elev_max_station_c=station_geom_mcd.Elev_.idxmax() # station of maximum elevation in analysis
    elev_max_c=station_geom_mcd.Elev_.max() # maximum elevaiton of stations in analysis

    # Extract list of station numbers for indexing. Alternative, you can set the list of stations manually!
    stations_list_c=station_geom_mcd.index.values.astype('int64') # default is stations in desired elevation band: stations_info_c.station.values.astype('int64')
    for name in forcing_list:
        for var in variables:
            globals()[var+'_'+name+'_ds_c']=eval(name+'_ds')[var].sel(time=slice(start_date_c, end_date_c), station=station_geom_mcd.index)
            print(var+'_'+name+'_ds_c')    

#%% PLOTS
# Monthly Tmax
fig, ax=plt.subplots(1,1,figsize=(6, 3))
plt.plot(wy_index, tmax_liv2013_ds_c.groupby('time.month').mean('time').mean('station').values[wy_numbers],'b-',linewidth=2, label='L13, L15')
plt.plot(wy_index, tmax_liv2015bc_ds_c.groupby('time.month').mean('time').mean('station').values[wy_numbers],'c-',linewidth=2, label='L15 BC')
plt.plot(wy_index, tmax_wrf2014_ds_c.groupby('time.month').mean('time').mean('station').values[wy_numbers],'r-', linewidth=2, label='WRF 16')
#plt.plot(wy_index, tmax_wrf2014bc_ds_c.groupby('time.month').mean('time').mean('station').values[wy_numbers],'b--', linewidth=2, label='WRFbc')
plt.xlabel ('Month',fontsize=14)
plt.ylabel('Temperature (deg C)',fontsize=14)
#plt.title(str(loc_name)+', Years:'+str(start_date_c.year)+'-'+str(end_date_c.year)+'\n All Met Stations (Elevation: '+str(min_elev_c)+'m -'+str(max_elev_c)+'m)', fontsize=16)
plt.title('Average Monthly Maximum Temperature\n'+plot_title, fontsize=16)
plt.legend()
plt.tick_params(labelsize=12)
plt.grid(which='both')
plt.xlim(1,12);
plt.xticks(wy_index, month_strings);

# Monthly Tmin
fig, ax=plt.subplots(1,1,figsize=(6, 3))
plt.plot(wy_index,tmin_liv2013_ds_c.groupby('time.month').mean('time').mean('station').values[wy_numbers],'b-',linewidth=2, label='L13, L15')
plt.plot(wy_index,tmin_liv2015bc_ds_c.groupby('time.month').mean('time').mean('station').values[wy_numbers],'c-',linewidth=2, label='L15 BC')
plt.plot(wy_index,tmin_wrf2014_ds_c.groupby('time.month').mean('time').mean('station').values[wy_numbers],'r-', linewidth=2, label='WRF 16')
#plt.plot(wy_index,tmin_wrf2014bc_ds_c.groupby('time.month').mean('time').mean('station').values[wy_numbers],'b--', linewidth=2, label='WRFbc')
plt.xlabel ('Month',fontsize=14)
plt.ylabel('Temperature (deg C)',fontsize=14)
plt.title('Average Monthly Minimum Temperature\n'+plot_title, fontsize=16)
#plt.title(str(loc_name)+'\nAverage Monthly Min Temperature\n Years: '+str(start_date_c.year)+'-'+str(end_date_c.year)+'; Elevation: '+str(min_elev_c)+'m -'+str(max_elev_c)+'m', fontsize=16)
plt.legend()
plt.tick_params(labelsize=12)
plt.grid(which='both')
plt.xlim(1,12);
plt.xticks(wy_index, month_strings);

# Monthly Tavg
#fig, ax=plt.subplots(1,1,figsize=(6, 4))
#plt.plot(wy_index, (tmax_liv2013_ds_c.groupby('time.month').mean('time').mean('station').values[wy_numbers]+ tmin_liv2013_ds_c.groupby('time.month').mean('time').mean('station').values[wy_numbers])/2,'r-',linewidth=2, label='Liv')
#plt.plot(wy_index,(tmax_liv2015bc_ds_c.groupby('time.month').mean('time').mean('station').values[wy_numbers]+ tmin_liv2015bc_ds_c.groupby('time.month').mean('time').mean('station').values[wy_numbers])/2,'r--',linewidth=2, label='Livbc')
#plt.plot(wy_index,(tmax_wrf2014_ds_c.groupby('time.month').mean('time').mean('station').values[wy_numbers]+ tmin_wrf2014_ds_c.groupby('time.month').mean('time').mean('station').values[wy_numbers])/2,'b-', linewidth=2, label='WRF')
#plt.plot(wy_index, (tmax_wrf2014bc_ds_c.groupby('time.month').mean('time').mean('station').values[wy_numbers]+ tmin_wrf2014bc_ds_c.groupby('time.month').mean('time').mean('station').values[wy_numbers])/2,'b--', linewidth=2, label='WRFbc')
#plt.xlabel ('Month',fontsize=14)
#plt.ylabel('Temperature (deg C)',fontsize=14)
#plt.title(str(loc_name)+'\nAverage Monthly Avg Temperature\n Years: '+str(start_date_c.year)+'-'+str(end_date_c.year)+'; Elevation: '+str(min_elev_c)+'m -'+str(max_elev_c)+'m', fontsize=16)
#plt.legend()
#plt.tick_params(labelsize=12)
#plt.grid(which='both')
#plt.xlim(1,12);
#plt.xticks(wy_index, month_strings);

# Monthly Precipitation
fig, ax=plt.subplots(1,1,figsize=(6,3))
plt.plot(wy_index, days_per_month*precip_liv2013_ds_c.groupby('time.month').mean('time').mean('station').values[wy_numbers],'b-',linewidth=2, label='L13')
plt.plot(wy_index, days_per_month*precip_liv2015_ds_c.groupby('time.month').mean('time').mean('station').values[wy_numbers],'g-',linewidth=2, label='L15')
plt.plot(wy_index, days_per_month*precip_liv2015bc_ds_c.groupby('time.month').mean('time').mean('station').values[wy_numbers],'c-',linewidth=2, label='L15 BC')
plt.plot(wy_index, days_per_month*precip_wrf2014_ds_c.groupby('time.month').mean('time').mean('station').values[wy_numbers],'r-',linewidth=2, label='WRF 16')
#plt.plot(wy_index, days_per_month*precip_wrf2014bc_ds_c.groupby('time.month').mean('time').mean('station').values[wy_numbers],'b--',linewidth=2, label='WRF 2014bc')
#plt.plot([1, 12],[0, 0], 'k-',linewidth=2)
plt.xlabel ('Month',fontsize=14)
plt.ylabel('Precipitation (mm)',fontsize=14)
plt.title('Average Monthly Precipitation\n'+plot_title, fontsize=16)
#plt.title(str(loc_name)+':\nAverage Monthly Precipitation\n Years: '+str(start_date_c.year)+'-'+str(end_date_c.year)+'; Elevation: '+str(min_elev_c)+'m -'+str(max_elev_c)+'m', fontsize=16)
plt.legend(loc='best')
plt.tick_params(labelsize=12)
plt.grid(which='both')
plt.xlim(1,12);
plt.xticks(wy_index, month_strings);

#%% COMPARE GRIDDED DATASETS TO OBSERVATION STATIONS
# Elwha RS elevation=110 m,  Cell 43 elevation= 468 m
### START INPUT
obs_variable='waterhole' # elwha_rs, sutherland, port_angeles, buckinghorse, waterhole

# Extract observation station information
obs_elev=obs_stations_df.elevation[obs_variable]
obs_title=obs_stations_df.title[obs_variable]
start_date_c=obs_stations_df.start_date[obs_variable]
end_date_c=obs_stations_df.end_date[obs_variable]
station=obs_stations_df.met_station[obs_variable]
    
# Compute lapse rate
station_elev=stations_info_c.elevation[station]
if station_elev>obs_elev:
    temp_lapse=6.5*(station_elev-obs_elev)/1000
else: temp_lapse=-6.5*(station_elev-obs_elev)/1000

# Extract relevant obs info
obs_start_year=str(start_date_c.year)
obs_end_year=str(end_date_c.year)
obs_precip=eval(obs_variable)['PRCP'][obs_start_year:obs_end_year].value.astype('float')/10 # mm
if obs_variable is not 'sutherland':
    obs_tmax=eval(obs_variable)['TMAX'][obs_start_year:obs_end_year].value.astype('float')/10 # deg C
    obs_tmin=eval(obs_variable)['TMIN'][obs_start_year:obs_end_year].value.astype('float')/10 # deg C

# Extract relevant met info
station_elev=stations_info_c.loc[station,'elevation']
wy_numbers_pd=[x+1 for x in wy_numbers]
for name in forcing_list:
    for var in variables:
        globals()[var+'_'+name+'_ds_c']=eval(name+'_ds')[var].sel(time=slice(start_date_c, end_date_c), station=station)
        print(var+'_'+name+'_ds_c')

# Tmax
if obs_variable is not 'sutherland':
    fig, ax=plt.subplots(1,1,figsize=(6,4))
    plt.plot(wy_index, tmax_liv2013_ds_c.groupby('time.month').mean('time').values[wy_numbers],'r-',linewidth=2, label='Liv ('+str(station_elev)+' m)')
    plt.plot(wy_index, tmax_liv2015bc_ds_c.groupby('time.month').mean('time').values[wy_numbers],'r--',linewidth=2, label='Livbc ('+str(station_elev)+' m)')
    plt.plot(wy_index, tmax_wrf2014_ds_c.groupby('time.month').mean('time').values[wy_numbers],'b-', linewidth=2, label='WRF ('+str(station_elev)+' m)')
    plt.plot(wy_index, tmax_wrf2014bc_ds_c.groupby('time.month').mean('time').values[wy_numbers],'b--', linewidth=2, label='WRFbc ('+str(station_elev)+' m)')
    plt.plot(wy_index, temp_lapse+tmax_liv2013_ds_c.groupby('time.month').mean('time').values[wy_numbers],'y--',linewidth=3, label='Liv-6.5lapse ('+str(obs_elev)+' m)')
    plt.plot(wy_index, (obs_tmax.groupby(obs_tmax.index.month).mean())[wy_numbers_pd], 'k-', linewidth=2, label=obs_title+' ('+str(obs_elev)+' m)')
    
    plt.xlabel ('Month',fontsize=14)
    plt.ylabel('Temperature (deg C)',fontsize=14)
    plt.title(obs_title+' ('+str(obs_elev)+' m) vs Gridded Met ('+str(station_elev)+' m) \nAverage Monthly Max Temperature, Years: '+str(start_date_c.year)+'-'+str(end_date_c.year), fontsize=16)
    #plt.legend(bbox_to_anchor=(0.9, -0.15), ncol=3)
    plt.legend()
    plt.tick_params(labelsize=12)
    plt.grid(which='both')
    plt.xlim(1,12);
    plt.xticks(wy_index, month_strings);
    
    # Tmin
    fig, ax=plt.subplots(1,1,figsize=(6,4))
    plt.plot(wy_index, tmin_liv2013_ds_c.groupby('time.month').mean('time').values[wy_numbers],'r-',linewidth=2, label='Liv ('+str(station_elev)+' m)')
    plt.plot(wy_index, tmin_liv2015bc_ds_c.groupby('time.month').mean('time').values[wy_numbers],'r--',linewidth=2, label='Liv bc ('+str(station_elev)+' m)')
    plt.plot(wy_index, tmin_wrf2014_ds_c.groupby('time.month').mean('time').values[wy_numbers],'b-', linewidth=2, label='WRF ('+str(station_elev)+' m)')
    plt.plot(wy_index, tmin_wrf2014bc_ds_c.groupby('time.month').mean('time').values[wy_numbers],'b--', linewidth=2, label='WRFbc ('+str(station_elev)+' m)')
    plt.plot(wy_index, temp_lapse+tmin_liv2013_ds_c.groupby('time.month').mean('time').values[wy_numbers],'y--',linewidth=3, label='Liv -6.5lapse ('+str(obs_elev)+' m)')
    plt.plot(wy_index, obs_tmin.groupby(obs_tmin.index.month).mean()[wy_numbers_pd], 'k-', linewidth=2, label=obs_title+' ('+str(obs_elev)+' m)')
    
    plt.xlabel ('Month',fontsize=14)
    plt.ylabel('Temperature (deg C)',fontsize=14)
    plt.title(obs_title+' ('+str(obs_elev)+' m) vs Gridded Met ('+str(station_elev)+' m) \nAverage Monthly Min Temperature, Years: '+str(start_date_c.year)+'-'+str(end_date_c.year), fontsize=16)
    #plt.legend(bbox_to_anchor=(0.9, -0.15), ncol=3)
    plt.legend()
    plt.tick_params(labelsize=12)
    plt.grid(which='both')
    plt.xlim(1,12);
    plt.xticks(wy_index, month_strings);
    
    # Tavg
    fig, ax=plt.subplots(1,1,figsize=(6, 4))
    plt.plot(wy_index, (tmax_liv2013_ds_c.groupby('time.month').mean('time').values[wy_numbers]+ tmin_liv2013_ds_c.groupby('time.month').mean('time').values[wy_numbers])/2,'r-',linewidth=2, label='Liv 13 ('+str(station_elev)+' m)')
    plt.plot(wy_index,(tmax_liv2015bc_ds_c.groupby('time.month').mean('time').values[wy_numbers]+ tmin_liv2015bc_ds_c.groupby('time.month').mean('time').values[wy_numbers])/2,'r--',linewidth=2, label='Liv bc ('+str(station_elev)+' m)')
    plt.plot(wy_index,(tmax_wrf2014_ds_c.groupby('time.month').mean('time').values[wy_numbers]+ tmin_wrf2014_ds_c.groupby('time.month').mean('time').values[wy_numbers])/2,'b-', linewidth=2, label='WRF ('+str(station_elev)+' m)')
    plt.plot(wy_index, (tmax_wrf2014bc_ds_c.groupby('time.month').mean('time').values[wy_numbers]+ tmin_wrf2014bc_ds_c.groupby('time.month').mean('time').values[wy_numbers])/2,'b--', linewidth=2, label='WRFbc ('+str(station_elev)+' m)')
    plt.plot(wy_index, (temp_lapse+tmax_liv2013_ds_c.groupby('time.month').mean('time').values[wy_numbers]+ temp_lapse+tmin_liv2013_ds_c.groupby('time.month').mean('time').values[wy_numbers])/2,'y--',linewidth=3, label='Liv -6.5lapse ('+str(obs_elev)+' m)')
    plt.plot(wy_index, (obs_tmax.groupby(obs_tmax.index.month).mean()/2+obs_tmin.groupby(obs_tmin.index.month).mean()/2)[wy_numbers_pd], 'k-', linewidth=2, label=obs_title+' ('+str(obs_elev)+' m)')
    
    plt.xlabel ('Month',fontsize=14)
    plt.ylabel('Temperature (deg C)',fontsize=14)
    plt.title(obs_title+' ('+str(obs_elev)+' m) vs Gridded Met ('+str(station_elev)+' m) \nAverage Monthly Avg Temperature, Years: '+str(start_date_c.year)+'-'+str(end_date_c.year), fontsize=16)
    #plt.legend(bbox_to_anchor=(0.9, -0.15), ncol=3)
    plt.legend()
    plt.tick_params(labelsize=12)
    plt.grid(which='both')
    plt.xlim(1,12);
    plt.xticks(wy_index, month_strings);

# Precip
fig, ax=plt.subplots(1,1,figsize=(6,4))
plt.plot(wy_index, days_per_month*precip_liv2013_ds_c.groupby('time.month').mean('time').values[wy_numbers],'r-',linewidth=2, label='Liv 13 ('+str(station_elev)+' m)')
plt.plot(wy_index, days_per_month*precip_liv2015_ds_c.groupby('time.month').mean('time').values[wy_numbers],'g-',linewidth=2, label='Liv 15 ('+str(station_elev)+' m)')
plt.plot(wy_index, days_per_month*precip_liv2015bc_ds_c.groupby('time.month').mean('time').values[wy_numbers],'g--',linewidth=2, label='Liv 15 bc ('+str(station_elev)+' m)')
plt.plot(wy_index, days_per_month*precip_wrf2014_ds_c.groupby('time.month').mean('time').values[wy_numbers],'b-', linewidth=2, label='WRF ('+str(station_elev)+' m)')
plt.plot(wy_index, days_per_month*precip_wrf2014bc_ds_c.groupby('time.month').mean('time').values[wy_numbers],'b--', linewidth=2, label='WRFbc ('+str(obs_elev)+' m)')
plt.plot(wy_index, days_per_month*obs_precip.groupby(obs_precip.index.month).mean()[wy_numbers_pd], 'k-', linewidth=2, label=obs_title+' ('+str(obs_elev)+' m)')

plt.xlabel ('Month',fontsize=14)
plt.ylabel('Precipitation (mm)',fontsize=14)
plt.title(obs_title+' ('+str(obs_elev)+' m) vs Gridded Met('+str(station_elev)+' m) \nAverage Monthly Precipitation, Years:'+str(start_date_c.year)+'-'+str(end_date_c.year), fontsize=16)
#plt.legend(bbox_to_anchor=(0.9, -0.15), ncol=3)
plt.legend()
plt.tick_params(labelsize=12)
plt.grid(which='both')
plt.xlim(1,12);
plt.xticks(wy_index, month_strings);

#%% BIAS CORRECTION: Shift Liv 2013 monthly mean to match another dataset monthly mean
# Define correction functions
#Compute difference between monthly means for some data (e.g,. Temp and precip) for two different gridded datasets (e.g., Liv, WRF)
def compute_diff(closertotruth, needscorrection, vararray):
    diff_ds=xr.Dataset()
    for var in vararray:
        diff_ds[var] = eval('.'.join([closertotruth,var])).groupby('time.month').mean('time')-eval('.'.join([needscorrection,var])).groupby('time.month').mean('time')
    return diff_ds

def bc_diff(needscorrectionarray, needscorrectionvararray, correctionarray, newarray):
    for var in needscorrectionvararray:
        eval(newarray)[var]=eval('.'.join([needscorrectionarray,var])).groupby('time.month')+eval('.'.join([correctionarray,var]))

def compute_ratio(closertotruth, needscorrection, vararray):
    ratio_ds=xr.Dataset()
    for var in vararray:
        ratio_ds[var] = eval('.'.join([closertotruth,var])).groupby('time.month').mean('time')/eval('.'.join([needscorrection,var])).groupby('time.month').mean('time')
    return ratio_ds

def bc_ratio(needscorrectionarray, needscorrectionvararray, correctionarray, newarray):
    for var in needscorrectionvararray:
        eval(newarray)[var]=eval('.'.join([needscorrectionarray,var])).groupby('time.month')*eval('.'.join([correctionarray,var]))



#%%
### START INPUT
min_elev_1=met_df.elevation.min() # default is min elevation: met_df.elevation.min()
max_elev_1=met_df.elevation.max() # default is max elevation: met_df.elevation.max()
start_date_1=datetime.date(1950,1,1)
end_date_1=datetime.date(2010,12,31)
closertotruth_ratio_1='liv2015_ds_1'
needscorrection_ratio_1='liv2013_ds_1'
vararray_ratio_1=['precip']
new_array_1='liv2013corr_ds_1'
liv2013corr_ds_1=xr.Dataset()

min_elev_2=met_df.elevation.min() # default is min elevation: met_df.elevation.min()
max_elev_2=met_df.elevation.max() # default is max elevation: met_df.elevation.max()
start_date_2=datetime.date(1950,1,1)
end_date_2=datetime.date(2010,12,31)
closertotruth_diff_2='liv2015bc_ds_2'
needscorrection_diff_2='liv2013_ds_2'
vararray_diff_2=['tmax','tmin']
closertotruth_ratio_2='liv2015bc_ds_2'
needscorrection_ratio_2='liv2013_ds_2'
vararray_ratio_2=['precip']
new_array_2='liv2013corr_ds_2'
liv2013corr_ds_2=xr.Dataset()

min_elev_3=met_df.elevation.min() # default is min elevation: met_df.elevation.min()
max_elev_3=met_df.elevation.max() # default is max elevation: met_df.elevation.max()
start_date_3=datetime.date(1950,1,1)
end_date_3=datetime.date(2010,12,31)
closertotruth_diff_3='wrf2014_ds_3'
needscorrection_diff_3='liv2013_ds_3'
vararray_diff_3=['tmax','tmin']
closertotruth_ratio_3='wrf2014_ds_3'
needscorrection_ratio_3='liv2013_ds_3'
vararray_ratio_3=['precip']
new_array_3='liv2013corr_ds_3'
liv2013corr_ds_3=xr.Dataset()

min_elev_4=met_df.elevation.min() # default is min elevation: met_df.elevation.min()
max_elev_4=met_df.elevation.max() # default is max elevation: met_df.elevation.max()
start_date_4=datetime.date(1950,1,1)
end_date_4=datetime.date(2010,12,31)
closertotruth_diff_4='wrf2014bc_ds_4'
needscorrection_diff_4='liv2013_ds_4'
vararray_diff_4=['tmax','tmin']
closertotruth_ratio_4='wrf2014bc_ds_4'
needscorrection_ratio_4='liv2013_ds_4'
vararray_ratio_4=['precip']
new_array_4='liv2013corr_ds_4'
liv2013corr_ds_4=xr.Dataset()

min_elev_4=met_df.elevation.min() # default is min elevation: met_df.elevation.min()
max_elev_4=met_df.elevation.max() # default is max elevation: met_df.elevation.max()
start_date_4=datetime.date(1950,1,1)
end_date_4=datetime.date(2010,12,31)
closertotruth_diff_4='wrf2014bc_ds_4'
needscorrection_diff_4='liv2013_ds_4'
vararray_diff_4=['tmax','tmin']
closertotruth_ratio_4='wrf2014bc_ds_4'
needscorrection_ratio_4='liv2013_ds_4'
vararray_ratio_4=['precip']
new_array_4='liv2013corr_ds_4'
liv2013corr_ds_4=xr.Dataset()
### END INPUT

#%%
# Step 0: Correct liv 2013 ltm to liv 2015 ltm
# Extract data of interest for each step
stations_info_1=met_df[(met_df.elevation >= min_elev_1) & (met_df.elevation <= max_elev_1)] # Find points in elevation range of interest
elev_min_station_1=stations_info_1.elevation.idxmin() # station of minimum elevation in analysis
elev_min_1=stations_info_1.elevation.min() # minimum elevaiton of stations in analysis
elev_max_station_1=stations_info_1.elevation.idxmax() # station of maximum elevation in analysis
elev_max_1=stations_info_1.elevation.max() # maximum elevaiton of stations in analysis
stations_list_1=stations_info_1.station.values.astype('int64')
for name in forcing_list:
    ds=xr.Dataset()
    for var in ['precip']:
         ds[var]=eval(name+'_ds')[var].sel(time=slice(start_date_1, end_date_1), station=stations_list_1)
    globals()[name+'_ds_1']=ds
    print (name+'_ds_1')

stations_info_2=met_df[(met_df.elevation >= min_elev_2) & (met_df.elevation <= max_elev_2)] # Find points in elevation range of interest
elev_min_station_2=stations_info_2.elevation.idxmin() # station of minimum elevation in analysis
elev_min_2=stations_info_2.elevation.min() # minimum elevaiton of stations in analysis
elev_max_station_2=stations_info_2.elevation.idxmax() # station of maximum elevation in analysis
elev_max_2=stations_info_2.elevation.max() # maximum elevaiton of stations in analysis
stations_list_2=stations_info_2.station.values.astype('int64')
for name in forcing_list:
    ds=xr.Dataset()
    for var in  ['tmax', 'tmin', 'precip']:
         ds[var]=eval(name+'_ds')[var].sel(time=slice(start_date_2, end_date_2), station=stations_list_2)
    globals()[name+'_ds_2']=ds
    print (name+'_ds_2')

stations_info_3=met_df[(met_df.elevation >= min_elev_3) & (met_df.elevation <= max_elev_3)] # Find points in elevation range of interest
elev_min_station_3=stations_info_3.elevation.idxmin() # station of minimum elevation in analysis
elev_min_3=stations_info_3.elevation.min() # minimum elevaiton of stations in analysis
elev_max_station_3=stations_info_3.elevation.idxmax() # station of maximum elevation in analysis
elev_max_3=stations_info_3.elevation.max() # maximum elevaiton of stations in analysis
stations_list_3=stations_info_3.station.values.astype('int64')
for name in forcing_list:
    ds=xr.Dataset()
    for var in ['tmax', 'tmin', 'precip']:
         ds[var]=eval(name+'_ds')[var].sel(time=slice(start_date_3, end_date_3), station=stations_list_3)
    globals()[name+'_ds_3']=ds
    print (name+'_ds_3')
    
stations_info_4=met_df[(met_df.elevation >= min_elev_4) & (met_df.elevation <= max_elev_4)] # Find points in elevation range of interest
elev_min_station_4=stations_info_4.elevation.idxmin() # station of minimum elevation in analysis
elev_min_4=stations_info_4.elevation.min() # minimum elevaiton of stations in analysis
elev_max_station_4=stations_info_4.elevation.idxmax() # station of maximum elevation in analysis
elev_max_4=stations_info_4.elevation.max() # maximum elevaiton of stations in analysis
stations_list_4=stations_info_4.station.values.astype('int64')
for name in forcing_list:
    ds=xr.Dataset()
    for var in ['tmax', 'tmin', 'precip']:
         ds[var]=eval(name+'_ds')[var].sel(time=slice(start_date_4, end_date_4), station=stations_list_4)
    globals()[name+'_ds_4']=ds
    print (name+'_ds_4')

#%%
# Make corrections
liv2013corr_ds_1_ratio=compute_ratio(closertotruth_ratio_1, needscorrection_ratio_1, vararray_ratio_1)
bc_ratio('liv2013_ds', vararray_ratio_1, 'liv2013corr_ds_1_ratio', new_array_1)

liv2013corr_ds_2_diff=compute_diff(closertotruth_diff_2, needscorrection_diff_2, vararray_diff_2)
bc_diff('liv2013_ds', vararray_diff_2, 'liv2013corr_ds_2_diff', new_array_2)
liv2013corr_ds_2_ratio=compute_ratio(closertotruth_ratio_2, needscorrection_ratio_2, vararray_ratio_2)
bc_ratio('liv2013_ds', vararray_ratio_2, 'liv2013corr_ds_2_ratio', new_array_2)

liv2013corr_ds_3_diff=compute_diff(closertotruth_diff_3, needscorrection_diff_3, vararray_diff_3)
bc_diff('liv2013_ds', vararray_diff_3, 'liv2013corr_ds_3_diff', new_array_3)
liv2013corr_ds_3_ratio=compute_ratio(closertotruth_ratio_3, needscorrection_ratio_3, vararray_ratio_3)
bc_ratio('liv2013_ds', vararray_ratio_3, 'liv2013corr_ds_3_ratio', new_array_3)

liv2013corr_ds_4_diff=compute_diff(closertotruth_diff_4, needscorrection_diff_4, vararray_diff_4)
bc_diff('liv2013_ds', vararray_diff_4, 'liv2013corr_ds_4_diff', new_array_4)
liv2013corr_ds_4_ratio=compute_ratio(closertotruth_ratio_4, needscorrection_ratio_4, vararray_ratio_4)
bc_ratio('liv2013_ds', vararray_ratio_4, 'liv2013corr_ds_4_ratio', new_array_4)

    
# Tmax
fig, ax=plt.subplots(1,1,figsize=(6, 4))
plt.plot(wy_index, liv2013_ds.sel(time=slice(start_date_1, end_date_1)).tmax.groupby('time.month').mean('time').mean('station').values[wy_numbers],'k-',linewidth=2, label='Liv')
plt.plot(wy_index, liv2015bc_ds.tmax.groupby('time.month').mean('time').mean('station').values[wy_numbers],'b-',linewidth=2, label='Livbc')
plt.plot(wy_index, wrf2014_ds.tmax.groupby('time.month').mean('time').mean('station').values[wy_numbers],'g-', linewidth=2, label='WRF')
plt.plot(wy_index, wrf2014bc_ds.tmax.groupby('time.month').mean('time').mean('station').values[wy_numbers],'r-', linewidth=2, label='WRFbc')
plt.plot(wy_index, liv2013corr_ds_2.sel(time=slice(start_date_2, end_date_2)).tmax.groupby('time.month').mean('time').mean('station').values[wy_numbers],'c--', linewidth=2, label='Liv_Liv15bc_adj')
plt.plot(wy_index, liv2013corr_ds_3.sel(time=slice(start_date_3, end_date_3)).tmax.groupby('time.month').mean('time').mean('station').values[wy_numbers],'m--', linewidth=2, label='Liv_WRF_adj')
plt.plot(wy_index, liv2013corr_ds_4.sel(time=slice(start_date_4, end_date_4)).tmax.groupby('time.month').mean('time').mean('station').values[wy_numbers],'y--', linewidth=2, label='Liv_WRFbc_adj')
plt.xlabel ('Month',fontsize=14)
plt.ylabel('Temperature (deg C)',fontsize=14)
plt.title(str(loc_name)+'\nAverage Monthly Max Temperature\n Years: '+str(start_date_1.year)+'-'+str(end_date_1.year)+'; Elevation: '+str(min_elev_1)+'m -'+str(max_elev_1)+'m', fontsize=16)
#plt.legend(bbox_to_anchor=(0.9, -0.15), ncol=3)
plt.legend()
plt.tick_params(labelsize=12)
plt.grid(which='both')
plt.xlim(1,12);
plt.xticks(wy_index, month_strings);

# Tmin
fig, ax=plt.subplots(1,1,figsize=(6, 4))
plt.plot(wy_index,liv2013_ds.sel(time=slice(start_date_1, end_date_1)).tmin.groupby('time.month').mean('time').mean('station').values[wy_numbers],'k-',linewidth=2, label='Liv')
plt.plot(wy_index,liv2015bc_ds.tmin.groupby('time.month').mean('time').mean('station').values[wy_numbers],'b-',linewidth=2, label='Livbc')
plt.plot(wy_index,wrf2014_ds.tmin.groupby('time.month').mean('time').mean('station').values[wy_numbers],'g-', linewidth=2, label='WRF')
plt.plot(wy_index,wrf2014bc_ds.tmin.groupby('time.month').mean('time').mean('station').values[wy_numbers],'r-', linewidth=2, label='WRFbc')
plt.plot(wy_index, liv2013corr_ds_2.sel(time=slice(start_date_2, end_date_2)).tmin.groupby('time.month').mean('time').mean('station').values[wy_numbers],'c--', linewidth=2, label='Liv_Liv15bc_adj')
plt.plot(wy_index, liv2013corr_ds_3.sel(time=slice(start_date_3, end_date_3)).tmin.groupby('time.month').mean('time').mean('station').values[wy_numbers],'m--', linewidth=2, label='Liv_WRF_adj')
plt.plot(wy_index, liv2013corr_ds_4.sel(time=slice(start_date_4, end_date_4)).tmin.groupby('time.month').mean('time').mean('station').values[wy_numbers],'y--', linewidth=2, label='Liv_WRFbc_adj')

plt.xlabel ('Month',fontsize=14)
plt.ylabel('Temperature (deg C)',fontsize=14)
plt.title(str(loc_name)+'\nAverage Monthly Min Temperature\n Years: '+str(start_date_1.year)+'-'+str(end_date_1.year)+'; Elevation: '+str(min_elev_1)+'m -'+str(max_elev_1)+'m', fontsize=16)
#plt.legend(bbox_to_anchor=(0.9, -0.15), ncol=3)
plt.legend()
plt.tick_params(labelsize=12)
plt.grid(which='both')
plt.xlim(1,12);
plt.xticks(wy_index, month_strings);

# Precip
fig, ax=plt.subplots(1,1,figsize=(6,4))
plt.plot(wy_index, days_per_month*liv2013_ds_1.sel(time=slice(start_date_1, end_date_1)).precip.groupby('time.month').mean('time').mean('station').values[wy_numbers],'k-',linewidth=2, label='Liv13')
plt.plot(wy_index, days_per_month*liv2015_ds_1.precip.groupby('time.month').mean('time').mean('station').values[wy_numbers],'b-',linewidth=2, label='Liv15')
plt.plot(wy_index, days_per_month*liv2015bc_ds_1.precip.groupby('time.month').mean('time').mean('station').values[wy_numbers],'b*-',linewidth=2, label='Liv 15bc')
plt.plot(wy_index, days_per_month*wrf2014_ds_1.precip.groupby('time.month').mean('time').mean('station').values[wy_numbers],'g-',linewidth=2, label='WRF')
plt.plot(wy_index, days_per_month*wrf2014bc_ds_1.precip.groupby('time.month').mean('time').mean('station').values[wy_numbers],'g*-',linewidth=2, label='WRFbc')
plt.plot(wy_index, days_per_month*liv2013corr_ds_1.sel(time=slice(start_date_1, end_date_1)).precip.groupby('time.month').mean('time').mean('station').values[wy_numbers],'r--', linewidth=2, label='Liv_Liv15_adj')
plt.plot(wy_index, days_per_month*liv2013corr_ds_2.sel(time=slice(start_date_2, end_date_2)).precip.groupby('time.month').mean('time').mean('station').values[wy_numbers],'c--', linewidth=2, label='Liv_Liv15bc_adj')
plt.plot(wy_index, days_per_month*liv2013corr_ds_3.sel(time=slice(start_date_3, end_date_3)).precip.groupby('time.month').mean('time').mean('station').values[wy_numbers],'m--', linewidth=2, label='Liv_WRF_adj')
plt.plot(wy_index, days_per_month*liv2013corr_ds_4.sel(time=slice(start_date_4, end_date_4)).precip.groupby('time.month').mean('time').mean('station').values[wy_numbers],'y--', linewidth=2, label='Liv_WRFbc_adj')
plt.plot([1, 12],[0, 0], 'k-',linewidth=2)

plt.xlabel ('Month',fontsize=14)
plt.ylabel('Precipitation (mm)',fontsize=14)
plt.title(str(loc_name)+':\nAverage Monthly Precipitation\n Years: '+str(start_date_1.year)+'-'+str(end_date_1.year)+'; Elevation: '+str(min_elev_1)+'m -'+str(max_elev_1)+'m', fontsize=16)
plt.legend(loc='best')
plt.tick_params(labelsize=12)
plt.grid(which='both')
plt.xlim(1,12);
plt.xticks(wy_index, month_strings);

#%% COMPARE GRIDDED DATASETS TO OBSERVATION STATIONS
# Elwha RS elevation=110 m,  Cell 43 elevation= 468 m
### START INPUT
obs_variable='buckinghorse' # elwha_rs, sutherland, port_angeles, buckinghorse, waterhole

### END INPUT
# Extract observation station information
obs_elev=obs_stations_df.elevation[obs_variable]
obs_title=obs_stations_df.title[obs_variable]
start_date_c=obs_stations_df.start_date[obs_variable] #obs_stations_df.start_date[obs_variable]   datetime.datetime (1950,1,1)
end_date_c=obs_stations_df.end_date[obs_variable] #obs_stations_df.end_date[obs_variable]   datetime.datetime (2010,12,31)
station=obs_stations_df.met_station[obs_variable]

# Extract relevant obs info
obs_start_year=str(start_date_c.year)
obs_end_year=str(end_date_c.year)
obs_precip=eval(obs_variable)['PRCP'][obs_start_year:obs_end_year].value.astype('float')/10 # mm
if obs_variable is not 'sutherland':
    obs_tmax=eval(obs_variable)['TMAX'][obs_start_year:obs_end_year].value.astype('float')/10 # deg C
    obs_tmin=eval(obs_variable)['TMIN'][obs_start_year:obs_end_year].value.astype('float')/10 # deg C

# Extract relevant met info
station_elev=stations_info_c.loc[station,'elevation']
wy_numbers_pd=[x+1 for x in wy_numbers]

# Tmax
if obs_variable is not 'sutherland':
    fig, ax=plt.subplots(1,1,figsize=(6,4))
    plt.plot(wy_index, liv2013_ds.tmax.sel(time=slice(start_date_c, end_date_c), station=station).groupby('time.month').mean('time').values[wy_numbers],'r-',linewidth=2, label='Liv')
    plt.plot(wy_index, liv2013corr_ds_2.tmax.sel(time=slice(start_date_c, end_date_c), station=station).groupby('time.month').mean('time').values[wy_numbers],'c--', linewidth=2, label='Liv_Liv15bc_adj')
    plt.plot(wy_index, liv2013corr_ds_3.tmax.sel(time=slice(start_date_c, end_date_c), station=station).groupby('time.month').mean('time').values[wy_numbers],'m--', linewidth=2, label='Liv_WRF_adj')
    plt.plot(wy_index, liv2013corr_ds_4.tmax.sel(time=slice(start_date_c, end_date_c), station=station).groupby('time.month').mean('time').values[wy_numbers],'y--', linewidth=2, label='Liv_WRFbc_adj')
    plt.plot(wy_index, obs_tmax.groupby(obs_tmax.index.month).mean()[wy_numbers_pd], 'k-', linewidth=2, label=obs_title+' ('+str(obs_elev)+' m)')
    
    plt.xlabel ('Month',fontsize=14)
    plt.ylabel('Temperature (deg C)',fontsize=14)
    plt.title(obs_title+' ('+str(obs_elev)+' m) vs Gridded Met ('+str(station_elev)+' m) \nAverage Monthly Max Temperature, Years: '+str(start_date_c.year)+'-'+str(end_date_c.year), fontsize=16)
    #plt.legend(bbox_to_anchor=(0.9, -0.15), ncol=3)
    plt.legend()
    plt.tick_params(labelsize=12)
    plt.grid(which='both')
    plt.xlim(1,12);
    plt.xticks(wy_index, month_strings);
   
    # Tmin
    fig, ax=plt.subplots(1,1,figsize=(6,4))
    plt.plot(wy_index, liv2013_ds.tmin.sel(time=slice(start_date_c, end_date_c), station=station).groupby('time.month').mean('time').values[wy_numbers],'r-',linewidth=2, label='Liv')
    plt.plot(wy_index, liv2013corr_ds_2.tmin.sel(time=slice(start_date_c, end_date_c), station=station).groupby('time.month').mean('time').values[wy_numbers],'c--', linewidth=2, label='Liv_Liv15bc_adj')
    plt.plot(wy_index, liv2013corr_ds_3.tmin.sel(time=slice(start_date_c, end_date_c), station=station).groupby('time.month').mean('time').values[wy_numbers],'m--', linewidth=2, label='Liv_WRF_adj')
    plt.plot(wy_index, liv2013corr_ds_4.tmin.sel(time=slice(start_date_c, end_date_c), station=station).groupby('time.month').mean('time').values[wy_numbers],'y--', linewidth=2, label='Liv_WRFbc_adj')
    plt.plot(wy_index, obs_tmin.groupby(obs_tmin.index.month).mean()[wy_numbers_pd], 'k-', linewidth=2, label=obs_title+' ('+str(obs_elev)+' m)')
    
    plt.xlabel ('Month',fontsize=14)
    plt.ylabel('Temperature (deg C)',fontsize=14)
    plt.title(obs_title+' ('+str(obs_elev)+' m) vs Gridded Met ('+str(station_elev)+' m) \nAverage Monthly Min Temperature, Years: '+str(start_date_c.year)+'-'+str(end_date_c.year), fontsize=16)
    #plt.legend(bbox_to_anchor=(0.9, -0.15), ncol=3)
    plt.legend()
    plt.tick_params(labelsize=12)
    plt.grid(which='both')
    plt.xlim(1,12);
    plt.xticks(wy_index, month_strings);

# Precip
fig, ax=plt.subplots(1,1,figsize=(6,4))
plt.plot(wy_index, days_per_month*liv2013_ds.precip.sel(time=slice(start_date_c, end_date_c), station=station).groupby('time.month').mean('time').values[wy_numbers],'r-',linewidth=2, label='Liv13')
plt.plot(wy_index, days_per_month*liv2013corr_ds_1.precip.sel(time=slice(start_date_c, end_date_c), station=station).groupby('time.month').mean('time').values[wy_numbers],'g--', linewidth=2, label='Liv_Liv15_adj')
plt.plot(wy_index, days_per_month*liv2013corr_ds_2.precip.sel(time=slice(start_date_c, end_date_c), station=station).groupby('time.month').mean('time').values[wy_numbers],'c--', linewidth=2, label='Liv_Liv15bc_adj')
plt.plot(wy_index, days_per_month*liv2013corr_ds_3.precip.sel(time=slice(start_date_c, end_date_c), station=station).groupby('time.month').mean('time').values[wy_numbers],'m--', linewidth=2, label='Liv_WRF_adj')
plt.plot(wy_index, days_per_month*liv2013corr_ds_4.precip.sel(time=slice(start_date_c, end_date_c), station=station).groupby('time.month').mean('time').values[wy_numbers],'y--', linewidth=2, label='Liv_WRFbc_adj')
plt.plot(wy_index, days_per_month*obs_precip.groupby(obs_precip.index.month).mean()[wy_numbers_pd], 'k-', linewidth=2, label=obs_title+' ('+str(obs_elev)+' m)')

plt.xlabel ('Month',fontsize=14)
plt.ylabel('Precipitation (mm)',fontsize=14)
plt.title(obs_title+' ('+str(obs_elev)+' m) vs Gridded Met('+str(station_elev)+' m) \nAverage Monthly Precipitation, Years:'+str(start_date_c.year)+'-'+str(end_date_c.year), fontsize=16)
#plt.legend(bbox_to_anchor=(0.9, -0.15), ncol=3)
plt.legend()
plt.tick_params(labelsize=12)
plt.grid(which='both')
plt.xlim(1,12);
plt.xticks(wy_index, month_strings);

#%%
test_station=15
plot_index=liv2013_ds.precip.time.values
plt.plot(plot_index, liv2013_ds.precip.sel(station=test_station).values,'k-',linewidth=2, label='Liv13')
plt.plot(plot_index, liv2013corr_ds_1.precip.sel(station=test_station).values,'r-', linewidth=2, label='Liv_Liv15_adj')
plt.plot(plot_index, liv2013corr_ds_2.precip.sel(station=test_station).values,'b-', linewidth=2, label='Liv_Liv15bc_adj')
plt.plot(plot_index, liv2013corr_ds_3.precip.sel(station=test_station).values,'g-', linewidth=2, label='Liv_WRF_adj')
plt.plot(plot_index, liv2013corr_ds_4.precip.sel(station=test_station).values,'y--', linewidth=2, label='Liv_WRFbc_adj')

plt.xlabel ('Date',fontsize=14)
plt.ylabel('Precipitation (mm)',fontsize=14)
plt.title('Station '+str(test_station)+':Daily Precipitation', fontsize=16)
plt.legend(loc='best')
plt.tick_params(labelsize=12)
plt.grid(which='both')


#%% Cumulative Precipitation versus Runoff at Lake Mills Gage
t1=np.datetime64('1994-10-01','D')
t2=np.datetime64('1997-10-01','D')
t3=np.datetime64('2004-10-01','D')
t4=np.datetime64('2011-10-01','D')
q_mills_wy_dates=np.concatenate((np.arange(t1, t2), np.arange(t3, t4)),axis=0)

station_area_weights_mills_ds=xr.DataArray(station_geom_mills.Area_m2/drainage_area_mills, coords=[station_geom_mills.index], dims=['station'])

cum_precip_liv2013_mills_wy_df=(liv2013_ds.sel(time=q_mills_wy_dates,station=station_geom_mills.index).precip*station_area_weights_mills_ds).cumsum('time').sum('station')
cum_precip_corr1_mills_wy_df=(liv2013corr_ds_1.sel(time=q_mills_wy_dates,station=station_geom_mills.index).precip*station_area_weights_mills_ds).cumsum('time').sum('station')
cum_precip_corr2_mills_wy_df=(liv2013corr_ds_2.sel(time=q_mills_wy_dates,station=station_geom_mills.index).precip*station_area_weights_mills_ds).cumsum('time').sum('station')
cum_precip_corr3_mills_wy_df=(liv2013corr_ds_3.sel(time=q_mills_wy_dates,station=station_geom_mills.index).precip*station_area_weights_mills_ds).cumsum('time').sum('station')
cum_precip_corr4_mills_wy_df=(liv2013corr_ds_4.sel(time=q_mills_wy_dates,station=station_geom_mills.index).precip*station_area_weights_mills_ds).cumsum('time').sum('station')
cum_q_liv2013_mills_wy_df=q_mills_df.loc[q_mills_wy_dates,'flow_mmday'].cumsum()

RR_Liv2013=cum_q_liv2013_mills_wy_df[-1]/cum_precip_liv2013_mills_wy_df[-1]
RR_corr1=cum_q_liv2013_mills_wy_df[-1]/cum_precip_corr1_mills_wy_df[-1]
RR_corr2=cum_q_liv2013_mills_wy_df[-1]/cum_precip_corr2_mills_wy_df[-1]
RR_corr3=cum_q_liv2013_mills_wy_df[-1]/cum_precip_corr3_mills_wy_df[-1]
RR_corr4=cum_q_liv2013_mills_wy_df[-1]/cum_precip_corr4_mills_wy_df[-1]
print('RR Liv 2013, Mills gage=',RR_Liv2013.values)
print('RR Liv 2015, Mills gage=',RR_corr1.values)
print('RR Liv 2015bc, Mills gage=',RR_corr2.values)
print('RR WRF, Mills gage=',RR_corr3.values)
print('RR WRFbc, Mills gage=',RR_corr4.values)

fig1, ax1=plt.subplots(1,1,figsize=(6,4))
plt.xticks(rotation=40)
lns1=ax1.plot(cum_precip_liv2013_mills_wy_df.time, cum_precip_liv2013_mills_wy_df,'r-', label='P- L13, RR='+str(RR_Liv2013.values.round(2)),linewidth=2)
lns3=ax1.plot(cum_precip_liv2013_mills_wy_df.time, cum_precip_corr3_mills_wy_df,'c-', label='P- WRF 16, RR='+str(RR_corr3.values.round(2)),linewidth=2)
lns2=ax1.plot(cum_precip_liv2013_mills_wy_df.time, cum_precip_corr1_mills_wy_df,'b-', label='P- L15, RR='+str(RR_corr1.values.round(2)),linewidth=2)
#lns5=ax1.plot(cum_precip_liv2013_mills_wy_df.time, cum_precip_corr4_mills_wy_df,'y-', label='P- WRF-PRISM, RR='+str(RR_corr4.values.round(2)),linewidth=2)
lns4=ax1.plot(cum_precip_liv2013_mills_wy_df.time, cum_precip_corr2_mills_wy_df,'g--', label='P-L15 BC, RR='+str(RR_corr2.values.round(2)),linewidth=2)
lns6=ax1.plot(cum_precip_liv2013_mills_wy_df.time, cum_q_liv2013_mills_wy_df,'k-', label='Q', linewidth=2)
plt.grid(which='both')
ax1.set_xlabel('Date',fontsize=16)
ax1.set_ylabel('Streamflow and\nPrecipitation (mm)',fontsize=16) 
plt.title('Lake Mills Gage Tributary\nCumulative Precipitation and Streamflow\n WY 1995-1997, 2005-2011 (8 yrs total)',fontsize=15)

lns=lns1+lns2+lns3+lns4+lns6
labs = [l.get_label() for l in lns]
by_label = OrderedDict(zip(labs, lns))
plt.legend(by_label.values(), by_label.keys(), loc='best')
#plt.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(0.9, -0.4), ncol=2)

#%% Cumulative Precipitation versus Runoff at McDonald Gage
q_mcd_wy_dates=slice(datetime.datetime(1918, 10, 1 ), datetime.datetime(2010, 9,30))
station_area_weights_mcd_ds=xr.DataArray(station_geom_mcd.Area_m2/drainage_area_mcd, coords=[station_geom_mcd.index], dims=['station'])

cum_precip_liv2013_mcd_wy_df=(liv2013_ds.sel(time=q_mcd_wy_dates,station=station_geom_mcd.index).precip*station_area_weights_mcd_ds).cumsum('time').sum('station')
cum_precip_corr1_mcd_wy_df=(liv2013corr_ds_1.sel(time=q_mcd_wy_dates,station=station_geom_mcd.index).precip*station_area_weights_mcd_ds).cumsum('time').sum('station')
cum_precip_corr2_mcd_wy_df=(liv2013corr_ds_2.sel(time=q_mcd_wy_dates,station=station_geom_mcd.index).precip*station_area_weights_mcd_ds).cumsum('time').sum('station')
cum_precip_corr3_mcd_wy_df=(liv2013corr_ds_3.sel(time=q_mcd_wy_dates,station=station_geom_mcd.index).precip*station_area_weights_mcd_ds).cumsum('time').sum('station')
cum_precip_corr4_mcd_wy_df=(liv2013corr_ds_4.sel(time=q_mcd_wy_dates,station=station_geom_mcd.index).precip*station_area_weights_mcd_ds).cumsum('time').sum('station')
cum_q_liv2013_mcd_wy_df=q_mcd_df.loc[q_mcd_wy_dates,'flow_mmday'].cumsum()

RR_Liv2013=cum_q_liv2013_mcd_wy_df[-1]/cum_precip_liv2013_mcd_wy_df[-1]
RR_corr1=cum_q_liv2013_mcd_wy_df[-1]/cum_precip_corr1_mcd_wy_df[-1]
RR_corr2=cum_q_liv2013_mcd_wy_df[-1]/cum_precip_corr2_mcd_wy_df[-1]
RR_corr3=cum_q_liv2013_mcd_wy_df[-1]/cum_precip_corr3_mcd_wy_df[-1]
RR_corr4=cum_q_liv2013_mcd_wy_df[-1]/cum_precip_corr4_mcd_wy_df[-1]
print('RR Liv 2013, McD gage=',RR_Liv2013.values)
print('RR Liv 2015, McD gage=',RR_corr1.values)
print('RR Liv 2015bc, McD gage=',RR_corr2.values)
print('RR WRF, McD gage=',RR_corr3.values)
print('RR WRFbc, McD gage=',RR_corr4.values)

fig1, ax1=plt.subplots(1,1,figsize=(6,4))
plt.xticks(rotation=40)
lns1=ax1.plot(cum_precip_liv2013_mcd_wy_df.time, cum_precip_liv2013_mcd_wy_df,'r-', label='P- L13, RR='+str(RR_Liv2013.values.round(2)),linewidth=2)
lns2=ax1.plot(cum_precip_liv2013_mcd_wy_df.time, cum_precip_corr1_mcd_wy_df,'b-', label='P- L15, RR='+str(RR_corr1.values.round(2)),linewidth=2)
lns3=ax1.plot(cum_precip_liv2013_mcd_wy_df.time, cum_precip_corr3_mcd_wy_df,'c-', label='P- WRF 16, RR='+str(RR_corr3.values.round(2)),linewidth=2)
#lns4=ax1.plot(cum_precip_liv2013_mcd_wy_df.time, cum_precip_corr4_mcd_wy_df,'y-', label='P- WRF-PRISM, RR='+str(RR_corr4.values.round(2)),linewidth=2)
lns5=ax1.plot(cum_precip_liv2013_mcd_wy_df.time, cum_precip_corr2_mcd_wy_df,'g--', label='P- L15 BC, RR='+str(RR_corr2.values.round(2)),linewidth=2)
lns6=ax1.plot(cum_precip_liv2013_mcd_wy_df.time, cum_q_liv2013_mcd_wy_df,'k-', label='Q', linewidth=2)
plt.grid(which='both')
ax1.set_xlabel('Date',fontsize=16)
ax1.set_ylabel('Streamflow and\nPrecipitation (mm)',fontsize=16) 
plt.title('McDonald Gage Tributary\nCumulative Precipitation and Streamflow\n WY 1919-2010 (91 years total)',fontsize=15)

lns=lns1+lns2+lns3+lns5+lns6
labs = [l.get_label() for l in lns]
by_label = OrderedDict(zip(labs, lns))
plt.legend(by_label.values(), by_label.keys(), loc='best')
#plt.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(0.9, -0.4), ncol=2)

#%% Save to CSV
# START INPUT
ds_name='liv2015adj' # Note that this has precip only, so tmax and tmin come from liv 2013
ds_var=liv2013corr_ds_1
# END INPUT

dest_dir=homedir+'/forcings/data/'+ds_name
os.chdir(dest_dir)

for i in range(0, len(met_df.index)):
    file_name='_'.join([ds_name,met_df.loc[i, 'latitude'],met_df.loc[i, 'longitude']])
    array=pd.DataFrame()
    array.loc[:,0]=ds_var.precip[i].values.T
    array.loc[:,1]=liv2013_ds.tmax[i].values.T
    array.loc[:,2]=liv2013_ds.tmin[i].values.T
    array.loc[:,3]=liv2013_ds.wind[i].values.T
    array.to_csv(os.path.join(dest_dir, os.path.basename(file_name)), sep='\t', header=None, index=False, float_format='%.4f')
#%%
# START INPUT    
ds_name='liv2015bcadj'
ds_var=liv2013corr_ds_2
# END INPUT

dest_dir=homedir+'/forcings/data/'+ds_name
os.chdir(dest_dir)

for i in range(0, len(met_df.index)):
    file_name='_'.join([ds_name,met_df.loc[i, 'latitude'],met_df.loc[i, 'longitude']])
    array=pd.DataFrame()
    array.loc[:,0]=ds_var.precip[i].values.T
    array.loc[:,1]=ds_var.tmax[i].values.T
    array.loc[:,2]=ds_var.tmin[i].values.T
    array.loc[:,3]=liv2013_ds.wind[i].values.T
    array.to_csv(os.path.join(dest_dir, os.path.basename(file_name)), sep='\t', header=None, index=False, float_format='%.4f')    
#%%
# START INPUT   
ds_name='WRFadj'
ds_var=liv2013corr_ds_3
# END INPUT

dest_dir=homedir+'/forcings/data/'+ds_name
os.chdir(dest_dir)

for i in range(0, len(met_df.index)):
    file_name='_'.join([ds_name,met_df.loc[i, 'latitude'],met_df.loc[i, 'longitude']])
    array=pd.DataFrame()
    array.loc[:,0]=ds_var.precip[i].values.T
    array.loc[:,1]=ds_var.tmax[i].values.T
    array.loc[:,2]=ds_var.tmin[i].values.T
    array.loc[:,3]=liv2013_ds.wind[i].values.T
    array.to_csv(os.path.join(dest_dir, os.path.basename(file_name)), sep='\t', header=None, index=False, float_format='%.4f')
#%%
# START INPUT  
ds_name='WRFbcadj'
ds_var=liv2013corr_ds_4
# END INPUT

dest_dir=homedir+'/forcings/data/'+ds_name
os.chdir(dest_dir)

for i in range(0, len(met_df.index)):
    file_name='_'.join([ds_name,met_df.loc[i, 'latitude'],met_df.loc[i, 'longitude']])
    array=pd.DataFrame()
    array.loc[:,0]=ds_var.precip[i].values.T
    array.loc[:,1]=ds_var.tmax[i].values.T
    array.loc[:,2]=ds_var.tmin[i].values.T
    array.loc[:,3]=liv2013_ds.wind[i].values.T
    array.to_csv(os.path.join(dest_dir, os.path.basename(file_name)), sep='\t', header=None, index=False, float_format='%.4f')
#%% Low elevation Bias Correction- Livneh, 2013

# dates for 1950 bias corrections
start_date_bc50=datetime.date(1950,1,1)
end_date_bc50=datetime.date(2010,12,31)
# dates for 1915 bias corrections
start_date_bc15=datetime.date(1915,1,1)
end_date_bc15=datetime.date(2011,12,31)

min_elev_le_bc1=met_df.elevation.min() # default is min elevation: met_df.elevation.min()
max_elev_le_bc1=500 # default is max elevation: met_df.elevation.max()

min_elev_all_bc1=met_df.elevation.min() # default is min elevation: met_df.elevation.min()
max_elev_all_bc1=met_df.elevation.max() # default is max elevation: met_df.elevation.max()

# Extract data of interest for each step
stations_info_le_bc1=met_df[(met_df.elevation >= min_elev_le_bc1) & (met_df.elevation <= max_elev_le_bc1)] # Find points in elevation range of interest
elev_min_station_le_bc1=stations_info_le_bc1.elevation.idxmin() # station of minimum elevation in analysis
elev_min_le_bc1=stations_info_le_bc1.elevation.min() # minimum elevaiton of stations in analysis
elev_max_station_le_bc1=stations_info_le_bc1.elevation.idxmax() # station of maximum elevation in analysis
elev_max_le_bc1=stations_info_le_bc1.elevation.max() # maximum elevaiton of stations in analysis
stations_list_le_bc1=stations_info_le_bc1.station.values.astype('int64')
for name in ['liv2013', 'wrf2014']:
    ds=xr.Dataset()
    for var in  ['tmax', 'tmin', 'precip']:
         ds[var]=eval(name+'_ds')[var].sel(time=slice(start_date_bc50, end_date_bc50), station=stations_list_le_bc1)
    globals()[name+'_ds_bc1']=ds
    print (name+'_ds_bc1')


# Make corrections
# Temperature
corr_le_bc1_diff=compute_diff('liv2013_ds_bc1', 'wrf2014_ds_bc1', ['tmax','tmin'])

mean_corr_le_bc1_tmax=corr_le_bc1_diff.tmax.groupby('month').mean()
mean_corr_le_bc1_tmin=corr_le_bc1_diff.tmin.groupby('month').mean()

wrf2014_lecorr_bc1_tmax=wrf2014_ds.tmax.groupby('time.month')+mean_corr_le_bc1_tmax
wrf2014_lecorr_bc1_tmin=wrf2014_ds.tmin.groupby('time.month')+mean_corr_le_bc1_tmin

corr_all_bc1_diff_tmax=wrf2014_lecorr_bc1_tmax.groupby('time.month').mean('time')-liv2013_ds.tmax.sel(time=slice(start_date_bc50, end_date_bc50)).groupby('time.month').mean('time')
corr_all_bc1_diff_tmin=wrf2014_lecorr_bc1_tmin.groupby('time.month').mean('time')-liv2013_ds.tmin.sel(time=slice(start_date_bc50, end_date_bc50)).groupby('time.month').mean('time')

liv2013_ds_bc1_tmax=liv2013_ds.tmax.sel(time=slice(start_date_bc15, end_date_bc15)).groupby('time.month')+corr_all_bc1_diff_tmax
liv2013_ds_bc1_tmin=liv2013_ds.tmin.sel(time=slice(start_date_bc15, end_date_bc15)).groupby('time.month')+corr_all_bc1_diff_tmin

# Precip
corr_le_bc1_ratio=compute_ratio('liv2013_ds_bc1', 'wrf2014_ds_bc1', ['precip']) # Compute ratio between low elevation means
mean_corr_le_bc1_precip=corr_le_bc1_ratio.precip.groupby('month').mean() # compute mean of ratio for low elevation stations
wrf2014_lecorr_bc1_precip=wrf2014_ds.precip.groupby('time.month')+wrf2014_ds.precip.groupby('time.month').mean('time')*(mean_corr_le_bc1_precip-1)  # replace mean of all stations
corr_all_bc1_diff_precip=wrf2014_lecorr_bc1_precip.groupby('time.month').mean('time')/liv2013_ds.precip.sel(time=slice(start_date_bc50, end_date_bc50)).groupby('time.month').mean('time')
# Compute mean ratio for all stations
liv2013_ds_bc1_precip=liv2013_ds.precip.sel(time=slice(start_date_bc15, end_date_bc15)).groupby('time.month')+liv2013_ds.precip.sel(time=slice(start_date_bc15, end_date_bc15)).groupby('time.month').mean('time')*(corr_all_bc1_diff_precip-1)


stations_plot=station_geom_mills.index # station_geom_mills.index OR stations_list_le_bc1
# Tmax
fig, ax=plt.subplots(1,1,figsize=(6, 4))
plt.plot(wy_index, liv2013_ds.sel(time=slice(start_date_bc50, end_date_bc50),station=stations_plot).tmax.groupby('time.month').mean('time').mean('station').values[wy_numbers],'b-',linewidth=2, label='Liv')
plt.plot(wy_index, wrf2014_ds.sel(station=stations_plot).tmax.groupby('time.month').mean('time').mean('station').values[wy_numbers],'r-', linewidth=2, label='WRF')
plt.plot(wy_index, liv2013_ds_bc1_tmax.sel(time=slice(start_date_bc50, end_date_bc50),station=stations_plot).groupby('time.month').mean('time').mean('station').values[wy_numbers],'k-', linewidth=2, label='BC')
plt.xlabel ('Month',fontsize=14)
plt.ylabel('Temperature (deg C)',fontsize=14)
plt.title('Bias-Corrected Average Monthly Max Temperature', fontsize=16)
#plt.legend(bbox_to_anchor=(0.9, -0.15), ncol=3)
plt.legend()
plt.tick_params(labelsize=12)
plt.grid(which='both')
plt.xlim(1,12);
plt.xticks(wy_index, month_strings);

# Tmin
fig, ax=plt.subplots(1,1,figsize=(6, 4))
plt.plot(wy_index, liv2013_ds.sel(time=slice(start_date_bc50, end_date_bc50),station=stations_plot).tmin.groupby('time.month').mean('time').mean('station').values[wy_numbers],'b-',linewidth=2, label='Liv')
plt.plot(wy_index, wrf2014_ds.sel(station=stations_plot).tmin.groupby('time.month').mean('time').mean('station').values[wy_numbers],'r-', linewidth=2, label='WRF')
plt.plot(wy_index, liv2013_ds_bc1_tmin.sel(time=slice(start_date_bc50, end_date_bc50),station=stations_plot).groupby('time.month').mean('time').mean('station').values[wy_numbers],'k-', linewidth=2, label='BC')
plt.xlabel ('Month',fontsize=14)
plt.ylabel('Temperature (deg C)',fontsize=14)
plt.title('Bias-Corrected Average Monthly Min Temperature', fontsize=16)
#plt.legend(bbox_to_anchor=(0.9, -0.15), ncol=3)
plt.legend()
plt.tick_params(labelsize=12)
plt.grid(which='both')
plt.xlim(1,12);
plt.xticks(wy_index, month_strings);

# Precip
fig, ax=plt.subplots(1,1,figsize=(6,4))
plt.plot(wy_index, days_per_month*liv2013_ds.sel(time=slice(start_date_bc50, end_date_bc50),station=stations_plot).precip.groupby('time.month').mean('time').mean('station').values[wy_numbers],'b-',linewidth=2, label='Liv')
plt.plot(wy_index, days_per_month*wrf2014_ds.sel(station=stations_plot).precip.groupby('time.month').mean('time').mean('station').values[wy_numbers],'r-', linewidth=2, label='WRF')
plt.plot(wy_index, days_per_month*liv2013_ds_bc1_precip.sel(time=slice(start_date_bc50, end_date_bc50),station=stations_plot).groupby('time.month').mean('time').mean('station').values[wy_numbers],'k-', linewidth=2, label='BC')
plt.plot([1, 12],[0, 0], 'k--',linewidth=1)

plt.xlabel ('Month',fontsize=14)
plt.ylabel('Precipitation (mm)',fontsize=14)
plt.title('Bias-Corrected Average Monthly Precipitation', fontsize=16)
plt.legend(loc='best')
plt.tick_params(labelsize=12)
plt.grid(which='both')
plt.xlim(1,12);
plt.xticks(wy_index, month_strings);


#%% Low elevation Bias Correction- Livneh, 2015

# dates for 1950 bias corrections
start_date_bc50=datetime.date(1950,1,1)
end_date_bc50=datetime.date(2010,12,31)
# dates for 1915 bias corrections
start_date_bc15=datetime.date(1915,1,1)
end_date_bc15=datetime.date(2011,12,31)

min_elev_le_bc1=met_df.elevation.min() # default is min elevation: met_df.elevation.min()
max_elev_le_bc1=500 # default is max elevation: met_df.elevation.max()

min_elev_all_bc1=met_df.elevation.min() # default is min elevation: met_df.elevation.min()
max_elev_all_bc1=met_df.elevation.max() # default is max elevation: met_df.elevation.max()

# Extract data of interest for each step
stations_info_le_bc1=met_df[(met_df.elevation >= min_elev_le_bc1) & (met_df.elevation <= max_elev_le_bc1)] # Find points in elevation range of interest
elev_min_station_le_bc1=stations_info_le_bc1.elevation.idxmin() # station of minimum elevation in analysis
elev_min_le_bc1=stations_info_le_bc1.elevation.min() # minimum elevaiton of stations in analysis
elev_max_station_le_bc1=stations_info_le_bc1.elevation.idxmax() # station of maximum elevation in analysis
elev_max_le_bc1=stations_info_le_bc1.elevation.max() # maximum elevaiton of stations in analysis
stations_list_le_bc1=stations_info_le_bc1.station.values.astype('int64')
for name in ['liv2015', 'wrf2014']:
    ds=xr.Dataset()
    for var in  ['tmax', 'tmin', 'precip']:
         ds[var]=eval(name+'_ds')[var].sel(time=slice(start_date_bc50, end_date_bc50), station=stations_list_le_bc1)
    globals()[name+'_ds_bc1']=ds
    print (name+'_ds_bc1')


# Make corrections
# Temperature
corr_le_bc1_diff=compute_diff('liv2015_ds_bc1', 'wrf2014_ds_bc1', ['tmax','tmin'])

mean_corr_le_bc1_tmax=corr_le_bc1_diff.tmax.groupby('month').mean()
mean_corr_le_bc1_tmin=corr_le_bc1_diff.tmin.groupby('month').mean()

wrf2014_lecorr_bc1_tmax=wrf2014_ds.tmax.groupby('time.month')+mean_corr_le_bc1_tmax
wrf2014_lecorr_bc1_tmin=wrf2014_ds.tmin.groupby('time.month')+mean_corr_le_bc1_tmin

corr_all_bc1_diff_tmax=wrf2014_lecorr_bc1_tmax.groupby('time.month').mean('time')-liv2013_ds.tmax.sel(time=slice(start_date_bc50, end_date_bc50)).groupby('time.month').mean('time')
corr_all_bc1_diff_tmin=wrf2014_lecorr_bc1_tmin.groupby('time.month').mean('time')-liv2013_ds.tmin.sel(time=slice(start_date_bc50, end_date_bc50)).groupby('time.month').mean('time')

liv2013_ds_bc1_tmax=liv2013_ds.tmax.sel(time=slice(start_date_bc15, end_date_bc15)).groupby('time.month')+corr_all_bc1_diff_tmax
liv2013_ds_bc1_tmin=liv2013_ds.tmin.sel(time=slice(start_date_bc15, end_date_bc15)).groupby('time.month')+corr_all_bc1_diff_tmin

#%%
# Precip
corr_le_bc1_ratio=compute_ratio('liv2015_ds_bc1', 'wrf2014_ds_bc1', ['precip']) # Compute ratio between low elevation means
mean_corr_le_bc1_precip=corr_le_bc1_ratio.precip.groupby('month').mean() # compute mean of ratio for low elevation stations


wrf2014_lecorr_bc1_precip=wrf2014_ds.precip.groupby('time.month')+wrf2014_ds.precip.groupby('time.month').mean('time')*(mean_corr_le_bc1_precip-1)  # replace mean of all stations

corr_all_bc1_diff_precip=wrf2014_lecorr_bc1_precip.groupby('time.month').mean('time')/liv2013_ds.precip.sel(time=slice(start_date_bc50, end_date_bc50)).groupby('time.month').mean('time')

# Compute mean ratio for all stations
liv2013_ds_bc1_precip=liv2013_ds.precip.sel(time=slice(start_date_bc15, end_date_bc15)).groupby('time.month')+liv2013_ds.precip.sel(time=slice(start_date_bc15, end_date_bc15)).groupby('time.month').mean('time')*(corr_all_bc1_diff_precip-1)

# PLOTS
stations_plot=station_geom_mills.index # stations_list_le_bc1 OR station_geom_mills.index
# Tmax
fig, ax=plt.subplots(1,1,figsize=(6, 4))
plt.plot(wy_index, liv2015_ds.sel(time=slice(start_date_bc50, end_date_bc50),station=stations_plot).tmax.groupby('time.month').mean('time').mean('station').values[wy_numbers],'b-',linewidth=2, label='Liv')
plt.plot(wy_index, wrf2014_ds.sel(station=stations_plot).tmax.groupby('time.month').mean('time').mean('station').values[wy_numbers],'r-', linewidth=2, label='WRF')
plt.plot(wy_index, liv2013_ds_bc1_tmax.sel(time=slice(start_date_bc50, end_date_bc50),station=stations_plot).groupby('time.month').mean('time').mean('station').values[wy_numbers],'k-', linewidth=2, label='BC')
plt.xlabel ('Month',fontsize=14)
plt.ylabel('Temperature (deg C)',fontsize=14)
plt.title('Bias-Corrected Average Monthly Max Temperature', fontsize=16)
#plt.legend(bbox_to_anchor=(0.9, -0.15), ncol=3)
plt.legend()
plt.tick_params(labelsize=12)
plt.grid(which='both')
plt.xlim(1,12);
plt.xticks(wy_index, month_strings);

# Tmin
fig, ax=plt.subplots(1,1,figsize=(6, 4))
plt.plot(wy_index, liv2015_ds.sel(time=slice(start_date_bc50, end_date_bc50),station=stations_plot).tmin.groupby('time.month').mean('time').mean('station').values[wy_numbers],'b-',linewidth=2, label='Liv')
plt.plot(wy_index, wrf2014_ds.sel(station=stations_plot).tmin.groupby('time.month').mean('time').mean('station').values[wy_numbers],'r-', linewidth=2, label='WRF')
plt.plot(wy_index, liv2013_ds_bc1_tmin.sel(time=slice(start_date_bc50, end_date_bc50),station=stations_plot).groupby('time.month').mean('time').mean('station').values[wy_numbers],'k-', linewidth=2, label='BC')
plt.xlabel ('Month',fontsize=14)
plt.ylabel('Temperature (deg C)',fontsize=14)
plt.title('Bias-Corrected Average Monthly Min Temperature', fontsize=16)
#plt.legend(bbox_to_anchor=(0.9, -0.15), ncol=3)
plt.legend()
plt.tick_params(labelsize=12)
plt.grid(which='both')
plt.xlim(1,12);
plt.xticks(wy_index, month_strings);

# Precip
fig, ax=plt.subplots(1,1,figsize=(6,4))
plt.plot(wy_index,  days_per_month*liv2015_ds.sel(time=slice(start_date_bc50, end_date_bc50),station=stations_plot).precip.groupby('time.month').mean('time').mean('station').values[wy_numbers],'b-',linewidth=2, label='Liv')
plt.plot(wy_index,  days_per_month*wrf2014_ds.sel(station=stations_plot).precip.groupby('time.month').mean('time').mean('station').values[wy_numbers],'r-', linewidth=2, label='WRF')
plt.plot(wy_index,  days_per_month*wrf2014_lecorr_bc1_precip.sel(time=slice(start_date_bc50, end_date_bc50),station=stations_plot).groupby('time.month').mean('time').mean('station').values[wy_numbers],'k-', linewidth=2, label='BC')
plt.plot([1, 12],[0, 0], 'k--',linewidth=1)

plt.xlabel ('Month',fontsize=14)
plt.ylabel('Precipitation (mm)',fontsize=14)
plt.title('Bias-Corrected Average Monthly Precipitation', fontsize=16)
plt.legend(loc='best')
plt.tick_params(labelsize=12)
plt.grid(which='both')
plt.xlim(1,12);
plt.xticks(wy_index, month_strings);

#%% Low elevation Bias Correction- PRISM, 2015

# dates for 1950 bias corrections
start_date_bc50=datetime.date(1950,1,1)
end_date_bc50=datetime.date(2010,12,31)
# dates for 1915 bias corrections
start_date_bc15=datetime.date(1915,1,1)
end_date_bc15=datetime.date(2011,12,31)

min_elev_le_bc1=met_df.elevation.min() # default is min elevation: met_df.elevation.min()
max_elev_le_bc1=500 # default is max elevation: met_df.elevation.max()

min_elev_all_bc1=met_df.elevation.min() # default is min elevation: met_df.elevation.min()
max_elev_all_bc1=met_df.elevation.max() # default is max elevation: met_df.elevation.max()

# Extract data of interest for each step
stations_info_le_bc1=met_df[(met_df.elevation >= min_elev_le_bc1) & (met_df.elevation <= max_elev_le_bc1)] # Find points in elevation range of interest
elev_min_station_le_bc1=stations_info_le_bc1.elevation.idxmin() # station of minimum elevation in analysis
elev_min_le_bc1=stations_info_le_bc1.elevation.min() # minimum elevaiton of stations in analysis
elev_max_station_le_bc1=stations_info_le_bc1.elevation.idxmax() # station of maximum elevation in analysis
elev_max_le_bc1=stations_info_le_bc1.elevation.max() # maximum elevaiton of stations in analysis
stations_list_le_bc1=stations_info_le_bc1.station.values.astype('int64')
for name in ['liv2015', 'wrf2014']:
    ds=xr.Dataset()
    for var in  ['tmax', 'tmin', 'precip']:
         ds[var]=eval(name+'_ds')[var].sel(time=slice(start_date_bc50, end_date_bc50), station=stations_list_le_bc1)
    globals()[name+'_ds_bc1']=ds
    print (name+'_ds_bc1')


# Make corrections
# Temperature
corr_le_bc1_diff=compute_diff('liv2015_ds_bc1', 'wrf2014_ds_bc1', ['tmax','tmin'])

mean_corr_le_bc1_tmax=corr_le_bc1_diff.tmax.groupby('month').mean()
mean_corr_le_bc1_tmin=corr_le_bc1_diff.tmin.groupby('month').mean()

wrf2014_lecorr_bc1_tmax=wrf2014_ds.tmax.groupby('time.month')+mean_corr_le_bc1_tmax
wrf2014_lecorr_bc1_tmin=wrf2014_ds.tmin.groupby('time.month')+mean_corr_le_bc1_tmin

corr_all_bc1_diff_tmax=wrf2014_lecorr_bc1_tmax.groupby('time.month').mean('time')-liv2013_ds.tmax.sel(time=slice(start_date_bc50, end_date_bc50)).groupby('time.month').mean('time')
corr_all_bc1_diff_tmin=wrf2014_lecorr_bc1_tmin.groupby('time.month').mean('time')-liv2013_ds.tmin.sel(time=slice(start_date_bc50, end_date_bc50)).groupby('time.month').mean('time')

liv2013_ds_bc1_tmax=liv2013_ds.tmax.sel(time=slice(start_date_bc15, end_date_bc15)).groupby('time.month')+corr_all_bc1_diff_tmax
liv2013_ds_bc1_tmin=liv2013_ds.tmin.sel(time=slice(start_date_bc15, end_date_bc15)).groupby('time.month')+corr_all_bc1_diff_tmin
#%%
# Precip
days_per_month_ds=xr.DataArray(days_per_month,dims=['month'])
corr_le_bc1_ratio=prism_ds.sel(station=stations_list_le_bc1)/(days_per_month_ds*wrf2014_ds_bc1.groupby('time.month').mean('time').precip)
mean_corr_le_bc1_precip=corr_le_bc1_ratio.groupby('month').mean()
mean_corr_le_bc1_precip
#%%

wrf2014_lecorr_bc1_precip=wrf2014_ds.precip.groupby('time.month')+wrf2014_ds.precip.groupby('time.month').mean('time')*(mean_corr_le_bc1_precip-1)  # replace mean of all stations
corr_all_bc1_diff_precip=wrf2014_lecorr_bc1_precip.groupby('time.month').mean('time')/liv2013_ds.precip.sel(time=slice(start_date_bc50, end_date_bc50)).groupby('time.month').mean('time')

# Compute mean ratio for all stations
liv2013_ds_bc1_precip=liv2013_ds.precip.sel(time=slice(start_date_bc15, end_date_bc15)).groupby('time.month')+liv2013_ds.precip.sel(time=slice(start_date_bc15, end_date_bc15)).groupby('time.month').mean('time')*(corr_all_bc1_diff_precip-1)

#%%
# PLOTS
stations_plot=station_geom_mills.index # stations_list_le_bc1 OR station_geom_mills.index
# Tmax
fig, ax=plt.subplots(1,1,figsize=(6, 4))
plt.plot(wy_index, liv2015_ds.sel(time=slice(start_date_bc50, end_date_bc50),station=stations_plot).tmax.groupby('time.month').mean('time').mean('station').values[wy_numbers],'b-',linewidth=2, label='Liv')
plt.plot(wy_index, wrf2014_ds.sel(station=stations_plot).tmax.groupby('time.month').mean('time').mean('station').values[wy_numbers],'r-', linewidth=2, label='WRF')
plt.plot(wy_index, liv2013_ds_bc1_tmax.sel(time=slice(start_date_bc50, end_date_bc50),station=stations_plot).groupby('time.month').mean('time').mean('station').values[wy_numbers],'k-', linewidth=2, label='BC')
plt.xlabel ('Month',fontsize=14)
plt.ylabel('Temperature (deg C)',fontsize=14)
plt.title('Bias-Corrected Average Monthly Max Temperature', fontsize=16)
#plt.legend(bbox_to_anchor=(0.9, -0.15), ncol=3)
plt.legend()
plt.tick_params(labelsize=12)
plt.grid(which='both')
plt.xlim(1,12);
plt.xticks(wy_index, month_strings);

# Tmin
fig, ax=plt.subplots(1,1,figsize=(6, 4))
plt.plot(wy_index, liv2015_ds.sel(time=slice(start_date_bc50, end_date_bc50),station=stations_plot).tmin.groupby('time.month').mean('time').mean('station').values[wy_numbers],'b-',linewidth=2, label='Liv')
plt.plot(wy_index, wrf2014_ds.sel(station=stations_plot).tmin.groupby('time.month').mean('time').mean('station').values[wy_numbers],'r-', linewidth=2, label='WRF')
plt.plot(wy_index, liv2013_ds_bc1_tmin.sel(time=slice(start_date_bc50, end_date_bc50),station=stations_plot).groupby('time.month').mean('time').mean('station').values[wy_numbers],'k-', linewidth=2, label='BC')
plt.xlabel ('Month',fontsize=14)
plt.ylabel('Temperature (deg C)',fontsize=14)
plt.title('Bias-Corrected Average Monthly Min Temperature', fontsize=16)
#plt.legend(bbox_to_anchor=(0.9, -0.15), ncol=3)
plt.legend()
plt.tick_params(labelsize=12)
plt.grid(which='both')
plt.xlim(1,12);
plt.xticks(wy_index, month_strings);

# Precip
fig, ax=plt.subplots(1,1,figsize=(6,4))
plt.plot(wy_index, days_per_month*liv2015_ds.sel(time=slice(start_date_bc50, end_date_bc50),station=stations_plot).precip.groupby('time.month').mean('time').mean('station').values[wy_numbers],'b-',linewidth=2, label='Liv')
plt.plot(wy_index, days_per_month*wrf2014_ds.sel(station=stations_plot).precip.groupby('time.month').mean('time').mean('station').values[wy_numbers],'r-', linewidth=2, label='WRF')
plt.plot(wy_index, days_per_month*wrf2014_lecorr_bc1_precip.sel(time=slice(start_date_bc50, end_date_bc50),station=stations_plot).groupby('time.month').mean('time').mean('station').values[wy_numbers],'k-', linewidth=2, label='BC')
plt.plot([1, 12],[0, 0], 'k--',linewidth=1)

plt.xlabel ('Month',fontsize=14)
plt.ylabel('Precipitation (mm)',fontsize=14)
plt.title('Bias-Corrected Average Monthly Precipitation', fontsize=16)
plt.legend(loc='best')
plt.tick_params(labelsize=12)
plt.grid(which='both')
plt.xlim(1,12);
plt.xticks(wy_index, month_strings);
#%% Switch up VIC soils
soil_dir='D:\GitHub\Elwha_Landlab\dhsvm\met\VIC\soil'
soil_input_file='soil_elwha'
soil_output_file='soil_elwha'

os.chdir(soil_dir)

#Read in table of VIC soil inputs -- assumes all Lat/Long set to zero
soil_base = pd.read_table(soil_input_file,header=None)

#Read in mappingfile from TreatGeoSelf()
maptable = pd.read_table(mappingfile,sep=",")

#Make a list of all lat/long values
latlong=soil_base.apply(lambda x:tuple([x[2],x[3]]), axis=1)

#Make a list Lat/Long files that need to switched up 
latlong_1=maptable.apply(lambda x:tuple([x[1],x[2]]), axis=1)

#Switch up from 0 to 1 so VIC will run for this Lat/Long point - print new output file (VIC model input file)
soil_base[0] = latlong.apply(lambda x: 1 if x in set(latlong_1) else 0)        
soil_base.to_csv(soil_output_file, header=False,index=False,sep="\t")
print(str(soil_base[0].sum()) +' VIC grid cells have successfully been switched up.') 
print('Check your home directory for your new VIC soil model input set to your list of Lat/Long grid centroids.')
#%%
os.chdir(homedir)
