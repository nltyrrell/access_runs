import numpy as np
import numpy.ma as ma
import iris as iris
import iris.plot as iplt
import iris.quickplot as qplt
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import mycmaps as mc
import scipy.stats as stats
import sys
import troposave as ta


# Load in Tsfc, Ttropos (or T300hPa?)
# calculate the regression between the timeseries
# Tsfc = alpha * Ttropo + epsilon

ncfile_path = '/home/nicholat/project/mit_tcm/access_runs/ncfiles/'

# import the ACCESS data using iris
temp = iris.load_cube(ncfile_path + 'temp.sfc.4ysl.ym.nc')
temp_plv = iris.load_cube(ncfile_path + 'temp.plv.4ysl.ym.nc')
rh_plv = iris.load_cube(ncfile_path + 'rhum.plv.4ysl.ym.nc')

# Calcualate the anomalies
temp_mean = temp[:,0,::].collapsed('time',iris.analysis.MEAN)
temp_anom = temp[:,0,::]-temp_mean

temp_plv_mean = temp_plv[:,:,::].collapsed('time',iris.analysis.MEAN)
temp_plv_anom = temp_plv-temp_plv_mean
rh_plv_mean = rh_plv[:,:,::].collapsed('time',iris.analysis.MEAN)
rh_plv_anom = rh_plv-rh_plv_mean

plev = 300
temp_tropo    = temp_plv_anom.extract(iris.Constraint(air_pressure=plev))
temp_tropo_mean = temp_tropo.collapsed('time',iris.analysis.MEAN)
temp_tropo_anom = temp_tropo-temp_tropo_mean

# Weighted average
uptrop = iris.Constraint(air_pressure = lambda p: 200 <= p <= 500)
temp_uptrop = temp_plv_anom.extract(uptrop)

trop_weights = ta.troposweights(startlev=4,endlev=10)
bweight = iris.util.broadcast_to_shape(trop_weights,temp_uptrop.shape,(1,))

t_weights = temp_uptrop.collapsed('air_pressure',
                           iris.analysis.MEAN,
                           weights=bweight)
t_noweights = temp_uptrop.collapsed('air_pressure',
                           iris.analysis.MEAN)

tdata = temp_plv_anom.data.copy()
t_myweights = ta.troposave(tdata,startlev=4,endlev=10)

linreg_map = np.zeros(temp_mean.shape)
linreg_map_nw = np.zeros(temp_mean.shape)
linreg_map_mw = np.zeros(temp_mean.shape)

for nlat, lat in enumerate(temp_anom.coord('latitude')):
    for nlon, lon in enumerate(temp_anom.coord('longitude')):
        # get the sfc temp timeseries at each lat, lon.
        tsurf  = temp_anom.extract(iris.Constraint(latitude=lat.points[0],longitude=lon.points[0])).data
        ttropo = t_weights.extract(iris.Constraint(latitude=lat.points[0],longitude=lon.points[0])).data
        ttropo_nw = t_noweights.extract(iris.Constraint(latitude=lat.points[0],longitude=lon.points[0])).data
        ttropo_mw = t_myweights[:,nlat,nlon]

        # Use ev_idx and lat_idx to get relevent phi value in RCM array
        linreg_map[nlat,nlon] = stats.linregress(tsurf,ttropo)[0]
        linreg_map_nw[nlat,nlon] = stats.linregress(tsurf,ttropo_nw)[0]
        linreg_map_mw[nlat,nlon] = stats.linregress(tsurf,ttropo_mw)[0]
        #print '-------------------'

# Get a lsmask, mask out sea values in phi_array
lsmask = iris.load_cube(ncfile_path + 'lsmask.nc')[0,0,::]
landmask = ~(ma.make_mask(lsmask.data.copy()) + np.zeros(temp_mean.shape)).astype(bool) # mask sea, show land
linreg_mask = ma.array(linreg_map, mask = landmask)
linreg_mask_nw = ma.array(linreg_map_nw, mask = landmask)
linreg_mask_mw = ma.array(linreg_map_mw, mask = landmask)

#Copy the soil moisture cube then change everything, clean it up, save phi as an netcdf
reg_cube = temp_mean.copy()
reg_cube.data[:] = linreg_mask
reg_cube.long_name = 'Lin Regression'
reg_cube.units = 'no_unit'
reg_cube.attributes['title'] = 'Lin Regression'
reg_cube.attributes['name'] = 'reg'
reg_cube.remove_coord('surface')
reg_cube.remove_coord('time')
iris.save(reg_cube,ncfile_path+'lreg.4ysl.tsfc.ttropo.nc')

reg_cube_nw = temp_mean.copy()
reg_cube_nw.data[:] = linreg_mask_nw
reg_cube_nw.long_name = 'Lin Regression'
reg_cube_nw.units = 'no_unit'
reg_cube_nw.attributes['title'] = 'Lin Regression'
reg_cube_nw.attributes['name'] = 'reg'
reg_cube_nw.remove_coord('surface')
reg_cube_nw.remove_coord('time')

reg_cube_mw = temp_mean.copy()
reg_cube_mw.data[:] = linreg_mask_mw
reg_cube_mw.long_name = 'Lin Regression'
reg_cube_mw.units = 'no_unit'
reg_cube_mw.attributes['title'] = 'Lin Regression'
reg_cube_mw.attributes['name'] = 'reg'
reg_cube_mw.remove_coord('surface')
reg_cube_mw.remove_coord('time')

qplt.pcmeshclf(reg_cube,vmin=-1,vmax=1,cmap=mc.jetwhite())

