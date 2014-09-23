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
cld_thlev = iris.load_cube(ncfile_path + 'cld.thlev.4ysl.ym.nc')
# change name of height coord to a standard name
cld_thlev.coord('Hybrid height').standard_name = 'atmosphere_hybrid_height_coordinate'

# Calcualate the anomalies
temp_mean = temp[:,0,::].collapsed('time',iris.analysis.MEAN)
temp_anom = temp[:,0,::]-temp_mean

temp_plv_mean = temp_plv[:,:,::].collapsed('time',iris.analysis.MEAN)
temp_plv_anom = temp_plv-temp_plv_mean
rh_plv_mean = rh_plv[:,:,::].collapsed('time',iris.analysis.MEAN)
rh_plv_anom = rh_plv-rh_plv_mean
rh_1000    = rh_plv_anom.extract(iris.Constraint(air_pressure=1000))
rh_925    = rh_plv_anom.extract(iris.Constraint(air_pressure=925))
rh_850    = rh_plv_anom.extract(iris.Constraint(air_pressure=850))

#plev = 300
#temp_300    = temp_plv_anom.extract(iris.Constraint(air_pressure=plev))
#temp_300_mean = temp_300.collapsed('time',iris.analysis.MEAN)
#temp_300_anom = temp_300-temp_300_mean

# Cloud are fractions
cld_mean = cld_thlev[:,:,::].collapsed('time',iris.analysis.MEAN)
cld_anom = cld_thlev - cld_mean
high = iris.Constraint(atmosphere_hybrid_height_coordinate = lambda h: 5000 <= h <= 15000)
low = iris.Constraint(atmosphere_hybrid_height_coordinate = lambda h: 0 <= h <= 5000)
high_cld = cld_anom.extract(high).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN)
low_cld = cld_anom.extract(low).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN)

# Weighted average
uptrop = iris.Constraint(air_pressure = lambda p: 200 <= p <= 500)
temp_uptrop = temp_plv_anom.extract(uptrop)

trop_weights = ta.troposweights(startlev=4,endlev=10)
bweight = iris.util.broadcast_to_shape(trop_weights,temp_uptrop.shape,(1,))

t_weights = temp_uptrop.collapsed('air_pressure',
                           iris.analysis.MEAN,
                           weights=bweight)

linreg_tstt = np.zeros(temp_mean.shape)
linreg_tsrh = np.zeros(temp_mean.shape)
linreg_ttrh = np.zeros(temp_mean.shape)

linreg_tshc = np.zeros(temp_mean.shape)
linreg_tslc = np.zeros(temp_mean.shape)

for nlat, lat in enumerate(temp_anom.coord('latitude')):
    for nlon, lon in enumerate(temp_anom.coord('longitude')):
        # get the sfc temp timeseries at each lat, lon.
        tsurf  = temp_anom.extract(iris.Constraint(latitude=lat.points[0],longitude=lon.points[0])).data
        ttropo = t_weights.extract(iris.Constraint(latitude=lat.points[0],longitude=lon.points[0])).data

        rh_1 = rh_1000.extract(iris.Constraint(latitude=lat.points[0],longitude=lon.points[0])).data
        hcld = high_cld.extract(iris.Constraint(latitude=lat.points[0],longitude=lon.points[0])).data
        lcld = low_cld.extract(iris.Constraint(latitude=lat.points[0],longitude=lon.points[0])).data
        #rh_9 = rh_925.extract(iris.Constraint(latitude=lat.points[0],longitude=lon.points[0])).data
        #rh_8 = rh_850.extract(iris.Constraint(latitude=lat.points[0],longitude=lon.points[0])).data

        # Use ev_idx and lat_idx to get relevent phi value in RCM array
        linreg = stats.linregress(tsurf,ttropo)
        linreg_map[nlat,nlon] = linreg[0]
        cor_map[nlat,nlon] = linreg[2]
        linreg = stats.linregress(tsurf,rh_1)
        linreg_tsrh[nlat,nlon] = linreg[0]
        cor_tsrh[nlat,nlon] = linreg[2]
        linreg = stats.linregress(ttropo,rh_1)
        linreg_ttrh[nlat,nlon] = linreg[0]
        cor_ttrh[nlat,nlon] = linreg[2]
        linreg = stats.linregress(tsurf,hcld)
        linreg_tshc[nlat,nlon] = linreg[0]
        cor_tshc[nlat,nlon] = linreg[2]
        linreg = stats.linregress(tsurf,lcld)
        linreg_tslc[nlat,nlon] = linreg[0]
        cor_tslc[nlat,nlon] = linreg[2]
        #linreg_rh9[nlat,nlon] = stats.linregress(tsurf,rh_9)[0]
        #linreg_rh8[nlat,nlon] = stats.linregress(tsurf,rh_8)[0]
        #print '-------------------'

# Get a lsmask, mask out sea values in phi_array
lsmask = iris.load_cube(ncfile_path + 'lsmask.nc')[0,0,::]
landmask = ~(ma.make_mask(lsmask.data.copy()) + np.zeros(temp_mean.shape)).astype(bool) # mask sea, show land
linreg_mask = ma.array(linreg_map, mask = landmask)

#Copy the soil moisture cube then change everything, clean it up, save phi as an netcdf
reg_cube = temp_mean.copy()
reg_cube.data[:] = linreg_mask
reg_cube.long_name = 'Lin Regression Tsurf Ttropos'
reg_cube.units = 'no_unit'
reg_cube.attributes['title'] = 'Lin Regression Tsurf Ttropos'
reg_cube.attributes['name'] = 'reg'
reg_cube.remove_coord('surface')
reg_cube.remove_coord('time')
iris.save(reg_cube,ncfile_path+'lreg.4ysl.tsfc.ttropo.nc')

reg_cube_tsrh = temp_mean.copy()
reg_cube_tsrh.data[:] = linreg_tsrh
reg_cube_tsrh.long_name = 'Lin Regression Tsurf RHsurf'
reg_cube_tsrh.units = 'no_unit'
reg_cube_tsrh.attributes['title'] = 'Lin Regression Tsurf RHsurf'
reg_cube_tsrh.attributes['name'] = 'reg'
reg_cube_tsrh.remove_coord('surface')
reg_cube_tsrh.remove_coord('time')
#iris.save(reg_cube,ncfile_path+'lreg.4ysl.tsfc.ttropo.nc')

reg_cube_ttrh = temp_mean.copy()
reg_cube_ttrh.data[:] = linreg_ttrh
reg_cube_ttrh.long_name = 'Lin Regression Ttropos RHsurf '
reg_cube_ttrh.units = 'no_unit'
reg_cube_ttrh.attributes['title'] = 'Lin Regression Ttropos RHsurf '
reg_cube_ttrh.attributes['name'] = 'reg'
reg_cube_ttrh.remove_coord('surface')
reg_cube_ttrh.remove_coord('time')
#iris.save(reg_cube,ncfile_path+'lreg.4ysl.tsfc.ttropo.nc')

reg_cube_tshc = temp_mean.copy()
reg_cube_tshc.data[:] = linreg_tshc
reg_cube_tshc.long_name = 'Lin Regression Tsurf High Cloud'
reg_cube_tshc.units = 'no_unit'
reg_cube_tshc.attributes['title'] = 'Lin Regression Tsurf High Cloud'
reg_cube_tshc.attributes['name'] = 'reg'
reg_cube_tshc.remove_coord('surface')
reg_cube_tshc.remove_coord('time')
#iris.save(reg_cube,ncfile_path+'lreg.4ysl.tsfc.ttropo.nc')

reg_cube_tslc = temp_mean.copy()
reg_cube_tslc.data[:] = linreg_tslc
reg_cube_tslc.long_name = 'Lin Regression Tsurf Low Cloud'
reg_cube_tslc.units = 'no_unit'
reg_cube_tslc.attributes['title'] = 'Lin Regression Tsurf Low Cloud'
reg_cube_tslc.attributes['name'] = 'reg'
reg_cube_tslc.remove_coord('surface')
reg_cube_tslc.remove_coord('time')
#iris.save(reg_cube,ncfile_path+'lreg.4ysl.tsfc.ttropo.nc')

plt.figure(1)
qplt.pcmeshclf(reg_cube,vmin=-1,vmax=1,cmap=mc.jetwhite())
plt.figure(2)
qplt.pcmeshclf(reg_cube_tsrh,vmin=-10,vmax=10,cmap=mc.jetwhite())
plt.figure(3)
qplt.pcmeshclf(reg_cube_ttrh,vmin=-1,vmax=1,cmap=mc.jetwhite())
plt.figure(4)
qplt.pcmeshclf(reg_cube_tshc,vmin=-1,vmax=1,cmap=mc.jetwhite())
plt.figure(5)
qplt.pcmeshclf(reg_cube_tslc,vmin=-1,vmax=1,cmap=mc.jetwhite())

