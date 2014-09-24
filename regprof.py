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
temp = iris.load_cube(ncfile_path + 'temp.sfc.4ysl.ym.nc')
temp_mean = temp[:,0,::].collapsed('time',iris.analysis.MEAN)

def linregcube(cube1,cube2,name1,name2,ncfile_path='/home/nicholat/project/mit_tcm/access_runs/ncfiles/',copy_cube=temp_mean):
    linreg_map = np.zeros(copy_cube.shape)
    cor_map = np.zeros(copy_cube.shape)
    print('Linreg/Cor map for '+name1+' and '+name2)
    
    for nlat, lat in enumerate(copy_cube.coord('latitude')):
        for nlon, lon in enumerate(copy_cube.coord('longitude')):
            # get the sfc temp timeseries at each lat, lon.
            var1 = cube1.extract(iris.Constraint(latitude=lat.points[0],longitude=lon.points[0])).data
            var2 = cube2.extract(iris.Constraint(latitude=lat.points[0],longitude=lon.points[0])).data

            linreg = stats.linregress(var1,var2)
            linreg_map[nlat,nlon] = linreg[0]
            cor_map[nlat,nlon] = linreg[2]

    #Copy the soil moisture cube then change everything, clean it up, save phi as an netcdf
    reg_cube = copy_cube.copy()
    reg_cube.data[:] = linreg_map
    reg_cube.long_name = 'Lin Regression '+ name1 +' '+ name2
    reg_cube.units = 'no_unit'
    reg_cube.attributes['title'] = 'Lin Regression '+ name1 +' '+ name2
    reg_cube.attributes['name'] = 'reg'
    reg_cube.remove_coord('surface')
    reg_cube.remove_coord('time')
    iris.save(reg_cube,ncfile_path+'lreg.4ysl.'+name1+'.'+name2+'.nc')
    #Copy the soil moisture cube then change everything, clean it up, save phi as an netcdf
    cor_cube = copy_cube.copy()
    cor_cube.data[:] = cor_map
    cor_cube.long_name = 'Correlation '+ name1 +' '+ name2
    cor_cube.units = 'no_unit'
    cor_cube.attributes['title'] = 'Correlation '+ name1 +' '+ name2
    cor_cube.attributes['name'] = 'r_val'
    cor_cube.remove_coord('surface')
    cor_cube.remove_coord('time')
    iris.save(cor_cube,ncfile_path+'cor.4ysl.'+name1+'.'+name2+'.nc')

    return reg_cube, cor_cube

# import the ACCESS data using iris
temp_plv = iris.load_cube(ncfile_path + 'temp.plv.4ysl.ym.nc')
rh_plv = iris.load_cube(ncfile_path + 'rhum.plv.4ysl.ym.nc')
cld_thlev = iris.load_cube(ncfile_path + 'cld.thlev.4ysl.ym.nc')
# change name of height coord to a standard name
cld_thlev.coord('Hybrid height').standard_name = 'atmosphere_hybrid_height_coordinate'

# Calcualate the anomalies
temp_anom = temp[:,0,::]-temp_mean
T_sfc = temp_anom

temp_plv_mean = temp_plv[:,:,::].collapsed('time',iris.analysis.MEAN)
temp_plv_anom = temp_plv-temp_plv_mean
rh_plv_mean = rh_plv[:,:,::].collapsed('time',iris.analysis.MEAN)
rh_plv_anom = rh_plv-rh_plv_mean
rh_1000_anom    = rh_plv_anom.extract(iris.Constraint(air_pressure=1000))
RH_sfc = rh_1000_anom

# Cloud are fractions
cld_mean = cld_thlev[:,:,::].collapsed('time',iris.analysis.MEAN)
cld_anom = cld_thlev - cld_mean
high = iris.Constraint(atmosphere_hybrid_height_coordinate = lambda h: 5000 <= h <= 15000)
low = iris.Constraint(atmosphere_hybrid_height_coordinate = lambda h: 0 <= h <= 5000)
high_cld = cld_anom.extract(high).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN)
low_cld = cld_anom.extract(low).collapsed('atmosphere_hybrid_height_coordinate',iris.analysis.MEAN)
Cld_High = high_cld
Cld_Low = low_cld

# Weighted average
uptrop = iris.Constraint(air_pressure = lambda p: 200 <= p <= 500)
temp_uptrop = temp_plv_anom.extract(uptrop)
trop_weights = ta.troposweights(startlev=4,endlev=10)
br_weight = iris.util.broadcast_to_shape(trop_weights,temp_uptrop.shape,(1,))

t_weights = temp_uptrop.collapsed('air_pressure',
                           iris.analysis.MEAN,
                           weights=br_weight)

T_tropo = t_weights

rh_uptrop = rh_plv_anom.extract(uptrop)
rh_weights = rh_uptrop.collapsed('air_pressure',
                           iris.analysis.MEAN,
                           weights=br_weight)
RH_tropo = rh_weights

precip   = iris.load_cube(ncfile_path + 'precip.4ysl.ym.nc')
precip_mean= precip[:,:,::].collapsed('time',iris.analysis.MEAN)
precip_anom = precip-precip_mean
Precip = precip_anom[:,0,::]

linregcube(T_sfc,T_tropo,'T_sfc','T_tropo')
linregcube(T_sfc,RH_sfc,'T_sfc','RH_sfc')
linregcube(T_sfc,RH_tropo,'T_sfc','RH_tropo')
linregcube(T_sfc,Cld_High,'T_sfc','Cld_High')
linregcube(T_sfc,Cld_Low,'T_sfc','Cld_Low')
linregcube(T_sfc,Precip,'T_sfc','Precip')
linregcube(T_tropo,RH_sfc,'T_tropo','RH_sfc')
linregcube(T_tropo,RH_tropo,'T_tropo','RH_tropo')






