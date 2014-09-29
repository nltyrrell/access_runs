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

temp_plv_mean = temp_plv[:,:,::].collapsed('time',iris.analysis.MEAN)
temp_plv_anom = temp_plv-temp_plv_mean

# Define regions
temp_plv_anom.coord('latitude').guess_bounds()
temp_plv_anom.coord('longitude').guess_bounds()

# ---------- define sst/tland----------
lsmask = iris.load_cube(ncfile_path + 'lsmask.nc')[0,0,::]
landmask = ~(ma.make_mask(lsmask.data.copy()) + np.zeros(temp_plv.shape)).astype(bool) # mask sea, show land
seamask = (ma.make_mask(lsmask.data.copy()) + np.zeros(temp_plv.shape)).astype(bool) # mask land, show sea

Tocean = temp_plv_anom.copy()
Tland = temp_plv_anom.copy()
Tocean.data = ma.array(Tocean.data, mask=seamask)
Tland.data = ma.array(Tland.data, mask=landmask)
# --------------

tropics = iris.Constraint(latitude = lambda v: -30 <= v <= 30)
Tocean_trop = Tocean.extract(tropics)
grid_areas_trop = iris.analysis.cartography.area_weights(Tocean_trop)

Tocean_trop_mean = Tocean_trop.collapsed(['latitude', 'longitude'],
                           iris.analysis.MEAN,
                           weights=grid_areas_trop)


SthAm = iris.Constraint(longitude=lambda l: (290 <= l <= 315), latitude = lambda l: (-23 <= l <= 0))
Tplv_SthAm = temp_plv_anom.extract(SthAm)
grid_areas_SthAm = iris.analysis.cartography.area_weights(Tplv_SthAm)

Tplv_SthAm_armean = Tplv_SthAm.collapsed(['latitude', 'longitude'],
                           iris.analysis.MEAN,
                           weights=grid_areas_SthAm)


hold = {}
# temp[height (100), evfr (4), ssts (11), l/o col (0-Oc,1-Ln)
# col which is used for regresion, '0' for ocean, '1' for land
lo_col = 0
n_bands = 20
reghold = np.zeros((2,n_bands,temp_plv.coord('air_pressure').shape[0]))
TrTochold = np.zeros(temp_plv.coord('air_pressure').shape[0])
pressure = temp_plv.coord('air_pressure').points

def jet(i,nloops):
    N_jet = plt.cm.jet.N
    colornum = plt.cm.jet((N_jet*i/(nloops-1)))
    return colornum

plt.ion()
plt.clf()
il = 0
for pn, p in enumerate(temp_plv.coord('air_pressure')): 
    Tplv = Tocean_trop_mean.extract(iris.Constraint(air_pressure=p.points[0])).data
    TOsfc = Tocean_trop_mean.extract(iris.Constraint(air_pressure=1000)).data
    linreg_ln = stats.linregress(TOsfc,Tplv)
    TrTochold[pn] = linreg_ln[0]
lat_range = np.linspace(0,45,n_bands).astype(int)
for lat in xrange(0,40,40/n_bands):
    Tplv_NH = Tland.extract(iris.Constraint(longitude = lambda z: (0 <= z <= 360), latitude = lambda l: (lat <= l <= lat+5)))
    Tplv_SH = Tland.extract(iris.Constraint(latitude = lambda l: (-lat-5 <= l <= -lat)))
    Tplv_NHmean = Tplv_NH.collapsed(['latitude', 'longitude'],
                           iris.analysis.MEAN,
                           weights=iris.analysis.cartography.area_weights(Tplv_NH))
    Tplv_SHmean = Tplv_SH.collapsed(['latitude', 'longitude'],
                           iris.analysis.MEAN,
                           weights=iris.analysis.cartography.area_weights(Tplv_SH))
    
    for pn, p in enumerate(temp_plv.coord('air_pressure')): 

        TplvNH = Tplv_NHmean.extract(iris.Constraint(air_pressure=p.points[0])).data
        TplvSH = Tplv_SHmean.extract(iris.Constraint(air_pressure=p.points[0])).data
        TsfcNH = Tplv_NHmean.extract(iris.Constraint(air_pressure=1000)).data
        TsfcSH = Tplv_SHmean.extract(iris.Constraint(air_pressure=1000)).data
        TOsfc = Tocean_trop_mean.extract(iris.Constraint(air_pressure=1000)).data
        linreg_lnNH = stats.linregress(TsfcNH,TplvNH)
        linreg_lnSH = stats.linregress(TsfcSH,TplvSH)
        reghold[0,il,pn] = linreg_lnNH[0]
        reghold[1,il,pn] = linreg_lnSH[0]
        #linreg_ln = stats.linregress(Tsfc,Tplv)
        #SthAmhold[1,pn] = linreg_ln[0]

    plt.plot(reghold[0,il],pressure,color=jet(il,n_bands))
    #plt.plot(reghold[1,il],pressure,'--',color=jet(il,n_bands))
    il = il + 1
plt.legend((lat_range.astype(str)),loc=0)
#plt.title('Regression with surface temp and temp profile.') 
plt.xlim(-0.5,1.5)
plt.ylabel('z [hPa]')
plt.xlabel('Regresion coeff.')
plt.gca().invert_yaxis()
#plt.savefig('./figures/'+run_name+'reg_prof_ocensfc_e0.05.eps')




