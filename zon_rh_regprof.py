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

# import the ACCESS data using iris
temp_plv = iris.load_cube(ncfile_path + 'temp.plv.4ysl.ym.nc')
rh_plv = iris.load_cube(ncfile_path + 'rhum.plv.4ysl.ym.nc')

rh_plv_mean = rh_plv[:,:,::].collapsed('time',iris.analysis.MEAN)
rh_plv_anom = rh_plv-rh_plv_mean
temp_plv_mean = temp_plv[:,:,::].collapsed('time',iris.analysis.MEAN)
temp_plv_anom = temp_plv-temp_plv_mean

# Define regions
rh_plv_anom.coord('latitude').guess_bounds()
rh_plv_anom.coord('longitude').guess_bounds()
temp_plv_anom.coord('latitude').guess_bounds()
temp_plv_anom.coord('longitude').guess_bounds()

# ---------- define sst/tland----------
lsmask = iris.load_cube(ncfile_path + 'lsmask.nc')[0,0,::]
landmask = ~(ma.make_mask(lsmask.data.copy()) + np.zeros(rh_plv.shape)).astype(bool) # mask sea, show land
seamask = (ma.make_mask(lsmask.data.copy()) + np.zeros(rh_plv.shape)).astype(bool) # mask land, show sea

RHocean = rh_plv_anom.copy()
Tocean = temp_plv_anom.copy()
RHland = rh_plv_anom.copy()
Tland = temp_plv_anom.copy()
RHocean.data = ma.array(RHocean.data, mask=seamask)
Tocean.data = ma.array(Tocean.data, mask=seamask)
RHland.data = ma.array(RHland.data, mask=landmask)
Tland.data = ma.array(Tland.data, mask=landmask)
# --------------

tropics = iris.Constraint(latitude = lambda v: -30 <= v <= 30)
Tocean_trop = Tocean.extract(tropics)
grid_areas_trop = iris.analysis.cartography.area_weights(Tocean_trop)

Tocean_trop_mean = Tocean_trop.collapsed(['latitude', 'longitude'],
                           iris.analysis.MEAN,
                           weights=grid_areas_trop)


SthAm = iris.Constraint(longitude=lambda l: (290 <= l <= 315), latitude = lambda l: (-23 <= l <= 0))
RHplv_SthAm = rh_plv_anom.extract(SthAm)
grid_areas_SthAm = iris.analysis.cartography.area_weights(RHplv_SthAm)

RHplv_SthAm_armean = RHplv_SthAm.collapsed(['latitude', 'longitude'],
                           iris.analysis.MEAN,
                           weights=grid_areas_SthAm)


hold = {}
# rh[height (100), evfr (4), ssts (11), l/o col (0-Oc,1-Ln)
# col which is used for regresion, '0' for ocean, '1' for land
lo_col = 0
n_bands = 13
reghold = np.zeros((2,n_bands,rh_plv.coord('air_pressure').shape[0]))
TrRHochold = np.zeros(rh_plv.coord('air_pressure').shape[0])
pressure = rh_plv.coord('air_pressure').points

def jet(i,nloops):
    N_jet = plt.cm.jet.N
    colornum = plt.cm.jet_r((N_jet*i/(nloops-1)))
    return colornum
handle = {} #np.zeros((n_bands))
plt.ion()
plt.clf()
il = 0
# for pn, p in enumerate(rh_plv.coord('air_pressure')): 
#    RHplv = RHocean_trop_mean.extract(iris.Constraint(air_pressure=p.points[0])).data
#    TOsfc = Tocean_trop_mean.extract(iris.Constraint(air_pressure=1000)).data
#    linreg_ln = stats.linregress(TOsfc,RHplv)
#    TrRHochold[pn] = linreg_ln[0]
lat_range = np.linspace(0,65,n_bands).astype(int)
for lat in xrange(0,65,65/n_bands):
    RHplv_NH = RHland.extract(iris.Constraint(longitude = lambda z: (0 <= z <= 360), latitude = lambda l: (lat <= l <= lat+5)))
    RHplv_SH = RHland.extract(iris.Constraint(latitude = lambda l: (-lat-5 <= l <= -lat)))
    RHplv_NHmean = RHplv_NH.collapsed(['latitude', 'longitude'],
                           iris.analysis.MEAN,
                           weights=iris.analysis.cartography.area_weights(RHplv_NH))
    RHplv_SHmean = RHplv_SH.collapsed(['latitude', 'longitude'],
                           iris.analysis.MEAN,
                           weights=iris.analysis.cartography.area_weights(RHplv_SH))
    
    Tplv_NH = Tland.extract(iris.Constraint(longitude = lambda z: (0 <= z <= 360), latitude = lambda l: (lat <= l <= lat+5)))
    Tplv_SH = Tland.extract(iris.Constraint(latitude = lambda l: (-lat-5 <= l <= -lat)))
    Tplv_NHmean = Tplv_NH.collapsed(['latitude', 'longitude'],
                           iris.analysis.MEAN,
                           weights=iris.analysis.cartography.area_weights(Tplv_NH))
    Tplv_SHmean = Tplv_SH.collapsed(['latitude', 'longitude'],
                           iris.analysis.MEAN,
                           weights=iris.analysis.cartography.area_weights(Tplv_SH))
    
    for pn, p in enumerate(rh_plv.coord('air_pressure')): 

        RHplvNH = RHplv_NHmean.extract(iris.Constraint(air_pressure=p.points[0])).data
        RHplvSH = RHplv_SHmean.extract(iris.Constraint(air_pressure=p.points[0])).data
        RHsfcNH = RHplv_NHmean.extract(iris.Constraint(air_pressure=1000)).data
        RHsfcSH = RHplv_SHmean.extract(iris.Constraint(air_pressure=1000)).data
        TsfcNH = Tplv_NHmean.extract(iris.Constraint(air_pressure=1000)).data
        TsfcSH = Tplv_SHmean.extract(iris.Constraint(air_pressure=1000)).data
        TOsfc = Tocean_trop_mean.extract(iris.Constraint(air_pressure=1000)).data
        linreg_lnNH = stats.linregress(TOsfc,RHplvNH)
        linreg_lnSH = stats.linregress(TOsfc,RHplvSH)
        reghold[0,il,pn] = linreg_lnNH[0]
        reghold[1,il,pn] = linreg_lnSH[0]
        #linreg_ln = stats.linregress(RHsfc,RHplv)
        #SthAmhold[1,pn] = linreg_ln[0]

    plt.plot(reghold[0,il],pressure,color=jet(il,n_bands))
    plt.plot(reghold[1,il],pressure,'--',color=jet(il,n_bands),label="_nolegend_")
    il = il + 1
plt.legend((lat_range.astype(str)),loc=0)
#plt.title('Regression with surface rh and rh profile.') 
plt.xlim(-9.0,7.0)
plt.ylabel('z [hPa]')
plt.xlabel('Regresion coeff.')
plt.gca().invert_yaxis()
#plt.savefig('./figures/'+run_name+'reg_prof_ocensfc_e0.05.eps')




