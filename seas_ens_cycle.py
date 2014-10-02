import iris
import numpy as np
import iris.coord_categorisation
import iris.plot as iplt
import iris.quickplot as qplt
import cartopy.crs as ccrs
import matplotlib.pyplot as plt


temp = iris.load_cube('ncfiles/temp.sfc.4ysl.nc')
iris.coord_categorisation.add_month_number(temp, 't', 'month_number')
temp_mean = temp[:,0,::].collapsed('t',iris.analysis.MEAN)
temp_anom = temp-temp_mean
mon_mean = temp_anom.aggregated_by('month_number', iris.analysis.MEAN)

temp_plv = iris.load_cube('ncfiles/temp.plv.4ysl.nc')
iris.coord_categorisation.add_month_number(temp_plv, 't', 'month_number')
temp_plv_mean = temp_plv[:,0,::].collapsed('t',iris.analysis.MEAN)
temp_plv_anom = temp_plv-temp_plv_mean
mon_mean_plv = temp_plv_anom.aggregated_by('month_number', iris.analysis.MEAN)



def bc(aux_cube, dim_cube, comparable_coord):
    if type(aux_cube.coord(comparable_coord)) is not iris.coords.AuxCoord:
        raise TypeError
    if type(dim_cube.coord(comparable_coord)) is not iris.coords.AuxCoord:
        raise TypeError
    
    aux_cube_coord = aux_cube.coord(comparable_coord)
    dim_cube_coord = dim_cube.coord(comparable_coord)
    aux_cube_dim, = aux_cube.coord_dims(aux_cube_coord)
    dim_cube_dim, = aux_cube.coord_dims(dim_cube_coord)
    
    s_aux_cube = [slice(None)]*len(aux_cube.shape)
    s_aux_cube[aux_cube_dim] = 0
    s_dim_cube = [slice(None)]*len(dim_cube.shape)
    s_dim_cube[dim_cube_dim] = 0
    a = aux_cube[tuple(s_aux_cube)]
    a.attributes = None
    a.cell_methods = None
    b = dim_cube[tuple(s_dim_cube)]
    b.attributes = None
    b.cell_methods = None
    if not a.is_compatible(b):
        iris.util.describe_diff(a, b)
        raise RuntimeError("Cubes are not compatible")
    
    ind = []
    for p in aux_cube.coord(comparable_coord).points:
        i = np.where(dim_cube.coord(comparable_coord).points == p)
        ind.append(i[0][0])

    s = [slice(None)]*len(dim_cube.shape)
    s[dim_cube_dim] = ind
    new_data = dim_cube.data[tuple(s)]
    new_cube = aux_cube.copy()
    new_cube.data = new_data
    new_cube.history = "%s comparable to %s in terms of %s" % (dim_cube.name(),
                                                               aux_cube.name(),
                                                               comparable_coord)
    
    return new_cube

# <codecell>
def remove_seascyc(cube, time_name='t'):
    iris.coord_categorisation.add_month_number(cube, 't', 'month_number')
    cube_mean = cube[:,0,::].collapsed('t',iris.analysis.MEAN)
    cube_anom = cube-cube_mean
    cube_mon_mean = cube_anom.aggregated_by('month_number', iris.analysis.MEAN)
    seasonal_cycle = bc(cube, cube_mon_mean, 'month_number')
    cube_rsc = cube_anom - seasonal_cycle
    return cube_rsc

def enscyc_ag(cube):
    ens = np.tile(np.linspace(1,48,48),24)
    #trim cube
    cube = cube[0:ens.shape[0]]
    cube.coord('month_number').long_name = '48_months'
    cube.coord('48_months').points = ens

    m48 = temp_rsc.aggregated_by('48_months',iris.analysis.MEAN)
    
seasonal_cycle = bc(temp_anom, mon_mean, 'month_number')
temp_rsc = temp_anom - seasonal_cycle

seasonal_cycle_plv = bc(temp_plv_anom, mon_mean_plv, 'month_number')
temp_plv_rsc = temp_plv_anom - seasonal_cycle_plv

# alter the auxilliary coordinate to aggregate over ens period
ens = np.tile(np.linspace(1,48,48),24)

temp_rsc = temp_rsc[0:ens.shape[0]]
temp_rsc.coord('month_number').long_name = '48_months'
temp_rsc.coord('48_months').points = ens

m48 = temp_rsc.aggregated_by('48_months',iris.analysis.MEAN)

temp_plv_rsc = temp_plv_rsc[0:ens.shape[0]]
temp_plv_rsc.coord('month_number').long_name = '48_months'
temp_plv_rsc.coord('48_months').points = ens

m48_plv = temp_plv_rsc.aggregated_by('48_months',iris.analysis.MEAN)


