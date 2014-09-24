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


ncfile_path = '/home/nicholat/project/mit_tcm/access_runs/ncfiles/'


reg_T_sfc_T_tropo = iris.load_cube(ncfile_path + 'lreg.4ysl.T_sfc.T_tropo.nc')
reg_T_sfc_RH_sfc = iris.load_cube(ncfile_path + 'lreg.4ysl.T_sfc.RH_sfc.nc')
reg_T_sfc_RH_tropo = iris.load_cube(ncfile_path + 'lreg.4ysl.T_sfc.RH_tropo.nc')
reg_T_sfc_Cld_High = iris.load_cube(ncfile_path + 'lreg.4ysl.T_sfc.Cld_High.nc')
reg_T_sfc_Cld_Low = iris.load_cube(ncfile_path + 'lreg.4ysl.T_sfc.Cld_Low.nc')
reg_T_tropo_RH_sfc = iris.load_cube(ncfile_path + 'lreg.4ysl.T_tropo.RH_sfc.nc')
reg_T_tropo_RH_tropo = iris.load_cube(ncfile_path + 'lreg.4ysl.T_tropo.RH_tropo.nc')
reg_T_sfc_Precip = iris.load_cube(ncfile_path + 'lreg.4ysl.T_sfc.Precip.nc')

cor_T_sfc_T_tropo = iris.load_cube(ncfile_path + 'cor.4ysl.T_sfc.T_tropo.nc')
cor_T_sfc_RH_sfc = iris.load_cube(ncfile_path + 'cor.4ysl.T_sfc.RH_sfc.nc')
cor_T_sfc_RH_tropo = iris.load_cube(ncfile_path + 'cor.4ysl.T_sfc.RH_tropo.nc')
cor_T_sfc_Cld_High = iris.load_cube(ncfile_path + 'cor.4ysl.T_sfc.Cld_High.nc')
cor_T_sfc_Cld_Low = iris.load_cube(ncfile_path + 'cor.4ysl.T_sfc.Cld_Low.nc')
cor_T_tropo_RH_sfc = iris.load_cube(ncfile_path + 'cor.4ysl.T_tropo.RH_sfc.nc')
cor_T_tropo_RH_tropo = iris.load_cube(ncfile_path + 'cor.4ysl.T_tropo.RH_tropo.nc')
cor_T_sfc_Precip = iris.load_cube(ncfile_path + 'cor.4ysl.T_sfc.Precip.nc')

plt.close('all')
plt.figure(1)
qplt.pcmeshclf(reg_T_sfc_T_tropo,vmin=-1,vmax=1,cmap=mc.jetwhite())
plt.savefig('figures/reg_T_sfc_T_tropo.png')
plt.figure(2)
qplt.pcmeshclf(reg_T_sfc_RH_sfc,vmin=-6,vmax=6,cmap=mc.jetwhite())
plt.savefig('figures/reg_T_sfc_RH_sfc.png')
plt.figure(3)
qplt.pcmeshclf(reg_T_sfc_RH_tropo,vmin=-6,vmax=6,cmap=mc.jetwhite())
plt.savefig('figures/reg_T_sfc_RH_tropo.png')
plt.figure(4)
qplt.pcmeshclf(reg_T_sfc_Cld_High,vmin=-5e-2,vmax=5e-2,cmap=mc.jetwhite())
plt.savefig('figures/reg_T_sfc_Cld_High.png')
plt.figure(5)
qplt.pcmeshclf(reg_T_sfc_Cld_Low,vmin=-5e-2,vmax=5e-2,cmap=mc.jetwhite())
plt.savefig('figures/reg_T_sfc_Cld_Low.png')
plt.figure(6)
qplt.pcmeshclf(reg_T_tropo_RH_sfc,vmin=-6,vmax=6,cmap=mc.jetwhite())
plt.savefig('figures/reg_T_tropo_RH_sfc.png')
plt.figure(7)
qplt.pcmeshclf(reg_T_tropo_RH_tropo,vmin=-6,vmax=6,cmap=mc.jetwhite())
plt.savefig('figures/reg_T_tropo_RH_tropo.png')
plt.figure(8)
qplt.pcmeshclf(reg_T_sfc_Precip,vmin=-1e-5,vmax=1e-5,cmap=mc.jetwhite())
plt.savefig('figures/reg_T_sfc_Precip.png')
# Plot correlations
plt.close('all')
plt.figure(1)
qplt.pcmeshclf(cor_T_sfc_T_tropo,vmin=-1,vmax=1,cmap=mc.jetwhite())
plt.savefig('figures/cor_T_sfc_T_tropo.png')
plt.figure(2)
qplt.pcmeshclf(cor_T_sfc_RH_sfc,vmin=-1,vmax=1,cmap=mc.jetwhite())
plt.savefig('figures/cor_T_sfc_RH_sfc.png')
plt.figure(3)
qplt.pcmeshclf(cor_T_sfc_RH_tropo,vmin=-1,vmax=1,cmap=mc.jetwhite())
plt.savefig('figures/cor_T_sfc_RH_tropo.png')
plt.figure(4)
qplt.pcmeshclf(cor_T_sfc_Cld_High,vmin=-1,vmax=1,cmap=mc.jetwhite())
plt.savefig('figures/cor_T_sfc_Cld_High.png')
plt.figure(5)
qplt.pcmeshclf(cor_T_sfc_Cld_Low,vmin=-1,vmax=1,cmap=mc.jetwhite())
plt.savefig('figures/cor_T_sfc_Cld_Low.png')
plt.figure(6)
qplt.pcmeshclf(cor_T_tropo_RH_sfc,vmin=-1,vmax=1,cmap=mc.jetwhite())
plt.savefig('figures/cor_T_tropo_RH_sfc.png')
plt.figure(7)
qplt.pcmeshclf(cor_T_tropo_RH_tropo,vmin=-1,vmax=1,cmap=mc.jetwhite())
plt.savefig('figures/cor_T_tropo_RH_tropo.png')
plt.figure(8)
qplt.pcmeshclf(cor_T_sfc_Precip,vmin=-1,vmax=1,cmap=mc.jetwhite())
plt.savefig('figures/cor_T_sfc_Precip.png')





