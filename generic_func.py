import numpy as np
from netCDF4 import Dataset
import matplotlib
import matplotlib.pylab as plt

matplotlib.rcParams.update({'font.size': 9})

def find_ls_idx(ls_arr, ls):
    idx = (np.abs(ls_arr-ls)).argmin() # finding index corresponding to wanted solar long
    assert np.min(np.abs(ls_arr-ls)) < 1, 'Solar Longitude {} not in the file'.format(ls)
    return idx

def temperature_level(filename, ls1):
    ### ls1, ls2 - solar longitude 1, 2 ###
    
    nc_file = filename
    data = Dataset(nc_file, mode='r')
    
    t = data.variables['T'][:] # perturbation potential temp
    p = data.variables['P'][:] # perturbation pressure
    pb = data.variables['PB'][:] # base state pressure
    ls = data.variables['L_S'][:] # solar longitude
    
    idx = find_ls_idx(ls, ls1)

    t0 = 300. #data.variables['T00'][:] # base temperature
    p0 = 610. #data.variables['P00'][:] # base pressure
    r_d = 191.8366
    cp = 767.3466
    g = 3.727
    gamma = r_d/cp
    
    pot_temp = t + t0 # potential temperature
    press = p + pb # pressure
    
    temp = pot_temp*(press/p0)**gamma # temperature
    
    temp_avg = temp[idx] # avg in long
    press_avg = press.mean(axis = 3).mean(axis = 2)[idx] # avg in lat, long
    
    return temp_avg, press_avg

def zonal_temperature(filename, ls1, ls2):
    ### ls1, ls2 - solar longitude 1, 2 ###
    
    nc_file = filename
    data = Dataset(nc_file, mode='r')
    
    t = data.variables['T'][:] # perturbation potential temp
    p = data.variables['P'][:] # perturbation pressure
    pb = data.variables['PB'][:] # base state pressure
    ls = data.variables['L_S'][:] # solar longitude
    
    idx1 = find_ls_idx(ls, ls1)
    idx2 = find_ls_idx(ls, ls2)

    t0 = 300. #data.variables['T00'][:] # base temperature
    p0 = 610. #data.variables['P00'][:] # base pressure
    r_d = 191.8366
    cp = 767.3466
    g = 3.727
    gamma = r_d/cp
    
    pot_temp = t + t0 # potential temperature
    press = p + pb # pressure
    
    temp = pot_temp*(press/p0)**gamma # temperature
    
    temp_avg = temp.mean(axis = 3)[idx1:idx2].mean(axis = 0) # avg in long and solar long
    press_avg = press.mean(axis = 3).mean(axis = 2)[idx1:idx2].mean(axis = 0) # avg in lat, long and solar long
    
    temp_total_avg = temp_avg.mean(axis=1).mean() # total avg temp
    h = (r_d*(temp_total_avg/g))*(np.log(p0/press_avg)/np.log(np.e)) # hypsometric equation?
    
    return temp_avg, press_avg, h/1000.

def zonal_wind(filename, ls1, ls2):
    nc_file = filename
    data = Dataset(nc_file, mode='r')
    
    u = data.variables['U'][:] # zonal wind speed in m/s
#    p = data.variables['P'][:] # perturbation pressure
#    pb = data.variables['PB'][:] # base state pressure
    ls = data.variables['L_S'][:] # solar longitude
    
    idx1 = find_ls_idx(ls, ls1)
    idx2 = find_ls_idx(ls, ls2)
    
    u_avg = u.mean(axis = 3)[idx1:idx2].mean(axis = 0) # avg in long and solar long
    
    return u_avg

def net_hr_aer(filename, ls1, ls2):
    nc_file = filename
    data = Dataset(nc_file, mode='r')
    
    hrvis = data.variables['HRAERVIS'][:] # heating rate in visible
    hrir = data.variables['HRAERIR'][:] # heating rate in infrared
    p = data.variables['P'][:] # perturbation pressure
    pb = data.variables['PB'][:] # base state pressure
    ls = data.variables['L_S'][:] # solar longitude
    
    idx1 = find_ls_idx(ls, ls1)
    idx2 = find_ls_idx(ls, ls2)
    
    press = p + pb # pressure

    hrvis_avg = hrvis.mean(axis = 3)[idx1:idx2].mean(axis = 0) # avg in long and solar long 
    hrir_avg = hrir.mean(axis = 3)[idx1:idx2].mean(axis = 0) # avg in long and solar long
    net_ir = hrvis_avg + hrir_avg # net heating rate due to aerosols
    
    press_avg = press.mean(axis = 3).mean(axis = 2)[idx1:idx2].mean(axis = 0) # avg in lat, long and solar long
    
    return net_ir, press_avg
    
    
#ls1, ls2 = 315,335 # range of the two solar longitudes we want to look at
#temp_avg, press_avg, h = zonal_temperature('./test_data/wrfout_d01_0002_00582', ls1, ls2)
#u_avg = zonal_wind('./test_data/wrfout_d01_0002_00582', ls1, ls2)
#
#lat = np.linspace(-90,90,36) # latitude
#xlong = np.linspace(0,360,72)
#lat, press = np.meshgrid(lat, press_avg)
#
#
#fig = plt.subplots(figsize=(14, 4))
#plt.subplot(1,2,1)
#ax1 = plt.contourf(lat, press, temp_avg, 11, cmap = 'viridis', origin='upper')
##plt.clabel(ax1, color='k')
#ax2 = plt.contour(lat, press, temp_avg, 11, colors='black')
#plt.clabel(ax2, color='k', inline = 1, fmt='%1i')
#plt.xlabel('Latitude (degree)')
#plt.ylabel('Pressure (Pa)')
#plt.title('Avg Zonal Temperature LS {}-{}'.format(ls1, ls2))
#plt.gca().invert_yaxis()
#plt.yscale('log')
##
##ax2 = plt.twinx()
##ax2.set_yticks(h)
##ax2.set_yscale('log')
##ax2.set_ylabel('Height (km)')
#plt.show()


#plt.subplot(1,2,2)
#ax1 = plt.contourf(lat, press, u_avg, 11, cmap = 'magma', origin='upper')
#ax2 = plt.contour(lat, press, u_avg, 11, colors='black')
#plt.clabel(ax2, color='k', inline = 1, fmt='%1i')
##plt.colorbar(ax1)
#plt.xlabel('Latitude (degree)')
#plt.ylabel('Pressure (Pa)')
#plt.title('Avg Zonal Wind LS {}-{}'.format(ls1, ls2))
#plt.gca().invert_yaxis()
#plt.yscale('log')
#plt.show()





