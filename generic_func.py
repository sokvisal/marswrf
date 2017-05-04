import numpy as np
from netCDF4 import Dataset
import matplotlib
import matplotlib.pylab as plt
import os
import glob

matplotlib.rcParams.update({'font.size': 9})

def tp(filename, data): # temperature and pressure
    ### ls1, ls2 - solar longitude 1, 2 ###
    
    t = data.variables['T'][:] # perturbation potential temp
    
    p = data.variables['P'][:] # perturbation pressure
    pb = data.variables['PB'][:] # base state pressure

    t0 = 300. #data.variables['T00'][:] # base temperature
    p0 = 610. #data.variables['P00'][:] # base pressure
    r_d = 191.8366
    cp = 767.3466
    g = 3.727
    gamma = r_d/cp
    
    pot_temp = t + t0 # potential temperature
    press = p + pb # pressure
    
    temp = pot_temp*(press/p0)**gamma # temperature
    
    temp = temp.mean(axis = 3) # avg in long
    press = press.mean(axis = 3) # avg in long
    
    filename1 = filename + '_t'
    np.save(filename1, temp)
    
    filename2 = filename + '_p'
    np.save(filename2, press)
    
def misc(filename, data): # miscellaneous
    
    ls = data.variables['L_S'][:] # solar long
    u = data.variables['U'][:] # zonal wind
    ph = data.variables['PH'][:] # perturbation geopot height
    phb = data.variables['PHB'][:] # base state geopot height
    
    geo_height = ph + phb
    geo_height = geo_height.mean(axis = 3)

    u = u.mean(axis = 3)

    filename1 = filename + '_ls'
    np.save(filename1, ls)
    
    filename2 = filename +'_u'
    np.save(filename2, u) 

    filename3 = filename + '_p'
    np.save(filename3, geo_height)
    
def init_reduction(filedir):
    
    filepath = filedir + '/wrfout_d01_0001*'
    print filepath
    if not os.path.exists(filedir+'/reduction'): os.mkdir(filedir+'/reduction')
#    print filedir
    
    for i in sorted(glob.glob(filepath)):
        print i
        nc_file = i
        data = Dataset(nc_file, mode='r')
        
        filedir = i.replace('wrfout','reduction/wrfout')
        
        tp(filedir, data)
        misc(filedir, data)
        
        
init_reduction('./../mars/planetWRF-dev-new/WRFV3/test/em_global_mars/')
 
def find_ls_idx(ls_arr, ls):
    idx = (np.abs(ls_arr-ls)).argmin() # finding index corresponding to wanted solar long
    assert np.min(np.abs(ls_arr-ls)) < 1, 'Solar Longitude {} not in the file'.format(ls)
    return idx

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
    
