import numpy as np
from netCDF4 import Dataset
import scipy.ndimage 
import os
import glob
import sys

def load_zm(filename, data, varlist):
    t0 = 300. #data.variables['T00'][:] # base temperature
    p0 = 610. #data.variables['P00'][:] # base pressure
    r_d = 191.8366
    cp = 767.3466
    g = 3.727
    gamma = r_d/cp
    
    ls = data.variables['L_S'][:] # solar long
    
    tmp = []
    for var in varlist:
        if var == 'T':
            t = data.variables['T'][:][:,25:30] # perturbation potential temp
            
            p = data.variables['P'][:][:,25:30] # perturbation pressure
            pb = data.variables['PB'][:][:,25:30] # base state pressure
            
            pot_temp = t + t0 # potential temperature
            p = p + pb # pressure
            del t, pb
            
            tmp.append((pot_temp*(p/p0)**gamma).mean(axis=1)) # temperature
            tmp.append(p.mean(axis=1))
        elif var == 'PH':
            ph = data.variables['PH'][:]
            phb = data.variables['PHB'][:]
            
            tmp.append((ph+phb).mean(axis = 3))
        elif var == 'U':
            u = data.variables['U'][:]
            tmp.append(u[:,:,16:20,61:65])
        elif var == 'PSFC' or var == 'TSK':
            u = data.variables[var][:]
            tmp.append(u)
        else: 
            tmp2 = data.variables[var][:]
            tmp.append(tmp2[:,:,16:20,60:64])
    return ls, tmp


filedir = './../model_run/dustL45/'
filepath = filedir + '/wrfout_d01*'
print (filepath)

if not os.path.exists(filedir+'/reduction'): 
    os.mkdir(filedir+'/reduction')

varlist = np.array(['T', 'PSFC', 'TSK'])
for num, i in enumerate(sorted(glob.glob(filepath)[:])):
    print (i)
    if num == 0:
        nc_file = i
        data = Dataset(nc_file, mode='r')
        
        lsd, tmp = load_zm(nc_file, data, varlist)
        print (np.array(tmp[0]).shape)
    else:
        nc_file = i
        data = Dataset(nc_file, mode='r')
		    
        lsd2, tmp2 = load_zm(nc_file, data, varlist)
        
        lsd = np.concatenate((lsd, lsd2))
        tmp = tmp + tmp2

def create_var(varnameList, units, data):
    tmp2 = dataset.createVariable(varnameList[0], np.float32, (varnameList[1],), zlib=True)
    tmp2.units = (units)
    tmp2[:] = data
def create2D_var(varnameList, units, data):
    tmp2 = dataset.createVariable(varnameList[0], np.float32, (varnameList[1], varnameList[2],), zlib=True)
    tmp2.units = (units)
    tmp2[:] = data
def create3D_var(varnameList, units, data):   
    tmp2 = dataset.createVariable(varnameList[0], np.float32, (varnameList[1], varnameList[2], varnameList[3],), zlib=True)
    tmp2.units = (units)
    tmp2[:] = data
def create4D_var(varnameList, units, data):   
    tmp2 = dataset.createVariable(varnameList[0], np.float32, (varnameList[1], varnameList[2], varnameList[3], varnameList[4],), zlib=True)
    tmp2.units = (units)
    tmp2[:] = data

dataset = Dataset('./dustL45_test.nc', 'w')

varlen = varlist.size + 1
time_dim = np.vstack(tmp[0::varlen]).shape[0]

solar_long = dataset.createDimension('time', time_dim)
pressure = dataset.createDimension('bottom_top', None)
latitude = dataset.createDimension('south_north', None)
longitude = dataset.createDimension('west_east', None)

ls = dataset.createVariable('LS', np.float32, ('time',))
lat = dataset.createVariable('LAT', np.float32, ('south_north',))
long = dataset.createVariable('LONG', np.float32, ('west_east',))

ls.units = ('Solar Longtitude')
ls[:] = lsd

create3D_var(['T50Pa', 'time', 'south_north', 'west_east'], 'K (measured @equator @50 Pa)', np.vstack(tmp[0::varlen]))
create3D_var(['P', 'time', 'south_north', 'west_east'], 'K (measured @equator @50 Pa)', np.vstack(tmp[1::varlen]))
create3D_var(['PSFC', 'time', 'south_north', 'west_east'], 'Pa (Surface Pressure)', np.vstack(tmp[2::varlen]))
create3D_var(['TSK', 'time', 'south_north', 'west_east'], 'K (Surface Temperature)', np.vstack(tmp[3::varlen]))

#lat.units = ('Degree')
#b = np.linspace(-180,180,72)
#a = np.linspace(-90,90,36)
#lat[:] = a[16:]
#long[:] = b

dataset.close()
