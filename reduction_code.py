import numpy as np
from netCDF4 import Dataset
import scipy.ndimage 
import os
import glob
import sys
import tarfile
import gzip

def name(**variable): 
    return [x for x in variable][0]

def load_temp(filename, data):
    t0 = 300. #data.variables['T00'][:] # base temperature
    p0 = 610. #data.variables['P00'][:] # base pressure
    r_d = 191.8366
    cp = 767.3466
    g = 3.727
    gamma = r_d/cp
    
    t = data.variables['T'][:] # perturbation potential temp
    
    p = data.variables['P'][:] # perturbation pressure
    pb = data.variables['PB'][:] # base state pressure
    
    ls = data.variables['L_S'][:] # solar long
    ph = data.variables['PH'][:] # perturbation geopot height
    phb = data.variables['PHB'][:] # base state geopot height
    
    u = data.variables['U'][:] # zonal wind
    
    pot_temp = t + t0 # potential temperature
    press = p + pb # pressure
    geo_height = ph + phb
    geo_height = geo_height.mean(axis = 3)
    
    temp = pot_temp*(press/p0)**gamma # temperature
    
    temp = temp.mean(axis = 3) # avg in long
    press = press.mean(axis = 3) # avg in long
    u = u.mean(axis = 3) # avg in long
    
    return ls, temp, press, geo_height, u

def load_misc(filename, data, var_name):
    var = var_name + '_PHY_AM'
    t_am = data.variables[var][:] # solar long
    
    var2 = var_name + '_PHY_PM'
    t_pm = data.variables[var2][:] # solar long

    t_d = t_pm.mean(axis=3) - t_am.mean(axis=3) 
    t_d_2pa = t_pm[:,42,:,:] - t_am[:,42,:,:]

    return (t_d)/2., t_d_2pa/2.

def load_misc3D(filename, data, var_name):
    temp = data.variables[var_name][:]
    ls = data.variables['L_S'][:]
    return temp, ls
    
def init_reduction(filedir):
    arg = sys.argv[1]
    tar_cond = int(sys.argv[2])
    def wrfout():
        filepath = filedir + '/wrfout_d01*'
        print (filepath)
        
        if not os.path.exists(filedir+'/reduction'): 
            os.mkdir(filedir+'/reduction')
        	#    print filedir
        	    
        for num, i in enumerate(sorted(glob.glob(filepath))):
            print(i)
            if num == 0:
                nc_file = i
                data = Dataset(nc_file, mode='r')
                
                ls, temp, press, geoH, u = load_temp(nc_file, data)
            else:
                nc_file = i
                data = Dataset(nc_file, mode='r')
        		    
                ls2, temp2, press2, geoH2, u2 = load_temp(nc_file, data)
        		    
                ls = np.concatenate((ls, ls2),axis=0)
                temp = np.concatenate((temp, temp2),axis=0)
                press = np.concatenate((press, press2),axis=0)
                geoH = np.concatenate((geoH, geoH2),axis=0)
                u = np.concatenate((u, u2),axis=0)
        	    
        filedir2 = i.replace('wrfout','reduction/wrfout')
        var_list = ['_ls','_temp','_press','_geoH', '_u']
        for num, i in enumerate([ls,temp,press,geoH,u]):
            print('Saving', var_list[num])
            np.save(filedir2 + var_list[num], i)
    
    def auxhist9():
        print ('hi')  
        print ('cool')
        filepath2 = filedir + '/auxhist9*'
        print (filepath2)
        
        for num, i in enumerate(sorted(glob.glob(filepath2))):
            print (i)
            if num == 0:
                nc_file = i
                data = Dataset(nc_file, mode='r')
		
                t_d, t_d_2Pa = load_misc(nc_file, data, 'T' )
            else: 
                nc_file = i
                data = Dataset(nc_file, mode='r')
		
                t_d2, t_d2_2Pa = load_misc(nc_file, data, 'T')
                t_d = np.concatenate((t_d, t_d2),axis=0)
                t_d_2Pa = np.concatenate((t_d_2Pa, t_d2_2Pa),axis=0)
	
        filedir3 = i.replace('auxhist9','reduction/wrfout')        
        var_list = ['_t_d','_t_d_2Pa']
        for num, i in enumerate([t_d,t_d_2Pa]):
            np.save(filedir3 + var_list[num], i)
 
    def auxhist5():    
        filepath2 = filedir+'/auxhist5*'
        
        for num, i in enumerate(sorted(glob.glob(filepath2))):
            print (i)
            if num == 0:
                nc_file = i
                data = Dataset(nc_file, mode='r')

                psfc, ls_psfc = load_misc3D(nc_file, data, 'PSFC' )
            else:
                nc_file = i
                data = Dataset(nc_file, mode='r')

                psfc2, ls_psfc2 = load_misc3D(nc_file, data, 'PSFC')
                psfc = np.concatenate((psfc, psfc2),axis=0)
                ls_psfc = np.concatenate((ls_psfc, ls_psfc2),axis=0)

        filedir3 = i.replace('auxhist5','reduction/wrfout')
        var_list = ['_psfc','_ls_psfc']
        for num, i in enumerate([psfc,ls_psfc]):
            print(psfc.shape)
            np.save(filedir3 + var_list[num], i)

    if arg == 'wrfout': wrfout()
    if arg == 'auxhist9': auxhist9()
    if arg == 'auxhist5': auxhist5()
    
    print (tar_cond)
    if tar_cond == 1:
	tar = tarfile.open(filedir+'/reduction.tar.gz', 'w:gz')
        for i in glob.glob(filedir+'/reduction/wrfout*'):
	    print ('Tarring file,', i)            
	    tar.add(i)
        tar.close()
init_reduction('./../planetWRF/WRFV3/run/new_wbm')
 
