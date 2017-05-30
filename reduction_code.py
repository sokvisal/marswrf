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
            t = data.variables['T'][:] # perturbation potential temp
            
            p = data.variables['P'][:] # perturbation pressure
            pb = data.variables['PB'][:] # base state pressure
            
            pot_temp = t + t0 # potential temperature
            p = p + pb # pressure
            
            tmp.append((pot_temp*(p/p0)**gamma).mean(axis = 3)) # temperature
            tmp.append(p.mean(axis = 3))
        else: 
            tmp2 = data.variables[var][:]
            tmp.append(tmp2.mean(axis = 3))
    return ls, np.array(tmp)

def load_misc(filename, data, var_name):
    var = var_name + '_PHY_AM'
    t_am = data.variables[var][:] # solar long
    
    var2 = var_name + '_PHY_PM'
    t_pm = data.variables[var2][:] # solar long
    
    ls = data.variables['L_S'][:]

    t_d = t_pm.mean(axis=3) - t_am.mean(axis=3)
    t_a = t_pm.mean(axis=3) + t_am.mean(axis=3) 
    
    t_d_2pa = t_pm[:,42,:,:] - t_am[:,42,:,:]
    t_a_2pa = t_pm[:,42,:,:] + t_am[:,42,:,:]

    return t_a/2., t_d/2.,  t_a_2pa/2., t_d_2pa/2., ls

def load_misc_zm(filename, data, var_name):
    temp = data.variables[var_name][:].mean(axis=3)
    return temp

def load_misc3D(filename, data, var_name):
    temp = data.variables[var_name][:]
    ls = data.variables['L_S'][:]
    return temp, ls

def load_misc4D(filename, data, var_name):
    temp = data.variables[var_name][:]
    return temp
    
def init_reduction(filedir):
    arg = sys.argv[1]
    tar_cond = int(sys.argv[2])
    def wrfout():
        filepath = filedir + '/wrfout_d01*'
        print (filepath)
        
        if not os.path.exists(filedir+'/reduction'): 
            os.mkdir(filedir+'/reduction')
        	#    print filedir
        varlist = ['T', 'U']
        for num, i in enumerate(sorted(glob.glob(filepath))):
            print (i )
            if num == 0:
                nc_file = i
                data = Dataset(nc_file, mode='r')
                
                ls, tmp = load_zm(nc_file, data, varlist)
            else:
                nc_file = i
                data = Dataset(nc_file, mode='r')
        		    
                ls2, tmp2 = load_zm(nc_file, data, varlist)
                
                ls = np.concatenate((ls, ls2))
                tmp = np.concatenate((tmp, tmp2),axis=1)
        
        filedir2 = i.replace(i,'{}/reduction/wrfout_{}'.format(filedir, 'ls'))
        np.save(filedir2, ls)
        del ls
        
        varlist.insert(1, 'P')
        for j, var in enumerate(varlist):
            filedir2 = i.replace(i,'{}/reduction/wrfout_ {}'.format(filedir, var))
            print('Saving', var)
            np.save(filedir2, tmp[j])
            
    def wrfout_ext():
        filepath = filedir + '/wrfout_d01*'
        print (filepath)
        
        if not os.path.exists(filedir+'/reduction'): 
            os.mkdir(filedir+'/reduction')
        	    
        for num, i in enumerate(sorted(glob.glob(filepath))):
            print(i)
            if num == 0:
                nc_file = i
                data = Dataset(nc_file, mode='r')
                
                t = load_misc4D(nc_file, data, 'T' )[:,16] # t at 2 km
                psfc, ls_psfc = load_misc3D(nc_file, data, 'PSFC' )
                v = load_misc_zm(nc_file, data, 'V' )
            else:
                nc_file = i
                data = Dataset(nc_file, mode='r')
        		    
                t2 = load_misc4D(nc_file, data, 'T')[:,16]
                psfc2, ls_psfc2 = load_misc3D(nc_file, data, 'PSFC')
                v2 = load_misc_zm(nc_file, data, 'V' )
        		    
                t = np.concatenate((t, t2),axis=0)
                psfc = np.concatenate((psfc, psfc2),axis=0)
                v = np.concatenate((v, v2),axis=0)
                ls_psfc = np.concatenate((ls_psfc, ls_psfc2),axis=0)
        	    
        filedir2 = i.replace(i,'{}/reduction/wrfout'.format(filedir))
        varlist = ['_T2KM', '_PSFC', '_ls_PSFC', '_V']
        for num, i in enumerate([t, psfc, ls_psfc, v]):
            print('Saving', varlist[num])
            np.save(filedir2 + varlist[num], i)
    
    def auxhist9():
        print ('hi')  
        print ('cool')
        filepath = filedir + '/auxhist9*'
        print (filepath)
        
        for num, i in enumerate(sorted(glob.glob(filepath))):
            print (i)
            if num == 0:
                nc_file = i
                data = Dataset(nc_file, mode='r')
		
                t_a, t_d, t_a_2Pa, t_d_2Pa, ls = load_misc(nc_file, data, 'T' )
            else: 
                nc_file = i
                data = Dataset(nc_file, mode='r')
		
                t_a2, t_d2, t_a2_2Pa, t_d2_2Pa, ls2 = load_misc(nc_file, data, 'T')
                
                t_a = np.concatenate((t_a, t_a2),axis=0)
                t_d = np.concatenate((t_d, t_d2),axis=0)
                
                t_a_2Pa = np.concatenate((t_a_2Pa, t_a2_2Pa),axis=0)
                t_d_2Pa = np.concatenate((t_d_2Pa, t_d2_2Pa),axis=0)
                
                ls = np.concatenate((ls, ls2), axis=0)
	
        filedir2 = i.replace(i,'{}/reduction/wrfout'.format(filedir))        
        var_list = ['_TAVG', '_TDIFF', '_TAVG2PA', '_TDIFF2PA', '_ls_AUX9']
        for num, i in enumerate([t_a, t_d, t_a_2Pa, t_d_2Pa, ls]):
            np.save(filedir2 + var_list[num], i)
 
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
                print (psfc2.shape, ls_psfc2.shape)
                psfc = np.concatenate((psfc, psfc2),axis=0)
                ls_psfc = np.concatenate((ls_psfc, ls_psfc2),axis=0)

        filedir3 = i.replace('auxhist5','reduction/wrfout')
        var_list = ['_psfc','_ls_psfc']
        for num, i in enumerate([psfc,ls_psfc]):
            print(psfc.shape)
            np.save(filedir3 + var_list[num], i)

    if arg == 'wrfout': wrfout()
    if arg == 'wrfout_ext': wrfout_ext()
    if arg == 'auxhist9': auxhist9()
    if arg == 'auxhist5': auxhist5()
    
    print (tar_cond)
    if tar_cond == 1:
        tar = tarfile.open(filedir+'/reduction.tar.gz', 'w:gz')
        for i in glob.glob(filedir+'/reduction/wrfout*'):
            print ('Tarring file,', i)            
            tar.add(i, arcname = i.replace(filedir, ''))
            tar.close()
init_reduction('./../data_marswrf/diag.r14p1dustL40/data')
 
