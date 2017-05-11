import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use('Agg')
import matplotlib
import matplotlib.pylab as plt
import os
import glob
import sys

matplotlib.rcParams.update({'font.size': 9})

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
    
    pot_temp = t + t0 # potential temperature
    press = p + pb # pressure
    geo_height = ph + phb
    geo_height = geo_height.mean(axis = 3)
    
    temp = pot_temp*(press/p0)**gamma # temperature
    
    temp = temp.mean(axis = 3) # avg in long
    press = press.mean(axis = 3) # avg in long
    
    return ls, temp, press, geo_height

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
    arg = sys.argv[1]
    def wrfout():
	    filepath = filedir + '/wrfout_d01*'
	    print filepath
	    if not os.path.exists(filedir+'/reduction'): os.mkdir(filedir+'/reduction')
	#    print filedir
	    
	    for num, i in enumerate(sorted(glob.glob(filepath))):
		print i
		if num == 0:
		    nc_file = i
		    data = Dataset(nc_file, mode='r')
		    
		    ls, temp, press, geoH = load_temp(nc_file, data)
		else:
		    nc_file = i
		    data = Dataset(nc_file, mode='r')
		    
		    ls2, temp2, press2, geoH2 = load_temp(nc_file, data)
		    
		    ls = np.concatenate((ls, ls2),axis=0)
		    temp = np.concatenate((temp, temp2),axis=0)
		    press = np.concatenate((press, press2),axis=0)
		    geoH = np.concatenate((geoH, geoH2),axis=0)
	    
	    filedir2 = i.replace('wrfout','reduction/wrfout')
	    var_list = ['_ls','_temp','_press','_geoH']
	    for num, i in enumerate([ls,temp,press,geoH]):
        	np.save(filedir2 + var_list[num], i)
    
    def auxhist9():
	print 'hi'  
	print 'cool'
	filepath2 = filedir + '/auxhist9*'
	print filepath2
	for num, i in enumerate(sorted(glob.glob(filepath2))):
	    print i
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
	    print i
	    if num == 0:
		nc_file = i
		data = Dataset(nc_file, mode='r')

		psfc, ls_psfc = load_misc3D(nc_file, data, 'PSFC' )
	    else:
		nc_file = i
		data = Dataset(nc_file, mode='r')

		psfc2, ls_psfc2 = load_misc3D(nc_file, data, 'PSFC')
		psfc = np.concatenate((psfc, psfc2),axis=0)
		ls_psfc2 = np.concatenate((ls_psfc, ls_psfc2),axis=0)

	filedir3 = i.replace('auxhist5','reduction/wrfout')
	var_list = ['_psfc','_ls_psfc']
	for num, i in enumerate([psfc,ls_psfc]):
	    np.save(filedir3 + var_list[num], i) 

    if arg == 'wrfout': return wrfout()
    if arg == 'auxhist9': return auxhist9()
    if arg == 'auxhist5': return auxhist5()
init_reduction('./../planetWRF/WRFV3/run/new_wbm')
	 
def find_ls_idx(ls_arr, ls):
    idx = (np.abs(ls_arr-ls)).argmin() # finding index corresponding to wanted solar long
    assert np.min(np.abs(ls_arr-ls)) < 1, 'Solar Longitude {} not in the file'.format(ls)
    return idx


def zonal_temperature(filedir, month, ls2):
    ### ls1, ls2 - solar longitude 1, 2 ###
    
    filepath = filedir + '*_t.npy'
    print filepath
    temp = np.empty([1,52,36])
    for i in glob.glob(filepath):
        print temp.shape
        temp = np.concatenate((temp, np.load(i)), axis=0)
    temp = temp[1:]
    
    filepath = filedir + '*_u.npy'
    print filepath
    u = np.empty([1,52,36])
    for i in glob.glob(filepath):
        u = np.concatenate((u, np.load(i)), axis=0)
    u = u[1:]
    
    filepath = filedir + '*_p.npy'
    print filepath
    p = np.empty([1,52,36])
    for i in glob.glob(filepath):
        print i
        p = np.concatenate((p, np.load(i)), axis=0)
    p = p[1:]
    print filedir+'*_am'
    
    if sorted(glob.glob(filedir+'*_am.npy')) or sorted(glob.glob(filedir+'*_am.npy')):
        print 'Looking at thermal tides'
        
        filepath = filedir + '*_t_am.npy'
        t_am = np.empty([1,52,36])
        for i in sorted(glob.glob(filepath)):
            t_am = np.concatenate((t_am, np.load(i)), axis=0)
        t_am = t_am[1:][::2]
        
        filepath = filedir + '*_t_pm.npy'
        t_pm = np.empty([1,52,36])
        for i in sorted(glob.glob(filepath)):
            t_pm = np.concatenate((t_pm, np.load(i)), axis=0)
        t_pm = t_pm[1:][2::2]
        
    filepath = filedir + '*_ls.npy'
    ls = []
    for i in glob.glob(filepath):
        ls = np.concatenate((ls, np.load(i)), axis=0)
    
    print np.where(ls ==360)
    idx = np.where(ls ==360)[0][1] # only looking at the second year
    ls = ls[idx:]
    temp = temp[idx:]
    u = u[idx:]
    p = p[idx:]
    
    zonal_t = []
    zonal_u = []
    zonal_p = []
    zonal_t_am, zonal_t_pm = [], []
    for i in np.arange(0, 12):
        idx = np.where((ls>i*30)&(ls<(i+1)*30))
        zonal_t.append(temp[idx].mean(axis=0))
        zonal_u.append(u[idx].mean(axis=0))
        zonal_p.append(p[idx].mean(axis=0))
        zonal_t_am.append(t_am[idx].mean(axis=0))
        zonal_t_pm.append(t_pm[idx].mean(axis=0))
    del temp,u,p
    zonal_t_am = np.array(zonal_t_am)
    zonal_t_pm = np.array(zonal_t_pm)
    zonal_t = np.array(zonal_t)
    zonal_u = np.array(zonal_u)
    zonal_p = np.array(zonal_p)
    print zonal_t.shape, zonal_p.shape
#    
    lat = np.linspace(-90,90,36) # latitude
    press = zonal_p.mean(axis=2).mean(axis=0)
    lat, press = np.meshgrid(lat, press)
    print lat.shape, press.shape
    
    print 'Plotting zonal temp' 
    fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(20,14))
    for i, ax in enumerate(axes.flat):
#        plt.subplot(4,3,i)
        temp = (zonal_t_pm[i-1] - zonal_t_am[i-1])/2.
        #temp = zonal_t[i-1]
        im = ax.contourf(lat, press, temp, cmap='viridis')
        ax.contour(lat, press, temp, linewidths=0.5, colors='black')
        ax.set_title(r'Avg Zonal T$_d$ LS {}-{}'.format((i)*30, (i+1)*30))
        ax.invert_yaxis()
        ax.set_yscale('log')
        
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.3, orientation='horizontal', pad=0.03)
    plt.savefig('test.png',bbox_inches='tight', dpi=600)
    
#    print 'Plotting zonal wind' 
#    fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(20,14))
#    for i, ax in enumerate(axes.flat):
#        u = zonal_u[i-1]
#        im = ax.contourf(lat, press, u, cmap='inferno')
#        ax.contour(lat, press, u, linewidths=0.5, colors='black')
#        ax.set_title('Avg Zonal Wind LS {}-{}'.format((i)*30, (i+1)*30))
#        ax.invert_yaxis()
#        ax.set_yscale('log')
        
#    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.3, orientation='horizontal', pad=0.03)
#    plt.savefig('test2.png',bbox_inches='tight', dpi=600)
    
    
#zonal_temperature('./test_data/reduction_mcd_kdm/',12,2)

def zonal_temperature2(filename, ls1, ls2):
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
    
