import numpy as np
from netCDF4 import Dataset
import dwell.fft as fft
import matplotlib
matplotlib.use('Agg')
#import matplotlib
import matplotlib.pylab as plt
import matplotlib.colors as colors
#import matplotlib.tri as tri
import scipy
import os
import glob
import sys

matplotlib.rcParams.update({'font.size': 9})

def find_ls_idx(ls_arr, ls):
    idx = (np.abs(ls_arr-ls)).argmin() # finding index corresponding to wanted solar long
    assert np.min(np.abs(ls_arr-ls)) < 1, 'Solar Longitude {} not in the file'.format(ls)
    return idx

def martians_month(ls, data):
    temp = []
    for i in np.arange(0, 12):
        idx = np.where((ls>i*30)&(ls<(i+1)*30))
        temp.append(data[idx].mean(axis=0))
    temp = np.array(temp)
    return temp

def redefine_latField(v):
    temp = np.zeros((52,36))
    for i in np.arange(v.shape[0]):
        for j in np.arange(v.shape[1]-1):
            temp[i,j] = np.mean(v[i,j]+v[i,j+1])
    return temp

def zonal_avg(filedir, month, ls2):
    print ('Looking at zonal avg')
    ### ls1, ls2 - solar longitude 1, 2 ###
    
    filepath = glob.glob(filedir + '*_temp.npy')[0]
    print (filepath)
    temp = np.load(filepath)
    
    filepath = glob.glob(filedir + '*_u.npy')[0]
    print (filepath)
    u = np.load(filepath)
    
    filepath = glob.glob(filedir + '*_press.npy')[0]
    print (filepath)
    p = np.load(filepath)
        
    filepath = glob.glob(filedir + '*_ls.npy')[0]
    ls = np.load(filepath)[::2]

    idx1 = np.where(ls == 360)[0][1] # only looking at the second year
    idx2 = np.where(ls == 360)[0][2]
    
    ls = ls[idx1:idx2]
    temp = temp[idx1:idx2]
    u = u[idx1:idx2]
    p = p[idx1:idx2]
    print (ls)
    
    zonal_t = []
    zonal_u = []
    zonal_p = []
    for i in np.arange(0, 12):
        print ('month', i+1)
        idx = np.where((ls>i*30)&(ls<(i+1)*30))
        zonal_t.append(temp[idx].mean(axis=0))
        zonal_u.append(u[idx].mean(axis=0))
        zonal_p.append(p[idx].mean(axis=0))
    del temp, p
    zonal_t = np.array(zonal_t)
    zonal_u = np.array(zonal_u)
    zonal_p = np.array(zonal_p)
    
    print ('Plotting some cool shit')
    fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(20,14))
    for i, ax in enumerate(axes.flat):
        press = zonal_p[i-1][4:]
        
        #press2 = zonal_p[i-1].mean(axis=1)
        lat = np.linspace(-90, 90, 36) 
        temp_press = np.linspace(1e-2, 900, 52)[4:]
        
        lat, temp_press = np.meshgrid(lat, temp_press)
        
        temp = zonal_t[i-1][4:]
        print (temp_press.shape, lat.shape, temp.shape)
#        temp = temp[4:]
        
        im = ax.contourf(lat, press, temp, 12, cmap='viridis')
        ax.contour(lat, press, temp, 12, colors='k')
        
        ax.set_title(r'Avg Zonal Temp LS {}-{}'.format((i)*30, (i+1)*30))
        ax.invert_yaxis()
        ax.set_yscale('log')
        ax.set_ylim([900, 1e-2])
    
    print ('Saving 1st cool shit')
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.3, orientation='horizontal', pad=0.03)
    plt.savefig('zonal_temp.png',bbox_inches='tight', dpi=600)
    
    print ('Saving 2nd cool shit')
    fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(20,14))
    for i, ax in enumerate(axes.flat):
        press = zonal_p[i-1][4:]
        
        #press2 = zonal_p[i-1].mean(axis=1)
        lat = np.linspace(-90, 90, 36) 
        temp_press = np.linspace(1e-2, 900, 52)[4:]
        
        lat, temp_press = np.meshgrid(lat, temp_press)
        
        u = zonal_u[i-1][4:]
        
        im = ax.contourf(lat, press, u, 12, cmap='inferno')
        ax.contour(lat, press, u, 12, colors='k')
        
        ax.set_title(r'Avg Zonal Wind LS {}-{}'.format((i)*30, (i+1)*30))
        ax.invert_yaxis()
        ax.set_yscale('log')
        ax.set_ylim([900, 1e-2])
    
    print ('Saving 2nd cool shit')
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.3, orientation='horizontal', pad=0.03)
    plt.savefig('zonal_wind.png',bbox_inches='tight', dpi=600)

def zonal_diff(filedir, month, ls2):
    print ('Looking at thermal tides')
    
    filepath = glob.glob(filedir + '*_press.npy')[0]
    print (filepath)
    p = np.load(filepath)
        
    filepath = glob.glob(filedir + '*_t_d.npy')[0]
    t_d = np.load(filepath)
    t_d = t_d
    
    filepath = glob.glob(filedir + '*_t_d_2Pa.npy')[0]
    t_d_2Pa = np.load(filepath)
    t_d_2Pa = t_d_2Pa 
    
    filepath = glob.glob(filedir + '*_ls.npy')[0]
    ls = np.load(filepath)[::2]

    idx1 = np.where(ls == 360)[0][1] # only looking at the second year
    idx2 = np.where(ls == 360)[0][2]
    
    ls = ls[idx1:idx2]
    t_d = t_d[idx1:idx2]
    t_d_2Pa = t_d_2Pa[idx1:idx2]
    p = p[idx1:idx2]
    
    print (t_d_2Pa.shape)
    
#    test_diff = t_d_2Pa.reshape((223,3,36,72)).mean(axis=1)
    ampl, phase, axis = fft.spec1d(t_d_2Pa, 1/72., use_axes = 2)
    ampl = ampl
    phase = phase
    print (np.min(ampl))
#    ampl.hstack(axis=0)
#    phase.hstack(axis=0)
    
    test = np.zeros((4,669,36)) # amplitude
    test2 = np.zeros((4,669,36)) # phase
    for j in np.arange(0,4):
        for i in np.arange(0,669):
            test[j,i] = ampl[i,:,j]
            test2[j,i] = phase[i,:,j]
    
    print ('Plotting some cool shit')
    print ('Saving 1st cool shit')
    lat = np.linspace(-90, 90, 36) 
    ls = np.linspace(0,360, 669)
    t_lat, t_ls = np.meshgrid(lat, ls)
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(8,8))
    for i, ax in enumerate(axes.flat):
        amplitude = test[i]
        phase = test2[i]
        im = ax.contourf(t_ls, t_lat, amplitude, cmap='viridis')
#        ax.contour(t_ls, t_lat, phase, 20)
        ax.set_title(r'Wavenumber {} Tide'.format(i))
#    plt.contour(t_ls, t_lat, np.tan(test2), colors='k')
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.3, orientation='horizontal', pad=0.03)
    
    zonal_t_d = martians_month(ls, t_d)
    zonal_t_2P = martians_month(ls, t_d_2Pa)
    zonal_p = martians_month(ls, p)
    
    fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(20,14))
    for i, ax in enumerate(axes.flat):
        press = zonal_p[i-1][4:]
        
        #press2 = zonal_p[i-1].mean(axis=1)
        lat = np.linspace(-90, 90, 36) 
        temp_press = np.linspace(1e-2, 900, 52)[4:]
        
        #lat, press = np.meshgrid(lat, press2[5:])
        lat, temp_press = np.meshgrid(lat, temp_press)
        
        temp = zonal_t_d[i-1][4:]
        
        im = ax.contourf(lat, press, temp, 12, cmap='viridis')
        ax.contour(lat, press, temp, 12, linewidths=0.5, colors='black')
        
        ax.set_title(r'Avg Zonal T$_d$ LS {0}-{1}'.format((i)*30, (i+1)*30))
        ax.invert_yaxis()
        ax.set_yscale('log')
        ax.set_ylim([900, 1e-2])
        
    print ('Saving 2nd cool shit')
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.3, orientation='horizontal', pad=0.03)
    plt.savefig('zonal_tdiff.png',bbox_inches='tight', dpi=600)
    
    fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(20,14))
    for i, ax in enumerate(axes.flat):
        t_2P = zonal_t_2P[i-1]
        
        #press2 = zonal_p[i-1].mean(axis=1)
        lat = np.linspace(-90, 90, 36) 
        lon = np.linspace(0, 360, 72)
        
        lon, lat = np.meshgrid(lon, lat)

        im = ax.contourf(lon, lat, t_2P, cmap='viridis')
        ax.contour(lon, lat, t_2P, linewidths=0.5, colors='black')
        ax.set_title(r'Avg T$_d$ LS {0}-{1} at 2 Pa'.format((i)*30, (i+1)*30))
        
    print ('Saving 3rd cool shit')
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.3, orientation='horizontal', pad=0.03)
    plt.savefig('zonal_tdiff_2Pa.png',bbox_inches='tight', dpi=600)

#zonal_avg('./test_data/reduction_mcd_wbm/',12,2)
#zonal_diff('./test_data/reduction_mcd_kdm/',12,2)

def fft_hovmoller(filedir):
    print ('Looking at spectral of surface presssure')
    
    filepath = glob.glob(filedir + '*0_psfc.npy')[0]
    print (filepath)
    psfc = np.load(filepath)
    
    filepath = glob.glob(filedir + '*_ls_psfc.npy')[0]
    print (filepath)
    ls = np.load(filepath)
    
    idx1 = np.where(ls == 360)[0][1] # only looking at the second year
    idx2 = np.where(ls == 360)[0][2]
    ls = ls[idx1:idx2]
    psfc = psfc[idx1:idx2]
    
    #avging for specific latitude
    lat = np.linspace(-90,90,36)
    idx = np.where((lat>-10)&(lat<10))[0]
    
    psfc = psfc[:,18,:]#.mean(axis=1)
    
    idx = np.where((ls>270)&(ls<301))[0]
    ls = ls[idx]
    psfc = psfc[idx] 
    psfc = psfc - psfc.mean(axis=0)
    print (psfc.shape)
#    
    lon = np.linspace(0,360,72)
    lon, ls = np.meshgrid(lon, ls)
    
    plt.figure(figsize=(10,4))
    plt.subplot(1,2,1)
    plt.contourf(lon, ls, psfc)
    plt.ylabel('Solar Longitude')
    plt.xlabel('Longitude')
    plt.title(r'Hovm$\mathrm{\"{o}}$ller Diagram of Surface Pressure')

    lat = np.linspace(-90,90,36)
    
    plt.subplot(1,2,2)
    ampl, phase, cycle, waven = fft.spec(psfc, 1./8, 1./72, axes=[0,1])
    
    padd = np.zeros((psfc.shape[0]*2, psfc.shape[1]*2))
    padd[:psfc.shape[0], :psfc.shape[1]] = psfc
    amplt = np.fft.fftshift(np.fft.fft2(padd))
    amplt = np.abs(amplt)**2
    c = np.fft.fftshift(np.fft.fftfreq(amplt.shape[0], .125))
    wavenu = np.fft.fftshift(np.fft.fftfreq(amplt.shape[1],2.51748252))[72:83]*360
#    cycle, wavenu = np.meshgrid(cycle, wavenu)
    amplt = np.log10(amplt).T[72:83]
    
#    plt.contourf(cycle, waven[:5], np.log10(ampl).T[:5])
    plt.contourf(c, wavenu, amplt)
    plt.colorbar()
    plt.ylabel('Wavenumber')
    plt.xlabel('Cycles/sol')
    plt.title(r'Amplitude of FFT of Hovm$\mathrm{\"{o}}$ller Diagram')
    
#fft_hovmoller('./test_data/reduction_mcd_wbm/')

def mer_stream_func(filedir):
    print ('Looking at zonal mean meridional function')
    
    filepath = glob.glob(filedir + '*_v.npy')[0]
    print (filepath)
    v = np.load(filepath)
    
    filepath = glob.glob(filedir + '*_press.npy')[0]
    print (filepath)
    p = np.load(filepath)
    
    filepath = glob.glob(filedir + '*_ls.npy')[0]
    print (filepath)
    ls = np.load(filepath)
    
    idx1 = np.where(ls == 360)[0][1] # only looking at the second year
    idx2 = np.where(ls == 360)[0][2]
    
    p = p[idx1:idx2]
    v = v[idx1:idx2]
    ls = ls[idx1:idx2]
    
    msf = np.zeros((3,52,36))
    for k in range(0,3):
    
        zonal_v = martians_month(ls, v)[k]
        zonal_p = martians_month(ls, p)[k]
        zonal_v = redefine_latField(zonal_v)
        
        lat = np.linspace(-90,90,36)
        
        a = 3389920.
        g = 3.727
        temp = np.zeros((52,36))
        for i in np.arange(1,temp.shape[0]):
            for j in np.arange(temp.shape[1]):
                temp[i,j] = 2*np.pi*(a/g)*np.cos(np.deg2rad(lat[j]))*np.trapz(zonal_v[:i,j][::-1],zonal_p[:i,j][::-1])
        msf[k] = temp
    
    msf = msf.mean(axis=0)
    press = np.arange(0,52)
    lat, press = np.meshgrid(lat, press)
    print (zonal_p.shape, lat.shape, msf.shape)
            
    plt.contourf(lat, zonal_p, msf)
    plt.gca().invert_yaxis()
    plt.yscale('log')
    plt.colorbar()
    
mer_stream_func('./test_data/reduction_mcd_kdm/')

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