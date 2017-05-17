import numpy as np
from netCDF4 import Dataset
import dwell.fft as fft
import matplotlib
matplotlib.use('Agg')
#import matplotlib
import matplotlib.pylab as plt
import matplotlib.colors as colors
from matplotlib.backends.backend_pdf import PdfPages
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

def martians_month(ls, data, avg):
    temp = []
    for i in np.arange(0, 12):
        idx = np.where((ls>i*30)&(ls<(i+1)*30))
        if avg:
#            temp = np.dstack((temp, data[idx].mean(axis=0)))
            temp.append(data[idx].mean(axis=0))
        else:
#            temp = np.dstack((temp, data[idx]))
            temp.append(data[idx])
            #temp.append(data[idx])
    temp = np.array(temp)
    return temp

def redefine_latField(v):
    temp = np.zeros((52,36))
    for i in np.arange(v.shape[0]):
        for j in np.arange(v.shape[1]-1):
            temp[i,j] = np.mean(v[i,j]+v[i,j+1])
    return temp

def zonal_plt_monthly(ydata, ls, data, title,  cmap):
    fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(20,14))
    for i, ax in enumerate(axes.flat):
        y = ydata[i-1][4:]
        
        #press2 = zonal_p[i-1].mean(axis=1)
        lat = np.linspace(-90, 90, 36) 
        temp_press = np.linspace(1e-2, 900, 52)[4:]
        
        lat, temp_press = np.meshgrid(lat, temp_press)
        
        d = data[i-1][4:]
        
        im = ax.contourf(lat, y, d, 12, cmap=cmap)
        ax.contour(lat, y, d, 12, linewidths=0.5, colors='k')
        
        ax.set_title(r'Avg Zonal {} LS {}-{} at 2 Pa'.format((title), (i)*30, (i+1)*30))
        ax.invert_yaxis()
        ax.set_yscale('log')
        ax.set_ylim([900, 1e-2])
        
#        print ('Saving 1st cool shit')
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.3, orientation='horizontal', pad=0.03)
    plt.savefig(pdFfigures, format='pdf', bbox_inches='tight', dpi=600)

def basemap_plt_monthly(data, ls, title, cmap):
    fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(20,14))
    for i, ax in enumerate(axes.flat):
        d = data[i-1]
        
        #press2 = zonal_p[i-1].mean(axis=1)
        lat = np.linspace(-90, 90, 36) 
        lon = np.linspace(0, 360, 72)
        
        lon, lat = np.meshgrid(lon, lat)

        im = ax.contourf(lon, lat, d, cmap='viridis')
        ax.contour(lon, lat, d, linewidths=0.5, colors='black')
        ax.set_title(r'Avg {} LS {}-{} at 2 Pa'.format((title), (i)*30, (i+1)*30))
        
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.3, orientation='horizontal', pad=0.03)
    plt.savefig(pdFfigures, format='pdf', bbox_inches='tight', dpi=600)

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

    zonal_t = martians_month(ls, temp, True)
    zonal_u = martians_month(ls, u, True)
    zonal_p = martians_month(ls, p, True)
    
    print ('Plotting some cool shit')
    print ('Saving 1st shit')
    zonal_plt_monthly(zonal_p, ls, zonal_t, 'Temp', 'viridis')
    print ('Saving 2nd shit')
    zonal_plt_monthly(zonal_p, ls, zonal_u, 'Wind', 'inferno')

def zonal_diff(filedir, var1, var2):
    if var1 == '*_t_d.npy': name = 'diff'
    else: name = 'avg'
    print ('Looking at thermal tides')
    
    filepath = glob.glob(filedir + '*_press.npy')[0]
    print (filepath)
    p = np.load(filepath)

    filepath = glob.glob(filedir + var1)[0]
    t_d = np.load(filepath)
    t_d = t_d
    
    filepath = glob.glob(filedir + var2)[0]
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
    
#    test_diff = t_d_2Pa.reshape((223,3,36,72)).mean(axis=1)
    ampl, phase, axis = fft.spec1d(t_d_2Pa, 1/72., use_axes = 2)
    ampl = ampl
    phase = phase
    
    # stacking the array with the respectful wavenumber
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
        ax.set_title(r'm = {}'.format(i))
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.3, orientation='horizontal', pad=0.06)
    plt.savefig('wavenumber_timeseries_tdiff.png')
    
    zonal_t_d = martians_month(ls, t_d, True)
    zonal_t_2P = martians_month(ls, t_d_2Pa, True)
    zonal_p = martians_month(ls, p, True)
    
    print ('Saving 2nd shit')
    zonal_plt_monthly(zonal_p, ls, zonal_t_d, 'T$_{\mathrm{'+name+'}}$', 'viridis')
        
    print ('Saving 3rd cool shit')
    basemap_plt_monthly(zonal_t_2P, ls, 'T$_{\mathrm{'+name+'}}$', 'viridis')

pdFfigures = PdfPages('test.pdf')
directory = './test_data/reduction_no_dust_wbm/'
with PdfPages(directory+'figures.pdf') as pdFfigures:
    zonal_avg(directory,12,2)
    zonal_diff(directory, '*_t_d.npy', '*_t_d_2Pa.npy')
    zonal_diff(directory, '*_t_a.npy', '*_t_a_2Pa.npy')

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
    
#    idx = np.where((ls>300)&(ls<331))[0]
#    ls = ls[idx]
#    psfc = psfc[idx]
    
    psfc = martians_month(ls, psfc, False)
    print (psfc.shape)
    for i in range(0,12):
        psfc[i] = psfc[i] - psfc[i].mean(axis=2).mean(axis=1)
    ampl, phase, cycle, waven = fft.spec(psfc, 1./8, 1./72, axes=[1,2])
    
#    lon = np.linspace(0,360,72)
#    lon, ls = np.meshgrid(lon, ls)
#    plt.figure(figsize=(10,4))
#    plt.subplot(1,2,1)
#    plt.contourf(lon, ls, psfc)
#    plt.ylabel('Solar Longitude')
#    plt.xlabel('Longitude')
#    plt.title(r'Hovm$\mathrm{\"{o}}$ller Diagram PSFC')
#
#    lat = np.linspace(-90,90,36)
    
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
    plt.contourf(cycle[9], waven[9], ampl[9], cmap='inferno')
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
    
    start_month = 3
    end_month = 4
    msf = np.zeros((3,52,36))
    for k in range(start_month,end_month):
    
        zonal_v = martians_month(ls, v)[k]
        zonal_p = martians_month(ls, p)[k]
        zonal_v = redefine_latField(zonal_v)
        
        lat = np.linspace(-90,90,36)
        print (zonal_p[:,0])
        a = 3389920.
        g = 3.727
        temp = np.zeros((52,36))
        for i in np.arange(1,temp.shape[0]):
            for j in np.arange(temp.shape[1]):
                temp[i,j] = 2*np.pi*(a/g)*np.cos(np.deg2rad(lat[j]))*np.trapz(zonal_v[:i,j][::-1],zonal_p[:i,j][::-1])
        msf[k-start_month] = temp
    
    msf = msf.mean(axis=0)
    press = np.arange(0,52)
    lat, press = np.meshgrid(lat, press)
    print (zonal_p.shape, lat.shape, msf.shape)
            
    plt.contourf(lat, zonal_p, msf)
    plt.gca().invert_yaxis()
    plt.yscale('log')
    plt.colorbar()
    
#mer_stream_func('./test_data/reduction_no_dust_wbm/')

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