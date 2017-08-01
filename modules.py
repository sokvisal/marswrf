from netCDF4 import Dataset
from pylab import *
import numpy as np
import glob
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from spharm import Spharmt
from tqdm import tqdm
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from scipy.ndimage.filters import uniform_filter1d
from scipy.signal import butter, lfilter

from matplotlib.colors import SymLogNorm

from matplotlib.ticker import AutoMinorLocator

# cmap = sns.cubehelix_palette(light=1, as_cmap=True, reverse=True)
cmap = ListedColormap(sns.color_palette("coolwarm", 15).as_hex())
cubehelix = sns.cubehelix_palette(dark=0, light=1, as_cmap=True)
sns.reset_orig()

matplotlib.rcParams['lines.linewidth'] = 0.6

matplotlib.rcParams['axes.linewidth'] = 1
matplotlib.rcParams['font.size'] = 7
matplotlib.rcParams['legend.frameon'] = False

matplotlib.rcParams['figure.dpi'] = 200
matplotlib.rcParams['axes.facecolor'] = '#F8F8FF'
matplotlib.rcParams['axes.grid'] = True
matplotlib.rcParams['axes.axisbelow'] = True
matplotlib.rcParams['axes.labelsize'] = 7

matplotlib.rcParams['grid.linestyle'] = '-.'
matplotlib.rcParams['grid.linewidth'] = 0.4

matplotlib.rcParams['xtick.major.size'] = 5
matplotlib.rcParams['xtick.major.width'] = 1
matplotlib.rcParams['xtick.minor.size'] = 3
matplotlib.rcParams['xtick.minor.width'] = 0.7

matplotlib.rcParams['ytick.major.size'] = 5
matplotlib.rcParams['ytick.major.width'] = 1
matplotlib.rcParams['ytick.minor.size'] = 3
matplotlib.rcParams['ytick.minor.width'] = 0.7
                   
matplotlib.rcParams['figure.figsize'] = (7, 3)
                   
                   
import cubehelix 
cb = cubehelix.cmap(startHue=220,endHue=-300,minSat=1,maxSat=2.5,minLight=.3,maxLight=.8,gamma=1.2)
cb3 = cubehelix.cmap(startHue=240,endHue=-300,minSat=1,maxSat=2.5,minLight=.3,maxLight=.8,gamma=.9)
cb2 = cubehelix.cmap(reverse=True, start=0., rot=0.5)
#cb2 = cubehelix.cmap(rot=1, reverse=True)

def martians_year(ls, data):
    #### only looking at "second year"
    idx = np.where(ls==360)[0]
    if idx[0] != 0 and idx.size > 1:
        idx0 = idx[0]
        idx1 = idx[1]
        return data[idx0:idx1]
    elif idx[0] != 0 and idx.size == 1:
        return data[idx[0]:]
    elif idx[0] == 0:
        idx1 = idx[1]
        idx2 = idx[2]
        return data[idx1:idx2]

def martians_month(ls, data):
    temp = []
    for i in np.arange(0, 12):
        idx = np.where((ls>i*30)&(ls<(i+1)*30))[0]
        temp.append(data[idx].mean(axis=0))
    temp = np.array(temp)
    return temp

def yearly_ls(ls):
    idx = np.where(ls==360)[0]
    counter = 1
    for i in np.arange(idx.size-1):
        ls[idx[i]:idx[i+1]] += counter*360
        ls[idx[i]] -= 360
        counter += 1
    if idx.size == 1:
        i = 0
        ls[idx[i]:] += counter*360
        ls[idx[i]] -= 360
    else:
        ls[idx[i+1]:] += counter*360
        ls[idx[i+1]] -= 360
    return ls

def window_stdev(arr, radius):
    #array = arr, radius = half width of window in bins
    #windows mean
    c1 = uniform_filter1d(arr, 2*radius, mode='constant', origin=-radius,axis=0)
    #windowed square mean
    c2 = uniform_filter1d(arr*arr, 2*radius, mode='constant', origin=-radius,axis=0)
    #std = windowed square mean - square(windowed mean)
    return ((c2 - c1*c1)**.5)[:-2*radius+1]

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5,axis=0):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data,axis=axis)
    return y

def spect_v(ls, data, tstep, lonstep, lowcut, highcut, wave):
    lowcut = 1./lowcut
    highcut = 1./highcut
        
    fftdata = np.fft.fftn(data, axes=[0,2])
    freq = np.fft.fftfreq(fftdata.shape[0], tstep)
    waven = np.fft.fftfreq(fftdata.shape[2], lonstep)*360
    
    idx1 = np.where((abs(freq)>lowcut)|(abs(freq)<highcut))[0]
#     print (idx1, waven)
    
#     set everything that satisfy condition as zero, ie only filtering storm system
    fftdata[idx1] = 0
    if wave == 0:
        wave1_idx= np.where(((waven)!=wave))[0]
        fftdata[:,:,wave1_idx] = 0
    elif wave == None:
        pass
    elif wave > 0:
        wave1_idx= np.where(((waven)!=wave))[0]
        fftdata[:,:,wave1_idx] = 0
        fftdata = fftdata*2

    filtered = np.fft.ifftn(fftdata, axes=[0,2])
    filtered = np.abs(filtered)**2

    temp = filtered.mean(axis=2)
    return np.sqrt(temp)

def T2km_filter_waven(directory):
    
    filedir = directory+'_auxhist5.nc'
    data = Dataset(filedir,'r')
    ls = data.variables['LS'][:]
    psfc = data.variables['TSK'][:]
    data.close()

    wave = [None,1,2,3]
    psfc = martians_year(ls, psfc)
    ls = martians_year(ls, ls)

    title = directory.replace('./../marswrf/test_data/reduction_diag.','')
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10,6), sharex=True)
    for i, ax in tqdm(enumerate(axes.flat)):
        tst = spect_v(ls, psfc,  1/8., 5., 1.5, 10., wave[i])

        im = ax.contourf(np.linspace(0,360, 223), np.linspace(-90,90,36), tst.T.reshape((36,223,24)).mean(axis=2), np.linspace(0,7,8), extend='both', cmap = cb)
        for c in im.collections:
                c.set_edgecolor("face")
                
        ax.set_title('{} TSK Amplitude [1.5-10 period] [s={}]'.format(title, wave[i]))
        if i+1 in [1,3]:
            ax.set_ylabel('Latitude [$^\circ$]')
        if i in [2,3]:
            ax.set_xlabel('Solar Longitude [$^\circ$]')
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    fig.tight_layout()
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.6, orientation='vertical', pad=0.01)
    
def zonal_plt_monthly(ydata, ls, data, title, level=9, norm=False, cmap=None):
    cmap=cmap or "viridis"
    
    fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(14,12))
    for i, ax in enumerate(axes.flat):
        y = ydata[i][6:]
        
        #press2 = zonal_p[i-1].mean(axis=1)
        lat = np.linspace(-90, 90, 36) 
        temp_press = np.linspace(1e-2, 900, ydata[i].shape[0])[6:]
        
        lat, temp_press = np.meshgrid(lat, temp_press)
        
        d = data[i][6:]

        if norm:
            im = ax.contourf(lat, y, d, levels=level, cmap=cmap, extend='both', norm=SymLogNorm(linthresh=np.abs(np.min(level)),vmin=np.min(d), vmax=np.max(d)))
#            if not np.isnan(d).any():
#                ax.contour(lat, y, d, levels=level, linewidths=0.5, colors='k', norm=SymLogNorm(linthresh=np.min(level),vmin=np.min(d), vmax=np.max(d)))
                
            ax.xaxis.set_minor_locator(AutoMinorLocator(4))
        else: 
            im = ax.contourf(lat, y, d, level, cmap=cmap, extend='both')
            for c in im.collections:
                c.set_edgecolor("face")
#             if not np.isnan(d).any():
#                 ax.contour(lat, y, d, level, linewidths=0.5, colors='k', extend='both')
        
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        ax.set_title(r'{} LS {}-{}'.format((title), (i)*10+180, (i+1)*10+180))
        if i in [0,4,8]: ax.set_ylabel('Pressure [Pa]')
        if i in [8,9,10,11]: ax.set_xlabel('Latitude [$^\circ$]')
        ax.set_yscale('log')
        ax.set_ylim([900, 1e-2])
        ax.set_xlim([-95, 95])
        ax.grid(True, which='both')
        
#        print ('Saving 1st cool shit')
    fig.tight_layout()
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.3, orientation='horizontal', pad=0.05)#, format='%.1e')

def bandpass_filter(filedir):
    
    data = Dataset(filedir, 'r')
    t = data.variables['T'][:]
    ls = data.variables['LS'][:]
    p = data.variables['P'][:].mean(axis=3)
    qice = data.variables['QICE'][:].mean(axis=3)
    data.close()
    
    t = martians_year(ls, t)
    qice = martians_year(ls, qice)
    p = martians_year(ls, p)
    ls = martians_year(ls, ls)
    
    fs = 8.
    lowcut = 1./10
    highcut = 1./1.5

    wlen = 1+20/0.25
    filtered = []
    for i in tqdm(np.arange(52)):
        tmp = t[:,i] - t[:,i].mean(axis=0)
        y = butter_bandpass_filter(tmp , lowcut, highcut, fs, order=4, axis=0)
        dy = (window_stdev(y, int(wlen/2)))
        filtered.append(dy)
    filtered = np.array(filtered).mean(axis=3)
    print (filtered.shape)
    
    dls = np.linspace(0, 360, dy.shape[0])
    #ls, press = np.meshgrid(dls, press)
    
    period = [220, 260, 300, 340]
#    period = [40, 80, 120, 160]
    ampl = []
    T = []
    u = []
    P = []
    dls = np.linspace(0, 360, dy.shape[0])
    ls = np.linspace(0,360,t.shape[0])
    for i in period:
        idx = np.where((dls>i)&(dls<(i+10)))[0]
        tmp = filtered[:,idx].mean(axis=1)
        ampl.append(tmp)
        
        idx2 = np.where((ls>i)&(ls<(i+10)))[0]
        T.append(t[idx2].mean(axis=0).mean(axis=2))
        P.append(p[idx2].mean(axis=0))
        u.append(qice[idx2].mean(axis=0))
        
    lat = np.linspace(-90,90,36)[22:]
    tmp = np.linspace(1e-2, 900,52)
    lat, tmp = np.meshgrid(lat, tmp)
    
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(8,10))
    for i, ax in enumerate(axes.flat):
        
        im2 = ax.contour(lat, P[i], ampl[i], levels = np.linspace(0,8,9), colors='k')
#        for c in im.collections:
#                c.set_edgecolor("face")
        im = ax.contourf(lat, P[i], u[i], np.logspace(-10,-4,7), cmap=cb, extend='both', norm=SymLogNorm(linthresh=1e-10,vmin=np.min(u[i]), vmax=np.max(u[i]))) #np.logspace(-9,-4,6)
        ax.contour(lat, P[i], T[i], np.linspace(110,240,14), colors='w')
        ax.set_title('Solar Longitude {}-{}{}'.format(period[i], period[i]+10, '$^\circ$'))
        
        plt.clabel(im2, fontsize=8, inline=True, use_clabeltext=True, fmt='%.3g')
        ax.set_yscale('log')
        ax.set_ylim([1000, 1e-2])
        if i in [0,2]: ax.set_ylabel('Pressure [Pa]')
        if i in [3,2]: ax.set_xlabel('Latitude [$^\circ$]')
        
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))

#    ax.set_title('Bandpass Filter Temperature at 2.5 km')
    fig.tight_layout()
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.4, orientation='vertical', pad=0.025)
#    plt.savefig('r14p1dustL45_butterworth_t_u_superimposed.pdf',  bbox_inches='tight', dpi=400)
    plt.show()