import numpy as np
from netCDF4 import Dataset
import os
import glob
import sys
from tqdm import tqdm


class createNC:
    def __init__(self, name, solarLong):
#        dataset.close()
        self.name = name
        self.dataset = Dataset(self.name, 'w')
        
        solar_long = self.dataset.createDimension('time', None)
        latitude = self.dataset.createDimension('south_north', None)
        longitude = self.dataset.createDimension('west_east', None)
        
        ls = self.dataset.createVariable('LS', np.float32, ('time',))
        lat = self.dataset.createVariable('LAT', np.float32, ('south_north',))
        long = self.dataset.createVariable('LONG', np.float32, ('west_east',))
        
        ls[:] = solarLong
            
    def checkDim(self, data, dataname, unitname):

        self.dataname = dataname
        self.unitname = unitname
        
        def create3D_var(varnameList, units, data):   
            tmp2 =  self.dataset.createVariable(varnameList[0], np.float32, (varnameList[1], varnameList[2], varnameList[3],), zlib=True)
            tmp2.units = (units)
            tmp2[:] = data
        def create4D_var(varnameList, units, data):   
            tmp2 =  self.dataset.createVariable(varnameList[0], np.float32, (varnameList[1], varnameList[2], varnameList[3], varnameList[4],), zlib=True)
            tmp2.units = (units)
            tmp2[:] = data
                
        if data.ndim == 3:
            create3D_var([self.dataname, 'time', 'south_north', 'west_east'], self.unitname, data)
        
    def close(self, close = True):
        if close:
            self.dataset.close()

class reduction:
    def __init__(self, directory):
        self.dir = directory
        
    def checkUnit(self, file):
        Data = Dataset(file, 'r')
         
        tmp = []
        for var in self.varList:
           tmp.append( Data.variables[var].description )
        return tmp    
        
    def loadData(self, file):
        t0 = 300. #data.variables['T00'][:] # base temperature
        p0 = 610. #data.variables['P00'][:] # base pressure
        r_d = 191.8366
        cp = 767.3466
        g = 3.727
        gamma = r_d/cp
        
        Data = Dataset(file, 'r')
        ls = Data.variables['L_S'][:] # solar long
        
        tmp = []
        for var in self.varList:
            temporary =  Data.variables[var][:]
            dim = temporary.ndim
    
            lat0 = self.latRange[0]
            lat1 = self.latRange[1]
            
            lon0 = self.lonRange[0]
            lon1 = self.lonRange[1]
            
            if var == 'T':
                t = Data.variables['T'][:][:,:,lat0:lat1,lon0:lon1] # perturbation potential temp
                
                p = Data.variables['P'][:][:,:,lat0:lat1,lon0:lon1] # perturbation pressure
                pb = Data.variables['PB'][:][:,:,lat0:lat1,lon0:lon1] # base state pressure
                
                pot_temp = t + t0 # potential temperature
                p = p + pb # pressure
                del t, pb
                
                tmp.append((pot_temp*(p/p0)**gamma)) # temperature
                tmp.append(p)
            elif var == 'PH':
                ph = Data.variables['PH'][:]
                phb = Data.variables['PHB'][:]
                
                tmp.append((ph+phb).mean(axis = 3))
            elif dim == 3:
                u = Data.variables[var][:][:,lat0:lat1,lon0:lon1]
                tmp.append(u)
            else: 
                tmp2 = Data.variables[var][:]
                tmp.append(tmp2[:])
        return ls, tmp
            
    def galeCrater(self, varList):
        self.latRange = [16,20]
        self.longRange = [60,64]
        self.varList = varList
        
        print(self.dir)
        
        ls = []
        for i, file in enumerate(sorted(glob.glob(self.dir+'/wrfout*')[:1])):
            print (file)
            if not i:
                ls, tmp = self.loadData(file)
            else:
                ls2, tmp2 = self.loadData(file)
                ls = ls + ls2
                tmp = tmp + tmp2
                
    def auxhist9(self, varList):
        self.latRange = [0,-1]
        self.lonRange = [0,-1]
        self.varList = varList
        
        print(self.dir)
        
        ls = []
        for i, file in enumerate(sorted(glob.glob(self.dir+'/auxhist9*'))):
            print (file)
            if not i:
                ls, tmp = self.loadData(file)
            else:
                ls2, tmp2 = self.loadData(file)
                ls = np.concatenate((ls, ls2))
                tmp = tmp + tmp2
        
        varlen = self.varList.size
        unitList = self.checkUnit(file)
        
        ncfile = createNC('test8.nc', ls)
        for i, var in enumerate(self.varList):
            reshapedData = np.vstack(tmp[i::varlen])
            ncfile.checkDim(reshapedData, self.varList[i], unitList[i])
        ncfile.close(True)     
        
a = reduction('./../data_marswrf/diag.r14p1dustL40/data/')
a.auxhist9(np.array(['TAU_OD2D_PM', 'TAU_CL2D_PM']))