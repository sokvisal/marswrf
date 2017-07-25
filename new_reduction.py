import numpy as np
from netCDF4 import Dataset
import os
import glob
import sys
from tqdm import tqdm


class createNC:
    def __init__(self, name, solarLon):
#        dataset.close()
        self.name = name
        self.dataset = Dataset(self.name, 'w')
        
        
        solar_long = self.dataset.createDimension('time', solarLon.size)
        pressure = self.dataset.createDimension('bottom_top', 52)
        latitude = self.dataset.createDimension('south_north', 36)
        longitude = self.dataset.createDimension('west_east', 72)
        
        ls = self.dataset.createVariable('LS', np.float32, ('time',))
        lat = self.dataset.createVariable('LAT', np.float32, ('south_north',))
        long = self.dataset.createVariable('LONG', np.float32, ('west_east',))
        
        ls[:] = solarLon
            
    def saveVar(self, data, dataname, unitname):

        self.dataname = dataname
        self.unitname = unitname
        self.dim = data.ndim
        
        def create3D_var(varnameList, units, data):   
            tmp2 =  self.dataset.createVariable(varnameList[0], np.float32, (varnameList[1], varnameList[2], varnameList[3],), zlib=True)
            tmp2.units = (units)
            tmp2[:] = data
        def create4D_var(varnameList, units, data):   
            tmp2 =  self.dataset.createVariable(varnameList[0], np.float32, (varnameList[1], varnameList[2], varnameList[3], varnameList[4],), zlib=True)
            tmp2.units = (units)
            tmp2[:] = data
                
        if self.dim == 3 and data.shape[1] == 36:
            create3D_var([self.dataname, 'time', 'south_north', 'west_east'], self.unitname, data)
        if self.dim == 3 and data.shape[1] == 52:
            create3D_var([self.dataname, 'time', 'bottom_top', 'south_north'], self.unitname, data)
        
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
    
    def checkZM(self, data):
        if self.zonalmean and self.dim==4:
            return data.mean(axis=3)
        if self.zonalmean and self.dim==3:
            return data
        if not self.zonalmean:
            return data
        
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
            self.dim = temporary.ndim
    
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
                
                tmp.append( self.checkZM(pot_temp*(p/p0)**gamma) ) # temperature
                tmp.append( self.checkZM(p) )
            elif var == 'PH':
                ph = Data.variables['PH'][:]
                phb = Data.variables['PHB'][:]
                
                tmp.append( self.checkZM(ph+phb) )
            elif self.dim == 3:
                data = Data.variables[var][:][:,lat0:lat1,lon0:lon1]
                tmp.append( self.checkZM(data) )
            elif self.dim == 4: 
                data = Data.variables[var][:][:,:,lat0:lat1,lon0:lon1]
                tmp.append( self.checkZM(data) )
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
                
    def wrfout(self, varList):
        self.latRange = [0,None]
        self.lonRange = [0,None]
        self.varList = varList
        self.zonalmean = True
        
        ls = []
        for i, file in enumerate(sorted(glob.glob(self.dir+'/wrfout*')[:1])):
            print (file)
            if not i:
                ls, tmp = self.loadData(file)
            else:
                ls2, tmp2 = self.loadData(file)
                ls = np.concatenate((ls, ls2))
                tmp = tmp + tmp2
        
        print (self.varList)
        self.varList = np.concatenate(([varList[0]],['P'], varList[1:]))
#        self.varList = np.array(self.varList)
        print (self.varList)
        varlen = self.varList.size
        unitList = self.checkUnit(file)
        
        ncfile = createNC('r14p1dustL45_wrfout_5.nc', ls)
        for i, var in enumerate(self.varList):
            reshapedData = np.vstack(tmp[i::varlen])
            print ( 'Saving {} ...'.format(var), reshapedData.shape )
            ncfile.saveVar(reshapedData, self.varList[i], unitList[i])
        ncfile.close(True)
                
    def auxhist9(self, varList):
        self.latRange = [0,None]
        self.lonRange = [0,None]
        self.varList = varList
        self.zonalmean = True
        
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
        
        ncfile = createNC('r14p1dustL45_auxhist9_7.nc', ls)
        for i, var in enumerate(self.varList):
            reshapedData = np.vstack(tmp[i::varlen])
            print ( 'Saving {} ...'.format(var) )
            ncfile.saveVar(reshapedData, self.varList[i], unitList[i])
        ncfile.close(True)     
    
    def auxhist5(self, varList):
        self.latRange = [0,None]
        self.lonRange = [0,None]
        self.varList = varList
        self.zonalmean = True
        
        print(self.dir)
        
        ls = []
        for i, file in enumerate(sorted(glob.glob(self.dir+'/auxhist5*'))):
            print (file)
            if not i:
                ls, tmp = self.loadData(file)
            else:
                ls2, tmp2 = self.loadData(file)
                ls = np.concatenate((ls, ls2))
                tmp = tmp + tmp2
        
        varlen = self.varList.size
        unitList = self.checkUnit(file)
        
        ncfile = createNC('r14p1dustL45_auxhist5.nc', ls)
        for i, var in enumerate(self.varList):
            reshapedData = np.vstack(tmp[i::varlen])
            print ( 'Saving {} ...'.format(var) )
            ncfile.saveVar(reshapedData, self.varList[i], unitList[i])
        ncfile.close(True)   
        
a = reduction('./../diag.r14p1dustL45/')
#a.wrfout(np.array(['T', 'U', 'V']))
a.auxhist5(np.array(['PSFC', 'TSK', 'HGT']))
#a.auxhist9(np.array(['T_PHY_AM', 'T_PHY_PM', 'TAU_OD2D_AM', 'TAU_OD2D_PM', 'TAU_CL2D_AM', 'TAU_CL2D_PM']))