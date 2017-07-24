import numpy as np
from netCDF4 import Dataset
import os
import glob
import sys
from tqdm import tqdm

class reduction:
    def __init__(self, directory):
        self.dir = directory
        
    def fullfield(self):
        for i in sorted(glob.glob(self.dir)):
            print (i)
            
    def loadData(self, file):
        t0 = 300. #data.variables['T00'][:] # base temperature
        p0 = 610. #data.variables['P00'][:] # base pressure
        r_d = 191.8366
        cp = 767.3466
        g = 3.727
        gamma = r_d/cp
        
        Data = Dataset(file, 'r')
        ls = Data.variables['L_S'][:] # solar long
        
        Data.close()
        
        for var in self.varList:
            print (var)
#            dim =  Data.variables[var][:].ndim
#    
#            lat0 = self.latRange[0]
#            lat1 = self.latRange[1]
#            
#            lon0 = self.lonRange[0]
#            lon1 = self.lonRange[1]
#            
#            if var == 'T':
#                t = Data.variables['T'][:][:,:,lat0:lat1,lon0:lon1] # perturbation potential temp
#                
#                p = Data.variables['P'][:][:,:,lat0:lat1,lon0:lon1] # perturbation pressure
#                pb = Data.variables['PB'][:][:,:,lat0:lat1,lon0:lon1] # base state pressure
#                
#                pot_temp = t + t0 # potential temperature
#                p = p + pb # pressure
#                del t, pb
#                
#                tmp.append((pot_temp*(p/p0)**gamma)) # temperature
#                tmp.append(p)
#            elif var == 'PH':
#                ph = data.variables['PH'][:]
#                phb = data.variables['PHB'][:]
#                
#                tmp.append((ph+phb).mean(axis = 3))
#            elif dim == 3:
#                u = data.variables[var][:][:,16:20,60:64]
#                tmp.append(u)
#            else: 
#                tmp2 = data.variables[var][:]
#                tmp.append(tmp2[:])
#            
#            return ls, tmp
            
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
                ls = ls + lsp2
                tmp = tmp + tmp2
                
        
a = reduction('./../diag.r14p1dustL45/')
a.galeCrater(np.array(['T','PH', 'TAU_CL', 'TAU_OD']))