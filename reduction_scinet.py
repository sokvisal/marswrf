import numpy as np
from netCDF4 import Dataset
import os
import glob
import sys

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
	    print (varnameList, type(data))
        def create4D_var(varnameList, units, data):   
            tmp2 =  self.dataset.createVariable(varnameList[0], np.float32, (varnameList[1], varnameList[2], varnameList[3], varnameList[4],), zlib=True)
            tmp2.units = (units)
            tmp2[:] = data
                
        if self.dim == 3 and data.shape[1] == 36:
            create3D_var([self.dataname, 'time', 'south_north', 'west_east'], self.unitname, data)
        if self.dim == 3 and data.shape[1] == 52:
            create3D_var([self.dataname, 'time', 'bottom_top', 'south_north'], self.unitname, data)
        if self.dim == 4:
            create4D_var([self.dataname, 'time', 'bottom_top', 'south_north', 'west_east'], self.unitname, data)
 
    def close(self, close = True):
        if close:
            self.dataset.close()

class reduction:
    def __init__(self, directory):
        self.dir = directory
	self.filename = directory.replace('./../pw.v.wet/WRFV3/run/', '')        

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

    def reff_ice(self, nice, qice, qcore,mu=1.0, rhoi=1000., rhoc=2500.):
        import numpy as np
	print(np.max(nice))
        qtot = qice+qcore
        nlow=1e3
        rho = np.ones_like(qtot)#
        reff=np.ma.array(np.nan+np.zeros_like(qice),mask=nice<nlow)    
        m=(nice>0)&(qtot>0)
        rho[m] = (qice[m]*rhoi+qcore[m]*rhoc)/qtot[m]
        reff[m] = pow((qtot[m]/nice[m])*(3/(4*np.pi*rho[m]))*(mu+3)**2/((mu+2)*(mu+1)),1./3)
        print(type(reff))
	return reff*1e6
    
    def reff_dust(self, ndust, qdust,mu=1.0, rhoc=2500.):
        import numpy as np
        nlow=1e1
        rho = np.ones_like(qdust)#
        reff=np.ma.array(100+np.zeros_like(qdust),mask=ndust<nlow)    
        m=(ndust>0)&(qdust>0)
        reff[m] = pow((qdust[m]/ndust[m])*(3/(4*np.pi*rhoc))*(mu+3)**2/((mu+2)*(mu+1)),1./3)
        return reff*1e6

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
            elif var == 'REFF_ICE':
                qnice = Data.variables['QNICE'][:]
                qice = Data.variables['QICE'][:]
                trc_ic = Data.variables['TRC_IC'][:]
                
                r_ice = self.reff_ice(qnice, qice, trc_ic)
                del qnice, qice, trc_ic
                

		tmp2 = self.checkZM(r_ice)
		print type(tmp2)
                tmp.append( self.checkZM(r_ice) )
            elif var == 'REFF_DUST':
                qndust = Data.variables['NDUST'][:]
                qdust = Data.variables['TRC01'][:]
                
                r_dust = self.reff_dust(qndust, qdust)
                del qndust, qdust
                
                tmp.append( self.checkZM(r_dust) )
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
        for i, file in enumerate(sorted(glob.glob(self.dir+'/wrfout*'))):
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
        
        ncfile = createNC(self.filename+'_wrfout.nc', ls)
        for i, var in enumerate(self.varList):
            reshapedData = np.vstack(tmp[i::varlen])
            print ( 'Saving {} ...'.format(var), reshapedData.shape )
            ncfile.saveVar(reshapedData, self.varList[i], unitList[i])
        ncfile.close(True)
                
    def auxhist9(self, varList):
        self.latRange = [0,None]
        self.lonRange = [0,None]
        self.varList = varList
        self.zonalmean = False
        
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
        
        tdiff = (np.vstack(tmp[1::varlen]) - np.vstack(tmp[0::varlen]))*0.5        
        tavg = (np.vstack(tmp[1::varlen]) + np.vstack(tmp[0::varlen]))*0.5
        newData = [tdiff.mean(axis=3), tavg.mean(axis=3)]
        
        self.varList[:2] = ['T_PHY_DIFF', 'T_PHY_AVG']
        unitList[:2] = ['Temperature difference at 2pm and 2am', 'Temperature average at 2pm and 2am']
        
        ncfile = createNC(self.filename+'_auxhist9.nc', ls)
        for i, var in enumerate(self.varList):
            if i  in [0,1]:
                reshapedData = newData[i]
            else:
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
        print (unitList)
        ncfile = createNC(self.filename+'_auxhist5.nc', ls)
        for i, var in enumerate(self.varList):
            reshapedData = np.vstack(tmp[i::varlen])
            print ( 'Saving {} ...'.format(var) )
            ncfile.saveVar(reshapedData, self.varList[i], unitList[i])
        ncfile.close(True) 

    def auxhist8(self, varList):
        self.latRange = [0,None]
        self.lonRange = [0,None]
        self.varList = varList
        self.zonalmean = True

        print(self.dir)

        ls = []
        for i, file in enumerate(sorted(glob.glob(self.dir+'/auxhist8*'))):
            print (file)
            if not i:
                ls, tmp = self.loadData(file)
            else:
                ls2, tmp2 = self.loadData(file)
                ls = np.concatenate((ls, ls2))
                tmp = tmp + tmp2

        varlen = self.varList.size
        unitList = self.checkUnit(file)
        print (unitList)
        ncfile = createNC(self.filename+'_auxhist8.nc', ls)
        for i, var in enumerate(self.varList):
            reshapedData = np.vstack(tmp[i::varlen])
            print ( 'Saving {} ...'.format(var) )
            ncfile.saveVar(reshapedData, self.varList[i], unitList[i])
        ncfile.close(True)  
        
a = reduction('./../pw.v.wet/WRFV3/run/wetL50m8')
a.wrfout(np.array(['T', 'T1_5', 'NLIF1', 'QV_COLUMN', 'QI_COLUMN', 'REFF_ICE', 'QICE', 'QNICE', 'TRC_IC']))
#a.auxhist5(np.array(['PSFC', 'TSK', 'HGT']))
#a.auxhist9(np.array(['T_PHY_AM', 'T_PHY_PM', 'TAU_OD2D_AM', 'TAU_OD2D_PM', 'TAU_CL2D_AM', 'TAU_CL2D_PM']))
#a.auxhist8(np.array(['TRC_IC', 'DUSTN_SED', 'DUSTQ_SED', 'ICEQ_SED', 'ICEN_SED']))
