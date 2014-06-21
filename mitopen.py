#!/bin/python
# -*- coding: utf-8 -*-

"""Komod open module
Opens different file types, mainly those produced by MITgcm.

mitbincoord	 - opens MITgcm binary XC.data and YC.data files
rmeta		 - reads .meta file and return information about .data file
mitbin		 - opens MITgcm binary data file
phcascii	 - opens PHC (http://psc.apl.washington.edu/Climatology.html) ascii file that have no delimiter
nc2d		 - Convert 2d fields from adxx* and xx* fles to netCDF format with use of Nio module.
var_nc2d	 - Convert 2d fields produced by MITgcm to netCDF format with use of Nio module. 
nc3d		 - Convert 3d fields from adxx* and xx* fles to netCDF format with use of Nio module.
var_nc3d	 - Convert 2d fields produced by MITgcm to netCDF format with use of Nio module. 
gatrib 		 - Return attrubutes for known variables.
ncep2bin	 - Converts NCEP reanalysis data from netCDF to binary

Copyright (C) 2010 - 2014 Nikolay Koldunov <koldunovn@gmail.com> 

"""

# -------------------------------------------------

import numpy
import os
try:
	import Nio
except ImportError:
	pass
	#print('Nio is not avalible, some functions will not work')
import glob
from netCDF4 import Dataset


def mitbincoord(xdim, ydim, xcdata='./', ycdata='./', bswap=1):
	"""Opens MITgcm binary XC.data and YC.data files

    Usage: mitbincoord(xdim, ydim, [xcdata], [ycdata],  [endianness])

    Input:
	xdim        = x dimension
	ydim        = y dimension
        xcdata      = path to XC.data file [default ./]
        ycdata      = path to YC.data file [default ./]
	bswap       = do we need a byte swap? Yes (1) or no (0) [default 1]
        

    Output:
        arrays of longitude (lon) and latitude (lat) 

    """
        fd_xc = open(xcdata+'XC.data', 'rb')
	fd_yc = open(ycdata+'YC.data', 'rb') 
	size = xdim*ydim
	shape = (xdim,ydim)
	datatype = 'float32'
		
	lon = numpy.fromfile(file=fd_xc, dtype = datatype, count=size)
	lon = lon.reshape(shape)
	if bswap==1:
		lon = lon.byteswap()
		
	lat = numpy.fromfile(file=fd_yc, dtype = datatype, count=size)
	lat = lat.reshape(shape)
	if bswap==1:
		lat = lat.byteswap()
	
	return lon, lat

def rmeta(filename):
	''' Reads .meta file and return information about .data file.
	Usage: rmeta(filename)

	Input:
	    filename        = name of the .meta file (with .meta extention)
	            
	Output:
	    ndim		- Number of dimensions.
	    xdim		- Xdim
	    ydim		- Ydim
	    zdim		- Zdim
	    datatype		- Datatype
	    nrecords		- Number of records.
	    timeStepNumber	- Time step number

	'''
	
	ifile = open(filename, 'r')
	lines = ifile.readlines()  
	ifile.close()

	ndim = int(lines[0].split()[3])
	ydim = int(lines[2].split()[0][:-1])
	xdim = int(lines[3].split()[0][:-1])

	if ndim == 2:
		increm = 0
		zdim = 1
	elif ndim == 3:
		zdim = int(lines[4].split()[0][:-1])
		increm = 1
	else:
		print("unsupported number of dimensions")
	
	datatype	 = lines[5+increm].split()[3][1:-1]
	nrecords 	 = int(lines[6+increm].split()[3])
	
	if any("timeStepNumber" in s for s in lines):
		timeStepNumber   = int(lines[7+increm].split()[3])
	else:
		timeStepNumber = 1	
	
	return ndim, xdim, ydim, zdim, datatype, nrecords, timeStepNumber
	
def mitbin(filename, xdim, ydim, zdim=1, tdim=1, bswap=1, datatype = 'float64'):
	"""Opens MITgcm binary file that represents one time step.

    Usage: mitbin(filename, xdim, ydim, [zdim], [tdim], [bswap])

    Input:
	filename    = path to the file
	xdim        = x dimension
	ydim        = y dimension
	zdim        = z dimension [default 1]
	tdim        - time dimension [default 1]
	bswap       = do we need a byte swap? Yes (1) or no (0) [default 1]
	datatype    = data type, e.g., 'float32'. By default - float64
        
    Output:
       tdim*zdim*xdim*ydim array of data

    """
    	fd_data = open(filename, 'rb')
	
	size = tdim*zdim*xdim*ydim
	shape = (tdim,zdim,xdim,ydim)

        
	data = numpy.fromfile(file=fd_data, dtype = datatype, count=size)
	
        data = data.reshape(shape)
	if bswap==1:
		data = data.byteswap()
		
	fd_data.close()
	
	return data

def mitbin2(filename, bswap=1, meta=None):
	'''Uses rmeta to get inforamtion about the file and return field extracted from it.
	
	Usage: mitbin2(filename, [bswap], [meta])
	
	Input:
	    filename    - path to the file.
	    bswap       - do we need a byte swap? Yes (1) or no (0) [default 1]
	    meta	- None - flag to fix problem with wrong adxx*.meta files. 
			  If meta = 'xx', use .meta file from xx files 
	
	Output:
	    nrecords*zdim*xdim*ydim numpy array of data.
	'''
	
	fd_data = open(filename, 'rb')
	if meta == None:
		ndim, xdim, ydim, zdim, datatype, nrecords, timeStepNumber = rmeta(filename[:-4]+"meta")
	elif meta == 'xx':
		ndim, xdim, ydim, zdim, datatype, nrecords, timeStepNumber = \
		rmeta("".join(filename.split('adxx_')[:-1])+'xx_'+filename.split('adxx_')[-1][:-4]+"meta")
	
	size = nrecords*zdim*xdim*ydim
	shape = (nrecords,zdim,xdim,ydim)
	
	data = numpy.fromfile(file=fd_data, dtype = datatype, count=size)
	data = data.reshape(shape)
	if bswap==1:
		data = data.byteswap()
	
	fd_data.close()
	
	return data
		
	
		
			

def phcascii(filename, xdim, ydim,  nlines, nwords=10, lword=8,  zdim=1):
	"""open PHC (http://psc.apl.washington.edu/Climatology.html) ascii file that have no delimiter

    Usage: phcascii(filename, xdim, ydim, lword, nlines, [zdim])

    Input:
	filename    = path to the file
	xdim        = x dimension
	ydim        = y dimension
	nlines      = number of lines in input file
	nwords      = number of "words", chnks of the data that represent one number (df=10)
	lword       = length of every "word" (df=8)
	zdim        = z dimension [default 1]
	
        
    Output:
       zdim*xdim*ydim array of data

    """
    	in_file = open(filename,'r')
    
    	dd = numpy.zeros((nlines,nwords))
    
    	for i in range(nlines):
	
		for nn in range(nwords):
			dd[i,nn] = float(in_file.read(lword))
		
		in_file.read(1) # this one is needed to read "/n" at the end of each line

    	in_file.close()
    
    	data = dd.reshape(zdim,xdim,ydim)
	
	return data


def nc2d(parameters=['adxx_atemp'], ofile='adxx', iteration='0', bswap=1,
         sstart_date = "seconds since 2002-10-01 07:00", deltaT=1200, 
         xx_period=240000.0, FillValue=-1.0e+23, meta=None, dump='no'):
	'''
	Convert 2d fields from adxx* and xx* fles to netCDF format with use of Nio module.
	In order to convert variables (like T, S, AREA) use var_nc2d.  
	Names of the files should be defined in form of the list, even if we have only one variable.
		
	I assume that if file contain more than one record it is xx or adxx file.
	I put everything on the C grid!
	
	You have to have following files in the the directory where you run your code:

	XC.data
	XC.meta
	YC.data
	YC.meta
	maskCtrlC.data
	maskCtrlC.meta
	
	Input:
	    parameters		- list with names of the variables.
	    ofile 		- name of the output file.
	    iteration		- iteration of optimisation, should be STRING!
	    bswap       	- do we need a byte swap? Yes (1) or no (0) [default 1]
	    sstart_date		- should be "seconds since", [default "seconds since 2002-10-01 07:00"
	    deltaT		- time step in seconds
	    xx_period		- xx_*period
	    FillValue		- missing value
	    meta		- flag to fix problem with wrong adxx*.meta files. 
				  If meta = 'xx', use .meta file from xx files 
	    dump 		- if dump='yes' will return numpy array with data
	    	
	'''
	lon = mitbin2('XC.data',bswap)[0,0,:,:]
	lat = mitbin2('YC.data',bswap)[0,0,:,:]
	lsmask = mitbin2('maskCtrlC.data',bswap)[:,0,:,:]

	
	if os.path.exists(ofile+".nc") == True:
		os.system("rm "+ofile+".nc")
	
	if meta == None:
		ndim, xdim, ydim, zdim, datatype, nrecords, timeStepNumber = rmeta(parameters[0]+"."+iteration.zfill(10)+".meta")
	elif meta == 'xx':
		ndim, xdim, ydim, zdim, datatype, nrecords, timeStepNumber = rmeta(parameters[0][2:]+"."+iteration.zfill(10)+".meta")
	if nrecords == 1:
		ttime = numpy.zeros((nrecords))
		ttime[0] = timeStepNumber*deltaT
	elif nrecords > 1:
		ttime = numpy.zeros((nrecords))
		for i in range (nrecords):
	    		ttime[i] = xx_period*i
    
	opt = Nio.options()
	opt.PreFill = False
	opt.HeaderReserveSpace = 4000
	f = Nio.open_file(ofile+".nc","w",opt)

	f.title = "MITgcm variables in netCDF format"
	f.create_dimension('x',xdim)
	f.create_dimension('y',ydim)
	f.create_dimension('time',ttime.shape[0])

	f.create_variable('time','d',('time',))
	f.variables['time'].units        = sstart_date 
	f.variables['time'][:] = ttime


	f.create_variable('latitude','d',('x','y'))
	f.variables['latitude'].long_name        = "latitude"
	f.variables['latitude'].units            = "degrees_north"
	f.variables['latitude'].standard_name    = "grid_latitude"
	f.variables['latitude'][:] = lat[:]

	f.create_variable('longitude','d',('x','y'))
	f.variables['longitude'].long_name        = "longitude"
	f.variables['longitude'].units            = "degrees_east"
	f.variables['longitude'].standard_name    = "grid_longitude"
	f.variables['longitude'][:] = lon[:]

	#vvariables = ["atemp","aqh", "uwind", "vwind", ]
	#vvariables = ["atemp"]

	for parameter in parameters:
		adatemp = mitbin2(parameter+"."+iteration.zfill(10)+".data", bswap=bswap, meta=meta)[:,0,:,:]

	#	adatemp = numpy.where(adatemp[:] > 1.0e+12, 0, adatemp[:])
		adatemp = numpy.where(adatemp[:] < -1.0e+20, FillValue, adatemp[:])
		adatemp = numpy.where(lsmask[:]==0, FillValue, adatemp[:])
		
       		f.create_variable(parameter,'d',('time','x','y'))

        	f.variables[parameter].long_name    = parameter
        	f.variables[parameter].units        = "xz"
        	f.variables[parameter]._FillValue   = FillValue
		f.variables[parameter].missing_value = FillValue
		#print(adatemp.shape())
		
	        f.variables[parameter][:] = adatemp
		
	f.close()
	if dump == 'yes':
		return adatemp
	
	

def var_nc2d(parameters=['AREA','HEFF'], ofile='MIT_output_2d', bswap=1, sstart_date = "seconds since 2002-10-01 07:00", deltaT=1800, FillValue=-1.0e+23, dump='no'):
	'''
	Convert 2d fields produced by MITgcm to netCDF format with use of Nio module. 
	Names of the files should be defined in form of the list, even if we have only one variable.

	I put everything on the C grid!
	
	You have to have following files in the the directory where you run your code:

	XC.data
	XC.meta
	YC.data
	YC.meta
	maskCtrlC.data
	maskCtrlC.meta
	
	Input:
	    parameters		- list with names of the variables (like AREA or AREAtave).
	    ofile 		- name of the output file.
	    bswap       	- do we need a byte swap? Yes (1) or no (0) [default 1]
	    sstart_date		- should be "seconds since", [default "seconds since 2002-10-01 07:00"
	    deltaT		- time step in seconds
	    FillValue		- missing value
	    dump 		- if dump='yes' will return numpy array with data
	'''
	lon = mitbin2('XC.data',bswap)[0,0,:,:]
	lat = mitbin2('YC.data',bswap)[0,0,:,:]
	lsmask = mitbin2('maskCtrlC.data',bswap)[:,0,:,:]
	fileList = glob.glob(parameters[0]+"*.data")
	
	if os.path.exists(ofile+".nc") == True:
		os.system("rm "+ofile+".nc")
	
	ndim, xdim, ydim, zdim, datatype, nrecords, timeStepNumber = rmeta(fileList[0][:-4]+"meta")
	
	
	ttime = numpy.zeros((len(fileList)))
	#ttime[0] = timeStepNumber*deltaT
	
	opt = Nio.options()
	opt.PreFill = False
	opt.HeaderReserveSpace = 4000
	f = Nio.open_file(ofile+".nc","w",opt)

	f.title = "MITgcm variables in netCDF format"
	f.create_dimension('x',xdim)
	f.create_dimension('y',ydim)
	f.create_dimension('time',ttime.shape[0])

	f.create_variable('time','d',('time',))
	f.variables['time'].units        = sstart_date 
	


	f.create_variable('latitude','d',('x','y'))
	f.variables['latitude'].long_name        = "latitude"
	f.variables['latitude'].units            = "degrees_north"
	f.variables['latitude'].standard_name    = "grid_latitude"
	f.variables['latitude'][:] = lat[:]

	f.create_variable('longitude','d',('x','y'))
	f.variables['longitude'].long_name        = "longitude"
	f.variables['longitude'].units            = "degrees_east"
	f.variables['longitude'].standard_name    = "grid_longitude"
	f.variables['longitude'][:] = lon[:]

	for parameter in parameters:
		
		f.create_variable(parameter,'d',('time','x','y'))

        	f.variables[parameter].long_name    = gatrib(parameter)[0]
        	f.variables[parameter].units        = gatrib(parameter)[1]
        	f.variables[parameter]._FillValue   = FillValue
		f.variables[parameter].missing_value = FillValue

		adatemp_final = numpy.zeros((len(fileList), xdim, ydim))
	        
		iterator = 0
		for fileName in fileList:
			
			adatemp = mitbin2(parameter+fileName[-16:], bswap=bswap)[0,0,:,:]
			ndim, xdim, ydim, zdim, datatype, nrecords, timeStepNumber = rmeta(fileName[:-4]+"meta")
			adatemp = numpy.where(adatemp[:] < -1.0e+20, FillValue, adatemp[:])
			adatemp = numpy.where(lsmask[:]==0, FillValue, adatemp[:])
			adatemp_final[iterator,:,:] = adatemp
			ttime[iterator] = timeStepNumber*deltaT
			iterator = iterator + 1
			
		f.variables[parameter][:] = adatemp_final
 	
	f.variables['time'][:] = ttime
	f.close()
	if dump == 'yes':
		return adatemp
	
def nc3d(parameters=['adxx_atemp'], ofile='adxx', iteration='0', bswap=1, sstart_date = "seconds since 2002-10-01 07:00", deltaT=1200, xx_period=240000.0, FillValue=-1.0e+23, meta=None, dump="no"):
	'''
	Convert 3d fields from adxx* and xx* fles to netCDF format with use of Nio module.
	Names of the files should be defined in form of the list, even if we have only one variable.

	I put everything on the C grid!
	
	You have to have following files in the the directory where you run your code:

	XC.data
	XC.meta
	YC.data
	YC.meta
	DRC.data
	DRC.meta
	maskCtrlC.data
	maskCtrlC.meta
	
	Input:
	    parameters		- list with names of the variables.
	    ofile 		- name of the output file.
	    iteration		- iteration of optimisation, should be STRING!
	    bswap       	- do we need a byte swap? Yes (1) or no (0) [default 1]
	    sstart_date		- should be "seconds since", [default "seconds since 2002-10-01 07:00"
	    deltaT		- time step in seconds
	    xx_period		- xx_*period
	    FillValue		- missing value
	    meta		- flag to fix problem with wrong adxx*.meta files. 
				  If meta = 'xx', use .meta file from xx files 
	    dump 		- if dump='yes' will return numpy array with data
	'''
	lon = mitbin2('XC.data',bswap)[0,0,:,:]
	lat = mitbin2('YC.data',bswap)[0,0,:,:]
	lev = mitbin2('DRC.data',bswap)[0,:,0,0]
	lev = numpy.cumsum(lev)
	lsmask = mitbin2('maskCtrlC.data',bswap)[:,:,:,:]
	
	if os.path.exists(ofile+".nc") == True:
		os.system("rm "+ofile+".nc")
	
	if meta == None:
		ndim, xdim, ydim, zdim, datatype, nrecords, timeStepNumber = rmeta(parameters[0]+"."+iteration.zfill(10)+".meta")
	elif meta == 'xx':
		ndim, xdim, ydim, zdim, datatype, nrecords, timeStepNumber = rmeta(parameters[0][2:]+"."+iteration.zfill(10)+".meta")
	if nrecords == 1:
		ttime = numpy.zeros((nrecords))
		ttime[0] = timeStepNumber*deltaT
	elif nrecords > 1:
		ttime = numpy.zeros((nrecords))
		for i in range (nrecords):
	    		ttime[i] = xx_period*i
    
	opt = Nio.options()
	opt.PreFill = False
	opt.HeaderReserveSpace = 4000
	f = Nio.open_file(ofile+".nc","w",opt)

	f.title = "MITgcm variables in netCDF format"
	f.create_dimension('x',xdim)
	f.create_dimension('y',ydim)
	f.create_dimension('z',zdim)
	f.create_dimension('time',ttime.shape[0])

	f.create_variable('time','d',('time',))
	f.variables['time'].units        = sstart_date 
	f.variables['time'][:] = ttime

	f.create_variable('z','d',('z',))
	f.variables['z'].units        = "meters" 
	f.variables['z'][:] =  lev[:]

	
	f.create_variable('latitude','d',('x','y'))
	f.variables['latitude'].long_name        = "latitude"
	f.variables['latitude'].units            = "degrees_north"
	f.variables['latitude'].standard_name    = "grid_latitude"
	f.variables['latitude'][:] = lat[:]

	f.create_variable('longitude','d',('x','y'))
	f.variables['longitude'].long_name        = "longitude"
	f.variables['longitude'].units            = "degrees_east"
	f.variables['longitude'].standard_name    = "grid_longitude"
	f.variables['longitude'][:] = lon[:]

	#vvariables = ["atemp","aqh", "uwind", "vwind", ]
	#vvariables = ["atemp"]

	for parameter in parameters:
		adatemp = mitbin2(parameter+"."+iteration.zfill(10)+".data", bswap=bswap, meta=meta)[:,:,:,:]

	#	adatemp = numpy.where(adatemp[:] > 1.0e+12, 0, adatemp[:])
		adatemp = numpy.where(adatemp[:] < -1.0e+20, FillValue, adatemp[:])
		adatemp = numpy.where(lsmask[:]==0, FillValue, adatemp[:])
       		f.create_variable(parameter,'d',('time','z','x','y'))
		
		nname, unit, grid = gatrib(parameter)
		
        	f.variables[parameter].long_name    = nname
        	f.variables[parameter].units        = unit
		f.variables[parameter].grid         = grid
        	f.variables[parameter]._FillValue   = FillValue
		#print(adatemp.shape())

	        f.variables[parameter][:] = adatemp
		
	f.close()
	if dump == 'yes':
		return adatemp

		
def var_nc3d(parameters=['Ttave'], ofile='MIT_output_3d', bswap=1, sstart_date = "seconds since 2002-10-01 07:00", deltaT=1200, FillValue=-1.0e+23, dump="no"):
	'''
	Convert 3d fields produced by MITgcm to netCDF format with use of Nio module. 
	Names of the files should be defined in form of the list, even if we have only one variable.

	I put everything on the C grid!
	
	You have to have following files in the the directory where you run your code:

	XC.data
	XC.meta
	YC.data
	YC.meta
	DRC.data
	DRC.meta
	maskCtrlC.data
	maskCtrlC.meta
	
	Input:
	    parameters		- list with names of the variables.
	    ofile 		- name of the output file.
	    iteration		- iteration of optimisation, should be STRING!
	    bswap       	- do we need a byte swap? Yes (1) or no (0) [default 1]
	    sstart_date		- should be "seconds since", [default "seconds since 2002-10-01 07:00"
	    deltaT		- time step in seconds
	    xx_period		- xx_*period
	    FillValue		- missing value
	    meta		- flag to fix problem with wrong adxx*.meta files. 
				  If meta = 'xx', use .meta file from xx files 
	    dump 		- if dump='yes' will return numpy array with data
	'''
	lon = mitbin2('XC.data',bswap)[0,0,:,:]
	lat = mitbin2('YC.data',bswap)[0,0,:,:]
	lev = mitbin2('DRC.data',bswap)[0,:,0,0]
	lev = numpy.cumsum(lev)
	lsmask = mitbin2('maskCtrlC.data',bswap)[:,:,:,:]
	
	fileList = glob.glob(parameters[0]+"*.data")
	
	if os.path.exists(ofile+".nc") == True:
		os.system("rm "+ofile+".nc")
	
	ndim, xdim, ydim, zdim, datatype, nrecords, timeStepNumber = rmeta(fileList[0][:-4]+"meta")

	ttime = numpy.zeros((len(fileList)))

    
	opt = Nio.options()
	opt.PreFill = False
	opt.HeaderReserveSpace = 4000
	f = Nio.open_file(ofile+".nc","w",opt)

	f.title = "MITgcm variables in netCDF format"
	f.create_dimension('x',xdim)
	f.create_dimension('y',ydim)
	f.create_dimension('z',zdim)
	f.create_dimension('time',ttime.shape[0])

	f.create_variable('time','d',('time',))
	f.variables['time'].units        = sstart_date 


	f.create_variable('z','d',('z',))
	f.variables['z'].units        = "meters" 
	f.variables['z'][:] =  lev[:]

	
	f.create_variable('latitude','d',('x','y'))
	f.variables['latitude'].long_name        = "latitude"
	f.variables['latitude'].units            = "degrees_north"
	f.variables['latitude'].standard_name    = "grid_latitude"
	f.variables['latitude'][:] = lat[:]

	f.create_variable('longitude','d',('x','y'))
	f.variables['longitude'].long_name        = "longitude"
	f.variables['longitude'].units            = "degrees_east"
	f.variables['longitude'].standard_name    = "grid_longitude"
	f.variables['longitude'][:] = lon[:]

	#vvariables = ["atemp","aqh", "uwind", "vwind", ]
	#vvariables = ["atemp"]

	for parameter in parameters:
		f.create_variable(parameter,'d',('time','z','x','y'))
		
		f.variables[parameter].long_name    = gatrib(parameter)[0]
        	f.variables[parameter].units        = gatrib(parameter)[1]
        	f.variables[parameter]._FillValue   = FillValue
		f.variables[parameter].missing_value = FillValue
	  
		adatemp_final = numpy.zeros((len(fileList), zdim, xdim, ydim))
		
		for ind, fileName in enumerate(fileList):
		  adatemp = mitbin2(parameter+fileName[-16:], bswap=bswap)[:,:,:,:]
		  ndim, xdim, ydim, zdim, datatype, nrecords, timeStepNumber = rmeta(fileName[:-4]+"meta")

		#	adatemp = numpy.where(adatemp[:] > 1.0e+12, 0, adatemp[:])
		  adatemp = numpy.where(adatemp[:] < -1.0e+20, FillValue, adatemp[:])
		  adatemp = numpy.where(lsmask[:]==0, FillValue, adatemp[:])
		  adatemp_final[ind,:,:,:] = adatemp
		  ttime[ind] = timeStepNumber*deltaT
		  

	        f.variables[parameter][:] = adatemp
		
	f.close()
	if dump == 'yes':
		return adatemp 
		
def gatrib(parname):
	'''Return attrubutes for known variables'''
	
	if (parname == 'T') or (parname == 'Ttave'):
		name = 'Potential Temperature'
		unit = 'deg. C'
		grid = 'TS'
      
	elif (parname == 'S') or (parname == 'Stave'):
		name = 'Salinity'
		unit = 'PSU'
		grid = 'TS'
       
        elif (parname == 'U') or (parname == 'Utave'):
		name = 'Zonal Velocity'
		unit = 'm/s'
		grid = 'U'
       
        elif (parname == 'V') or (parname == 'Vtave'):
		name = 'Meridional Velocity'
		unit = 'm/s'
		grid = 'V'
       
	elif (parname == 'W') or (parname == 'Wtave'):
		name = 'Vertical Velocity'
		unit = 'm/s'
		grid = 'W'
        elif (parname == '"ETA') or (parname == 'ETAtave'):
		name = 'Sea Surface Height'
		unit = 'm'
		grid = 'TS'
      
        elif (parname == 'PHL') or (parname == 'PHLtave'):
		name = 'Bottom Dyn. Height Anom.'
		unit = 'm2/s2'
		grid = 'TS'
      
        elif (parname == 'UICE') or (parname == 'UICEtave'):
		name = 'Zonal Ice Velocity'
		unit = 'm/s'
		grid = 'U'
      
        elif (parname == 'VICE') or (parname == 'VICEtave'):
		name = 'Meridional Ice Velocity'
		unit = 'm/s'
		grid = 'V'
      
        elif (parname == 'HEFF') or (parname == 'HEFFtave'):
		name = 'Eff. Ice Thickness'
		unit = 'm'
		grid = 'TS'
     
	elif (parname == 'HSNOW') or (parname == 'HSNOW'):
		name = 'Eff. Snow Thickness'
		unit = 'm'
		grid = 'TS'

	elif (parname == 'AREA') or (parname == 'AREAtave'):
		name = 'Ice Concentration'
		unit = ' '
		grid = 'TS'

	elif (parname == 'HSALT') or (parname == 'HSALTtave'):
		name = 'Eff. Ice Salinity'
		unit = 'g/m2'
		grid = 'TS'

	elif (parname == 'QNET') or (parname == 'QNETtave'):
		name = 'Net Heat Flux'
		unit = 'W/m2'
		grid = 'TS'

	elif (parname == 'QSW') or (parname == 'QSWtave'):
		name = 'Net Shortwave Heat Flux'
		unit = 'W/m2'
		grid = 'TS'

	elif (parname == 'EmPmR') or (parname == 'EmPmRtave'):
		name = 'Net Freshwater Flux'
		unit = 'm/s'
		grid = 'TS'

	elif (parname == 'FU') or (parname == 'FUtave'):
		name = 'Zonal Wind Stress'
		unit = 'N/m2'
		grid = 'U'

	elif (parname == 'FV') or (parname == 'FVtave'):
		name = 'Meridional Wind Stress'
		unit = 'N/m2'
		grid = 'V'

	elif (parname == 'UWIND') or (parname == 'UWINDtave'):
		name = 'Zonal Wind'
		unit = 'm/s'
		grid = 'U'

	elif (parname == 'VWIND') or (parname == 'VWINDtave'):
		name = 'Meridional Wind'
		unit = 'm/s'
		grid = 'V'

	else:
		name = parname
		unit = 'some'
		grid = 'have no idea'

	return name, unit, grid



def ncep2bin(ifile, variable, bswap=1, coef=1.):
    '''Converts NCEP reanalysis data from netCDF file
    to binary file.

    Input:
    ifile - input file
    variable - name of the variable from the netCDF file
    bswap - do we need swipe bytes or not 
    coef - data will be multiplied by this coefficient
           (in case of prate we have to mult by 0.001)

    '''
    f = Dataset(ifile)
    #here ::-1 is for fliping data up side down.
    if coef != 1.:
    	a = numpy.float64(f.variables[variable][:,::-1,:])*coef
    	a = numpy.float32(a)
    else:
    	a = f.variables[variable][:,::-1,:]*coef
    
    if bswap==1:
        a = a.byteswap()
    a.tofile(ifile[:-2]+'bin')
    print('convert '+ifile+' to '+ifile[:-2]+'bin')

