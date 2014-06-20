#!/bin/python
# -*- coding: utf-8 -*-


"""Komod plot module 
Contain mostly set of wrapper functions for map plotting with [PyNGL](http://www.pyngl.ucar.edu/).
Can be used with any 2D data, not necessarily MITgcm.

coord2d 		- Convert 1d coordinates to 2d coordinates
colormap 		- Define custom colormaps
arctpl 			- plot contours of 2D, 3D or 4D field in the Arctic region from ndarray
arctpltnc 		- plot contours of 2D, 3D or 4D field in the Arctic region from netCDF file
globpltnc 		- plot contours of 2D, 3D or 4D field on the global map from netCDF file. You can also specify custom region.
pltgrd			- plot model grid from lat and lon arrays
pltgrdnc		- plot model grid from netCDF file
plt_vectors     	- plot vector data from 2D fields
plt_vectors_scalars 	- plot vectors and scalar filed (e.g. wind direction and speed)
plt_vectors_colors      - plot colored vectors
 
Copyright (C) 2010 - 2013 Nikolay Koldunov <koldunovn@gmail.com> 
"""

# -------------------------------------------------
try:
	import Ngl
except ImportError:
	pass
	#print('Ngl is not avalible, some functions will not work')

try:
	import Nio
except ImportError:
	pass
	#print('Nio is not avalible, some functions will not work')

import numpy
import os
from netcdftime import num2date
import matplotlib.pyplot as plt


def coord2d(lon, lat, dshape):
	""" Convert 1d coordinates to 2d coordinates.
	This function was created before I find out about numpy.meshgrid :)
	Keep it for backward compatibility.
	
	Usage; coord2d(lat, lon, dshape)
	Input:
	    lat - 1d latitudes
	    lon - 1d longitudes
	    dshape - shape of the data field (x,y)
	Output:
	    2D lat and lon fields
	"""
	lonn = numpy.zeros(dshape)
	latt = numpy.zeros(dshape)
	
		    
	if len(lon) == dshape[0]:
	  for i in range(dshape[1]):
	    lonn[:,i] = lon[:]
	elif len(lon) == dshape[1]:
	  for i in range(dshape[0]):
	    lonn[i,:] = lon[:]
	
	if len(lat) == dshape[0]:
	  for i in range(dshape[1]):
	    latt[:,i] = lat[:]
	elif len(lat) == dshape[1]:
	  for i in range(dshape[0]):
	    latt[i,:] = lat[:]
	
	return lonn, latt 


def colormap(colormapname):
	""" Define custom colormaps.
	
	Usage: colormap(colormapname)
	Input:
		colormapname - name of the custom color map
	Output:
		return array that represents custom color map, 
		or the name of the colormap if there is no custom colormap with specified name.
	"""
	if colormapname == 'ice_conc':
		cmap = numpy.array([[1.00, 1.00, 1.00], [0.00, 0.00, 0.00], \
		    		[1.00, 1.00, 1.00],[0.34,0.53,1.00], \
                    		[0.13,0.75,1.00],[0.16,1.00,1.00], \
		    		[0.33,1.00,1.00],[0.66,1.00,1.00], \
		    		[1.00,1.00,0.33],[1.00,0.75,0.00], \
		    		[1.00,0.54,0.00],[1.00,0.00,0.00]],'f')
				
	elif colormapname == 'spread':
		cmap = numpy.zeros((104,3),'f')
		cmap[0] = [1.,1.,1.]
		cmap[1] = [0.,0.,0.]
		cmap[2] = [.5,.5,.5]
		cmap[3] = [.8,.8,.8]

		iofc = 151
		iolc = 250
		for i in xrange(151,251):
 			p = (1.*iolc-i)/(1.*iolc-1.*iofc)
  			q = (i-1.*iofc)/(1.*iolc-1.*iofc)
  			cmap[i-147] = [0.,p,q]
			

	else:
		cmap = colormapname
				
	return cmap

def reg(region, minLon=0, maxLon=360 , minLat=-80 , maxLat=85):
  
	mapDict = {}
	
  	if region == 'Arctic':
		mapDict.update({"mpProjection":"LambertEqualArea"}) 
		mapDict.update({"mpDataBaseVersion":"LowRes"}) 
		mapDict.update({"mpGeophysicalLineThicknessF": 4}) 
		mapDict.update({"mpLimitMode":"Corners"}) 
		mapDict.update({"mpLeftCornerLatF":54}) 
		mapDict.update({"mpLeftCornerLonF":-50}) 
		mapDict.update({"mpRightCornerLatF":57}) 
		mapDict.update({"mpRightCornerLonF":140}) 
		mapDict.update({"mpCenterLonF":0.}) 
		mapDict.update({"mpCenterLatF":90})
		mapDict.update({"mpMinLonF":minLon})
		mapDict.update({"mpMaxLonF":maxLon})
		mapDict.update({"mpMinLatF":minLat})
		mapDict.update({"mpMaxLatF":maxLat})
		mapDict.update({"mpMaxLatF":maxLat})
		mapDict.update({"mpGridLonSpacingF":15})
		mapDict.update({"mpGridLatSpacingF":5})

	elif region == 'Global':
		
		mapDict.update({"mpProjection":"CylindricalEquidistant"}) 
		mapDict.update({"mpDataBaseVersion":"LowRes"}) 
		mapDict.update({"mpGeophysicalLineThicknessF": 2}) 
		mapDict.update({"mpLimitMode":"LatLon"}) 
		mapDict.update({"mpLeftCornerLatF":0}) 
		mapDict.update({"mpLeftCornerLonF":0}) 
		mapDict.update({"mpRightCornerLatF":0}) 
		mapDict.update({"mpRightCornerLonF":0}) 
		mapDict.update({"mpCenterLonF":180.}) 
		mapDict.update({"mpCenterLatF":0})
		mapDict.update({"mpMinLonF":minLon})
		mapDict.update({"mpMaxLonF":maxLon})
		mapDict.update({"mpMinLatF":minLat})
		mapDict.update({"mpMaxLatF":maxLat})
		mapDict.update({"mpMaxLatF":maxLat})
		mapDict.update({"mpGridLonSpacingF":15})
		mapDict.update({"mpGridLatSpacingF":5})
		
		
	elif region == 'NAtlantic':
		
		mapDict.update({"mpProjection":"LambertEqualArea"}) 
		mapDict.update({"mpDataBaseVersion":"MediumRes"}) 
		mapDict.update({"mpGeophysicalLineThicknessF": 4}) 
		mapDict.update({"mpLimitMode":"Corners"}) 
		mapDict.update({"mpLeftCornerLatF":65}) 
		mapDict.update({"mpLeftCornerLonF":-20}) 
		mapDict.update({"mpRightCornerLatF":70}) 
		mapDict.update({"mpRightCornerLonF":100}) 
		mapDict.update({"mpCenterLonF":0.}) 
		mapDict.update({"mpCenterLatF":90})
		mapDict.update({"mpMinLonF":minLon})
		mapDict.update({"mpMaxLonF":maxLon})
		mapDict.update({"mpMinLatF":minLat})
		mapDict.update({"mpMaxLatF":maxLat})
		mapDict.update({"mpMaxLatF":maxLat})
		mapDict.update({"mpGridLonSpacingF":15})
		mapDict.update({"mpGridLatSpacingF":5})
		
	elif region == 'NAtlanticOcean':
		
		mapDict.update({"mpProjection":"Mercator"}) 
		mapDict.update({"mpDataBaseVersion":"MediumRes"}) 
		mapDict.update({"mpGeophysicalLineThicknessF": 4}) 
		mapDict.update({"mpLimitMode":"LatLon"}) 
		mapDict.update({"mpLeftCornerLatF":-65}) 
		mapDict.update({"mpLeftCornerLonF":-20}) 
		mapDict.update({"mpRightCornerLatF":70}) 
		mapDict.update({"mpRightCornerLonF":40}) 
		mapDict.update({"mpCenterLonF":0.}) 
		mapDict.update({"mpCenterLatF":0})
		mapDict.update({"mpMinLonF":-100})
		mapDict.update({"mpMaxLonF":20})
		mapDict.update({"mpMinLatF":0})
		mapDict.update({"mpMaxLatF":75})
		mapDict.update({"mpGridLonSpacingF":15})
		mapDict.update({"mpGridLatSpacingF":5})
	
	elif region == 'stanna':
	  
	  	mapDict.update({"mpProjection":"LambertEqualArea"}) 
		mapDict.update({"mpDataBaseVersion":"MediumRes"}) 
		mapDict.update({"mpGeophysicalLineThicknessF": 4}) 
		mapDict.update({"mpLimitMode":"Corners"}) 
		mapDict.update({"mpLeftCornerLatF":80}) 
		mapDict.update({"mpLeftCornerLonF":25}) 
		mapDict.update({"mpRightCornerLatF":70}) 
		mapDict.update({"mpRightCornerLonF":100}) 
		mapDict.update({"mpCenterLonF":0.}) 
		mapDict.update({"mpCenterLatF":90})
		mapDict.update({"mpMinLonF":minLon})
		mapDict.update({"mpMaxLonF":maxLon})
		mapDict.update({"mpMinLatF":minLat})
		mapDict.update({"mpMaxLatF":maxLat})
		mapDict.update({"mpMaxLatF":maxLat})
		mapDict.update({"mpGridLonSpacingF":15})
		mapDict.update({"mpGridLatSpacingF":5})
	
	elif region == "Monarch":
		
	  	mapDict.update({"mpProjection":"LambertEqualArea"}) 
		mapDict.update({"mpDataBaseVersion":"LowRes"}) 
		mapDict.update({"mpGeophysicalLineThicknessF": 4}) 
		mapDict.update({"mpLimitMode":"Corners"}) 
		mapDict.update({"mpLeftCornerLatF":33}) 
		mapDict.update({"mpLeftCornerLonF":-70}) 
		mapDict.update({"mpRightCornerLatF":40}) 
		mapDict.update({"mpRightCornerLonF":110}) 
		mapDict.update({"mpCenterLonF":-34.}) 
		mapDict.update({"mpCenterLatF":90})
		mapDict.update({"mpMinLonF":minLon})
		mapDict.update({"mpMaxLonF":maxLon})
		mapDict.update({"mpMinLatF":minLat})
		mapDict.update({"mpMaxLatF":maxLat})
		mapDict.update({"mpMaxLatF":maxLat})
		mapDict.update({"mpGridLonSpacingF":15})
		mapDict.update({"mpGridLatSpacingF":5})
	  
	elif region == 'halo':
	      	
	      	mapDict.update({"mpProjection":"LambertEqualArea"}) 
		mapDict.update({"mpDataBaseVersion":"MediumRes"}) 
		mapDict.update({"mpGeophysicalLineThicknessF": 4}) 
		mapDict.update({"mpLimitMode":"Corners"}) 
		mapDict.update({"mpLeftCornerLatF":85}) 
		mapDict.update({"mpLeftCornerLonF":0}) 
		mapDict.update({"mpRightCornerLatF":65}) 
		mapDict.update({"mpRightCornerLonF":120}) 
		mapDict.update({"mpCenterLonF":0.}) 
		mapDict.update({"mpCenterLatF":90})
		mapDict.update({"mpMinLonF":minLon})
		mapDict.update({"mpMaxLonF":maxLon})
		mapDict.update({"mpMinLatF":minLat})
		mapDict.update({"mpMaxLatF":maxLat})
		mapDict.update({"mpMaxLatF":maxLat})
		mapDict.update({"mpGridLonSpacingF":15})
		mapDict.update({"mpGridLatSpacingF":5})
	  
	elif region == 'FramStAnna':
	  
	  	mapDict.update({"mpProjection":"LambertEqualArea"}) 
		mapDict.update({"mpDataBaseVersion":"MediumRes"}) 
		mapDict.update({"mpGeophysicalLineThicknessF": 4}) 
		mapDict.update({"mpLimitMode":"Corners"}) 
		mapDict.update({"mpLeftCornerLatF":75}) 
		mapDict.update({"mpLeftCornerLonF":-5}) 
		mapDict.update({"mpRightCornerLatF":75}) 
		mapDict.update({"mpRightCornerLonF":120}) 
		mapDict.update({"mpCenterLonF":0.}) 
		mapDict.update({"mpCenterLatF":90})
		mapDict.update({"mpMinLonF":minLon})
		mapDict.update({"mpMaxLonF":maxLon})
		mapDict.update({"mpMinLatF":minLat})
		mapDict.update({"mpMaxLatF":maxLat})
		mapDict.update({"mpMaxLatF":maxLat})
		mapDict.update({"mpGridLonSpacingF":15})
		mapDict.update({"mpGridLatSpacingF":5})

	elif region == 'Fram2':
	  
	  	mapDict.update({"mpProjection":"LambertEqualArea"}) 
		mapDict.update({"mpDataBaseVersion":"MediumRes"}) 
		mapDict.update({"mpGeophysicalLineThicknessF": 4}) 
		mapDict.update({"mpLimitMode":"Corners"}) 
		mapDict.update({"mpLeftCornerLatF":75}) 
		mapDict.update({"mpLeftCornerLonF":0}) 
		mapDict.update({"mpRightCornerLatF":77}) 
		mapDict.update({"mpRightCornerLonF":90}) 
		mapDict.update({"mpCenterLonF":0.}) 
		mapDict.update({"mpCenterLatF":90})
		mapDict.update({"mpMinLonF":minLon})
		mapDict.update({"mpMaxLonF":maxLon})
		mapDict.update({"mpMinLatF":minLat})
		mapDict.update({"mpMaxLatF":maxLat})
		mapDict.update({"mpMaxLatF":maxLat})
		mapDict.update({"mpGridLonSpacingF":15})
		mapDict.update({"mpGridLatSpacingF":5})

			
	return mapDict

				
def arctpl(lon, lat, data, datamin=None, datamax=None, datastep=None,\
    showfig=True, psname="output", colormap_name='testcmap', start_color=2,\
    end_color=-2, vtitle="values", raster_fill=False, miss=None,levon=False,\
    llevel=None, add_cyclic=False, region='Arctic', \
    minLon=0, maxLon=360 , minLat=-80 , maxLat=85):
	""" Plot variable from array in the Arctic region.
	
	Usage:
		arctpl(lon, lat, data, datamin=None, datamax=None, datastep=None,\
		 showfig=True, psname="output", colormap_name='posneg_1', \
		 start_color=2, end_color=-2, vtitle="values"):
	
	Input:
		lon 		- 2D array of longitudes
		lat 		- 2D array of latitudes
		data 		- 2D, 3D or 4D array of scalar data.
				  Two last dimensions should be lon/lat
		datamin 	- minimum value, if not specified will be 
				  calculated as data.min() for the first time step and level
		datamax 	- maximum value, if not specified will be 
				  calculated as data.max() for the first time step and level
		datastep 	- interval between isolines, if not specified will be
				  calculated as abs(datamax-datamin)/20
		showfig 	- if True display the figure with gv
		psname 		- name of the output .ps file, default - "output"
		colormap_name 	- name of the PyNGL of custom colormap. 
				  Custom colormaps should be defined in komod.colormap function.
				  Color table for PyNGL can be found here 
				  http://www.pyngl.ucar.edu/Graphics/color_table_gallery.shtml
		start_color 	- number of the color in the PyNGL 
				  or custom colormap to start from, default = 2
		end_color  	- number of the color in the PyNGL 
				  or custom colormap to end with, default = -2
		vtitle 		- units of the data
		raster_fill	- switch between AreaFill and RasterFill.
		miss            - missing value
		levon           - do we need to plot level value in header?
		llevel          - data vector with values of levels (levon should be True)
		add_cyclic      - add cyclic points
		region 		- one of predefined regions (for list of regions see the reg function)
		minLon      - if region="Global", you can select specific part of the globe 
					  by specifying min/max lon/lat. Lons should be in 0-360 format, and
					  plotting over 0 meredian is still a problem. 
		maxLon
		minLat
		maxLat
		
	Output:
		.ps file, output.ps by default
	"""

	if lon.shape.__len__() == 1:
	  lon,lat = numpy.meshgrid(lon,lat)
	
	if add_cyclic==True:
	  lon = Ngl.add_cyclic(lon)
	  lat = Ngl.add_cyclic(lat)
	
	  
	rlist            = Ngl.Resources()
	rlist.wkColorMap = colormap(colormap_name)
	wks_type = "ps"
	wks = Ngl.open_wks(wks_type,psname,rlist)
	
	resources = Ngl.Resources()
	
	resources.sfXArray        = lon[:]
	resources.sfYArray        = lat[:]
	
	mapDict = reg(region, minLon, maxLon , minLat , maxLat)
	
	resources.mpProjection          = mapDict['mpProjection']
		
	resources.mpDataBaseVersion     = mapDict['mpDataBaseVersion']
	resources.mpGeophysicalLineThicknessF = mapDict['mpGeophysicalLineThicknessF']
	resources.mpLimitMode           = mapDict['mpLimitMode']
	resources.mpLeftCornerLatF      = mapDict['mpLeftCornerLatF']
	resources.mpLeftCornerLonF      = mapDict['mpLeftCornerLonF']
	resources.mpRightCornerLatF     = mapDict['mpRightCornerLatF']
	resources.mpRightCornerLonF     = mapDict['mpRightCornerLonF']
	resources.mpCenterLonF          = mapDict['mpCenterLonF']
	resources.mpCenterLatF          = mapDict['mpCenterLatF']
	resources.mpMinLonF             = mapDict['mpMinLonF']
	resources.mpMaxLonF             = mapDict['mpMaxLonF']
	resources.mpMinLatF             = mapDict['mpMinLatF']
	resources.mpMaxLatF             = mapDict['mpMaxLatF']
	resources.mpGridLonSpacingF     = mapDict['mpGridLonSpacingF']
	resources.mpGridLatSpacingF     = mapDict['mpGridLatSpacingF']
		
		

		
	resources.mpShapeMode    = 'FixedAspectFitBB'
	resources.cnFillDotSizeF    = 1
	

	
	resources.mpFillOn     = True 
	igray = Ngl.new_color(wks,0.7,0.7,0.7)
	
	resources.mpFillColors = [0,-1,igray,-1]
	
	resources.cnLineDrawOrder       = "Predraw"
	
	resources.cnFillOn              = True
	resources.cnFillDrawOrder       = "Predraw"
	resources.cnLineLabelsOn        = False
	resources.nglSpreadColorStart   = start_color
	resources.nglSpreadColorEnd     = end_color
	
	
	resources.cnLevelSelectionMode  = "ExplicitLevels" # Define own levels.
	#resources.cnLevels              = numpy.arange(datamin, datamax, datastep)
	
	resources.lbTitleString             = vtitle
	resources.lbTitleFontHeightF        = 0.022
	resources.lbLabelFontHeightF        = 0.018
	resources.lbTitleOffsetF            = -0.30
	resources.lbBoxMinorExtentF         = 0.15
	resources.pmLabelBarOrthogonalPosF  = -0.0
	resources.lbOrientation             = "Horizontal"
	
	
	if raster_fill==True:
		resources.cnFillMode           = "RasterFill"
		
	resources.cnLinesOn             = False
	resources.cnMaxDataValueFormat  = ".4f"
	resources.mpGridAndLimbOn       = False
	#resources.mpLabelsOn           = True
	#resources.mpOutlineOn          = False
	resources.tmXTOn                = False
	resources.tmXBOn                = False
	resources.tmYLOn                = False
	resources.tmYROn                = False
	
	if miss is not None:
		resources.sfMissingValueV = miss
	elif hasattr(data, "fill_value"):
		resources.sfMissingValueV = float(data.fill_value)
	
	
	if data.shape.__len__() == 2:
		if datamin is None:
			datamin = data.min()
		if datamax is None:
			datamax = data.max()
		if datastep is None:
			datastep = abs(datamax-datamin)/20.
	
		resources.cnLevelSelectionMode  = "ExplicitLevels" # Define own levels.
		levels = numpy.arange(datamin, datamax, datastep)
		levels = numpy.where((levels<0.00000001)&(levels>-0.00000001),0,levels)
		resources.cnLevels              = levels
		
		map = Ngl.contour_map(wks,data[:,:],resources)
		
	
	if data.shape.__len__() == 3:
		for ttime in range(data.shape[0]):
			if datamin is None:
				datamin = data[ttime,:,:].min()

			if datamax is None:
				datamax = data[ttime,:,:].max()

			if datastep is None:
				datastep = abs(datamax-datamin)/20.

	
			resources.cnLevelSelectionMode  = "ExplicitLevels" # Define own levels.
			levels = numpy.arange(datamin, datamax, datastep)
			levels = numpy.where((levels<0.00000001)&(levels>-0.00000001),0,levels)
			resources.cnLevels              = levels
			
			if levon==True:
				level_text = llevel[ttime]
				resources.tiMainString = "Lev "+str(level_text)
				
			
			map = Ngl.contour_map(wks,data[ttime,:,:],resources)
	
	if data.shape.__len__() == 4:
		for ttime in range(data.shape[0]):
			for llev in range(data.shape[1]):
				if datamin is None:
					datamin = data[ttime,llev,:,:].min()

				if datamax is None:
					datamax = data[ttime,llev,:,:].max()

				if datastep is None:
					datastep = abs(datamax-datamin)/20.

	
				resources.cnLevelSelectionMode  = "ExplicitLevels" # Define own levels.
				levels = numpy.arange(datamin, datamax, datastep)
				levels = numpy.where((levels<0.00000001)&(levels>-0.00000001),0,levels)
				resources.cnLevels              = levels
				
				if levon==True:
					level_text = llevel[llev]
					resources.tiMainString = "Lev "+str(level_text)
				
				map = Ngl.contour_map(wks,data[ttime,llev,:,:],resources)

	if showfig==True:
			os.system("gv "+psname+".ps")
		
	Ngl.delete_wks(wks)

	
def arctpltnc( data_file, variable_name, lon="lon", lat="lat", \
	region = "Arctic", datamin=None, datamax=None, datastep=None, \
	showfig=True, psname="output", colormap_name='testcmap', start_color=2, \
	end_color=-2, timeon=True, raster_fill=False, add_cyclic=False, mpFill=True,\
	levon=False,llevel="level",cnLevels=None, lbBoxFractions=None, sscale = 1):
	""" Plot variable from netCDF file in the Arctic region
	Usage:
		arctpl(data_file, variable_name, lon="lon", lat="lat", datamin=None, 
			   datamax=None, datastep=None, showfig=True, psname="output", 
			   colormap_name='testcmap', start_color=2, end_color=-2, 
			   timeon=True, raster_fill=False, add_cyclic=False, mpFill=True):
	
	Input:
		lon 		- 2D array of longitudes
		lat 		- 2D array of latitudes
		data 		- 2D, 3D or 4D array of scalar data. Two last dimensions should be lon/lat.
				  Time dimension must be relative ("time-units since reference-time"), You can convert from absolute to relative time by -r option in cdo.
)
		datamin 	- minimum value, if not specified will be calculated as data.min() for the first time step and level
		datamax 	- maximum value, if not specified will be calculated as data.max() for the first time step and level
		datastep 	- interval between isolines, if not specified will be calculated as abs(datamax-datamin)/20
		showfig 	- if True display the figure with gv
		psname 		- name of the output .ps file, default - "output"
		colormap_name	- name of the PyNGL of custom colormap. Custom colormaps should be defined in komod.colormap function. 
				  Color table for PyNGL can be found here http://www.pyngl.ucar.edu/Graphics/color_table_gallery.shtml
		
		start_color 	- number of the color in the PyNGL or custom colormap to start from, default = 2
		end_color  	- number of the color in the PyNGL or custom colormap to end with, default = -2
		vtitle 		- units of the data
		raster_fill	- switch between AreaFill and RasterFill.
		add_cyclic 	- do we need to add_cyclic or not, Default - False.
		mpFill    	- do we need to fill land or not, Default - True
		levon		- do we need to show level or not (timeon should be True), Default - False
		llevel		- name for the level variable, default = "level"
	"""
	ffile = Nio.open_file(data_file)
	data  = ffile.variables[variable_name][:]
	lonn = ffile.variables[lon][:]
	latt = ffile.variables[lat][:]
	
	if lonn.shape.__len__() == 1:
	  lonn,latt = coord2d(lonn,latt,data.shape[-2:])
	elif lonn.shape.__len__() == 2:
          lonn = lonn
          latt = latt
	elif lonn.shape.__len__() == 3:
          lonn = lonn[0,:,:]
          latt = latt[0,:,:] 
	
	
	if add_cyclic==True:
	  lonnn = Ngl.add_cyclic(lonn)
	  lattt = Ngl.add_cyclic(latt)
	else:
	  lonnn = lonn
	  lattt = latt
	
	

	rlist            = Ngl.Resources()
	rlist.wkColorMap = colormap(colormap_name)
	wks_type = "ps"
	wks = Ngl.open_wks(wks_type,psname,rlist)
	
	resources = Ngl.Resources()
	
	if hasattr(ffile.variables[variable_name],"_FillValue"):
		resources.sfMissingValueV = float(ffile.variables[variable_name]._FillValue[0])
		
	if hasattr(ffile.variables[variable_name],"missing_value"):
		resources.sfMissingValueV = float(ffile.variables[variable_name]._FillValue[0])
	
	resources.sfXArray        = lonnn[:]
	resources.sfYArray        = lattt[:]
	
	mapDict = reg(region)
	
	resources.mpProjection          = mapDict['mpProjection']
		
	resources.mpDataBaseVersion     = mapDict['mpDataBaseVersion']
	resources.mpGeophysicalLineThicknessF = mapDict['mpGeophysicalLineThicknessF']
	resources.mpLimitMode           = mapDict['mpLimitMode']
	resources.mpLeftCornerLatF      = mapDict['mpLeftCornerLatF']
	resources.mpLeftCornerLonF      = mapDict['mpLeftCornerLonF']
	resources.mpRightCornerLatF     = mapDict['mpRightCornerLatF']
	resources.mpRightCornerLonF     = mapDict['mpRightCornerLonF']
	resources.mpCenterLonF          = mapDict['mpCenterLonF']
	resources.mpCenterLatF          = mapDict['mpCenterLatF']
	resources.mpMinLonF             = mapDict['mpMinLonF']
	resources.mpMaxLonF             = mapDict['mpMaxLonF']
	resources.mpMinLatF             = mapDict['mpMinLatF']
	resources.mpMaxLatF             = mapDict['mpMaxLatF']
	resources.mpGridLonSpacingF     = mapDict['mpGridLonSpacingF']
	resources.mpGridLatSpacingF     = mapDict['mpGridLatSpacingF']

	resources.mpShapeMode    = 'FixedAspectFitBB'
	resources.cnFillDotSizeF    = 1

	resources.mpFillOn     = mpFill 
	igray = Ngl.new_color(wks,0.7,0.7,0.7)
	
	resources.mpFillColors = [0,-1,igray,-1]
	
	resources.cnLineDrawOrder       = "Predraw"
	
	resources.cnFillOn              = True
	resources.cnFillDrawOrder       = "Predraw"
	resources.cnLineLabelsOn        = False
	resources.nglSpreadColorStart   = start_color
	resources.nglSpreadColorEnd     = end_color
	
	
	
	if hasattr(ffile.variables[variable_name],"units"):
		resources.lbTitleString = str(ffile.variables[variable_name].units)
	else:
		resources.lbTitleString             = "values"
			
	
		
			
	
	resources.lbTitleFontHeightF        = 0.022
	resources.lbLabelFontHeightF        = 0.018
	resources.lbTitleOffsetF            = -0.30
	resources.lbBoxMinorExtentF         = 0.15
	resources.pmLabelBarOrthogonalPosF  = -0.0
	resources.lbOrientation             = "Horizontal"
	
	if raster_fill==True:
		resources.cnFillMode           = "RasterFill"
	
	resources.cnLinesOn             = False
	resources.cnMaxDataValueFormat  = ".4f"
	resources.mpGridAndLimbOn       = False
	#resources.mpLabelsOn           = True
	#resources.mpOutlineOn          = False
	resources.tmXTOn                = False
	resources.tmXBOn                = False
	resources.tmYLOn                = False
	resources.tmYROn                = False
	
	

	
	if data.shape.__len__() == 2:
		
		data_plot = data[:,:]
		
		if add_cyclic==True:
		  data_plot = Ngl.add_cyclic(data_plot)
		  
		if hasattr(ffile.variables[variable_name],"scale_factor"):
		  data_plot = data_plot*float(ffile.variables[variable_name].scale_factor[0])
		  
		if hasattr(ffile.variables[variable_name],"add_offset"):
		  data_plot = data_plot+float(ffile.variables[variable_name].add_offset[0])
		  
		if datamin is None:
			datamin = data.min()
		if datamax is None:
			datamax = data.max()
		if datastep is None:
			datastep = abs(datamax-datamin)/20.
	
		resources.cnLevelSelectionMode  = "ExplicitLevels" # Define own levels.
		if cnLevels is None:
			resources.cnLevels              = numpy.arange(datamin, datamax, datastep)
		else:
			resources.cnLevels = cnLevels
		if lbBoxFractions is not None:
			resources.lbBoxSizing               = 'ExplicitSizing'
			resources.lbBoxFractions            = lbBoxFractions
		  
		map = Ngl.contour_map(wks,data_plot*sscale,resources) 
		
	
	if data.shape.__len__() == 3:
		for ttime in range(data.shape[0]):
			
			if timeon==True:
				if 'Climatology' in ffile.variables["time"].description:
					ddate_str = 'Climatology'
				else:
					ddate = num2date(ffile.variables["time"][ttime], ffile.variables["time"].units)
					ddate_str = ddate.ctime()[4:]

				if levon==True:
					level_text = ffile.variables[llevel][ttime]
					resources.tiMainString = ddate_str+" Lev "+str(level_text)
				else:
					resources.tiMainString = ddate_str

			data_plot = data[ttime,:,:]
			
			if add_cyclic==True:
			  data_plot = Ngl.add_cyclic(data_plot)
			
			if hasattr(ffile.variables[variable_name],"scale_factor"):
			  data_plot = data_plot*float(ffile.variables[variable_name].scale_factor[0])
		  
			if hasattr(ffile.variables[variable_name],"add_offset"):
			  data_plot = data_plot+float(ffile.variables[variable_name].add_offset[0])
			
			if datamin is None:
			  datamin = data_plot.min()

			if datamax is None:
			  datamax = data_plot.max()

			if datastep is None:
			  datastep = abs(datamax-datamin)/20.

			resources.cnLevelSelectionMode  = "ExplicitLevels" # Define own levels.
			if cnLevels is None:
				resources.cnLevels              = numpy.arange(datamin, datamax, datastep)
			else:
				resources.cnLevels = cnLevels
			
			if lbBoxFractions is not None:
				resources.lbBoxSizing               = 'ExplicitSizing'
				resources.lbBoxFractions            = lbBoxFractions
		  
			map = Ngl.contour_map(wks,data_plot*sscale,resources) 
		

	
	if data.shape.__len__() == 4:
		for ttime in range(data.shape[0]):
			for llev in range(data.shape[1]):


				if timeon==True:
					if 'Climatology' in ffile.variables["time"].description:
						ddate_str = 'Climatology'
					else:
						ddate = num2date(ffile.variables["time"][ttime], ffile.variables["time"].units)
						ddate_str = ddate.ctime()[4:]

					if levon==True:
						level_text = ffile.variables[llevel][llev]
						resources.tiMainString = ddate_str+" Lev "+str(level_text)
					else:
						resources.tiMainString = ddate_str
					
						 
					
				data_plot = data[ttime,llev,:,:]
		
				if add_cyclic==True:
				  data_plot = Ngl.add_cyclic(data_plot)
				  
				if hasattr(ffile.variables[variable_name],"scale_factor"):
				  data_plot = data_plot*float(ffile.variables[variable_name].scale_factor[0])
		  
				if hasattr(ffile.variables[variable_name],"add_offset"):
				  data_plot = data_plot+float(ffile.variables[variable_name].add_offset[0])
				
				if datamin is None:
				  datamin = data_plot.min()

				if datamax is None:
				   datamax = data_plot.max()

				if datastep is None:
				  datastep = abs(datamax-datamin)/20
				
				resources.cnLevelSelectionMode  = "ExplicitLevels" # Define own levels.
				
				if cnLevels is None:
					resources.cnLevels              = numpy.arange(datamin, datamax, datastep)
				else:
					resources.cnLevels = cnLevels
				
				if lbBoxFractions is not None:
					resources.lbBoxSizing               = 'ExplicitSizing'
					resources.lbBoxFractions            = lbBoxFractions
		  
				map = Ngl.contour_map(wks,data_plot*sscale,resources) 
				

	if showfig==True:
			os.system("gv "+psname+".ps")
					
					
def globplt(lon, lat, data, datamin=None, datamax=None, \
	datastep=None, showfig=True, psname="output", colormap_name='testcmap',\
	 start_color=2, end_color=-2, vtitle="values", raster_fill=False, \
	 miss=None,levon=False,llevel=None,  minLon=0, maxLon=360 , minLat=-80 , maxLat=85, mpFill=True):
	""" Plot variable from array in the Global region
	
	Usage:
		arctpl(lon, lat, data, datamin=None, datamax=None, datastep=None, showfig=True, psname="output", colormap_name='posneg_1', start_color=2, end_color=-2, vtitle="values"):
	
	Input:
		lon 		- 2D array of longitudes
		lat 		- 2D array of latitudes
		data 		- 2D, 3D or 4D array of scalar data.Two last dimensions should be lon/lat
		datamin 	- minimum value, if not specified will be calculated as data.min() for the first time step and level
		datamax 	- maximum value, if not specified will be calculated as data.max() for the first time step and level
		datastep 	- interval between isolines, if not specified will be calculated as abs(datamax-datamin)/20
		showfig 	- if True display the figure with gv
		psname 		- name of the output .ps file, default - "output"
		colormap_name 	- name of the PyNGL of custom colormap. Custom colormaps should be defined in komod.colormap function. Color table for PyNGL can be found here http://www.pyngl.ucar.edu/Graphics/color_table_gallery.shtml
		start_color 	- number of the color in the PyNGL or custom colormap to start from, default = 2
		end_color  	- number of the color in the PyNGL or custom colormap to end with, default = -2
		vtitle 		- units of the data
		raster_fill	- switch between AreaFill and RasterFill.
		miss            - missing value
		levon           - do we need to plot level value in header?
		llevel          - data vector with values of levels (levon should be True)
		
	Output:
		.ps file, output.ps by default
	"""
	
	rlist            = Ngl.Resources()
	rlist.wkColorMap = colormap(colormap_name)
	wks_type = "ps"
	wks = Ngl.open_wks(wks_type,psname,rlist)
	
	resources = Ngl.Resources()
	
	resources.sfXArray        = lon[:]
	resources.sfYArray        = lat[:]
	
	
	resources.mpProjection          = "CylindricalEquidistant"
	resources.mpDataBaseVersion     = "LowRes"
	resources.mpLimitMode           = "LatLon"
	resources.mpMinLonF             = minLon
	resources.mpMaxLonF             = maxLon
	resources.mpMinLatF             = minLat
	resources.mpMaxLatF             = maxLat
	
	resources.mpShapeMode    = 'FixedAspectFitBB'
	resources.cnFillDotSizeF    = 1
	
	resources.mpCenterLonF           = 0.
	resources.mpCenterLatF           = 0.
	
	resources.mpFillOn     = True 
	igray = Ngl.new_color(wks,0.7,0.7,0.7)
	
	resources.mpFillColors = [0,-1,igray,-1]
	
	resources.cnLineDrawOrder       = "Predraw"
	
	resources.cnFillOn              = True
	resources.cnFillDrawOrder       = "Predraw"
	resources.cnLineLabelsOn        = False
	resources.nglSpreadColorStart   = start_color
	resources.nglSpreadColorEnd     = end_color
	
	
	resources.cnLevelSelectionMode  = "ExplicitLevels" # Define own levels.
	resources.cnLevels              = numpy.arange(datamin, datamax, datastep)
	
	resources.lbTitleString             = vtitle
	resources.lbTitleFontHeightF        = 0.022
	resources.lbLabelFontHeightF        = 0.018
	resources.lbTitleOffsetF            = -0.40
	resources.lbBoxMinorExtentF         = 0.15
	resources.pmLabelBarOrthogonalPosF  = -0.06
	resources.lbOrientation             = "Horizontal"
	
	
	if raster_fill==True:
		resources.cnFillMode           = "RasterFill"
		
	resources.cnLinesOn             = False
	resources.cnMaxDataValueFormat  = ".4f"
	resources.mpGridAndLimbOn       = False
	#resources.mpLabelsOn           = True
	#resources.mpOutlineOn          = False
	resources.tmXTOn                = False
	resources.tmXBOn                = False
	resources.tmYLOn                = False
	resources.tmYROn                = False
	
	if miss is not None:
		resources.sfMissingValueV = miss
	
	
	if data.shape.__len__() == 2:
		if datamin is None:
			datamin = data.min()
		if datamax is None:
			datamax = data.max()
		if datastep is None:
			datastep = abs(datamax-datamin)/20.
	
		resources.cnLevelSelectionMode  = "ExplicitLevels" # Define own levels.
		resources.cnLevels              = numpy.arange(datamin, datamax, datastep)
		map = Ngl.contour_map(wks,data[:,:],resources)
		
	
	if data.shape.__len__() == 3:
		for ttime in range(data.shape[0]):
			if datamin is None:
				datamin = data[ttime,:,:].min()

			if datamax is None:
				datamax = data[ttime,:,:].max()

			if datastep is None:
				datastep = abs(datamax-datamin)/20.

	
			resources.cnLevelSelectionMode  = "ExplicitLevels" # Define own levels.
			resources.cnLevels              = numpy.arange(datamin, datamax, datastep)
			
			if levon==True:
				level_text = llevel[ttime]
				resources.tiMainString = "Lev "+str(level_text)
				
			
			map = Ngl.contour_map(wks,data[ttime,:,:],resources)
	
	if data.shape.__len__() == 4:
		for ttime in range(data.shape[0]):
			for llev in range(data.shape[1]):
				if datamin is None:
					datamin = data[ttime,llev,:,:].min()

				if datamax is None:
					datamax = data[ttime,llev,:,:].max()

				if datastep is None:
					datastep = abs(datamax-datamin)/20.

	
				resources.cnLevelSelectionMode  = "ExplicitLevels" # Define own levels.
				resources.cnLevels              = numpy.arange(datamin, datamax, datastep)
				
				if levon==True:
					level_text = llevel[llev]
					resources.tiMainString = "Lev "+str(level_text)
				
				map = Ngl.contour_map(wks,data[ttime,llev,:,:],resources)

	if showfig==True:
			os.system("gv "+psname+".ps")
		
	Ngl.delete_wks(wks)




def globpltnc( data_file, variable_name, lon="lon", lat="lat", datamin=None, datamax=None, datastep=None, showfig=True, psname="output", colormap_name='testcmap', start_color=2, end_color=-2, timeon=True, raster_fill=False, add_cyclic=False, minLon=0, maxLon=360 , minLat=-80 , maxLat=85, mpFill=True):
	""" Plot variable from netCDF file for the Globe. You can also choose specific region by specifying max/min lat/lon
	Usage:
	arctpl(data_file, variable_name, lon="lon", lat="lat", datamin=None, datamax=None, datastep=None, showfig=True, psname="output", colormap_name='testcmap', start_color=2, end_color=-2, timeon=True, raster_fill=False, add_cyclic=False, minLon=0, maxLon=360 , minLat=-80 , maxLat=85, mpFill=True):
	
	Input:
		lon 		- 2D array of longitudes
		lat		- 2D array of latitudes
		data 		- 2D, 3D or 4D array of scalar data. Two last dimensions should be lon/lat.
		    Time dimension must be relative ("time-units since reference-time"), You can convert from absolute to relative time by -r option in cdo.
)
		datamin 	- minimum value, if not specified will be calculated as data.min() for the first time step and level
		datamax 	- maximum value, if not specified will be calculated as data.max() for the first time step and level
		datastep 	- interval between isolines, if not specified will be calculated as abs(datamax-datamin)/20
		showfig 	- if True display the figure with gv
		psname 		- name of the output .ps file, default - "output"
		colormap_name 	- name of the PyNGL of custom colormap. Custom colormaps should be defined in komod.colormap function. Color table for PyNGL can be found here http://www.pyngl.ucar.edu/Graphics/color_table_gallery.shtml
		start_color 	- number of the color in the PyNGL or custom colormap to start from, default = 2
		end_color   	- number of the color in the PyNGL or custom colormap to end with, default = -2
		vtitle 		- units of the data
		raster_fill	- switch between AreaFill and RasterFill.
		add_cyclic 	- do we need to add_cyclic or not, Default - False.
		minLon 		- minimum longitude
		maxLon 		- maximum longitude
		minLat 		- minimum latitude
		maxLat 		- maximum latitude
		mpFill     	- do we need to fill land or not, Default - True
	"""
	ffile = Nio.open_file(data_file)
	data  = ffile.variables[variable_name][:]
	lonn = ffile.variables[lon][:]
	latt = ffile.variables[lat][:]
	
	if lonn.shape.__len__() == 1:
	  lonnn,lattt = coord2d(lonn,latt,data.shape[-2:])
	
	if add_cyclic==True:
	  lonnn = Ngl.add_cyclic(lonn)
	  lattt = Ngl.add_cyclic(latt)
	else:
	  lonnn = lonn
	  lattt = latt
	  
	  
	

	rlist            = Ngl.Resources()
	rlist.wkColorMap = colormap(colormap_name)
	wks_type = "ps"
	wks = Ngl.open_wks(wks_type,psname,rlist)
	
	resources = Ngl.Resources()
	
	if hasattr(ffile.variables[variable_name],"_FillValue"):
		resources.sfMissingValueV = float(ffile.variables[variable_name]._FillValue[0])
		
	if hasattr(ffile.variables[variable_name],"missing_value"):
		resources.sfMissingValueV = float(ffile.variables[variable_name]._FillValue[0])
	
	
	resources.sfXArray        = lonnn[:]
	resources.sfYArray        = lattt[:]
	
	resources.mpProjection          = "CylindricalEquidistant"
	resources.mpDataBaseVersion     = "LowRes"
	resources.mpLimitMode           = "LatLon"
	resources.mpMinLonF             = minLon
	resources.mpMaxLonF             = maxLon
	resources.mpMinLatF             = minLat
	resources.mpMaxLatF             = maxLat
	
	resources.mpShapeMode    = 'FixedAspectFitBB'
	resources.cnFillDotSizeF    = 1
	
	resources.mpCenterLonF           = 0.
	resources.mpCenterLatF           = 0.
	
	resources.mpFillOn     = mpFill 
	igray = Ngl.new_color(wks,0.7,0.7,0.7)
	
	resources.mpFillColors = [0,-1,igray,-1]
	
	resources.cnLineDrawOrder       = "Predraw"
	
	resources.cnFillOn              = True
	resources.cnFillDrawOrder       = "Predraw"
	resources.cnLineLabelsOn        = False
	resources.nglSpreadColorStart   = start_color
	resources.nglSpreadColorEnd     = end_color
	
	
	
	if hasattr(ffile.variables[variable_name],"units"):
		resources.lbTitleString = str(ffile.variables[variable_name].units)
	else:
		resources.lbTitleString             = "values"

	resources.lbLabelFontHeightF        = 0.012
	resources.tiMainFontHeightF         = 0.012
	resources.lbTitleFontHeightF        = 0.012
	resources.lbTitleOffsetF            = -0.40
	resources.lbBoxMinorExtentF         = 0.15
	resources.pmLabelBarOrthogonalPosF  = -0.06
	resources.lbOrientation             = "Horizontal"
	
	if raster_fill==True:
		resources.cnFillMode           = "RasterFill"
	
	resources.cnLinesOn             = False
	resources.cnMaxDataValueFormat  = ".4f"
	resources.mpGridAndLimbOn       = False
	#resources.mpLabelsOn           = True
	#resources.mpOutlineOn          = False
	resources.tmXTOn                = False
	resources.tmXBOn                = False
	resources.tmYLOn                = False
	resources.tmYROn                = False
	
	
	if data.shape.__len__() == 2:
		
		data_plot = data[:,:]
		
		if add_cyclic==True:
		  data_plot = Ngl.add_cyclic(data_plot)
		  
		if hasattr(ffile.variables[variable_name],"scale_factor"):
		  data_plot = data_plot*float(ffile.variables[variable_name].scale_factor[0])
		  
		if hasattr(ffile.variables[variable_name],"add_offset"):
		  data_plot = data_plot+float(ffile.variables[variable_name].add_offset[0])
		  
		if datamin is None:
			datamin = data.min()
		if datamax is None:
			datamax = data.max()
		if datastep is None:
			datastep = abs(datamax-datamin)/20
	
		resources.cnLevelSelectionMode  = "ExplicitLevels" # Define own levels.
		resources.cnLevels              = numpy.arange(datamin, datamax, datastep)
		  
		map = Ngl.contour_map(wks,data_plot,resources) 
		
	
	if data.shape.__len__() == 3:
		for ttime in range(data.shape[0]):
			
			if timeon==True:
				ddate = num2date(ffile.variables["time"][ttime], ffile.variables["time"].units)
				resources.tiMainString = ddate.ctime()[4:]

			data_plot = data[ttime,:,:]
			
			if add_cyclic==True:
			  data_plot = Ngl.add_cyclic(data_plot)
			
			if hasattr(ffile.variables[variable_name],"scale_factor"):
			  data_plot = data_plot*float(ffile.variables[variable_name].scale_factor[0])
		  
			if hasattr(ffile.variables[variable_name],"add_offset"):
			  data_plot = data_plot+float(ffile.variables[variable_name].add_offset[0])
			
			if datamin is None:
			  datamin = data_plot.min()

			if datamax is None:
			  datamax = data_plot.max()

			if datastep is None:
			  datastep = abs(datamax-datamin)/20

			resources.cnLevelSelectionMode  = "ExplicitLevels" # Define own levels.
			resources.cnLevels              = numpy.arange(datamin, datamax, datastep)
			
			
		  
			map = Ngl.contour_map(wks,data_plot,resources) 
		

	
	if data.shape.__len__() == 4:
		for ttime in range(data.shape[0]):
			for llev in range(data.shape[1]):


				if timeon==True:
					ddate = num2date(ffile.variables["time"][ttime], ffile.variables["time"].units)
					resources.tiMainString = ddate.ctime()[4:]
					
				data_plot = data[ttime,llev,:,:]
		
				if add_cyclic==True:
				  data_plot = Ngl.add_cyclic(data_plot)
				  
				if hasattr(ffile.variables[variable_name],"scale_factor"):
				  data_plot = data_plot*float(ffile.variables[variable_name].scale_factor[0])
		  
				if hasattr(ffile.variables[variable_name],"add_offset"):
				  data_plot = data_plot+float(ffile.variables[variable_name].add_offset[0])
				
				if datamin is None:
				  datamin = data_plot.min()

				if datamax is None:
				   datamax = data_plot.max()

				if datastep is None:
				  datastep = abs(datamax-datamin)/20
				
				resources.cnLevelSelectionMode  = "ExplicitLevels" # Define own levels.
				resources.cnLevels              = numpy.arange(datamin, datamax, datastep)
		  
				map = Ngl.contour_map(wks,data_plot,resources) 
				

	if showfig==True:
			os.system("gv "+psname+".ps")
			
			
def pltgrd(lon,lat, region="Global",every=1, minLon=0, maxLon=360 , minLat=-80 , maxLat=85, coastThick=4, add_cyclic=False):
  """ Plot model grid from lat and lon arrays
  Usage:
	pltgrd(lon,lat, arctic=False, minLon=0, maxLon=360 , minLat=-80 , maxLat=85, coastThick=4)
  Imput:
	lon 		- longitude (ndarray)
	lat 		- latitude (ndarray)
	region		- one of predefined regions (for list of regions see the reg function)
	every 		- grid spacing (int)
	minLon 		- minimum longitude
	maxLon 		- maximum longitude
	minLat 		- minimum latitude
	maxLat 		- maximum latitude
	coastThick      - thickness of the coastlines
  Output: grid.ps"""
  if add_cyclic==True:
	    lon = Ngl.add_cyclic(lon)
	    lat = Ngl.add_cyclic(lat)
	    
  data = numpy.random.random(lon.shape)
  
  rlist            = Ngl.Resources()
  rlist.wkColorMap = 'posneg_1'
  wks_type = "ps"
  wks = Ngl.open_wks(wks_type,'grid',rlist)

  resources = Ngl.Resources()

  resources.sfXArray        = lon[::every,::every]
  resources.sfYArray        = lat[::every,::every]

  mapDict = reg(region, minLon, maxLon , minLat , maxLat)
  
  resources.mpProjection          = mapDict['mpProjection']
		
  resources.mpDataBaseVersion     = mapDict['mpDataBaseVersion']
  resources.mpGeophysicalLineThicknessF = mapDict['mpGeophysicalLineThicknessF']
  resources.mpLimitMode           = mapDict['mpLimitMode']
  resources.mpLeftCornerLatF      = mapDict['mpLeftCornerLatF']
  resources.mpLeftCornerLonF      = mapDict['mpLeftCornerLonF']
  resources.mpRightCornerLatF     = mapDict['mpRightCornerLatF']
  resources.mpRightCornerLonF     = mapDict['mpRightCornerLonF']
  resources.mpCenterLonF          = mapDict['mpCenterLonF']
  resources.mpCenterLatF          = mapDict['mpCenterLatF']
  resources.mpMinLonF             = mapDict['mpMinLonF']
  resources.mpMaxLonF             = mapDict['mpMaxLonF']
  resources.mpMinLatF             = mapDict['mpMinLatF']
  resources.mpMaxLatF             = mapDict['mpMaxLatF']
  resources.mpGridLonSpacingF     = mapDict['mpGridLonSpacingF']
  resources.mpGridLatSpacingF     = mapDict['mpGridLatSpacingF']
	
  resources.mpShapeMode    = 'FixedAspectFitBB'
  resources.cnFillDotSizeF    = 1
	
  resources.mpCenterLonF           = 0.
  resources.mpCenterLatF           = 0.

  resources.mpShapeMode    = 'FixedAspectFitBB'

  resources.cnFillDotSizeF    = 1

  resources.mpFillOn     = False 
  igray = Ngl.new_color(wks,0.7,0.7,0.7)
  resources.mpFillColors = [0,-1,igray,-1]

  resources.cnLineDrawOrder      = "Predraw"

  resources.cnFillOn             = True
  resources.cnFillDrawOrder       = "Predraw"
  resources.cnLineLabelsOn        = False
  
  resources.cnLevelSelectionMode = "ExplicitLevels" # Define own levels.
  resources.cnLevels             = numpy.arange(1,1.1,0.1)

  resources.lbLabelBarOn              = False
  
  resources.cnFillMode           = "CellFill"
  resources.cnLinesOn            = False
  resources.mpGridAndLimbOn      = False
  resources.tmXTOn                = False
  resources.tmXBOn                = False
  resources.tmYLOn                = False
  resources.tmYROn                = False
  resources.cnCellFillEdgeColor   = 1
  resources.cnCellFillMissingValEdgeColor   = 1
  
  map = Ngl.contour_map(wks,(data[::every,::every]),resources)


def pltgrdnc(data_file, variable_name, lon="lon", lat="lat", every=1,\
		     coastThick=4, region="Global", minLon=0, maxLon=360 , minLat=-80 , maxLat=85, add_cyclic=False):
  """ Plot model grid from netCDF file
  Usage:
	arcgrd(data_file, variable_name, lon="lon", lat="lat", arctic=False, minLon=0, maxLon=360 , minLat=-80 , maxLat=85, add_cyclic=False):
  Imput:
	data_file       - path to the data file
	variable_name   - name of the variable in the data file
	lon 		- longitude (name of the lon variable in netCDF file)
	lat 		- latitude (name of the lat variable in netCDF file)
	lat 		- latitude (ndarray)
	region		- one of predefined regions (for list of regions see the reg function)
	every 		- grid spacing (int)
	coastThick  - coast line thickness
	minLon 		- minimum longitude
	maxLon 		- maximum longitude
	minLat 		- minimum latitude
	maxLat 		- maximum latitude
  Output: grid.ps"""
  
  ffile = Nio.open_file(data_file)
  dataSh  = ffile.variables[variable_name][:]
  lonn  = ffile.variables[lon][:]
  latt  = ffile.variables[lat][:]
  

  
  if lonn.shape.__len__() == 1:
    lonnn,lattt = coord2d(lonn,latt,dataSh.shape[-2:])
  elif lonn.shape.__len__() == 2:
    lonnn = lonn
    lattt = latt
  elif lonn.shape.__len__() == 3:
    lonnn = lonn[0,:,:]
    lattt = latt[0,:,:]  


  
  if add_cyclic==True:
    lonnn = Ngl.add_cyclic(lonnn)
    lattt = Ngl.add_cyclic(lattt)
  else:
    lonnn = lonnn
    lattt = lattt

  data = numpy.random.random(lonnn.shape)

  
  rlist            = Ngl.Resources()
  rlist.wkColorMap = 'posneg_1'
  wks_type = "ps"
  wks = Ngl.open_wks(wks_type,'grid',rlist)

  resources = Ngl.Resources()

  resources.sfXArray        = lonnn[::every,::every]
  resources.sfYArray        = lattt[::every,::every]
  
  mapDict = reg(region, minLon, maxLon , minLat , maxLat)
  
  resources.mpProjection          = mapDict['mpProjection']
		
  resources.mpDataBaseVersion     = mapDict['mpDataBaseVersion']
  resources.mpGeophysicalLineThicknessF = mapDict['mpGeophysicalLineThicknessF']
  resources.mpLimitMode           = mapDict['mpLimitMode']
  resources.mpLeftCornerLatF      = mapDict['mpLeftCornerLatF']
  resources.mpLeftCornerLonF      = mapDict['mpLeftCornerLonF']
  resources.mpRightCornerLatF     = mapDict['mpRightCornerLatF']
  resources.mpRightCornerLonF     = mapDict['mpRightCornerLonF']
  resources.mpCenterLonF          = mapDict['mpCenterLonF']
  resources.mpCenterLatF          = mapDict['mpCenterLatF']
  resources.mpMinLonF             = mapDict['mpMinLonF']
  resources.mpMaxLonF             = mapDict['mpMaxLonF']
  resources.mpMinLatF             = mapDict['mpMinLatF']
  resources.mpMaxLatF             = mapDict['mpMaxLatF']
  resources.mpGridLonSpacingF     = mapDict['mpGridLonSpacingF']
  resources.mpGridLatSpacingF     = mapDict['mpGridLatSpacingF']
		

  resources.mpShapeMode    = 'FixedAspectFitBB'

  resources.cnFillDotSizeF    = 1

  resources.mpFillOn     = False 
  igray = Ngl.new_color(wks,0.7,0.7,0.7)
  resources.mpFillColors = [0,-1,igray,-1]

  resources.cnLineDrawOrder      = "Predraw"

  resources.cnFillOn             = True
  resources.cnFillDrawOrder       = "Predraw"
  resources.cnLineLabelsOn        = False
  
  resources.cnLevelSelectionMode = "ExplicitLevels" # Define own levels.
  resources.cnLevels             = numpy.arange(1,1.1,0.1)

  resources.lbLabelBarOn              = False
  
  resources.cnFillMode           = "CellFill"
  resources.cnLinesOn            = False
  resources.mpGridAndLimbOn      = False
  resources.tmXTOn                = False
  resources.tmXBOn                = False
  resources.tmYLOn                = False
  resources.tmYROn                = False
  resources.cnCellFillEdgeColor   = 1
  resources.cnCellFillMissingValEdgeColor   = 1
  
  map = Ngl.contour_map(wks,(data[::every,::every]),resources)

  
def pltgrd_line(lon, lat, lon1, lat1, lon2, lat2, npoints=10, every=1, region="Global", minLon=0, maxLon=360 , minLat=-80 , maxLat=85, psname="grid_line"):	
  """ Plot model grid from lat and lon arrays and line defined by coordinates of points defining the coordinates of the polyline
  Usage:
	pltgrd_line(lon, lat, lon1, lat1, lon2, lat2, npoints=10, every=1, arctic=False, minLon=0, maxLon=360 , minLat=-80 , maxLat=85)
  Imput:
	lon 		- longitude ( 2d ndarray)
	lat 		- latitude ( 2d ndarray)
	region		- one of predefined regions (for list of regions see the reg function)
	every 		- grid spacing (int)
	npoints 	- number of points, that will be used to drow the line.
	minLon 		- minimum longitude
	maxLon 		- maximum longitude
	minLat 		- minimum latitude
	maxLat 		- maximum latitude
  Output: grid_line.ps"""
  plat, plon = Ngl.gc_interp(lat1, lon1, lat2, lon2, npoints)
  data = numpy.random.random(lon.shape)
	
  rlist = Ngl.Resources()
  rlist.wkColorMap = 'posneg_1'
  wks_type = "ps"
  wks = Ngl.open_wks (wks_type,psname,rlist)
	
	
  res            = Ngl.Resources()   # map resources
  res.nglFrame   = False         # don't advance frame
  res.vpWidthF   = 0.80          # make map bigger
  res.vpHeightF  = 0.80
	
  res.sfXArray        = lon[::every,::every]
  res.sfYArray        = lat[::every,::every]
	
  mapDict = reg(region, minLon, maxLon , minLat , maxLat)
  
  res.mpProjection          = mapDict['mpProjection']
		
  res.mpDataBaseVersion     = mapDict['mpDataBaseVersion']
  res.mpGeophysicalLineThicknessF = mapDict['mpGeophysicalLineThicknessF']
  res.mpLimitMode           = mapDict['mpLimitMode']
  res.mpLeftCornerLatF      = mapDict['mpLeftCornerLatF']
  res.mpLeftCornerLonF      = mapDict['mpLeftCornerLonF']
  res.mpRightCornerLatF     = mapDict['mpRightCornerLatF']
  res.mpRightCornerLonF     = mapDict['mpRightCornerLonF']
  res.mpCenterLonF          = mapDict['mpCenterLonF']
  res.mpCenterLatF          = mapDict['mpCenterLatF']
  res.mpMinLonF             = mapDict['mpMinLonF']
  res.mpMaxLonF             = mapDict['mpMaxLonF']
  res.mpMinLatF             = mapDict['mpMinLatF']
  res.mpMaxLatF             = mapDict['mpMaxLatF']
  res.mpGridLonSpacingF     = mapDict['mpGridLonSpacingF']
  res.mpGridLatSpacingF     = mapDict['mpGridLatSpacingF']
		
  res.mpShapeMode    = 'FixedAspectFitBB'
	
  res.cnFillDotSizeF    = 1
	
  res.mpFillOn     = False 
  igray = Ngl.new_color(wks,0.7,0.7,0.7)
  res.mpFillColors = [0,-1,igray,-1]
	
  res.cnLineDrawOrder      = "Predraw"
	
  res.cnFillOn             = True
  res.cnFillDrawOrder       = "Predraw"
  res.cnLineLabelsOn        = False
	
  res.cnLevelSelectionMode = "ExplicitLevels" # Define own levels.
  res.cnLevels             = numpy.arange(1,1.1,0.1)
	
  res.lbLabelBarOn              = False
	
  res.cnFillMode           = "CellFill"
  res.cnLinesOn            = False
  res.mpGridAndLimbOn      = False
  res.tmXTOn                = False
  res.tmXBOn                = False
  res.tmYLOn                = False
  res.tmYROn                = False
  res.cnCellFillEdgeColor   = 1
  res.cnCellFillMissingValEdgeColor   = 1
	
  map = Ngl.contour_map(wks,(data[::every,::every]),res)
	
	
  pres                  = Ngl.Resources()        # polyline resources
  pres.gsLineThicknessF = 3.0                # line thickness
	
  Ngl.polyline(wks,map,plon,plat,pres)
	
	
  Ngl.frame(wks)
  
  os.system("gv "+psname+".ps" )
	
  return(plon, plat)
  #Ngl.end()

def pltgrd_line2(lon, lat, plon, plat, every=1, region="Global", minLon=0, maxLon=360 , minLat=-80 , maxLat=85):	
  """ Plot model grid from lat and lon arrays and line defined by points defining the coordinates of the polyline
  Usage:
	pltgrd_line(lon, lat, lon1, lat1, lon2, lat2, npoints=10, every=1, arctic=False, minLon=0, maxLon=360 , minLat=-80 , maxLat=85)
  Imput:
	lon 		- longitude (ndarray)
	lat 		- latitude (ndarray)
	plon, plat  - 1-dimensional arrays defining the coordinates of the polyline
	region		- one of predefined regions (for list of regions see the reg function)
	every 		- grid spacing (int)
	minLon 		- minimum longitude
	maxLon 		- maximum longitude
	minLat 		- minimum latitude
	maxLat 		- maximum latitude
  Output: grid_line.ps"""
  #plat, plon = Ngl.gc_interp(lat1, lon1, lat2, lon2, npoints)
  data = numpy.random.random(lon.shape)
	
  rlist = Ngl.Resources()
  rlist.wkColorMap = 'posneg_1'
  wks_type = "ps"
  wks = Ngl.open_wks (wks_type,"grid_line",rlist)
	
	
  res            = Ngl.Resources()   # map resources
  res.nglFrame   = False         # don't advance frame
  res.vpWidthF   = 0.80          # make map bigger
  res.vpHeightF  = 0.80
	
  res.sfXArray        = lon[::every,::every]
  res.sfYArray        = lat[::every,::every]
	
  res.mpProjection          = mapDict['mpProjection']
		
  res.mpDataBaseVersion     = mapDict['mpDataBaseVersion']
  res.mpGeophysicalLineThicknessF = mapDict['mpGeophysicalLineThicknessF']
  res.mpLimitMode           = mapDict['mpLimitMode']
  res.mpLeftCornerLatF      = mapDict['mpLeftCornerLatF']
  res.mpLeftCornerLonF      = mapDict['mpLeftCornerLonF']
  res.mpRightCornerLatF     = mapDict['mpRightCornerLatF']
  res.mpRightCornerLonF     = mapDict['mpRightCornerLonF']
  res.mpCenterLonF          = mapDict['mpCenterLonF']
  res.mpCenterLatF          = mapDict['mpCenterLatF']
  res.mpMinLonF             = mapDict['mpMinLonF']
  res.mpMaxLonF             = mapDict['mpMaxLonF']
  res.mpMinLatF             = mapDict['mpMinLatF']
  res.mpMaxLatF             = mapDict['mpMaxLatF']
  res.mpGridLonSpacingF     = mapDict['mpGridLonSpacingF']
  res.mpGridLatSpacingF     = mapDict['mpGridLatSpacingF']
		
  res.mpShapeMode    = 'FixedAspectFitBB'
	
  res.cnFillDotSizeF    = 1
	
  res.mpFillOn     = False 
  igray = Ngl.new_color(wks,0.7,0.7,0.7)
  res.mpFillColors = [0,-1,igray,-1]
	
  res.cnLineDrawOrder      = "Predraw"
	
  res.cnFillOn             = True
  res.cnFillDrawOrder       = "Predraw"
  res.cnLineLabelsOn        = False
	
  res.cnLevelSelectionMode = "ExplicitLevels" # Define own levels.
  res.cnLevels             = numpy.arange(1,1.1,0.1)
	
  res.lbLabelBarOn              = False
	
  res.cnFillMode           = "CellFill"
  res.cnLinesOn            = False
  res.mpGridAndLimbOn      = False
  res.tmXTOn                = False
  res.tmXBOn                = False
  res.tmYLOn                = False
  res.tmYROn                = False
  res.cnCellFillEdgeColor   = 1
  res.cnCellFillMissingValEdgeColor   = 1
	
  map = Ngl.contour_map(wks,(data[::every,::every]),res)
	
	
  pres                  = Ngl.Resources()        # polyline resources
  pres.gsLineThicknessF = 3.0                # line thickness
	
  Ngl.polyline(wks,map,plon,plat,pres)
	
	
  Ngl.frame(wks)
	
  #Ngl.end()

def plt_transect(distances, levels, vvalues, start_color = 5, end_color = -2, \
				 datamin = -2., datamax = 2., datastep = 0.2, psname="output",\
				 miss=None, vpWidthF=0.6, vpHeightF=0.6, colormap_name="posneg_1",\
				 cnLinesOn = True,cnLineLabelsOn = True, showfig=False):

	""" Plot transect. Input data for this function can be calculated by the get_transect function.

  	Input:
  		distances		- x distances, 1D (x_kilometers from get_transect function)
  		levels 			- depth of vertical levels, 1D
  		vvalues			- data to be plotted, 2d (data_prof from get_transect function)
  		start_color		- the index of the first color in the color table that should be used for a color contour or vector plot
  		end_color 		- the index of the last color in the color table that should be used for a color contour or vector plot
  		datamin 		- data minimum
  		datamax 		- data maximum
  		datastep 		- interval between contours
  		psname 			- name of the output file
	  	miss 			- missing values
	  	vpWidthF		- specifies the width of View object's bounding box in NDC units. 
  		vpHeightF		- specifies the height of View object's bounding box in NDC units
  		colormap_name	- name of the colormap
  		cnLinesOn 		- show contour lines
  		cnLineLabelsOn  - show labels of the contour lines
  		showfig 		- if True display the figure with gv

    	Output: output.ps"""
	wkres = Ngl.Resources()
	wkres.wkColorMap = colormap(colormap_name)
	wks_type = "ps"
	wks = Ngl.open_wks(wks_type, psname ,wkres)
	
	resources = Ngl.Resources()
	if miss != None:
		resources.sfMissingValueV = miss
	resources.sfYArray        = levels
	
	resources.nglYAxisType = "LinearAxis"
	resources.nglXAxisType = "LinearAxis"
	
	resources.cnFillOn        = True
	resources.trYReverse      = True
	
	resources.cnLinesOn            = cnLinesOn
	resources.cnLineLabelsOn        = False
	
	resources.nglSpreadColorStart   = start_color
	resources.nglSpreadColorEnd     = end_color
	
	resources.cnLevelSelectionMode  = "ExplicitLevels" # Define own levels.
	resources.cnLevels              = numpy.arange(datamin, datamax, datastep)
			
	resources.sfXArray        = distances
	resources.vpWidthF  = vpWidthF
	resources.vpHeightF = vpHeightF
	print(distances)
	plot = Ngl.contour(wks,vvalues,resources)
	
	if showfig==True:
		os.system("gv "+psname+".ps" )

def plt_vectors(lon, lat, u_wind, v_wind, MinFracLengthF = 0.005, \
				RefMagnitudeF = 0.3, RefLengthF = 0.05, LineArrowHeadMaxSizeF= 0.005, \
				LineArrowHeadMinSizeF=0.005, RefAnnoString1 = 'vector', \
				psname="output", sstep=1, region = 'Global', showfig=False,\
				 minLon=0, maxLon=360 , minLat=-80 , maxLat=85, vfMissingUValueV = None):
	'''Plot vectors over map
	'''
	rlist            = Ngl.Resources()
	rlist.wkColorMap = ["White","Black","Tan1","SkyBlue","Red"]
	wks_type = "ps"
	wks = Ngl.open_wks(wks_type,psname,rlist)
	
	resources = Ngl.Resources()
	resources.nglDraw  = False
	resources.nglFrame = False
	
	xxx_reshape   = u_wind[::sstep,::sstep]
	yyy_reshape   = v_wind[::sstep,::sstep]
	
	lon_reshape = lon[::sstep,::sstep]
	lat_reshape = lat[::sstep,::sstep]
	
	
	#resources.vfPolarData = True
	resources.vfUDataArray = xxx_reshape
	resources.vfVDataArray = yyy_reshape
	if vfMissingUValueV != None:
		vcres.vfMissingUValueV = vfMissingUValueV
	
	resources.vfXArray     = lon_reshape
	resources.vfYArray     = lat_reshape
	
	mapDict = reg(region, minLon, maxLon , minLat , maxLat)
	
	resources.mpProjection          = mapDict['mpProjection']
		
	resources.mpDataBaseVersion     = mapDict['mpDataBaseVersion']
	resources.mpGeophysicalLineThicknessF = mapDict['mpGeophysicalLineThicknessF']
	resources.mpLimitMode           = mapDict['mpLimitMode']
	resources.mpLeftCornerLatF      = mapDict['mpLeftCornerLatF']
	resources.mpLeftCornerLonF      = mapDict['mpLeftCornerLonF']
	resources.mpRightCornerLatF     = mapDict['mpRightCornerLatF']
	resources.mpRightCornerLonF     = mapDict['mpRightCornerLonF']
	resources.mpCenterLonF          = mapDict['mpCenterLonF']
	resources.mpCenterLatF          = mapDict['mpCenterLatF']
	resources.mpMinLonF             = mapDict['mpMinLonF']
	resources.mpMaxLonF             = mapDict['mpMaxLonF']
	resources.mpMinLatF             = mapDict['mpMinLatF']
	resources.mpMaxLatF             = mapDict['mpMaxLatF']
	resources.mpGridLonSpacingF     = mapDict['mpGridLonSpacingF']
	resources.mpGridLatSpacingF     = mapDict['mpGridLatSpacingF']
		
	
	
	resources.mpFillOn     = True 
	igray = Ngl.new_color(wks,0.7,0.7,0.7)
	resources.mpFillColors = [0,-1,igray,-1]
	resources.mpGridAndLimbOn      = False
	resources.tmXTOn                = False
	resources.tmXBOn                = False
	resources.tmYLOn                = False
	resources.tmYROn                = False
	
	resources.vcMinFracLengthF   = MinFracLengthF
	resources.vcRefMagnitudeF    = RefMagnitudeF
	resources.mpFillDrawOrder    = 'PostDraw'
	resources.vcRefLengthF       = RefLengthF 
	#resources.vcGlyphStyle       = "CurlyVector"
	#resources.vcFillArrowWidthF   = 5
	resources.vcLineArrowHeadMaxSizeF = LineArrowHeadMaxSizeF
	resources.vcLineArrowHeadMinSizeF = LineArrowHeadMinSizeF
	resources.vcPositionMode          = 'ArrowTail'
	resources.vcRefAnnoString1        = RefAnnoString1
	resources.vcRefAnnoFontHeightF    = 0.020
	resources.vcRefAnnoString2On      = False 
	#resources.vcLineArrowThicknessF   = 2
	resources.vcLineArrowThicknessF   = 3
	
	resources.vcRefAnnoOrthogonalPosF = -0.126      # Move ref anno up into plot.
	

	
	plot = Ngl.vector_map(wks,xxx_reshape,yyy_reshape,resources)
	
#	Ngl.maximize_plot(wks,plot)    # Maximize size of plot in frame.
	Ngl.draw(plot)
	Ngl.frame(wks)
	
	if showfig==True:
		os.system("gv "+psname+".ps")

def plt_vectors_scalars(lon, lat, u_wind, v_wind, scalar_data, datamin = -2., datamax = 2., datastep = 0.2, vfMissingUValueV = None, sfMissingValueV = None, MinFracLengthF = 0.005, RefMagnitudeF = 0.3, RefLengthF = 0.05, LineArrowHeadMaxSizeF= 0.005, LineArrowHeadMinSizeF=0.005, RefAnnoString1 = 'vector', psname="output", sstep=1, region = 'Global', showfig=False, TitleString = 'variable', tiMainString = None , minLon=0, maxLon=360 , minLat=-80 , maxLat=85, colormap_name='posneg_1', vcLineArrowThicknessF=1):
	'''Plot vectors and scalar value over map
  Usage:
	plt_vectors(lon, lat, u_wind, v_wind, MinFracLengthF = 0.005, RefMagnitudeF = 0.3, RefLengthF = 0.05, LineArrowHeadMaxSizeF= 0.005, LineArrowHeadMinSizeF=0.005, RefAnnoString1 = 'vector', psname="output", sstep=1, region = 'Global', minLon=0, maxLon=360 , minLat=-80 , maxLat=85)'''
	
	rlist            = Ngl.Resources()
	rlist.wkColorMap = colormap(colormap_name)
	wks_type = "ps"
	wks = Ngl.open_wks(wks_type, psname, rlist)
	
	mapres = Ngl.Resources()
	cnres  = Ngl.Resources()
	vcres  = Ngl.Resources()
	
	mapres.nglDraw  = False
	mapres.nglFrame = False
	cnres.nglDraw  = False
	cnres.nglFrame = False
	vcres.nglDraw  = False
	vcres.nglFrame = False
	
	xxx_reshape   = u_wind[::sstep,::sstep]
	yyy_reshape   = v_wind[::sstep,::sstep]
	
	lon_reshape = lon[::sstep,::sstep]
	lat_reshape = lat[::sstep,::sstep]
	
	#resources.vfPolarData = True
	vcres.vfUDataArray = xxx_reshape
	vcres.vfVDataArray = yyy_reshape
	if vfMissingUValueV != None:
		vcres.vfMissingUValueV = vfMissingUValueV
	if tiMainString != None:
		cnres.tiMainString = tiMainString
	
	vcres.vfXArray     = lon_reshape
	vcres.vfYArray     = lat_reshape
	
	cnres.sfXArray        = lon
	cnres.sfYArray        = lat
	if sfMissingValueV != None:
		cnres.sfMissingValueV = sfMissingValueV
	
	mapDict = reg(region, minLon, maxLon , minLat , maxLat)
	
	mapres.mpProjection          = mapDict['mpProjection']
		
	mapres.mpDataBaseVersion     = mapDict['mpDataBaseVersion']
	mapres.mpGeophysicalLineThicknessF = mapDict['mpGeophysicalLineThicknessF']
	mapres.mpLimitMode           = mapDict['mpLimitMode']
	mapres.mpLeftCornerLatF      = mapDict['mpLeftCornerLatF']
	mapres.mpLeftCornerLonF      = mapDict['mpLeftCornerLonF']
	mapres.mpRightCornerLatF     = mapDict['mpRightCornerLatF']
	mapres.mpRightCornerLonF     = mapDict['mpRightCornerLonF']
	mapres.mpCenterLonF          = mapDict['mpCenterLonF']
	mapres.mpCenterLatF          = mapDict['mpCenterLatF']
	mapres.mpMinLonF             = mapDict['mpMinLonF']
	mapres.mpMaxLonF             = mapDict['mpMaxLonF']
	mapres.mpMinLatF             = mapDict['mpMinLatF']
	mapres.mpMaxLatF             = mapDict['mpMaxLatF']
	mapres.mpGridLonSpacingF     = mapDict['mpGridLonSpacingF']
	mapres.mpGridLatSpacingF     = mapDict['mpGridLatSpacingF']
	
	mapres.mpFillOn     = True 
	igray = Ngl.new_color(wks,0.7,0.7,0.7)
	mapres.mpFillColors = [0,-1,igray,-1]
	mapres.mpGridAndLimbOn      = False
	mapres.tmXTOn                = False
	mapres.tmXBOn                = False
	mapres.tmYLOn                = False
	mapres.tmYROn                = False
	#mapres.mpGridSpacingF	     =2.
	#mapres.mpGridLatSpacingF     = 2
	#mapres.mpGridLonSpacingF     = 20
	
	vcres.vcMinFracLengthF   = MinFracLengthF
	vcres.vcRefMagnitudeF    = RefMagnitudeF
	mapres.mpFillDrawOrder    = 'PostDraw'
	vcres.vcRefLengthF       = RefLengthF 
	vcres.vcGlyphStyle       = "CurlyVector"
	#resources.vcFillArrowWidthF   = 5
	vcres.vcLineArrowHeadMaxSizeF = LineArrowHeadMaxSizeF
	vcres.vcLineArrowHeadMinSizeF = LineArrowHeadMinSizeF
	vcres.vcPositionMode          = 'ArrowTail'
	#vcres.vcRefAnnoString1        = '0.1 km~S~2~N~/M s~S~-1~N~'
	vcres.vcRefAnnoFontHeightF    = 0.020
	vcres.vcRefAnnoString2On      = False 
	vcres.vcLineArrowThicknessF   = vcLineArrowThicknessF
	
	
	vcres.vcRefAnnoOrthogonalPosF = -0.155      # Move ref anno up into plot.
	
	cnres.cnFillDotSizeF    = 1
	cnres.cnLineDrawOrder      = "Predraw"
	
	cnres.cnFillOn             = True
	cnres.cnFillDrawOrder       = "Predraw"
	cnres.cnLineLabelsOn        = False
	cnres.nglSpreadColorStart  = 5
	cnres.nglSpreadColorEnd      = -2
	
	cnres.cnLevelSelectionMode  = "ExplicitLevels" # Define own levels.
	cnres.cnLevels              = numpy.arange(datamin, datamax, datastep)
	
	cnres.lbTitleString             = TitleString
	cnres.lbTitleFontHeightF        = 0.022
	cnres.lbLabelFontHeightF        = 0.018
	cnres.lbTitleOffsetF            = -0.40
	cnres.lbBoxMinorExtentF         = 0.15
	cnres.pmLabelBarOrthogonalPosF  = -0.06
	cnres.lbOrientation             = "Horizontal"
	
	
	#resources.cnFillMode           = "RasterFill"
	cnres.cnLinesOn            = False
	cnres.cnMaxDataValueFormat = ".4f"
	
	map_plot          = Ngl.map(wks,mapres)
	std_plot          = Ngl.contour(wks, scalar_data, cnres)
	adxx_uv           = Ngl.vector(wks, xxx_reshape, yyy_reshape, vcres)
	
	Ngl.overlay(map_plot,adxx_uv)
	Ngl.overlay(map_plot,std_plot)
	
	Ngl.maximize_plot(wks,map_plot)    # Maximize size of plot in frame.
	Ngl.draw(map_plot)
	Ngl.frame(wks)
	
	if showfig==True:
		os.system("gv "+psname+".ps")

def plt_vectors_colors(lon, lat, u_wind, v_wind, scalar_data, 
	datamin = -2., datamax = 2., datastep = 0.2, vfMissingUValueV = None, 
	sfMissingValueV = None, MinFracLengthF = 0.005, RefMagnitudeF = 0.3, 
	RefLengthF = 0.05, LineArrowHeadMaxSizeF= 0.005, LineArrowHeadMinSizeF=0.005, 
	RefAnnoString1 = 'vector', psname="output", sstep=1, region = 'Global', 
	showfig=False, TitleString = 'variable', minLon=0, maxLon=360 , minLat=-80 , maxLat=85, 
	colormap_name='posneg_1',timeon=False,levon=False, level=None, timstep=None, tunits=None, 
	vcLineArrowThicknessF=1, overfigure=True):
	'''Plot vectors and scalar value over map
  Usage:
	plt_vectors(lon, lat, u_wind, v_wind, MinFracLengthF = 0.005, RefMagnitudeF = 0.3, RefLengthF = 0.05, LineArrowHeadMaxSizeF= 0.005, LineArrowHeadMinSizeF=0.005, RefAnnoString1 = 'vector', psname="output", sstep=1, region = 'Global', minLon=0, maxLon=360 , minLat=-80 , maxLat=85)'''
	
	rlist            = Ngl.Resources()
	rlist.wkColorMap = colormap(colormap_name)
	wks_type = "ps"
	wks = Ngl.open_wks(wks_type, psname, rlist)
	
	
	vcres  = Ngl.Resources()
	

	xxx_reshape   = u_wind[::sstep,::sstep]
	yyy_reshape   = v_wind[::sstep,::sstep]
	
	lon_reshape = lon[::sstep,::sstep]
	lat_reshape = lat[::sstep,::sstep]
	
	#resources.vfPolarData = True
	vcres.vfUDataArray = xxx_reshape
	vcres.vfVDataArray = yyy_reshape
	if vfMissingUValueV != None:
		vcres.vfMissingUValueV = vfMissingUValueV
	
	vcres.vfXArray     = lon_reshape
	vcres.vfYArray     = lat_reshape
	
	if timeon==True:
		ddate = num2date(timstep, tunits)
		if levon==True:
			level_text = level
			vcres.tiMainString = ddate.ctime()[4:]+" Lev "+str(numpy.round(level_text))
		else:
			vcres.tiMainString = ddate.ctime()[4:]
	
	
	mapDict = reg(region, minLon, maxLon , minLat , maxLat)
	
	vcres.mpProjection          = mapDict['mpProjection']
		
	vcres.mpDataBaseVersion     = mapDict['mpDataBaseVersion']
	vcres.mpGeophysicalLineThicknessF = mapDict['mpGeophysicalLineThicknessF']
	vcres.mpLimitMode           = mapDict['mpLimitMode']
	vcres.mpLeftCornerLatF      = mapDict['mpLeftCornerLatF']
	vcres.mpLeftCornerLonF      = mapDict['mpLeftCornerLonF']
	vcres.mpRightCornerLatF     = mapDict['mpRightCornerLatF']
	vcres.mpRightCornerLonF     = mapDict['mpRightCornerLonF']
	vcres.mpCenterLonF          = mapDict['mpCenterLonF']
	vcres.mpCenterLatF          = mapDict['mpCenterLatF']
	vcres.mpMinLonF             = mapDict['mpMinLonF']
	vcres.mpMaxLonF             = mapDict['mpMaxLonF']
	vcres.mpMinLatF             = mapDict['mpMinLatF']
	vcres.mpMaxLatF             = mapDict['mpMaxLatF']
	vcres.mpGridLonSpacingF     = mapDict['mpGridLonSpacingF']
	vcres.mpGridLatSpacingF     = mapDict['mpGridLatSpacingF']
		
	vcres.mpFillOn     = True 
	igray = Ngl.new_color(wks,0.7,0.7,0.7)
	vcres.mpFillColors = [0,-1,igray,-1]
	vcres.mpGridAndLimbOn      = False
	vcres.tmXTOn                = False
	vcres.tmXBOn                = False
	vcres.tmYLOn                = False
	vcres.tmYROn                = False
	
	vcres.vcMinFracLengthF   = MinFracLengthF
	vcres.vcRefMagnitudeF    = RefMagnitudeF
	vcres.mpFillDrawOrder    = 'PostDraw'
	vcres.vcRefLengthF       = RefLengthF 
	#vcres.vcGlyphStyle       = "CurlyVector"
	#resources.vcFillArrowWidthF   = 5
	vcres.vcLineArrowHeadMaxSizeF = LineArrowHeadMaxSizeF
	vcres.vcLineArrowHeadMinSizeF = LineArrowHeadMinSizeF
	vcres.vcPositionMode          = 'ArrowTail'
	#vcres.vcRefAnnoString1        = '0.1 km~S~2~N~/M s~S~-1~N~'
	vcres.vcRefAnnoFontHeightF    = 0.020
	vcres.vcRefAnnoString2On      = False 
	vcres.vcLineArrowThicknessF   = vcLineArrowThicknessF
	
	vcres.vcRefAnnoOn = False # Move ref anno up into plot.
	
	vcres.vcRefAnnoArrowLineColor = "Black"
	vcres.vcMaxLevelCount      = 5
	vcres.vcMonoLineArrowColor = False
	vcres.pmLabelBarDisplayMode         = "Always"
	vcres.lbTitleString          = TitleString
	
	vcres.vcLevelSelectionMode   = "ExplicitLevels"
	vcres.vcLevels             = numpy.arange(datamin, datamax, datastep)
	
	vcres.lbTitleFontHeightF        = 0.022
	vcres.lbLabelFontHeightF        = 0.018
	vcres.lbTitleOffsetF            = -0.40
	vcres.lbBoxMinorExtentF         = 0.15
	vcres.pmLabelBarOrthogonalPosF  = -0.06
	vcres.lbOrientation             = "Horizontal"
	vcres.nglSpreadColorStart  = 5
	vcres.nglSpreadColorEnd      = -5
	
	plot = Ngl.vector_map(wks, xxx_reshape, yyy_reshape, vcres)

		
	
	if showfig==True:
		os.system("gv "+psname+".ps")
	
	Ngl.delete_wks(wks)

def plot(data, data2=None, xvalues = None, xvalues2=None):
	"""simpe plot of one or two variables"""
	if xvalues is None:
		xvalues = numpy.arange(data.shape[0])
	
		
	if data2 is not None:		
		if xvalues2 is None:
			xvalues2 = numpy.arange(data2.shape[0])
		plt.plot(xvalues,data,xvalues2,data2)
	else:
		plt.plot(xvalues,data)	
	plt.show()

def plt_adxx_stat(varname,adxxdata,mmask):
	
    '''
    plot min max and mean of adxx fields.
    Input:
        varname  = name of the variable (e.g. "atemp"), string 
        adxxdata = 4D array with data
        mmask    = land sea mask 

    '''
    
    plt.subplot(2,3,1)
    plt.imshow(numpy.flipud(adxxdata.max(axis=0)[0,:,:]*mmask))
    plt.colorbar(orientation='horizontal')
    plt.title('adxx_'+varname+' max')
        
        
    plt.subplot(2,3,2)
    plt.imshow(numpy.flipud(adxxdata.min(axis=0)[0,:,:]*mmask))
    plt.colorbar(orientation='horizontal')
    plt.title('adxx_'+varname+' min')
        
        
    plt.subplot(2,3,3)
    plt.imshow(numpy.flipud(adxxdata.mean(axis=0)[0,:,:]*mmask))
    plt.colorbar(orientation='horizontal')
    plt.title('adxx_'+varname+' mean')
       
        
    plt.subplot(2,3,4)
    plt.plot(adxxdata.max(axis=2).max(axis=2))
        
    plt.subplot(2,3,5)
    plt.plot(adxxdata.min(axis=2).min(axis=2))
        
    plt.subplot(2,3,6)
    plt.plot(adxxdata.mean(axis=2).mean(axis=2))
        
    plt.tight_layout()
    plt.show