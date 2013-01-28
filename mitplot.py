#!/bin/python
# -*- coding: utf-8 -*-


"""Komod plot module 
coord2d 	- Convert 1d coordinates to 2d coordinates
colormap 	- Define custom colormaps
arctpl 		- plot contours of 2D, 3D or 4D field in the Arctic region from ndarray
arctpltnc 	- plot contours of 2D, 3D or 4D field in the Arctic region from netCDF file
globpltnc 	- plot contours of 2D, 3D or 4D field on the global map from netCDF file. You can also specify custom region.
pltgrd		- plot model grid from lat and lon arrays
pltgrdnc	- plot model grid from netCDF file

Nikolay Koldunov 23 February 2010
"""

# -------------------------------------------------
import Ngl
import Nio
import numpy
import os
from netcdftime import num2date
import matplotlib.pyplot as plt


def coord2d(lon, lat, dshape):
	""" Convert 1d coordinates to 2d coordinates.
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
	""" define custom colormaps
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


				
def arctpl(lon, lat, data, datamin=None, datamax=None, datastep=None, showfig=True, psname="output", colormap_name='testcmap', start_color=2, end_color=-2, vtitle="values", raster_fill=False, miss=None,levon=False,llevel=None, add_cyclic=False, region='Arctic', minLon=0, maxLon=360 , minLat=-80 , maxLat=85):
	""" Plot variable from array in the Arctic region
	
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
		region 		- one of predefined regions (Arctic(def), Global, NAtlantic)
		
	Output:
		.ps file, output.ps by default
	"""
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
	
	
	if region == 'Arctic':
		resources.mpProjection          = "LambertEqualArea"
		
		resources.mpDataBaseVersion     = "LowRes"
		resources.mpGeophysicalLineThicknessF = 4
		resources.mpLimitMode           = 'Corners'
		resources.mpLeftCornerLatF      = 54
		resources.mpLeftCornerLonF      = -50
		resources.mpRightCornerLatF     = 57
		resources.mpRightCornerLonF     = 140
		resources.mpCenterLonF           = 0.
		resources.mpCenterLatF           = 90.
	elif region == 'Global':
		resources.mpProjection          = "CylindricalEquidistant"
		resources.mpDataBaseVersion     = "MediumRes"
		resources.mpLimitMode           = "LatLon"
		resources.mpMinLonF             = minLon
		resources.mpMaxLonF             = maxLon
		resources.mpMinLatF             = minLat
		resources.mpMaxLatF             = maxLat
		resources.mpCenterLonF           = 180.
		resources.mpCenterLatF           = 0.
	elif region == 'NAtlantic':
		resources.mpProjection          = "LambertEqualArea"
		resources.mpDataBaseVersion     = "MediumRes"
		resources.mpLimitMode           = 'Corners'
		resources.mpLeftCornerLatF      = 65
		resources.mpLeftCornerLonF      = -20
		resources.mpRightCornerLatF     = 70
		resources.mpRightCornerLonF     = 100
		resources.mpCenterLonF           = 0.
		resources.mpCenterLatF           = 90.
	
	elif region == 'stanna':
		resources.mpProjection          = "LambertEqualArea"
		resources.mpDataBaseVersion     = "MediumRes"
		resources.mpLimitMode           = 'Corners'
		resources.mpLeftCornerLatF      = 80
		resources.mpLeftCornerLonF      = 25
		resources.mpRightCornerLatF     = 70
		resources.mpRightCornerLonF     = 100
		resources.mpCenterLonF           = 0.
		resources.mpCenterLatF           = 90.	
		resources.mpGridLonSpacingF    = 15
		resources.mpGridLatSpacingF    = 5
		
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
	resources.cnLevels              = numpy.arange(datamin, datamax, datastep)
	
	resources.lbTitleString             = vtitle
	resources.lbTitleFontHeightF        = 0.022
	resources.lbLabelFontHeightF        = 0.018
	resources.lbTitleOffsetF            = -0.40
	resources.lbBoxMinorExtentF         = 0.15
	resources.pmLabelBarOrthogonalPosF  = -0.03
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

	
def arctpltnc( data_file, variable_name, lon="lon", lat="lat", region = "Arctic", datamin=None, datamax=None, datastep=None, showfig=True, psname="output", colormap_name='testcmap', start_color=2, end_color=-2, timeon=True, raster_fill=False, add_cyclic=False, mpFill=True,levon=False,llevel="level",cnLevels=None, lbBoxFractions=None, sscale = 1):
	""" Plot variable from netCDF file in the Arctic region
	Usage:
		arctpl(data_file, variable_name, lon="lon", lat="lat", datamin=None, datamax=None, datastep=None, showfig=True, psname="output", colormap_name='testcmap', start_color=2, end_color=-2, timeon=True, raster_fill=False, add_cyclic=False, mpFill=True):
	
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
	if region == "Arctic":
		resources.mpProjection          = "LambertEqualArea"
		resources.mpDataBaseVersion     = "LowRes"
		resources.mpLimitMode           = 'Corners'
		resources.mpLeftCornerLatF      = 54
		resources.mpLeftCornerLonF      = -50
		resources.mpRightCornerLatF     = 57
		resources.mpRightCornerLonF     = 140
		resources.mpCenterLonF           = 0.
		resources.mpCenterLatF           = 90.

	if region == "Monarch":  
        	resources.mpProjection          = "LambertEqualArea"

 	   	resources.mpDataBaseVersion     = "LowRes"
   		#resources.mpGeophysicalLineThicknessF = coastThick
   		resources.mpLimitMode           = 'Corners'
   		resources.mpLeftCornerLatF      = 33
    		resources.mpLeftCornerLonF      = -70
    		resources.mpRightCornerLatF     = 40
    		resources.mpRightCornerLonF     = 110
    		resources.mpCenterLonF           = -34.
    		resources.mpCenterLatF           = 90.

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
				ddate = num2date(ffile.variables["time"][ttime], ffile.variables["time"].units)
				if levon==True:
					level_text = ffile.variables[llevel][ttime]
					resources.tiMainString = ddate.ctime()[4:]+" Lev "+str(level_text)
				else:
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
					ddate = num2date(ffile.variables["time"][ttime], ffile.variables["time"].units)
					if levon==True:
						level_text = ffile.variables[llevel][llev]
						resources.tiMainString = ddate.ctime()[4:]+" Lev "+str(level_text)
					else:
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
					
					
def globplt(lon, lat, data, datamin=None, datamax=None, datastep=None, showfig=True, psname="output", colormap_name='testcmap', start_color=2, end_color=-2, vtitle="values", raster_fill=False, miss=None,levon=False,llevel=None,  minLon=0, maxLon=360 , minLat=-80 , maxLat=85, mpFill=True):
	""" Plot variable from array in the Arctic region
	
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
			
			
def pltgrd(lon,lat, arctic=False, minLon=0, maxLon=360 , minLat=-80 , maxLat=85, coastThick=4, add_cyclic=False):
  """ Plot model grid from lat and lon arrays
  Usage:
	pltgrd(lon,lat, arctic=False, minLon=0, maxLon=360 , minLat=-80 , maxLat=85, coastThick=4)
  Imput:
	lon 		- longitude (ndarray)
	lat 		- latitude (ndarray)
	arctic		- if True we will plot predefined Arctic Region
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

  resources.sfXArray        = lon[:,:]
  resources.sfYArray        = lat[:,:]
  
  if arctic == True:
    resources.mpProjection          = "LambertEqualArea"

    resources.mpDataBaseVersion     = "LowRes"
    resources.mpGeophysicalLineThicknessF = coastThick
    resources.mpLimitMode           = 'Corners'
    resources.mpLeftCornerLatF      = 45
    resources.mpLeftCornerLonF      = -41
    resources.mpRightCornerLatF     = 48
    resources.mpRightCornerLonF     = 138
    resources.mpCenterLonF           = 0.
    resources.mpCenterLatF           = 90.
  else:
    resources.mpProjection          = "CylindricalEquidistant"
    resources.mpDataBaseVersion     = "LowRes"
    resources.mpLimitMode           = "LatLon"
    resources.mpMinLonF             = minLon
    resources.mpMaxLonF             = maxLon
    resources.mpMinLatF             = minLat
    resources.mpMaxLatF             = maxLat
    resources.mpGeophysicalLineThicknessF = coastThick
	
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
  
  map = Ngl.contour_map(wks,(data[:,:]),resources)


def pltgrdnc(data_file, variable_name, lon="lon", lat="lat", every=1, coastThick=4, region="Global", minLon=0, maxLon=360 , minLat=-80 , maxLat=85, add_cyclic=False):
  """ Plot model grid from netCDF file
  Usage:
	arcgrd(data_file, variable_name, lon="lon", lat="lat", arctic=False, minLon=0, maxLon=360 , minLat=-80 , maxLat=85, add_cyclic=False):
  Imput:
	data_file       - path to the data file
	variable_name   - name of the variable in the data file
	lon 		- longitude (name of the lon variable in netCDF file)
	lat 		- latitude (name of the lat variable in netCDF file)
	arctic		- if True we will plot predefined Arctic Region
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
  
  if region == 'Arctic':
	resources.mpProjection          = "LambertEqualArea"
		
	resources.mpDataBaseVersion     = "LowRes"
	resources.mpGeophysicalLineThicknessF = coastThick
	resources.mpLimitMode           = 'Corners'
	resources.mpLeftCornerLatF      = 54
	resources.mpLeftCornerLonF      = -50
	resources.mpRightCornerLatF     = 57
	resources.mpRightCornerLonF     = 140
	resources.mpCenterLonF           = 0.
	resources.mpCenterLatF           = 90.
  if region == "Monarch":  
        resources.mpProjection          = "LambertEqualArea"

    	resources.mpDataBaseVersion     = "LowRes"
   	resources.mpGeophysicalLineThicknessF = coastThick
   	resources.mpLimitMode           = 'Corners'
   	resources.mpLeftCornerLatF      = 20
    	resources.mpLeftCornerLonF      = -41
    	resources.mpRightCornerLatF     = 30
    	resources.mpRightCornerLonF     = 138
    	resources.mpCenterLonF           = 0.
    	resources.mpCenterLatF           = 90.

  elif region == 'Global':
	resources.mpProjection          = "LambertEqualArea"
	resources.mpDataBaseVersion     = "LowRes"
	resources.mpGeophysicalLineThicknessF = coastThick
	resources.mpLimitMode           = "LatLon"
	resources.mpMinLonF             = minLon
	resources.mpMaxLonF             = maxLon
	resources.mpMinLatF             = minLat
	resources.mpMaxLatF             = maxLat
	resources.mpCenterLonF           = 0.
	resources.mpCenterLatF           = 0.
  elif region == 'NAtlantic':
	resources.mpProjection          = "LambertEqualArea"
	resources.mpDataBaseVersion     = "MediumRes"
	resources.mpGeophysicalLineThicknessF = coastThick
	resources.mpLimitMode           = 'Corners'
	resources.mpLeftCornerLatF      = 65
	resources.mpLeftCornerLonF      = -20
	resources.mpRightCornerLatF     = 70
	resources.mpRightCornerLonF     = 100
	resources.mpCenterLonF           = 0.
	resources.mpCenterLatF           = 90.
  elif region == 'halo':
	resources.mpProjection          = "LambertEqualArea"
	resources.mpDataBaseVersion     = "MediumRes"
	resources.mpGeophysicalLineThicknessF = coastThick
	resources.mpLimitMode           = 'Corners'
	resources.mpLeftCornerLatF      = 85
	resources.mpLeftCornerLonF      = 0
	resources.mpRightCornerLatF     = 65
	resources.mpRightCornerLonF     = 120
	resources.mpCenterLonF           = 0.
	resources.mpCenterLatF           = 90.

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

  
def pltgrd_line(lon, lat, lon1, lat1, lon2, lat2, npoints=10, every=1, arctic=False, minLon=0, maxLon=360 , minLat=-80 , maxLat=85, psname="grid_line"):	
  """ Plot model grid from lat and lon arrays and line defined by coordinates of the first and the last point
  Usage:
	pltgrd_line(lon, lat, lon1, lat1, lon2, lat2, npoints=10, every=1, arctic=False, minLon=0, maxLon=360 , minLat=-80 , maxLat=85)
  Imput:
	lon 		- longitude (ndarray)
	lat 		- latitude (ndarray)
	arctic		- if True we will plot predefined Arctic Region
	minLon 		- minimum longitude
	maxLon 		- maximum longitude
	minLat 		- minimum latitude
	maxLat 		- maximum latitude
  Output: grid.ps"""
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
	
  if arctic == True:
    res.mpProjection          = "LambertEqualArea"

    res.mpDataBaseVersion     = "LowRes"
    res.mpGeophysicalLineThicknessF = 4
    res.mpLimitMode           = 'Corners'
    res.mpLeftCornerLatF      = 45
    res.mpLeftCornerLonF      = -41
    res.mpRightCornerLatF     = 48
    res.mpRightCornerLonF     = 138
    res.mpCenterLonF           = 0.
    res.mpCenterLatF           = 90.
  else:
    res.mpProjection          = "CylindricalEquidistant"
    res.mpDataBaseVersion     = "LowRes"
    res.mpLimitMode           = "LatLon"
    res.mpMinLonF             = minLon
    res.mpMaxLonF             = maxLon
    res.mpMinLatF             = minLat
    res.mpMaxLatF             = maxLat
	    	
    res.mpCenterLonF           = 0.
    res.mpCenterLatF           = 0.
		
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

def pltgrd_line2(lon, lat, plon, plat, every=1, arctic=False, minLon=0, maxLon=360 , minLat=-80 , maxLat=85):	
  """ Plot model grid from lat and lon arrays and line defined by coordinates of the first and the last point
  Usage:
	pltgrd_line(lon, lat, lon1, lat1, lon2, lat2, npoints=10, every=1, arctic=False, minLon=0, maxLon=360 , minLat=-80 , maxLat=85)
  Imput:
	lon 		- longitude (ndarray)
	lat 		- latitude (ndarray)
	arctic		- if True we will plot predefined Arctic Region
	minLon 		- minimum longitude
	maxLon 		- maximum longitude
	minLat 		- minimum latitude
	maxLat 		- maximum latitude
  Output: grid.ps"""
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
	
  if arctic == True:
    res.mpProjection          = "LambertEqualArea"

    res.mpDataBaseVersion     = "LowRes"
    res.mpGeophysicalLineThicknessF = 4
    res.mpLimitMode           = 'Corners'
    res.mpLeftCornerLatF      = 45
    res.mpLeftCornerLonF      = -41
    res.mpRightCornerLatF     = 48
    res.mpRightCornerLonF     = 138
    res.mpCenterLonF           = 0.
    res.mpCenterLatF           = 90.
  else:
    res.mpProjection          = "CylindricalEquidistant"
    res.mpDataBaseVersion     = "LowRes"
    res.mpLimitMode           = "LatLon"
    res.mpMinLonF             = minLon
    res.mpMaxLonF             = maxLon
    res.mpMinLatF             = minLat
    res.mpMaxLatF             = maxLat
	    	
    res.mpCenterLonF           = 0.
    res.mpCenterLatF           = 0.
		
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

def plt_transect(distances, levels, vvalues, start_color = 5, end_color = -2, datamin = -2., datamax = 2., datastep = 0.2, psname="output", miss=None, vpWidthF=0.6, vpHeightF=0.6, colormap_name="posneg_1",cnLinesOn = True,cnLineLabelsOn = True, showfig=False):

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

def plt_vectors(lon, lat, u_wind, v_wind, MinFracLengthF = 0.005, RefMagnitudeF = 0.3, RefLengthF = 0.05, LineArrowHeadMaxSizeF= 0.005, LineArrowHeadMinSizeF=0.005, RefAnnoString1 = 'vector', psname="output", sstep=1, region = 'Global', showfig=False, minLon=0, maxLon=360 , minLat=-80 , maxLat=85, vfMissingUValueV = None):
	'''Plot vectors over map
  Usage:
	plt_vectors(lon, lat, u_wind, v_wind, MinFracLengthF = 0.005, RefMagnitudeF = 0.3, RefLengthF = 0.05, LineArrowHeadMaxSizeF= 0.005, LineArrowHeadMinSizeF=0.005, RefAnnoString1 = 'vector', psname="output", sstep=1, region = 'Global', minLon=0, maxLon=360 , minLat=-80 , maxLat=85)
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
	
	
	
	if region == 'Arctic':
		resources.mpProjection          = "LambertEqualArea"
		
		resources.mpDataBaseVersion     = "LowRes"
		resources.mpGeophysicalLineThicknessF = 4
		resources.mpLimitMode           = 'Corners'
		resources.mpLeftCornerLatF      = 54
		resources.mpLeftCornerLonF      = -50
		resources.mpRightCornerLatF     = 57
		resources.mpRightCornerLonF     = 140
		resources.mpCenterLonF           = 0.
		resources.mpCenterLatF           = 90.
	elif region == 'Global':
		resources.mpProjection          = "LambertEqualArea"
		resources.mpDataBaseVersion     = "LowRes"
		resources.mpLimitMode           = "LatLon"
		resources.mpMinLonF             = minLon
		resources.mpMaxLonF             = maxLon
		resources.mpMinLatF             = minLat
		resources.mpMaxLatF             = maxLat
		resources.mpCenterLonF           = 0.
		resources.mpCenterLatF           = 0.
	elif region == 'NAtlantic':
		resources.mpProjection          = "LambertEqualArea"
		resources.mpDataBaseVersion     = "MediumRes"
		resources.mpLimitMode           = 'Corners'
		resources.mpLeftCornerLatF      = 65
		resources.mpLeftCornerLonF      = -20
		resources.mpRightCornerLatF     = 70
		resources.mpRightCornerLonF     = 100
		resources.mpCenterLonF           = 0.
		resources.mpCenterLatF           = 90.
	elif region == 'halo':
		resources.mpProjection          = "LambertEqualArea"
		resources.mpDataBaseVersion     = "MediumRes"
		resources.mpLimitMode           = 'Corners'
		resources.mpLeftCornerLatF      = 85
		resources.mpLeftCornerLonF      = 0
		resources.mpRightCornerLatF     = 65
		resources.mpRightCornerLonF     = 120
		resources.mpCenterLonF           = 0.
		resources.mpCenterLatF           = 90.
	
	
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
	
	Ngl.maximize_plot(wks,plot)    # Maximize size of plot in frame.
	Ngl.draw(plot)
	Ngl.frame(wks)
	
	if showfig==True:
		os.system("gv "+psname+".ps")

def plt_vectors_scalars(lon, lat, u_wind, v_wind, scalar_data, datamin = -2., datamax = 2., datastep = 0.2, vfMissingUValueV = None, sfMissingValueV = None, MinFracLengthF = 0.005, RefMagnitudeF = 0.3, RefLengthF = 0.05, LineArrowHeadMaxSizeF= 0.005, LineArrowHeadMinSizeF=0.005, RefAnnoString1 = 'vector', psname="output", sstep=1, region = 'Global', showfig=False, TitleString = 'variable', minLon=0, maxLon=360 , minLat=-80 , maxLat=85, colormap_name='posneg_1', vcLineArrowThicknessF=1):
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
	
	vcres.vfXArray     = lon_reshape
	vcres.vfYArray     = lat_reshape
	
	cnres.sfXArray        = lon
	cnres.sfYArray        = lat
	if sfMissingValueV != None:
		cnres.sfMissingValueV = sfMissingValueV
	
	
	if region == 'Arctic':
		mapres.mpProjection          = "LambertEqualArea"
		mapres.mpDataBaseVersion     = "LowRes"
		mapres.mpLimitMode           = 'Corners'
		mapres.mpLeftCornerLatF      = 54
		mapres.mpLeftCornerLonF      = -50
		mapres.mpRightCornerLatF     = 57
		mapres.mpRightCornerLonF     = 140
		mapres.mpCenterLonF           = 0.
		mapres.mpCenterLatF           = 90.
			
	elif region == 'Global':
		mapres.mpProjection          = "LambertEqualArea"
		mapres.mpDataBaseVersion     = "LowRes"
		mapres.mpLimitMode           = "LatLon"
		mapres.mpMinLonF             = minLon
		mapres.mpMaxLonF             = maxLon
		mapres.mpMinLatF             = minLat
		mapres.mpMaxLatF             = maxLat
		mapres.mpCenterLonF           = 0.
		mapres.mpCenterLatF           = 0.
	elif region == 'NAtlantic':
		mapres.mpProjection          = "LambertEqualArea"
		mapres.mpDataBaseVersion     = "MediumRes"
		mapres.mpLimitMode           = 'Corners'
		mapres.mpLeftCornerLatF      = 65
		mapres.mpLeftCornerLonF      = -20
		mapres.mpRightCornerLatF     = 70
		mapres.mpRightCornerLonF     = 100
		mapres.mpCenterLonF           = 0.
		mapres.mpCenterLatF           = 90.
	elif region == 'FramStAnna':
		mapres.mpProjection          = "LambertEqualArea"
		mapres.mpDataBaseVersion     = "MediumRes"
		mapres.mpLimitMode           = 'Corners'
		mapres.mpLeftCornerLatF      = 75
		mapres.mpLeftCornerLonF      = -5
		mapres.mpRightCornerLatF     = 75
		mapres.mpRightCornerLonF     = 120
		mapres.mpCenterLonF           = 0.
		mapres.mpCenterLatF           = 90.
	if region == 'Fram2':
    		mapres.mpProjection          = "LambertEqualArea"
    		mapres.mpDataBaseVersion     = "MediumRes"
    		mapres.mpLimitMode           = 'Corners'
    		mapres.mpLeftCornerLatF      = 75
    		mapres.mpLeftCornerLonF      = 0
    		mapres.mpRightCornerLatF     = 77
    		mapres.mpRightCornerLonF     = 90
    		mapres.mpCenterLonF           = 0.
    		mapres.mpCenterLatF           = 90.
    		mapres.mpGridLonSpacingF    = 15
    		mapres.mpGridLatSpacingF    = 5

	elif region == 'halo':
		mapres.mpProjection          = "LambertEqualArea"
		mapres.mpDataBaseVersion     = "MediumRes"
		mapres.mpLimitMode           = 'Corners'
		mapres.mpLeftCornerLatF      = 82
		mapres.mpLeftCornerLonF      = 50
		mapres.mpRightCornerLatF     = 75
		mapres.mpRightCornerLonF     = 120
		mapres.mpCenterLonF           = 0.
		mapres.mpCenterLatF           = 90.

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

def plt_vectors_colors(lon, lat, u_wind, v_wind, scalar_data, datamin = -2., datamax = 2., datastep = 0.2, vfMissingUValueV = None, sfMissingValueV = None, MinFracLengthF = 0.005, RefMagnitudeF = 0.3, RefLengthF = 0.05, LineArrowHeadMaxSizeF= 0.005, LineArrowHeadMinSizeF=0.005, RefAnnoString1 = 'vector', psname="output", sstep=1, region = 'Global', showfig=False, TitleString = 'variable', minLon=0, maxLon=360 , minLat=-80 , maxLat=85, colormap_name='posneg_1',timeon=False,levon=False, level=None, timstep=None, tunits=None, vcLineArrowThicknessF=1, overfigure=True):
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
	
	
	if region == 'Arctic':
		vcres.mpProjection          = "LambertEqualArea"
		vcres.mpDataBaseVersion     = "LowRes"
		vcres.mpLimitMode           = 'Corners'
		vcres.mpLeftCornerLatF      = 54
		vcres.mpLeftCornerLonF      = -50
		vcres.mpRightCornerLatF     = 57
		vcres.mpRightCornerLonF     = 140
		vcres.mpCenterLonF           = 0.
		vcres.mpCenterLatF           = 90.
			
	elif region == 'Global':
		vcres.mpProjection          = "LambertEqualArea"
		vcres.mpDataBaseVersion     = "LowRes"
		vcres.mpLimitMode           = "LatLon"
		vcres.mpMinLonF             = minLon
		vcres.mpMaxLonF             = maxLon
		vcres.mpMinLatF             = minLat
		vcres.mpMaxLatF             = maxLat
		vcres.mpCenterLonF           = 0.
		vcres.mpCenterLatF           = 0.
	elif region == 'NAtlantic':
		vcres.mpProjection          = "LambertEqualArea"
		vcres.mpDataBaseVersion     = "MediumRes"
		vcres.mpLimitMode           = 'Corners'
		vcres.mpLeftCornerLatF      = 65
		vcres.mpLeftCornerLonF      = -20
		vcres.mpRightCornerLatF     = 70
		vcres.mpRightCornerLonF     = 100
		vcres.mpCenterLonF           = 0.
		vcres.mpCenterLatF           = 90.
	elif region == 'FramStAnna':
		vcres.mpProjection          = "LambertEqualArea"
		vcres.mpDataBaseVersion     = "MediumRes"
		vcres.mpLimitMode           = 'Corners'
		vcres.mpLeftCornerLatF      = 75
		vcres.mpLeftCornerLonF      = -5
		vcres.mpRightCornerLatF     = 75
		vcres.mpRightCornerLonF     = 120
		vcres.mpCenterLonF           = 0.
		vcres.mpCenterLatF           = 90.
	elif region == 'Arctic2':
		vcres.mpProjection          = "LambertEqualArea"
		vcres.mpDataBaseVersion     = "MediumRes"
		vcres.mpLimitMode           = 'Corners'
		vcres.mpLeftCornerLatF      = 70
		vcres.mpLeftCornerLonF      = -60
		vcres.mpRightCornerLatF     = 60
		vcres.mpRightCornerLonF     = 140
		vcres.mpCenterLonF           = 0.
		vcres.mpCenterLatF           = 90.
	elif region == 'Anna':
		vcres.mpProjection          = "LambertEqualArea"
		vcres.mpDataBaseVersion     = "MediumRes"
		vcres.mpLimitMode           = 'Corners'
		vcres.mpLeftCornerLatF      = 83
		vcres.mpLeftCornerLonF      = 20
		vcres.mpRightCornerLatF     = 77
		vcres.mpRightCornerLonF     = 90
		vcres.mpCenterLonF           = 0.
		vcres.mpCenterLatF           = 90.
	elif region == 'Anna2':
		vcres.mpProjection          = "LambertEqualArea"
		vcres.mpDataBaseVersion     = "MediumRes"
		vcres.mpLimitMode           = 'Corners'
		vcres.mpLeftCornerLatF      = 82
		vcres.mpLeftCornerLonF      = 60
		vcres.mpRightCornerLatF     = 80
		vcres.mpRightCornerLonF     = 82
		vcres.mpCenterLonF           = 0.
		vcres.mpCenterLatF           = 90.
		
	elif region == 'Anna3':
		vcres.mpProjection          = "LambertEqualArea"
		vcres.mpDataBaseVersion     = "MediumRes"
		vcres.mpLimitMode           = 'Corners'
		vcres.mpLeftCornerLatF      = 70
		vcres.mpLeftCornerLonF      = 10
		vcres.mpRightCornerLatF     = 77
		vcres.mpRightCornerLonF     = 90
		vcres.mpCenterLonF           = 0.
		vcres.mpCenterLatF           = 90.
	elif region == 'halo':
		vcres.mpProjection          = "LambertEqualArea"
		vcres.mpDataBaseVersion     = "MediumRes"
		vcres.mpLimitMode           = 'Corners'
		vcres.mpLeftCornerLatF      = 85
		vcres.mpLeftCornerLonF      = 0
		vcres.mpRightCornerLatF     = 65
		vcres.mpRightCornerLonF     = 120
		vcres.mpCenterLonF           = 0.
		vcres.mpCenterLatF           = 90.
		
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
	