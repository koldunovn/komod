#!/bin/python
# -*- coding: utf-8 -*-


"""Komod module with different small utilits that makes life easier
 
pssplit 	- split multipage .ps file in to one-page files.
pssplit2eps	- split multipage .ps file in to one-page files and then convert them to eps
mshow           - show 2d matrix
epspair         - create .pdf file with images arranged in two columns
pointfind       - find grid point that is closest to some lat/lon point 


Nikolay Koldunov 18 May 2010
"""

# -------------------------------------------------
import os
import matplotlib.pyplot as plt
import numpy
np = numpy
try:
	import Ngl
except ImportError:
	pass
	#print('Ngl is not avalible, some functions will not work')
import smtplib

try:
	import pandas as pd
except ImportError:
	pass
import glob
from netCDF4 import Dataset, MFDataset

def pssplit(filename, npages):
	"""split multipage .ps file in to one-page files.
	It use unix psselect utilite
	    
	   Usage: pssplit(filename, npages)

    Input:
	filename    - file name
	npages      - number of pages

    Output:
        one-page ps files 

    """
	for page in range(npages):
		os.system("psselect -p"+str(page)+" "+filename+" "+filename[:-3]+"_"+str(page).zfill(3)+".ps") 
		

def pssplit2eps(filename, npages):
	"""split multipage .ps file in to one-page files and then convert them to eps 
	It use unix psselect and ps2eps (http://www.tm.uka.de/~bless/ps2eps) utilites
	    
	   Usage: pssplit2eps(filename, npages)

    Input:
	filename    - file name
	npages      - number of pages

    Output:
        one-page eps files 

    """
	for page in range(npages):
		os.system("psselect -p"+str(page)+" "+filename+" "+filename[:-3]+"_"+str(page).zfill(3)+".ps") 
		os.system("ps2eps -f "+filename[:-3]+"_"+str(page).zfill(3)+".ps")
		os.system("rm "+filename[:-3]+"_"+str(page).zfill(3)+".ps")


def mshow(matrix, cmap=None, norm=None, aspect=None, interpolation=None,
       alpha=1.0, vmin=None, vmax=None, origin=None, extent=None):
	""" show 2d matrix
	this is just pyplot.imshow that displays figure imidiately and with color bar
	Usage: 
		mshow(matrix,cmap=None, norm=None, aspect=None, interpolation=None,
       			   alpha=1.0, vmin=None, vmax=None, origin=None, extent=None)
	Input:
		matrix - 2d matrix
	   	description of other options can be found at: http://matplotlib.sourceforge.net/api/pyplot_api.html#matplotlib.pyplot.imshow
	Output:
		2D plot of matrix values
	"""
	plt.imshow(matrix,cmap=cmap, norm=norm, aspect=aspect, interpolation=interpolation,
       alpha=1.0, vmin=vmin, vmax=vmax, origin=origin, extent=extent)
	plt.colorbar()
	plt.show()

def eps2(filename1, filename2, npair, ofname="output",nstart=0, showfig=True):
	""" create .pdf file with images arranged in two columns 
	    (useful, for example, when you want to compare output of two different runs)
	Dependencies:
		python.sty	- have to be in your working directory or Latex path.
				  Can be obtained from here http://www.imada.sdu.dk/~ehmsen/python.sty
		latex, dvipdf, gv
	Usage: 
		epspair(filename1, filename2, npair, ofname="output",nstart=0,)
	Input:
		filename1	- constant part of .eps files
				  that will be displayed in the left column.
				  The rest of the name should be 3 char digit with zeros in front and .eps extension
				  E.g. my_file_007.eps (here filename1="my_file_"), run1_040.eps, model1_900.eps 
				  This files can be generated by komod.pssplit or komod.pssplit2eps
				  
		filename2	- the same as filename1, but for images of the right column.
		npair 		- number of pairs (rows in out two column table)
		ofname		- output file name, default = "output"
		nstart		- first number we will start with
		showfig		- if True display the figure with gv
	Output:
		ofname.pdf file with images arranged in two columns
		
		
	
	"""

	ofnametex = ofname+".tex"
		
	f = open(ofnametex, 'w')
	
	
	f.write('''\documentclass   {article}
	\usepackage{python}
	\usepackage{graphicx}
	\usepackage[top=1cm, bottom=1cm, left=1cm, right=2cm]{geometry}
	\\begin{document}
	
	\\begin{python}
	for i in range('''+str(nstart)+''', '''+str(npair)+'''):
	
	    print r'\\begin{center}'
	    print r'\parbox{0.40\linewidth}{\includegraphics[width=4.5cm, angle=0]{'''+filename1+''''+str(i).zfill(3)+'.eps}\\\\'
	    print r' \centering  }'
	    print r'\hspace{0.1\linewidth}'
	    print r'\parbox{0.40\linewidth}{\includegraphics[width=4.5cm, angle=0]{'''+filename2+''''+str(i).zfill(3)+'.eps}\\\\'
	    print r'\centering }'
	
	    print r'\end{center}'
	
	
	
	\end{python}
	
	\end{document}
	
	''')
	
	f.close()
	
	os.system("latex --shell-escape "+ofnametex)
	os.system("dvipdf "+ofnametex[:-3]+"dvi")
	if showfig==True:
		os.system("gv "+ofnametex[:-3]+"pdf")
		

def eps4(filename1, filename2, filename3, filename4, nquad, ofname="output",nstart=0, showfig=True):
	""" create .pdf file with images arranged in 4 columns 
	    (useful, for example, when you want to compare output of 4 different runs)
	Dependencies:
		python.sty	- have to be in your working directory or Latex path.
				  Can be obtained from here http://www.imada.sdu.dk/~ehmsen/python.sty
		latex, dvipdf, gv
	Usage: 
		epsquad(filename1, filename2, filename3, filename4, npair, ofname="output",nstart=0,)
	Input:
		filename1	- constant part of .eps files
				  that will be displayed in the left (1st) column.
				  The rest of the name should be 3 char digit with zeros in front and .eps extension
				  E.g. my_file_007.eps (here filename1="my_file_"), run1_040.eps, model1_900.eps 
				  This files can be generated by komod.pssplit or komod.pssplit2eps
				  
		filename2	- the same as filename1, but for images of the 2nd column.
		filename3	- the same as filename1, but for images of the 3rd column.
		filename4	- the same as filename1, but for images of the 4th column.
		nquad 		- number of quads (rows in the 4 column table)
		ofname		- output file name, default = "output"
		nstart		- first number we will start with
		showfig		- if True display the figure with gv
	Output:
		ofname.pdf file with images arranged in two columns
		
		
	
	"""

	ofnametex = ofname+".tex"
		
	f = open(ofnametex, 'w')
	
	
	f.write('''\documentclass{article}
	\usepackage{python}
	\usepackage{graphicx}
	\usepackage[top=1cm, bottom=1cm, left=1cm, right=2cm]{geometry}
	\\begin{document}
	
	\\begin{python}
	for i in range('''+str(nstart)+''', '''+str(nquad)+'''):
	
	    print r'\\begin{center}'
	    print r'\parbox{0.16\linewidth}{\includegraphics[width=4.5cm, angle=0]{'''+filename1+''''+str(i).zfill(3)+'.eps}\\\\'
	    print r' \centering  }'
	    print r'\hspace{0.1\linewidth}'
	    print r'\parbox{0.16\linewidth}{\includegraphics[width=4.5cm, angle=0]{'''+filename2+''''+str(i).zfill(3)+'.eps}\\\\'
	    print r'\centering }'
	    print r'\hspace{0.1\linewidth}'
	    print r'\parbox{0.16\linewidth}{\includegraphics[width=4.5cm, angle=0]{'''+filename3+''''+str(i).zfill(3)+'.eps}\\\\'
	    print r'\centering }'
	    print r'\hspace{0.1\linewidth}'	
	    print r'\parbox{0.16\linewidth}{\includegraphics[width=4.5cm, angle=0]{'''+filename4+''''+str(i).zfill(3)+'.eps}\\\\'
	    print r'\centering }'
	
	
	    print r'\end{center}'
	
	
	
	\end{python}
	
	\end{document}
	
	''')
	
	f.close()
	
	os.system("latex --shell-escape "+ofnametex)
	os.system("dvipdf "+ofnametex[:-3]+"dvi")
	if showfig==True:
		os.system("gv "+ofnametex[:-3]+"pdf")
		




def eps5(filename1, filename2, filename3, filename4, filename5, nquad, ofname="output",nstart=0, showfig=True):
	""" create .pdf file with images arranged in 5 columns 
	    (useful, for example, when you want to compare output of 5 different runs)
	Dependencies:
		python.sty	- have to be in your working directory or Latex path.
				  Can be obtained from here http://www.imada.sdu.dk/~ehmsen/python.sty
		latex, dvipdf, gv
	Usage: 
		epsquad(filename1, filename2, filename3, filename4, filename5, npair, ofname="output",nstart=0,)
	Input:
		filename1	- constant part of .eps files
				  that will be displayed in the left (1st) column.
				  The rest of the name should be 3 char digit with zeros in front and .eps extension
				  E.g. my_file_007.eps (here filename1="my_file_"), run1_040.eps, model1_900.eps 
				  This files can be generated by komod.pssplit or komod.pssplit2eps
				  
		filename2	- the same as filename1, but for images of the 2nd column.
		filename3	- the same as filename1, but for images of the 3rd column.
		filename4	- the same as filename1, but for images of the 4th column.
		filename5	- the same as filename1, but for images of the 5th column.
		nquad 		- number of quads (rows in the 4 column table)
		ofname		- output file name, default = "output"
		nstart		- first number we will start with
		showfig		- if True display the figure with gv
	Output:
		ofname.pdf file with images arranged in two columns
		
		
	
	"""

	ofnametex = ofname+".tex"
		
	f = open(ofnametex, 'w')
	
	
	f.write('''\documentclass{article}
	\usepackage{python}
	\usepackage{graphicx}
	\usepackage[top=1cm, bottom=1cm, left=1cm, right=2cm]{geometry}
	\\begin{document}
	
	\\begin{python}
	for i in range('''+str(nstart)+''', '''+str(nquad)+'''):
	
	    print r'\\begin{center}'
	    print r'\parbox{0.105\linewidth}{\includegraphics[width=3.7cm, angle=0]{'''+filename1+''''+str(i).zfill(3)+'.eps}\\\\'
	    print r' \centering  }'
	    print r'\hspace{0.1\linewidth}'
	    print r'\parbox{0.105\linewidth}{\includegraphics[width=3.7cm, angle=0]{'''+filename2+''''+str(i).zfill(3)+'.eps}\\\\'
	    print r'\centering }'
	    print r'\hspace{0.1\linewidth}'
	    print r'\parbox{0.105\linewidth}{\includegraphics[width=3.7cm, angle=0]{'''+filename3+''''+str(i).zfill(3)+'.eps}\\\\'
	    print r'\centering }'
	    print r'\hspace{0.1\linewidth}'	
	    print r'\parbox{0.105\linewidth}{\includegraphics[width=3.7cm, angle=0]{'''+filename4+''''+str(i).zfill(3)+'.eps}\\\\'
	    print r'\centering }'
	    print r'\hspace{0.1\linewidth}'	
	    print r'\parbox{0.105\linewidth}{\includegraphics[width=3.7cm, angle=0]{'''+filename5+''''+str(i).zfill(3)+'.eps}\\\\'
	    print r'\centering }'
	
	
	    print r'\end{center}'
	
	
	
	\end{python}
	
	\end{document}
	
	''')
	
	f.close()
	
	os.system("latex --shell-escape "+ofnametex)
	os.system("dvipdf "+ofnametex[:-3]+"dvi")
	if showfig==True:
		os.system("gv "+ofnametex[:-3]+"pdf")
		






def pointfind(plat, plon, lat, lon, pdif = 1):
	""" return indeces and values of the grid point closest to lat/lon point
	Usage: 
		pointfind(plat, plon, lat, lon, pdif = 0.5)
	Input:
		plat - latitude of the point
		plon - longitude of the point
		lat  - 2d array of latutudes
		lon  - 2d array of longitudes
		pdif  - initial +- window
	Output:
		indeces and values of the points that fulfil conditions 
		
		
	
	"""
	
	fff = 10
	while (fff > 1):
		
		#conditions for latitude (lat - 2d array of latitudes)
		c_lat=(lat>(plat-pdif))&(lat<(plat+pdif))
		#conditions for longiyude (lon - 2d array of longitudes)
		c_lon=(lon>(plon-pdif))&(lon<(plon+pdif))
		
		#combine both conditions together
		c_all=c_lat&c_lon
		
		#values of the points that fulfil conditions
		platf = lat[numpy.nonzero(c_all)]
		plonf = lon[numpy.nonzero(c_all)]
		
				
		#indeces of the poin that fulfil conditions 
		g = numpy.nonzero(c_all)
		
		
		#check if we have found uniq solution
		fff = platf.shape[0]
		# decrease window to reduce amount of solutions if we have more than one
		#print(pdif)
		pdif = pdif-0.001
	print("coordinates of the point that fulfil conditions: "+str(platf)+" "+str(plonf))
	print("indeces of the point that fulfil conditions: "+str(g[0])+" "+str(g[1]))
	
	return(g, platf, plonf)

def pointfind2(plat, plon, lat, lon, pdif=1):
	""" return indeces and values of the grid point closest to lat/lon point
	the same as pointfind but could be faster
	Usage: 
		pointfind(plat, plon, lat, lon, pdif = 0.5)
	Input:
		plat - latitude of the point
		plon - longitude of the point
		lat  - 2d array of latutudes
		lon  - 2d array of longitudes
		pdif  - we don't need it but leave it to have the same input as pointfind
	Output:
		indeces and values of the points that fulfil conditions 
		
		
	
	"""

	dist_min = 1000000.
	
	
	for i in range(lon.shape[0]):
		for j in range(lon.shape[1]):
			dist = Ngl.gc_dist(plat,plon,lat[i,j],lon[i,j])
			if dist_min > dist:
				dist_min = dist
				i_min = i
				j_min = j
				lat_min = lat[i,j]
				lon_min = lon[i,j]
	
	print(i_min,j_min,lat_min,lon_min)
	gg1 = i_min, j_min
	
	return(gg1, lat_min, lon_min)




def get_transect(lat, lon, data, lat1, lon1, lat2, lon2, npoints = 10, pdif = 1, norep=False):
	

	plat, plon = Ngl.gc_interp(lat1,lon1,lat2,lon2,npoints)
	
	m_grid = numpy.zeros((plon.shape[0]))
	n_grid = numpy.zeros((plat.shape[0]))
	grid_lats = numpy.zeros((plat.shape[0]))
	grid_lons = numpy.zeros((plat.shape[0]))
	
	for k in range(plon.shape[0]):
		coord, trash1, trash1 = pointfind2(plat[k], plon[k], lat, lon, pdif)
		m_grid[k] = coord[0]
		n_grid[k] = coord[1]
		grid_lats[k] = lat[m_grid[k],n_grid[k]]
		grid_lons[k] = lon[m_grid[k],n_grid[k]]
	
	
	if len(data.shape) == 4:
		data = data[0,:,:,:]
	elif len(data.shape) == 5:
		data = data[0,0,:,:,:]
	
		
	data_prof = numpy.zeros(([data.shape[0],npoints]))   
	
	
	
	for j in range(data_prof.shape[1]):
		data_prof[:,j] = data[:,m_grid[j],n_grid[j]]
	
	
	x_kilometers = numpy.zeros((grid_lats.shape[0]))
	
	for j in range(grid_lats.shape[0] - 1):
		angle = Ngl.gc_dist(grid_lats[j], grid_lons[j], grid_lats[j+1], grid_lons[j+1] )
		x_kilometers[j+1] = Ngl.gc_convert(angle, 2) + x_kilometers[j]
	
	
	numpy.savez("outfile.npz", data_prof, x_kilometers)
	
	
	
	if norep == True:
		data_prof_norep = numpy.zeros(([data_prof.shape[0],1]))
		x_kilometers_norep = numpy.array(([]))
		m_grid_norep       = numpy.array(([]))
		n_grid_norep       = numpy.array(([]))
		present_location = -999999.999

		for k in range(x_kilometers.shape[0]):
			if int(x_kilometers[k]) != int(present_location):
				present_location = x_kilometers[k]
				data_prof_norep = numpy.hstack((data_prof_norep, data_prof[:,k:k+1])) 
				x_kilometers_norep = numpy.append(x_kilometers_norep, x_kilometers[k:k+1])
				m_grid_norep = numpy.append(m_grid_norep, m_grid[k:k+1])
				n_grid_norep = numpy.append(n_grid_norep, n_grid[k:k+1])
			

		data_prof_norep = data_prof_norep[:,1:] 
		data_prof = data_prof_norep
		x_kilometers = x_kilometers_norep
		m_grid = m_grid_norep
		n_grid = n_grid_norep
		
	print(x_kilometers)
	print(m_grid)
	print(n_grid)
	
	return(data_prof, x_kilometers, m_grid, n_grid )
	
def get_cost(cost_list, cost_var):
    ''' Get cost value from costfunction???? files
    
    Input:
    cost_list - list of lines from costfunction???? file
    cost_var  - name of the variable
    
    Output:
    value of cost_var
    
    '''
    for line in cost_list:
        if line.startswith(cost_var):
            nn = float(line.split()[2].replace('D','E'))
            return nn

def get_2cost(modelYear, cost_var, iterations=[0,1], pathToCost='./CPCOST/'):
    '''For specific variable gets cost values from costfinction files
    and compare them (difference in %)
    
    Input:
    
    modelYear  - year of interest
    cost_var   - name of the variable
    iterations - list with numbers of iterations
    pathToCost - path to files with costfunction files, organised in YYYY/itX/ way
    
    Output:
    
    cost values for two variables and theit difference in %
    
    '''
    
    if cost_var in ['f_smrarea', 'f_smrsst', 'f_smrsss']:
        cost_var = ' '+cost_var
    
        f1 = open(pathToCost+"/"+str(modelYear)+'/costfunction_seaice000'+str(iterations[0]), 'r')
        f2 = open(pathToCost+"/"+str(modelYear)+'/costfunction_seaice000'+str(iterations[1]), 'r')
        
    else:
        f1 = open(pathToCost+"/"+str(modelYear)+'/costfunction000'+str(iterations[0]), 'r')
        f2 = open(pathToCost+"/"+str(modelYear)+'/costfunction000'+str(iterations[1]), 'r')
    

    f1d = f1.read().splitlines()
    f2d = f2.read().splitlines()
    
    val1 = get_cost(f1d, cost_var)
    val2 = get_cost(f2d, cost_var)
    
    if val1 != 0:
        perc = 100 - ((val2/(val1))*100)
    else:
        perc=0
    
    return val1, val2, perc
 
def sendemail(from_addr, to_addr_list, cc_addr_list,
              subject, message,
              login, password,
              smtpserver='smtp.gmail.com:587'):
    header  = 'From: %s\n' % from_addr
    header += 'To: %s\n' % ','.join(to_addr_list)
    header += 'Cc: %s\n' % ','.join(cc_addr_list)
    header += 'Subject: %s\n\n' % subject
    message = header + message
 
    server = smtplib.SMTP(smtpserver)
    server.starttls()
    server.login(login,password)
    problems = server.sendmail(from_addr, to_addr_list, message)
    server.quit()

def calc_area(data, area_size, lat=None, blat=None):
    ''' Calculates sea ice area
    
    Input:
    
    data      - 2d field with AREA values
    area_size - usually dxc*dyc
    lat       - 2d array of latitudes
    blat      - boundary latitude (southern border)
    
    Output:
    
    sea ice area in m^2
    '''

    if blat != None:
        ii, jj = np.where((lat<blat))
        data[ii,jj] = 0.
        
    data[:,:][data[:,:]<0.15]=0.
    data_masked      = np.ma.masked_equal(data[:,:], 0)
    data_masked_mult = data_masked.astype('float64')*area_size.astype('float64')
        
    return np.ma.sum(data_masked_mult)

def calc_extent(data, area_size, lat=None, blat=None):
    ''' Calculates sea ice area
    
    Input:
    
    data      - 2d field with AREA values
    area_size - usually dxc*dyc
    lat       - 2d array of latitudes
    blat      - boundary latitude (southern border)
    
    Output:
    
    sea ice extent in m^2
    '''

    if blat != None:
        ii, jj = np.where((lat<blat))
        data[ii,jj] = 0.
        
    data[:,:][data[:,:]<0.15]=0.
    data[:,:][data[:,:]>=0.15]=1.
    data_masked      = np.ma.masked_equal(data[:,:], 0)
    data_masked_mult = data_masked.astype('float64')*area_size.astype('float64')
    
    return np.ma.sum(data_masked_mult)  

def ice_comp_model_to_osi(pathToModel, modelYear, modelIteration, boundLat, pathToOSI, param = 'area'):
    '''
    Plot sea ice area from satellite data and several model iterations
    '''
    fsat = MFDataset(pathToOSI+'/OSI'+modelYear+'??.nc')
    osi_lat = fsat.variables['lat'][:]
    osi_lon = fsat.variables['lon'][:]

    osi_ice = fsat.variables['ice_conc'][0,:,:]
    area_osi = np.ones(osi_ice.shape)*100

    osi_area = []
    for mm in range(12):
        if param == 'area':
            osi_area.append(calc_area(fsat.variables['ice_conc'][mm,:,:]/100,\
                                      area_osi, osi_lat, blat=boundLat)/1e6)
        elif param == 'extent':
        	osi_area.append(calc_extent(fsat.variables['ice_conc'][mm,:,:]/100,\
        	                area_osi, osi_lat, blat=boundLat)/1e6)
        
    
    g  = Dataset('/scratch/local1/POL06/1_test/grid.cdf')
    dxc = g.variables['dxc'][0,:,:]
    dyc = g.variables['dyc'][0,:,:]
    lat = g.variables['yc'][0,:,:]
    dxcXdyc = dxc*dyc

    area_model=np.zeros((len(modelIteration), 12))

    for (it, iteration) in enumerate(modelIteration):
        fm = MFDataset(pathToModel+'/'+modelYear+'/'+'it'+str(iteration)+'/fw/*.cdf')
        for mm in range(12):
        	if param == 'area':
        	    area_model[it,mm] = calc_area(fm.variables['area'][mm,:,:],\
        	    	                          dxcXdyc, lat, blat=boundLat)/10e11
        	elif param == 'extent':
        		area_model[it,mm] = calc_extent(fm.variables['area'][mm,:,:],\
        			                          dxcXdyc, lat, blat=boundLat)/10e11
    
        fm.close()
    
    dates = pd.date_range(modelYear+'-01', str(int(modelYear)+1)+'-01', freq='M')
    dd = pd.DataFrame(index=dates)

    dd['Satellite']=osi_area

    for (it , iteration) in enumerate(modelIteration):
        dd['it'+str(iteration)]=area_model[it,:]

    return dd.plot(figsize=(10,5))


def ice_comp_model_to_osi_table(pathToModel, modelYears, modelIteration, boundLat, pathToOSI, param = 'area'):

    diff_array = numpy.zeros((len(modelYears), 12))

    for (nnum, yyear) in enumerate(modelYears):
    
        fsat = MFDataset(pathToOSI+'/OSI'+yyear+'??.nc')
        osi_lat = fsat.variables['lat'][:]
        osi_lon = fsat.variables['lon'][:]

        osi_ice = fsat.variables['ice_conc'][0,:,:]
        area_osi = np.ones(osi_ice.shape)*100

        osi_area = []
        for mm in range(12):
            if param == 'area':
                osi_area.append(calc_area(fsat.variables['ice_conc'][mm,:,:]/100,\
                                      area_osi, osi_lat, blat=boundLat)/1e6)
            elif param == 'extent':
                osi_area.append(calc_extent(fsat.variables['ice_conc'][mm,:,:]/100,\
                            area_osi, osi_lat, blat=boundLat)/1e6)
        
    
        g  = Dataset('/scratch/local1/POL06/1_test/grid.cdf')
        dxc = g.variables['dxc'][0,:,:]
        dyc = g.variables['dyc'][0,:,:]
        lat = g.variables['yc'][0,:,:]
        dxcXdyc = dxc*dyc

        area_model=np.zeros((len(modelIteration), 12))
    
        if modelIteration[0] == 'last':
            gg = glob.glob(pathToModel+'/'+yyear+'/'+'it*')
            gg.sort()
            lastit = [int(gg[-1].split('/')[-1].split('t')[-1])]
        else:
            lastit = modelIteration

        for (it, iteration) in enumerate(lastit):
            fm = MFDataset(pathToModel+'/'+yyear+'/'+'it'+str(iteration)+'/fw/*.cdf')
            for mm in range(12):
                if param == 'area':
                    area_model[it,mm] = calc_area(fm.variables['area'][mm,:,:],\
        	    	                          dxcXdyc, lat, blat=boundLat)/10e11
                elif param == 'extent':
                    area_model[it,mm] = calc_extent(fm.variables['area'][mm,:,:],\
        			                          dxcXdyc, lat, blat=boundLat)/10e11
    
            fm.close()
        diff_array[nnum,:] = area_model[0,:]-osi_area[:]
    return diff_array

    
