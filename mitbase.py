#!/bin/python
# -*- coding: utf-8 -*-


"""Komod plot module 
Contain mostly set of wrapper functions for map plotting with Basemap.
Can be used with any 2D data, not necessarily MITgcm."""

import numpy as np
try:
    from mpl_toolkits.basemap import Basemap
except:
    print('Basemap is not found, mitbase functions will not work')
    
def regbase(region):
    '''Takes name of the region and returns dictionary with
    information necessary for creation of the Basemap instance
    '''

    mapDict = {}

    if region == 'Arctic':
        mapDict['projection'] = 'npstere'
        mapDict['boundinglat'] = 60
        mapDict['lon_0'] = 0
        mapDict['resolution'] = 'l'
    if region == 'Glacindia':
        mapDict['projection'] = 'rotpole'
        mapDict['o_lon_p']=-123.34
        mapDict['o_lat_p']=79.95
        mapDict['lon_0'] = 180-123.34
        mapDict['lat_0']=90-79.95
        mapDict['llcrnrlon']=59
        mapDict['llcrnrlat']=4
        mapDict['urcrnrlon']=105
        mapDict['urcrnrlat']=38
        mapDict['resolution']="l"
    if region == "CORDEX_SA":
        mapDict['projection'] = 'rotpole'
        mapDict['o_lon_p']=-123.34
        mapDict['o_lat_p']=79.95
        mapDict['lon_0'] = 180-123.34
        mapDict['lat_0']=90-79.95
        mapDict['llcrnrlon']=21
        mapDict['llcrnrlat']=-19
        mapDict['urcrnrlon']=125
        mapDict['urcrnrlat']=45
        mapDict['resolution']="l"
                    
    return mapDict

def bp(lon, lat, data, region = 'Arctic', ptype = 'contourf', fillcont=False, **kwargs):
    
    '''Basic Basemap plot function. Use coordinates (1d or 2d), data and name of the region
     as an input and plot data. Region defines in the "regbase" function.

     You can also provide any argument for matplotlib plotting functions.

     Usage:
         bp(lon, lat, data, region = 'Arctic', ptype = 'contourf', **kwargs)
     
     Input:
        lon         - 2D or 1D array of longitudes
        lat         - 2D or 1D array of latitudes
        data        - 2D array of scalar data.
        region      - one of the predefined regions (for list of regions see the "regbase" function)
        ptype       - plot type (contour, contourf, pcolor, pcolormesh)
        **kwargs    - arguments for plotting functions

     Output:
        Basemap instance.
    '''
    
    mapDict = regbase(region)
   

    # Create Basemap instance
    if mapDict['projection'] == 'npstere':
        m = Basemap(projection=mapDict['projection'],boundinglat=mapDict['boundinglat'],\
                lon_0=mapDict['lon_0'],resolution=mapDict['resolution'])
    elif mapDict['projection'] == 'rotpole':
        m = Basemap(o_lon_p=mapDict['o_lon_p'],o_lat_p=mapDict['o_lat_p'],\
                    lon_0=mapDict['lon_0'],lat_0=mapDict['lat_0'],\
                    llcrnrlon=mapDict['llcrnrlon'],llcrnrlat=mapDict['llcrnrlat'],\
                    urcrnrlon=mapDict['urcrnrlon'],urcrnrlat=mapDict['urcrnrlat'],\
                    resolution=mapDict['resolution'],projection=mapDict['projection'])
    
    # Check if we have proper number of dimensions for lon (and hopefully lat as well)
    if lon.shape.__len__() == 1:
        lon, lat = np.meshgrid(lon, lat)
    elif lon.shape.__len__() > 2:
        raise Exception("Coordinate variables (lon) has too many dimensions")
    
    # Convert lat/lon to map coordinates
    x, y = m(lon, lat)

    # Make the map look better
    if fillcont == True:
        m.fillcontinents(color='gray',lake_color='gray')
    else:
        m.drawcoastlines()
    #m.drawparallels(np.arange(-80.,81.,20.))
    #m.drawmeridians(np.arange(-180.,181.,20.))
    m.drawmapboundary(fill_color='white')
    
    # Draw values on the map
    if ptype == 'contourf':
        cs = m.contourf(x,y,data, **kwargs)
    elif ptype == 'pcolormesh':
        cs = m.pcolormesh(x,y,data,**kwargs)
    elif ptype == 'contour':
        cs = m.contour(x,y,data,**kwargs)
    elif ptype == 'pcolor':
        cs = m.pcolor(x,y,data,**kwargs)
    else:
        raise Exception("Plot type not supported. Valid plot types are: contour, contourf, pcolor, pcolormesh ")
    
    return m


