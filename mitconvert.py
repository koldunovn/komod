#!/bin/python
# -*- coding: utf-8 -*-

"""Komod open module
Converts different files to netcdf

Nikolay Koldunov 10 July 2012
"""

import numpy
import Nio
import os

def adxxcnv(xdim, ydim, xcdata='./', ycdata='./', bswap=1):
	"""open MITgcm binary XC.data and YC.data files

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
