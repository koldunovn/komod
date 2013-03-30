komod
=====

Pyhton module that helps to postprocess and plot results from the ocean model (MITgcm).

This module is not intended for active public use, so the code probably will be dirty and comments rare :) 

Dependencies
=====
  - numpy, scipy, matplotlib
  - [PyNGL](http://www.pyngl.ucar.edu/)
  - psselect and [ps2eps](http://www.tm.uka.de/~bless/ps2eps) for some functions

Documentation:
=====
There are several collections of functions inside this module:

  - mitopen.py - Opens different file types, mainly those produced by MITgcm.
  - mitplot.py - Contain mostly set of wrapper functions for map plotting with [PyNGL](http://www.pyngl.ucar.edu/).
                 Can be used with any 2D data, not necessarily MITgcm.
  - mittime.py - Convert timesteps to time and vice versa for MITgcm.
  - ut.py      - Set of different small utilits that makes life easier.

Documentation in IPython notebook format located in doc folder, and also available through nbviewer:

  - [mitopen](http://nbviewer.ipython.org/urls/raw.github.com/koldunovn/komod/master/doc/mitopen.ipynb) 
  - [mitplot](http://nbviewer.ipython.org/urls/raw.github.com/koldunovn/komod/master/doc/mitplot.ipynb)



