#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Global tools for python
Created on Fri May  4 15:12:27 2018

@author: Anthony Schrapffer
@institution: CIMA, Universidad de Buenos Aires (Argentina)
"""
# To import : 
# import sys
# sys.path.append("/home/anthony/TotiTools/ORCH-Project/") # To Change
#################
### Libraries ###
#################
import sys 
import numpy as np
from datetime import date
from netCDF4 import Dataset as NetCDFFile
from netCDF4 import num2date
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.patches as mpatches
import numpy.ma as ma
import matplotlib.lines as mlines
#mpl.use('Agg')
from mpl_toolkits.basemap import Basemap
from matplotlib.ticker import AutoMinorLocator
import matplotlib.ticker as ticker
rc('text', usetex=True)
import matplotlib.dates as mdates
import matplotlib.cbook as cbook


# Resume function
# OptionPearser

### Extraction Variables ###
def get_var(ncdir, varname, dim):
    # dim = 0 get ncdf variable
    salid = NetCDFFile(ncdir, 'r') 
    if dim == 0:
        var = salid.variables[varname]
    elif dim == 1:    
        var = salid.variables[varname][:]
    elif dim == 2:    
        var = salid.variables[varname][:,:]
    elif dim == 3:    
        var = salid.variables[varname][:,:,:]
    elif dim == 4:    
        var = salid.variables[varname][:,:,:,:]
    else: 
        print "Not recognized dim"
        var = None
    return var

### Plot map simple ###



### Plot Data simple ###
# Data 1D

# Data 2D : choose a location

