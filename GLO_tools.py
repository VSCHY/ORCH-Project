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
def basicmap(indir, lonname, latname, varname, timestep, lonlattype="vect", save =False, name = "test"):
    
    olon = getvar(indir, lonname, 0)
    olat = getvar(indir, latname, 0)
    var = getvar(indir, varname, 0)

    nlon = np.min(olon)
    xlon = np.max(olon)
    nlat = np.min(olat)
    xlat = np.max(olat)

    fig=plt.figure(figsize=(4,3),dpi=300)

    m = Basemap(projection='cyl', resolution="l",llcrnrlon=nlon,
            llcrnrlat=nlat, urcrnrlon=xlon, urcrnrlat=xlat)
    if lonlattype == "vect":
        print "lonlat: vectors"
        lon, lat = np.meshgrid(olon[:], olat[:])
        xi, yi = m(lon, lat)
    else:
        xi, yi = m(olon[:,:], olat[:,:])

    #Draw
    cs = m.contourf(xi,yi,var[:,:,timestep],cmap=plt.get_cmap("rainbow"))
    # Map characterization
    m.drawparallels(np.arange(-90, 90,45), labels=[1,0,0,0], fontsize=2)
    m.drawmeridians(np.arange(-180., 180.,45), labels=[0,0,0,1], fontsize=2)
    m.drawcoastlines(linewidth=0.6)
    m.drawcountries(color='k',linewidth=0.6)
    # Colorbar
    cbar=m.colorbar(cs,location='right',pad=0.15)
    cbar.ax.tick_params(labelsize=10,direction="in",length = 2,pad=4)
    plt.setp(cbar.ax.get_xticklabels(), fontsize=18,weight='bold')
    
    #Adjust final result / Global title
    plt.subplots_adjust(wspace = 0.4, right=0.9,left=0.15,bottom=0.1, top= 0.95,hspace=0.32) 
    #Save
    if save: plt.savefig(chem_graph+name+".png",dpi=300)
    plt.show()
    plt.close()

### Plot Data simple ###
# Data 1D

# Data 2D : choose a location

