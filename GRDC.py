#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 18:04:34 2018

@author: anthony
"""

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


GRDC = "/home/anthony/Documents/Doctorat/ORCHIDEE/Z0-Originals/GRDC_Monthly_June14.nc"

salid = NetCDFFile(GRDC, 'r')
var = salid.variables["name"]
area = salid.variables["area"]
lon = salid.variables["lon"]
lat = salid.variables["lat"]
river = salid.variables["river"]
mhydro = salid.variables["mergedhydro"]
hydro = salid.variables["hydrographs"]
chydro = salid.variables["calculatedhydro"]
time = salid.variables["time"]
dtime = num2date(time[:], time.units)

name=var[:]

def convname(L):
    iend = lastletter(L)
    i=0
    name=""
    while i<iend:
        name=name+L[i]
        i=i+1
    return name
    
    
def lastletter(L):
    i=0 
    while L[59-i] == "":
        i=i+1
    return 59-i+1
        
def GRDCname(GRDC):
    salid = NetCDFFile(GRDC, 'r')
    var = salid.variables["name"][:]
    NAME = []
    i=0
    while i<len(var):
        NAME.append(convname(var[i]))
        if "Porto Murt" in NAME[i]: print i,NAME[i]
        i=i+1
    return NAME

NAME = GRDCname(GRDC)
i = 2917
print name[i]
print area[i]
print lon[i]
print lat[i]
h= hydro[:,i]
mh=mhydro[:,i]
ch=chydro[:,i]

th = np.where(ma.getmaskarray(h)==False)[0]
tmh = np.where(np.isnan(mh)==False)[0]
tch = np.where(ma.getmaskarray(ch)==False)[0]

th0 = th[0]; th1 = th[len(th)-1]
tmh0 = tmh[0]; tmh1 = tmh[len(tmh)-1]
tch0 = tch[0]; tch1 = tch[len(tch)-1]

print "hydrographs: ", dtime[th0],"   ", dtime[th1]
print "merged: ", dtime[tmh0],"   ", dtime[tmh1]
print "calculated: ", dtime[tch0],"   ", dtime[tch1]
