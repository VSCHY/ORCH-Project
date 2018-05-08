#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Extraction input ORCHIDEE (land coordinates)
Created on Wed May  2 12:47:43 2018

@author: Anthony Schrapffer
@institution: CIMA, Universidad de Buenos Aires (Argentina)
"""

# Libraries 
import os
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset as NetCDFFile
from netCDF4 import num2date

import matplotlib.pyplot as plt
from matplotlib import rc
from mpl_toolkits.basemap import Basemap
from optparse import OptionParser
import time

import sys
#sys.path.append("/home/anthony/TotiTools/ORCH-Project/")
sys.path.append("./")
import GLO_tools as GLO

ERROR = "*** --- *** ERROR *** --- ***"

######################
######## MAIN ########
######################

def get_DIM(indir):
    """
    Get the dimension (lon, lat, time) of the input file.
    indir: string, direction of the input file.
    """
    dimx = len(GLO.get_var(indir, "lon", 1))
    dimy = len(GLO.get_var(indir, "lat", 1))
    dimt = len(GLO.get_var(indir, "time", 1))
    return dimx, dimy, dimt

def set_ORCHForcEnv(indir, outdir, datan, dimx, dimy, dimt):
    """
    Create and set the dimension, variable and attribute of output file.
    indir: string, direction of the input file.
    outdir: string, direction of the output file.
    datan: string, name of the variable.
    dimx, dimy, dimt: int, dimension of (lon, lat, time)
    """
    # Create the file
    foo = NetCDFFile(outdir, 'w')
    
    # Creation dimensions
    foo.createDimension("lon", dimx)
    foo.createDimension("lat", dimy)
    foo.createDimension('time', dimt)

    # Environment variables
    for varn in ["lon", "lat", "time"]:
        ovar = GLO.get_var(indir, varn, 0)
        newvar = foo.createVariable(varn, ovar.dtype, ovar.dimensions)
        newvar[:] = ovar[:]
        for attrn in ovar.ncattrs():          
            attrv = getattr(ovar, attrn)
            newvar.setncattr(attrn, attrv)

    # Variable Environment
    odata = GLO.get_var(indir, datan, 0)
    newdata = foo.createVariable(datan, odata.dtype, ('lat', 'lon', 'time'))
    for attrn in odata.ncattrs():
        attrv = getattr(odata, attrn)
        newdata.setncattr(attrn, attrv)
    return newdata, foo

def convertlonlat_timstep(indir, datan, timestep, dimx, dimy, landref):
    """
    
    """
    print "Conversion", timestep
    ncvar = GLO.get_var(indir, datan, 0)
    out = ma.masked_all((dimy,dimx))
    var = ncvar[timestep,:]
    ind=0
    while ind<len(var):
        if ind % 100000. == 0 : print ind
        ref = landref[ind]
        k = int(ref/dimx)
        i = ref-k*dimx-1
        j = k-1
        out[j,i] = var[ind]
        ind=ind+1
    return out

def savedata(newdata, beg, Len, indir, outdir, datan, dimx, dimy):
    """
    """
    
    print beg, beg+Len
    landref = GLO.get_var(indir, "land", 0)
    for t in range(0,Len):
        print "loop",t, Len
        start = time.time()
        newdata[:,:,beg+t] = convertlonlat_timstep(indir, datan, beg+t, dimx, dimy, landref)
        print('Took', time.time()-start, 'seconds.')
    print "done"
    print "put data done"

def conversionlonlat(indir, outdir, datan, Lo):
    if os.path.isfile(outdir):
        print ERROR
        print "Output file already exists"
        return 
    convstart = time.time()
    print "Start CONVERSION at ", time.ctime() 
    print "***"
    print "Get dimension"
    dimx, dimy, dimt = get_DIM(indir)
    L=int(Lo)
    
    print "Set data environment"
    newdata, foo = set_ORCHForcEnv(indir, outdir, datan, dimx, dimy, dimt)
    N = int(dimt/float(L))
    R = dimt-N*L+1
    
    print "Start data conversion"
    for i in range(0,N):
        print ""
        print "Range : ",i,"/",N
        print "####"
        savedata(newdata, i*L, L, indir, outdir, datan, dimx, dimy)
    savedata(newdata, N*L, R, indir, outdir, datan, dimx, dimy)
    foo.sync()
    foo.close()
    print "Finished at ",time.time()-convstart

######################
####Option Parser ####
######################
operations = ["conversionlonlat"]

parser = OptionParser()
parser.add_option("-o", "--operation", type='choice', dest="operation", choices=operations, help="operation to make", metavar="OPER")
parser.add_option("-i", "--in", dest="indir", help="Input file")
parser.add_option("-w", "--out", dest="outdir", help="Output file")
parser.add_option("-v", "--var", dest="varname", help="Name of the variable")
parser.add_option("-l", "--len", dest="lengthloop", help="Length of the loop")

(opts, args) = parser.parse_args()

oper=opts.operation

if oper == 'conversionlonlat':
    conversionlonlat(opts.indir, opts.outdir, opts.varname, opts.lengthloop)

### Test ###
# python $DIRFORC -o conversionlonlat -i $indir -w $outdir -v $varn -l 500

