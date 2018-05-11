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



def set_ORCHForcEnv(indir, outdir, datan, dimx, dimy, dimt, dt=True, beg=True):
    """
    Create and set the dimension, variable and attribute of output file.

    indir: string, direction of the input file.
    outdir: string, direction of the output file.
    datan: string, name of the variable.
    dimx, dimy, dimt: int, dimension of (lon, lat, time)
    dt: True or int, True if full file, int for dtimestep (parallelization)
    beg: True or int, int for parallelization (begin step)
    """
    # Create the file
    foo = NetCDFFile(outdir, 'w')
    
    # Creation dimensions
    foo.createDimension("lon", dimx)
    foo.createDimension("lat", dimy)
    if dt:
        foo.createDimension('time', dimt)
    else:
        foo.createDimension('time', dt)

    # Environment variables
    for varn in ["lon", "lat", "time"]:
        ovar = GLO.get_var(indir, varn, 0)
        newvar = foo.createVariable(varn, ovar.dtype, ovar.dimensions)
        if varn="time"
            if dt:
                newvar[:] = ovar[:]
            else:
                newvar[:] = ovar[beg:beg+dt] 
        else:
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
    Convertion of lon / lat of ORCHIDEE forcing for a specific timestep.

    indir: str, Input NetCDF direction (ORCHIDEE Forcing).
    datan: str, name of the data in the Input file.
    timestep: int, Timestep index.
    dimx: int, dimension longitude.
    dimy: int, dimension latitude.
    landref: Land index (From netCDF).
    """
    print "Conversion", timestep
    ncvar = GLO.get_var(indir, datan, 0)
    out = ma.masked_all((dimy,dimx))
    var = ncvar[timestep,:]
    ind=0
    while ind<len(var):
        #if ind % 100000. == 0 : print ind
        ref = landref[ind]
        k = int(ref/dimx)
        i = ref-k*dimx-1
        j = k-1
        out[j,i] = var[ind]
        ind=ind+1
    return out



def savedata(beg, Len, indir, datan, dimx, dimy):
    """
    Convert lon / lat for ORCHIDEE Forcing for a LEN length loop beginning by beg.

    # newdata: SUPPRESS - Just if modify output data in this funcion
    beg: int, beginning timestep for operation.
    Len: int, length of the loop.
    indir: str, Input NetCDF direction (ORCHIDEE Forcing).
    datan: str, name of the data in the Input file.
    dimx: int, dimension longitude.
    dimy: int, dimension latitude.        
    """
    print beg, beg+Len
    landref = GLO.get_var(indir, "land", 0)
    out=ma.masked_all((dimy, dimx,Len))
    for t in range(0,Len):
        print "loop",t, Len
        start = time.time()
        out[:,:,t]= convertlonlat_timstep(indir, datan, beg+t, dimx, dimy, landref)
        #newdata[:,:,beg+t] = convertlonlat_timstep(indir, datan, beg+t, dimx, dimy, landref)
        print('Took', time.time()-start, 'seconds.')
    print "done"
    print "put data done"



def conversionlonlat(indir, outdir, datan, Lo):
    """
    Operation of conversion lon/lat ORCHIDEE Forcing.

    indir: str, Input NetCDF direction (ORCHIDEE Forcing).
    outdir: str, Output NetCDF direction (ORCHIDEE Forcing)
    datan: str, name of the data in the Input file.
    Lo: int, length of the loop.
    """
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

    print ""
    print "### Start data conversion ###"
    for i in range(0,N):
        timestartrange = time.time()
        print ""
        print "Range : ",i,"/",N
        print "####"
        newdata[:,:,i*L:(i+1)*L] = savedata(i*L, L, indir, datan, dimx, dimy)
    savedata(N*L, R, indir, datan, dimx, dimy)
    print ('Took', (time.time()-timestartrange)/60, 'Minutes.')
    foo.sync()
    foo.close()
    print "Finished at ",time.ctime()
    print "Took",(time-time()-convstart)/60, "minutes"




# get nombre de part puis lancer pour 0, 1, 2, 3 jusque nb part (voir avt)
# outdirpart  = outdirpart0.nc & parto = 0
def conversionlonlat_parallel(indir, outdirpart, datan, Lo, parto):
    """
    Operation of conversion lon/lat ORCHIDEE Forcing used for parallelization.

    indir: str, Input NetCDF direction (ORCHIDEE Forcing).
    outdirpart: str, Output NetCDF direction (ORCHIDEE Forcing)
    datan: str, name of the data in the Input file.
    Lo: int, length of the loop.
    parto: int, part of the parallelization (start at 0)
    """
    if os.path.isfile(outdirpart):
        print ERROR
        print "Output file already exists"
        return
    # Inuput parameters
    part = int(parto)
    L=int(Lo)
    print ""
    print "### Part",part,"###"
    convstart = time.time()
    print "### Start CONVERSION at ", time.ctime(), "###"
    print "***"

    print "Get dimension"
    dimx, dimy, dimt = get_DIM(indir)
    # Is necessary part ?
    if L*part>=dimt:
        print "Unnecessary part"
        return
    # Parameters
    N = int(dimt/float(L))
    R = dimt-N*L+1
    # Start operation
    print "Set data environment"
    if part<N:
        newdata, foo = set_ORCHForcEnv(indir, outdirpart, datan, dimx, dimy, dimt,500,part*L)
    else:
        newdata, foo = set_ORCHForcEnv(indir, outdirpart, datan, dimx, dimy, dimt,dimt-part*L,part*L)

    print "Start data conversion"
    if part<N:
        print ""
        print "Range : ",i,"/",N
        print "####"
        newdata[:,:,:] = savedata(i*L, L, indir, datan, dimx, dimy)
    else:
        print "Rest :", i*500, "/", dimt
    	newdata[:,:,:] = savedata(N*L, R, indir, datan, dimx, dimy)
    foo.sync()
    foo.close()
    print 
    print "Finished at ",time.ctime()
    print "Took",(time-time()-convstart)/60, "minutes"


######################
####Option Parser ####
######################
operations = ["conversionlonlat", "conversionlonlat_parallel"]

parser = OptionParser()
parser.add_option("-o", "--operation", type='choice', dest="operation", choices=operations, help="operation to make", metavar="OPER")
parser.add_option("-i", "--in", dest="indir", help="Input file")
parser.add_option("-w", "--out", dest="outdir", help="Output file")
parser.add_option("-v", "--var", dest="varname", help="Name of the variable")
parser.add_option("-l", "--len", dest="lengthloop", help="Length of the loop")
parser.add_option("-p", "--part", dest="part", help="Part of the loop")

(opts, args) = parser.parse_args()

oper=opts.operation

if oper == 'conversionlonlat':
    conversionlonlat(opts.indir, opts.outdir, opts.varname, opts.lengthloop)
elif oper == 'conversionlonlat_parallel":
    conversionlonlat_parallel(opts.indir, opts.outdir, opts.varname, opts.lengthloop, opts.part)

### Test ###
# python $DIRFORC -o conversionlonlat -i $indir -w $outdir -v $varn -l 500
# python $DIRFORC -o conversionlonlat_parallel -i $indir -w $outdir -v $varn -l 50 -p

