#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Compare E2OFD to others precipitation data
Created on Thu Apr 19 15:09:41 2018

@author: Anthony Schrapffer
@institution: CIMA, Universidad de Buenos Aires (Argentina)
"""
import sys
import os
sys.path.append(os.getcwd())
import ORCHIDEE_tools as OR
import sys 
import numpy as np
from datetime import date
from netCDF4 import Dataset as NetCDFFile
from netCDF4 import num2date
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.patches as mpatches
import numpy.ma as ma
import matplotlib.lines as mlines
from mpl_toolkits.basemap import Basemap
from matplotlib.ticker import AutoMinorLocator
import matplotlib.ticker as ticker
#rc('text', usetex=True)
import matplotlib.dates as mdates
import matplotlib.cbook as cbook



# Get data for upstream area of station
def varmask_stn_monthobs(stname, chem_file, chem_grdc_rd, chem_grdc, chem_grid, namegr, variablename, timename, y1, y2): 
    """
    Get the variable mean on the upstream area of the station.
    
    stname : str, name of the station.
    chem_file : str, direction of the data file.
    chem_grdc_rd : str, direction of GRDC_river_desc corresponding file.
    chem_grdc : str, direction of the GRDC file.
    chem_grid : str, direction of the cell area file corresponding with chem_file & grdc_rd.
    namegr : list, list of station name in the GRDC file.
    variablename : str, name of the studied variable in the file.
    timename : str, name of the time variable in the file
    y1, y2 : int, year of beginning and end (included)
    """
    # Get index
    index = OR.stgrdcindex(stname,namegr) # pour grdc_rd
    
    # Get mask and time
    mask = OR.getstn_grdc_rd(chem_grdc_rd, index)
    if mask is None:
        return None, None
    mask = ma.array(mask)
    dtime = OR.importTIME(chem_file, timename)

    # Get upstream area
    gridarea = OR.importgridarea(chem_grid)
    
    # Get time series of integrated value on upstream area
    M=ma.zeros((len(dtime)))
    i=0
    print "Dealing with subbasin integration : "
    while i<len(M):
        if i==(len(M)/4): print "25%"
        if i==(len(M)/2): print "50%"
        if i==(len(M)*3/4): print "75%"
        vari=OR.importTIMEvalue(chem_file, variablename, i)
        # mm/d = kg/m^2/d = m^3/1000/m^2/d = m/1000/s/86400 | * upstream en m^2
        stot = ma.sum(mask[:,:]*gridarea[:,:]) 
        M[i] = ma.sum(vari[:,:]*mask[:,:]*gridarea[:,:])/stot
        i=i+1
    return M, dtime


# Cut data for wanted period
def timcut(var, dtime, y1, y2):
    """
    Get the index of beginning and end to cut the dataon this period.
    
    var : list, data to cut on time period.
    dtime : list, time corresponding to data.
    y1, y2 : int, year of beginning and end (included)
    """
    ibeg = OR.finddatemonth(dtime,1,y1)
    iend = OR.finddatemonth(dtime,1,y2)
    varnew = var[ibeg:iend+1] # if (time lat lon - to check)
    dtimenew = dtime[ibeg:iend+1]
    return varnew, dtimenew




# Get annual cycle    
def annualcycle(M, dtime, y1, y2):
    """
    Get the annual cycle.
    
    M : array, data from which we want annual cycle.
    dtime : list, time corresponding to data.
    y1, y2 : int, year of beginning and end (included)
    """
    varnew, dtimenew = timcut(M, dtime, y1, y2)
    years = dtimenew[len(dtimenew)-1].year-dtimenew[0].year+1
    print years    
    # reorganize data
    M0=ma.zeros((12,years))
    i=0
    while i<len(varnew):
        M0[dtimenew[i].month-1,dtimenew[i].year-y1]=varnew[i]
        i=i+1
    # Calculate Annual cycle
    #print np.shape(M0)
    M1=ma.mean(M0,1) # 1 ou zero??
    return M1


    
# Get the data from obs & convert to annual cycle for station
def getdata_plotstn_obs_annualcycle(stname, chem_file, chem_grdc_rd, chem_grdc, chem_grid, variablename, timename, y1, y2):
    """
    Get & convert the data from obs to annual cycle from station upstream area.
    
    stname: str, name of the station.
    chem_file : str, direction of the data file.
    chem_grdc_rd : str, direction of GRDC_river_desc corresponding file.
    chem_grdc : str, direction of the GRDC file.
    chem_grid : str, direction of the cell area file corresponding with chem_file & grdc_rd.
    namegr : list, list of station name in the GRDC file.
    variablename : str, name of the studied variable in the file.
    timename : str, name of the time variable in the file
    y1, y2 : int, year of beginning and end (included)
    """
    namegr=OR.importGRDCname(chem_grdc)
    M, dtime = varmask_stn_monthobs(stname, chem_file, chem_grdc_rd, chem_grdc, chem_grid, namegr, variablename, timename, y1, y2)
    M1 = annualcycle(M, dtime, y1, y2)
    return M1



# Plot the obs annual cycle for station area
def plotstn_obs_annualcycle(stname, L, chem_grdc_rd, chem_grdc, chem_grid, dgraphs, y1, y2, style, basin):
    """
    Plot the annual cycle.
    
    stname: str, name of the station.
    L : array, array of details on the data file.
    chem_grdc_rd : str, direction of GRDC_river_desc corresponding file.
    chem_grdc : str, direction of the GRDC file.
    chem_grid : str, direction of the cell area file corresponding with chem_file & grdc_rd.
    dgraphs : str, direction where to save the result
    y1, y2 : int, year of beginning and end (included).
    style : array, style array for visual details.
    basin : str, name of the basin.
    
    L=[La, Lb]
    La = [chem_file (dir. file), "GPCC" (name), "precip" (name variable), "time" (name time), 
              k(coeff miltiplicatif - unite), True (True ou direction grdc rd si diff), chem_grid (autre)]
    
    style = [s1,s2...]
    s1 = ["-" (line style) ,"" (marker),"k"(color),0.6 (line size)]
    """
    # Extract the data
    OBS = np.zeros((12,len(L)))
    i=0
    while i<len(L):
        print L[i][1]
        if L[i][5] is True:
            OBS[:,i] = getdata_plotstn_obs_annualcycle(stname, L[i][0], chem_grdc_rd, chem_grdc, chem_grid, L[i][2], L[i][3], y1, y2)
        else:
            OBS[:,i] = getdata_plotstn_obs_annualcycle(stname, L[i][0], L[i][5], chem_grdc, L[i][6], L[i][2], L[i][3], y1, y2)
        # relation entre i et variablename et timename
        i=i+1
        
    # Prepare annual plot
    X=np.arange(1,13,1)
    LabMonths=["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","Jan"]
    
    ### Plot ###
    print "Start plotting"
    fig=plt.figure(figsize=(4.5,2.5),dpi=250)
    ax1 = plt.subplot2grid((1, 10), (0, 0), colspan=7)  
    
    # Create Legend
    LEG=[]
    i=0
    while i<len(L):
        LEG.append(mlines.Line2D([], [], color=style[i][2], marker=style[i][1],label=L[i][1],ls=style[i][0],ms=4))
        i=i+1
    # Plot data
    i=0
    a=0
    while i<len(L):
        print OBS[:,i]/L[i][4]
        if ma.max(OBS[:,i]/L[i][4])>a: a = ma.max(OBS[:,i]/L[i][4])
        ax1.plot(X, OBS[:,i]/L[i][4], color = style[i][2] , marker = style[i][1],ls=style[i][0], ms=2,lw=0.5) 
        # voir si /31 - mm/month - mm/day
        i=i+1

    print a
    plt.ylim( 0, a*1.2)
    ax1.set_ylabel('($mm/day$)',fontsize=9,labelpad=3,rotation=90)
    plt.setp(ax1.get_yticklabels(), fontsize=4)

    ax1.set_xticks(X)
    ax1.set_xticklabels(LabMonths, fontsize=6, rotation=-45)
    ax1.tick_params(axis='y', which='major',pad=1.0,labelsize=6)    
    
    print "Plot map"
    OR.addcardgrdcnew(stname, chem_grdc, basin, chem_grdc_rd, False)
    
    lg = ax1.legend(bbox_to_anchor=(1.05, 0.6, 0.2, 0.4),handles=LEG,fontsize=4,title=r'Legend',loc = 2)
    lg.get_title().set_fontsize(5)
    # Finalize    
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.1, top=0.93,wspace= 0.)
    
    fig.suptitle(r'Annual cycle '+stname.replace("\xd6","o") +" "+str(y1)+"-"+str(y2), fontsize=8,y=0.985, ha="left", x=0.1)
    fig.savefig(dgraphs+stname.replace(" ","-").replace("/","-").replace("\xd6","o")+"-annualcycle_OBS-"+str(y1)+str(y2)+".png",dpi=350)
    plt.close()
    return OBS



### Plot annual cycle for observation for list of stations ###
def plotallstn_obs_annualcycle(Lst, L, chem_grdc_rd, chem_grdc, chem_grid, dgraphs, y1, y2, style, basin):
    """
    Plot the annual cycle for a list of stations.
    
    Lst: list, list of stations' names.
    L : array, array of details on the data file.
    chem_grdc_rd : str, direction of GRDC_river_desc corresponding file.
    chem_grdc : str, direction of the GRDC file.
    chem_grid : str, direction of the cell area file corresponding with chem_file & grdc_rd.
    dgraphs : str, direction where to save the result
    y1, y2 : int, year of beginning and end (included).
    style : array, style array for visual details.
    basin : str, name of the basin.
    
    L=[La, Lb]
    La = [chem_file (dir. file), "GPCC" (name), "precip" (name variable), "time" (name time), True (True ou direction grdc rd si diff)]
    
    style = [s1,s2...]
    s1 = ["-" (line style) ,"" (marker),"k"(color),0.6 (line size)]
    """
    i=0
    while i<len(Lst):
        print "####"
        print i+1,"/",len(Lst), " :",Lst[i]
        plotstn_obs_annualcycle(Lst[i], L, chem_grdc_rd, chem_grdc, chem_grid, dgraphs, y1, y2, style, basin)
        i=i+1
    print "Finished"
    
    
### Plot annual cycle for observation for basins stations ###    
def plotallbasin_obs_annualcycle(L, chem_grdc_rd, chem_grdc, chem_grid, dgraphs, y1, y2, style, basin):
    """
    Plot the annual cycle for all the available stations of a basin.
    
    L : array, array of details on the data file.
    chem_grdc_rd : str, direction of GRDC_river_desc corresponding file.
    chem_grdc : str, direction of the GRDC file.
    chem_grid : str, direction of the cell area file corresponding with chem_file & grdc_rd.
    dgraphs : str, direction where to save the result
    y1, y2 : int, year of beginning and end (included).
    style : array, style array for visual details.
    basin : str, name of the basin.
    
    L=[La, Lb]
    La = [chem_file (dir. file), "GPCC" (name), "precip" (name variable), "time" (name time), True (True ou direction grdc rd si diff)]
    
    style = [s1,s2...]
    s1 = ["-" (line style) ,"" (marker),"k"(color),0.6 (line size)]
    """
    print "###",basin,"###"
    doc=open(dgraphs+basin+"stn.txt","w")
    doc.write("### "+basin+" ###")
    doc.close()
    Lavst = OR.AvailableStn(chem_grdc_rd, chem_grdc, basin, AR=False, BR=False)
    
    if len(Lavst)==0: 
        print "No Available Station"
        return
    Lst=[]
    i=0
    while i<len(Lavst):
        Lst.append(Lavst[i][0])
        i=i+1
    plotallstn_obs_annualcycle(Lst, L, chem_grdc_rd, chem_grdc, chem_grid, dgraphs, y1, y2, style, basin)
    return







# Plot the obs timeserie for station area
def plotstn_obs_timeserie(stname, L, chem_grdc_rd, chem_grdc, chem_grid, dgraphs, y1, y2, style, basin):
    """
    Plot the annual cycle.
    
    stname: str, name of the station.
    L : array, array of details on the data file.
    chem_grdc_rd : str, direction of GRDC_river_desc corresponding file.
    chem_grdc : str, direction of the GRDC file.
    chem_grid : str, direction of the cell area file corresponding with chem_file & grdc_rd.
    dgraphs : str, direction where to save the result
    y1, y2 : int, year of beginning and end (included).
    style : array, style array for visual details.
    basin : str, name of the basin.
    
    L=[La, Lb]
    La = [chem_file (dir. file), "GPCC" (name), "precip" (name variable), "time" (name time), 
              k(coeff miltiplicatif - unite), True (True ou direction grdc rd si diff), chem_grid (autre)]
    
    style = [s1,s2...]
    s1 = ["-" (line style) ,"" (marker),"k"(color),0.6 (line size)]
    """
    # Extract the data
    OBS = np.zeros((y2-y1+1)*12,len(L)))
    i=0
    while i<len(L):
        print L[i][1]
        if L[i][5] is True:
            OBS[:,i] = getdata_plotstn_obs_timeserie(stname, L[i][0], chem_grdc_rd, chem_grdc, chem_grid, L[i][2], L[i][3], y1, y2)
        else:
            OBS[:,i] = getdata_plotstn_obs_timeserie(stname, L[i][0], L[i][5], chem_grdc, L[i][6], L[i][2], L[i][3], y1, y2)
        # relation entre i et variablename et timename
        i=i+1
    
    ### Plot ###
    print "Start plotting"
    fig=plt.figure(figsize=(4.5,2.5),dpi=250)
    ax1 = plt.subplot2grid((1, 10), (0, 0), colspan=7)  
    
    # Create Legend
    LEG=[]
    i=0
    while i<len(L):
        LEG.append(mlines.Line2D([], [], color=style[i][2], marker=style[i][1],label=L[i][1],ls=style[i][0],ms=4))
        i=i+1
    # Plot data
    i=0
    a=0
    while i<len(L):
        print OBS[:,i]/float(L[i][4])
        if ma.max(OBS[:,i]/float(L[i][4]))>a: a = ma.max(OBS[:,i]/float(L[i][4]))
        ax1.plot(X, OBS[:,i]/float(L[i][4]), color = style[i][2] , marker = style[i][1],ls=style[i][0], ms=2,lw=0.5) 
        # voir si /31 - mm/month - mm/day
        i=i+1

    print a
    plt.ylim( 0, a*1.2)
    ax1.set_ylabel('($mm/day$)',fontsize=6,labelpad=3,rotation=90)
    plt.setp(ax1.get_yticklabels(), fontsize=4)

    # xtick
    xtickstimeMonth(y1, y2 , ax1)   
    
    print "Plot map"
    OR.addcardgrdcnew(stname, chem_grdc, basin, chem_grdc_rd, False)
    
    lg = ax1.legend(bbox_to_anchor=(1.05, 0.6, 0.2, 0.4),handles=LEG,fontsize=4,title=r'Legend',loc = 2)
    lg.get_title().set_fontsize(5)
    # Finalize    
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.1, top=0.93,wspace= 0.)
    
    fig.suptitle(r'Time Serie '+stname.replace("\xd6","o") +" "+str(y1)+"-"+str(y2), fontsize=8,y=0.985)#loc="left"
    fig.savefig(dgraphs+stname.replace(" ","-").replace("/","-").replace("\xd6","o")+"-timeserie_OBS-"+str(y1)+str(y2)+".png",dpi=350)
    plt.close()
    return 


### Plot annual cycle for observation for list of stations ###
def plotallstn_obs_timeserie(Lst, L, chem_grdc_rd, chem_grdc, chem_grid, dgraphs, y1, y2, style, basin):
    """
    Plot the annual cycle for a list of stations.
    
    Lst: list, list of stations' names.
    L : array, array of details on the data file.
    chem_grdc_rd : str, direction of GRDC_river_desc corresponding file.
    chem_grdc : str, direction of the GRDC file.
    chem_grid : str, direction of the cell area file corresponding with chem_file & grdc_rd.
    dgraphs : str, direction where to save the result
    y1, y2 : int, year of beginning and end (included).
    style : array, style array for visual details.
    basin : str, name of the basin.
    
    L=[La, Lb]
    La = [chem_file (dir. file), "GPCC" (name), "precip" (name variable), "time" (name time), True (True ou direction grdc rd si diff)]
    
    style = [s1,s2...]
    s1 = ["-" (line style) ,"" (marker),"k"(color),0.6 (line size)]
    """
    i=0
    while i<len(Lst):
        print "####"
        print i+1,"/",len(Lst), " :",Lst[i]
        plotstn_obs_timeserie(Lst[i], L, chem_grdc_rd, chem_grdc, chem_grid, dgraphs, y1, y2, style, basin)
        i=i+1
    print "Finished"

    
### Plot annual cycle for observation for basins stations ###    
def plotallbasin_obs_timeserie(L, chem_grdc_rd, chem_grdc, chem_grid, dgraphs, y1, y2, style, basin):
    """
    Plot the annual cycle for all the available stations of a basin.
    
    L : array, array of details on the data file.
    chem_grdc_rd : str, direction of GRDC_river_desc corresponding file.
    chem_grdc : str, direction of the GRDC file.
    chem_grid : str, direction of the cell area file corresponding with chem_file & grdc_rd.
    dgraphs : str, direction where to save the result
    y1, y2 : int, year of beginning and end (included).
    style : array, style array for visual details.
    basin : str, name of the basin.
    
    L=[La, Lb]
    La = [chem_file (dir. file), "GPCC" (name), "precip" (name variable), "time" (name time), True (True ou direction grdc rd si diff)]
    
    style = [s1,s2...]
    s1 = ["-" (line style) ,"" (marker),"k"(color),0.6 (line size)]
    """
    print "###",basin,"###"
    doc=open(dgraphs+basin+"stn.txt","w")
    doc.write("### "+basin+" ###")
    doc.close()
    Lavst = OR.AvailableStn(chem_grdc_rd, chem_grdc, basin, AR=False, BR=False)
    
    if len(Lavst)==0: 
        print "No Available Station"
        return
    Lst=[]
    i=0
    while i<len(Lavst):
        Lst.append(Lavst[i][0])
        i=i+1
    plotallstn_obs_timeserie(Lst, L, chem_grdc_rd, chem_grdc, chem_grid, dgraphs, y1, y2, style, basin)
    return



"""
################
#### Prueba #### 

style=[["-","","k",0.6],
       ["-","x","r",0.6],
       ["--","o","#dc6900",0.4],
       ["-.","","b",0.4],
       ["--","","g",0.6]]

d = "/home/anthony/Documents/Doctorat/ORCHIDEE/8-Comparison_Resolutions/8.3-FullBasinSAM/Originals/"
chem_grdc_rd = d +"grdc_river_desc.nc"
chem_grdc = "/home/anthony/Documents/Doctorat/ORCHIDEE/Z0-Originals/GRDC_Monthly_June14.nc"
chem_grid = d + "SAM_E2OFD_S0_orchidee_cellarea.nc"
dgraphs = "/home/anthony/Documents/Doctorat/ORCHIDEE/8-Comparison_Resolutions/8.4-CompPR-E2OFD/"
y1= 1980
y2= 1990
basin = "Parana"

chem_file= "/home/anthony/Documents/Doctorat/ORCHIDEE/8-Comparison_Resolutions/8.4-CompPR-E2OFD/precip.1x1_GPCC.mon.meanout2.nc"
La=[chem_file, "GPCC", "precip", "time", True]
# [direction, name, variablename, timename, k-coeff multiplication =1 daily =31 monthly, 3CN - altgrdc (True ou direction grdc rd), chemgrid si alt]
# Erreur creer par monthly, besoin de [31,28.25,31,30,31,30,31,31,30,31,30,31]
L=[La]

M0, M1= plotstn_obs_annualcycle("Corrientes", L, chem_grdc_rd, chem_grdc, chem_grid, dgraphs, y1, y2, style, basin)


Convertir to monthly 
GPCC : 1º, precip(time, lat, lon) en mm/mn, time // (X.5º)
GPCP : 0.5º, precip(time, lat, lon) en mm/day, time // 
3CN : 0.5 º, rr(time, latitude, longitude) (mm), time, 
CRU : 0.5º, pre(time, lat, lon) mm/month
CPC_UNI : 0.5º, precip(time, level001, latitude, longitude) mm/day (meme si mensuel)
E2OFD : 0.25º, rain, a mettre sur le serveur (demander si possible) daily mm/d
"""


#### TO DO ####

# Carte avec precipitation (toutes)
# Carte avec biases
# Other input file cut on smallest
