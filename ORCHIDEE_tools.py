#!/usr/bin/env python2
# -*- coding: iso-8859-15 -*-
# Tools for Python / ORCHIDEE
# From Anthony Schrapffer
# Working at CIMA (Argentina)
#
# This work is licendes under a Creative Commons 
#   Attribution-ShareAlike 4.0 International License (http://creativecommons.org/licenses/by-sa/4.0)
#
# To import : 
# import sys
# sys.path.append("/home/anthony/TotiTools") # To change
#
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

#################
### Variables ###
#################
error="ERROR - error - ERROR"
### For version 0.5° ###
# List of name - validated station on La Plata basin
Lstold=["Fecho Dos Morros","La Punilla","Andira","Barbosa Ferraz","Zanja Del Tigre",
    "Balsa Santa Maria","Salto Carlos Botelho","Porto Felicio (Jus. Us. Jaguara)","Santa Cruz Do Timbo","Villa Montes","Jupia - Jusante","Miraflores",
    "Porto Do Alegre","Caceres (Dnpvn)","Novo Porto Taquara","El Tunal",
    "Paso Lucero","Porto Murtinho (Fb/Dnos)","Salto Del Guaira","Salto Osorio Jusante","Fazenda Santa Fe","La Paz",
    "Anta Muerta","Abaixo Barra Do Rio Verde","Cabra Corral (1967: La Puerta)","Corrientes","Posadas",
    "Aguas Do Vere","Algarrobito (1971: San Telmo)","Jataizinho",
    "Porto Esperanca (Dnos)","Uniao Da Vitoria","Porto Paraiso Do Norte",
    "Acima Do Corrego Grande","Timbues"]

# List of correspondant Longitud in the model grid
LonMod=[-57.75,-66.25,-50.25,-52.25,-64.25,-53.25,-51.25,-47.75,-50.75,-63.75,-51.75,
        -65.25,-56.75,-57.75,-52.75,-64.75,-58.25,-57.75,-54.25,-52.75,-50.25,-62.75,
        -64.75,-50.25,-65.75,-58.75,-56.25,-52.75,-64.75,-50.75,-57.25,-50.75,-52.25,
        -54.75,-60.75]

# List of correspondant Latitud in the model grid
LatMod=[-21.25,-25.75,-23.25,-24.25,-22.75,-24.25,-21.25,-20.25,-26.75,-21.25,-21.25,
        -25.25,-17.25,-16.25,-23.25,-25.25,-28.75,-21.25,-23.75,-25.75,-18.75,-22.25,
        -23.25,-17.75,-25.25,-28.25,-27.25,-26.25,-22.25,-23.75,-19.75,-26.25,-23.25,
        -16.25,-32.25]

# Style to plot
style=[["-","","k",0.6],
       ["-","x","r",0.6],
       ["--","o","#dc6900",0.6],
       ["-.","","b",0.4],
       ["--","","g",0.6]]

# List of basins with geographical limits to plot
Basins=np.array([["Parana",-67.,-36.,-40.,-15.],["Amazon",-1,-1,-1,-1]])

### EXPLIQUER DEFINITION FILE TYPE ETC
### Plus pour vieux modele selection station et longitud latitude correspondante 
###  L2, LonMod, LatMod
#DIR, filetype (decrire), NAME)
#La = [dir_E2OFD_GRDC, "GRDCnew", "GRDC"]
#L1=[La, Lb, Lc, Lf]
### WHAT YOU NEED TO DEFINE IN MAIN
#style=[["-","","k",0.6],
#       ["-","x","r",0.6],
#       ["--","o","#dc6900",0.6],
#       ["-.","","b",0.4],
#       ["--","","g",0.6]]
### DGRAPHS plus redefinir ds annual cycle 
#DIR GRDC original file
# DIR grdc rd file 
# Dir cell area grid file made with CDO


### ADD FILE WITH BDHI EXTRACTION
### VOIR COMMENT FAIRE POUR QUE CE FICHIER SE RETROUVE AUTOMATIQUEMENT
### CF QUIL LA CONSIDERE COMME UNE LIBRAIRIE !!!
### METTRE SUR GIT YOLO
### AJOUTER FICHIRE P FAIRE TOURNER ORCHIDEE ET FICHIER EXTRACTION CONVERSION DEPUIS ARGENTINE ET CICLAD (PARTIE GRDC)



###############
### Content ###
###############
### Importation ###
# importGRDCname: import list of stations name from GRDC file
# importOLDSIM: import hydrograph / lon / lat / time from ORCHIDEE
# importGRDCnew: import GRDC hydrograph / time / index correspondance from ORCHIDEE output for GRDC stations
#                         Fuxing GRDC Module
# importNEWSIM: import hydrograph / time / index correspondance for GRDC stations from ORCHIDEE output
# importTIME: import time on date format from a NetCDF file 
# importTIMEvalue: import the 2D value of a NetCDF file for a specific timestep
# importvariable: import a 1D/2D/3D variable from a netCDF file
# import gridarea: import the cell area array from a netCDF file - made with cdo cellarea
#
### Preparation ###
#



########################################################################
############################## FUNCTIONS ###############################
########################################################################
###################
### Importation ###
###################
### Import list of stations name from GRDC file
def importGRDCname(chem_GRDC):
    GRDC = NetCDFFile(chem_GRDC, 'r')
    timv=GRDC.variables['time']
    time=timv[:]
    #dtimegr = num2date(time,timv.units)
    #hydrog = GRDC.variables['hydrographs'][:,:] #( time,stations) Original Hydrographs
    Name = GRDC.variables['name'][:]
    return Name #Seul util pour le moment, data du new
### Import hydrograph / lon / lat / time from ORCHIDEE
def importOLDSIM(chem_orchold):
    salid = NetCDFFile(chem_orchold, 'r')
    timv1=salid.variables['time_counter']
    #old sim from hydra ?
    #timv1=salid.variables['time']
    time1=timv1[:]  
    dtime1 = num2date(time1,timv1.units)
    lon1= salid.variables['lon'][:]
    lat1= salid.variables['lat'][:]
    hydro1 = salid.variables['hydrographs'][:,:,:] # (time, lat, lon)
    return  dtime1, lon1, lat1, hydro1
### Import GRDC hydrograph / time / index correspondance from ORCHIDEE output for GRDC stations Fuxing GRDC Module    
def importGRDCnew(chem_GRDC):
    salid = NetCDFFile(chem_GRDC, 'r')
    timv2=salid.variables['time_counter']
    time2=timv2[:]
    dtime2 = num2date(time2,timv2.units)
    hydro2 = salid.variables['Dis_Stn_GRDC'][:,:] # (Time , cell domain - index station)
    index2=salid.variables['Index_Stn_GRDC']
    hydro2
    return dtime2, hydro2, index2 
### Import hydrograph / time / index correspondance for GRDC stations from ORCHIDEE output
def importNEWSIM(chem_orchnew):
    salid = NetCDFFile(chem_orchnew, 'r')
    timv3=salid.variables['time_counter']
    time3=timv3[:]
    dtime3 = num2date(time3,timv3.units)
    hydro3 = salid.variables['Dis_Stn_Model'][:,:] # (Time , cell domain - index station)
    index3=salid.variables["Index_Stn_GRDC"][:]
    hydro3
    return dtime3, hydro3, index3
### Import time on date format from a NetCDF file
def importTIME(chem_file, variable):
    """
    Importation of a 3D variable from NetCDF.
    chem_file: string, direction of the file.
    variable: string, name of the variable.
    """
    salid = NetCDFFile(chem_file, 'r')
    timv=salid.variables[variable]
    time=timv[:]
    dtime = num2date(time,timv.units)
    return dtime
### Import the 2D value of a NetCDF file for a specific timestep
def importTIMEvalue(chem_file, variable, i):
    salid = NetCDFFile(chem_file, 'r')
    vari=salid.variables[variable][i,:,:] # time, lat, lon
    return vari
### Import a 1D/2D/3D variable from a netCDF file
def importvariable(chem_file, varname, dim):
    salid = NetCDFFile(chem_file, 'r')
    if dim == 1:
        var = salid.variables[varname][:]
    if dim == 2:
        var = salid.variables[varname][:,:]
    if dim == 3:
        var = salid.variables[varname][:,:,:]
    return var
### Import the cell area array (made from cdo cellarea)
def importgridarea(chem_file):
    gridarea = importvariable(chem_file, "cell_area", 2) # (lat,lon) in km^2
    return gridarea

#############
### Tools ###
#############
### TIME###
### Get index of the beginning a a year in a time array
def datebeg(dtime,y):
    """
    Get the index of the begining of the year y in the dtime file 
    dtime: netcdf time list convert with num2date
    y: int, year
    """
    i=0
    while dtime[i].year != y:
        if i==len(dtime)-1:
            print "Year not included"
            return None
        else:
            i=i+1
    if dtime[i].month!=1:
        print error
    return i
### Get the index of the end of a year in a time array
def dateend(dtime,y):
    """
    Get the index of the end of the year y in the dtime file 
    dtime: netcdf time list convert with num2date
    y: int, year
    """
    i=0
    while (dtime[i].year != y or dtime[i].month!= 12 or dtime[i].day!=31):
        if i==len(dtime)-1:
            print "Year not included"
            break
        else:
            i=i+1
    return i
### Find the index of a specific date in a time array
def finddate(dtime,d,m,y):
    """
    Get the index of an exact date of the year.
    dtime: netcdf time list convert with num2date.
    d: int, day.
    m: int, month.
    y: int, year.
    """
    i=0
    while (dtime[i].year != y or dtime[i].month != m or dtime[i].day != d) and i<len(dtime):
        i=i+1
    return i
### Find the index of a specific month (case of monthly data)
def finddatemonth(dtime,m,y):
    """
    Get the index of an exact date of the year.
    dtime: netcdf time list convert with num2date.
    m: int, month.
    y: int, year.
    """
    i=0
    while (dtime[i].year != y or dtime[i].month != m) and i<len(dtime):
        i=i+1
    return i
### Get the monthly mean data from daily data between y1 and y2 included
def monthmeantot(data, dtime, y1, y2):
    M=ma.zeros((y2-y1+1)*12)
    y=y1
    i=0
    while y<y2+1: 
        m=0
        while m<12:
            ii=finddate(dtime,1,m+1,y)
            d1=date(y,m+1,1)
            if m==11:
                m2=0; yn=y+1
            else:
                m2=m+1; yn=y
            d2=date(yn,m2+1,1)
            monlen=(d2-d1).days
            M[i]=ma.mean(data[ii:(ii+monlen)])
            i=i+1
            m=m+1
        y=y+1
    return M
### Get the monthly mean for a specific month

def monthmean(H, dtime, mon, y):
    """
    Return de mean of H for the month mon of y year, indexed with dtime
    H: list, list of value from which we want to extract the mean value.
    dtime: netcdf time list convert with num2date.
    mon: int, month. (indexed 1-12) ATTENTION !!!
    y: int, year.
    """
    debug=False
    i=finddate(dtime,1,mon,y)
    if debug: print i 
    if debug: print "lendtime" ,len(dtime)," lenH ",len(H)
    d1=date(y,mon,1)
    if debug: print "d1",d1
    # Get date of next month
    if mon==12: # if we have to change of year
        m2=1
        y2=y+1
    else:
        m2=mon+1
        y2=y
    d2=date(y2,m2,1)
    if debug: print "d2",d2
    monlen = (d2-d1).days
    monmean = ma.mean(H[i:(i+monlen+1)])        
    return monmean
### Others ###
### Get the index of a specific longitud latitud in a lon/lat file
def lonlatij(lon,lat,ilon,ilat):
    """
    Get the index of the latitude & longitude in the netCDF file
    lon: longitud array of the netCDF file
    lat: latitud array of the netCDF file
    ilon: longitud we are interested in 
    ilat: latitude we are interested in
    !!!CAUTION!!! ilon and ilat must be part of the grid
    """
    nlon=np.where(lon==ilon)[0][0]
    nlat=np.where(lat==ilat)[0][0]
    return nlon, nlat
### Get the GRDC Monthly data between y1 and y2 (included)
def getGRDCdata(stn, chem_GRDC, y1, y2):
    GRDC = NetCDFFile(chem_GRDC, 'r')
    Name = GRDC.variables['name'][:]
    timv = GRDC.variables['time']
    time = timv[:]
    dtimegr = num2date(time,timv.units)
    i = stgrdcindex(stn,Name)-1 #Because indexation python from 0 
    hydrographs = GRDC.variables['hydrographs'][:,i]
    
    t1 = datebeg(dtimegr,y1)
    t2 = datebeg(dtimegr,y2)
    print Name[i]
    H=hydrographs[t1:t2+12]
    return H, dtimegr[t1:t2+12]
#Convert GRDC name
def convertgrdcname(Namelist):
    i=len(Namelist)-1
    while Namelist[i]==" ":
        i=i-1
    Name=""
    j=0
    while j<=i:
        Name=Name+Namelist[j]
        j=j+1
    return Name

###################
### Preparation ###
###################

### Get the index of a station on the original GRDC file
def stgrdcindex(stn,namegr):
    """ Find the index number of a station from its name in the GRDC output.
    (for finding it in the new orchidee version - Fuxing module - GRDC)
    stname, string : name of the station.
    namgr, string : name of the GRDC file, from which was made the simulation !!!!!
    """
    lenst=len(stn)
    i=0
    while i<len(namegr):
       stname = ""
       for ic in range(len(namegr[i,:])):
           stname  = stname + namegr[i,ic]
       if stname[0:lenst]==stn: # Name from GRDC = stn (station wanted)
           return i+1 # Because index is from 1-190, python starts from 0 not 1
       i=i+1
    print error
    print "Not present in GRDC file."
### Get the index of a station on the GRDC ORCHIDEE output (from Fuxing module)
def stoutputindex(stn, namegr, index):
    """
    Get the GRDC new index for the station stn.
    stname: string, name of the station.
    namgr: array, list of name from the GRDC original file, from which was made the simulation !!!!! 
    index: array, list of index from the GRDC output file.
    """
    ind=stgrdcindex(stn, namegr)
    i=0
    while index[i]!=ind and i<len(index)-1:
        i=i+1
    if index[i]==ind:
        return i
    else:
        print error
        print "Not an output referenced station"
        return None
### Get the mask of the subbasin of a station ###
def getstn_grdc_rd(chem_grdc_rd, index):
    """
    Get the the variable in the grdc_river_desc.nc that correspond
    to the index number.
    Return None if not present ! 
    chem_grdc_rd: string, direction of the grdc_river_desc file.
    index: int, index of the station, in relation to original GRDC file.
    /!\ index start at 1 !!
    """
    salid = NetCDFFile(chem_grdc_rd, 'r')
    for m in salid.variables.keys():
        if salid.variables[m].Index_of_GRDC_Station==index:
            a=salid.variables[m][:,:] # (lat, lon)
            return a
    return None
### Get the location of a station in the model ###
def get_lonlat_stn_model(chem_GRDC, chem_GRDCnew, stname):
    """
    Get the corresponding longitud and latitud for the station stname in the model.
    """
    Namegr = importGRDCname(chem_GRDC)
    index = importvariable(chem_GRDCnew, 'Index_Stn_GRDC',1)
    ind = stoutputindex(stname, Namegr, index) 
    
    salid = NetCDFFile(chem_GRDCnew, 'r')
    lon=salid.variables["Lon_Stn_Model"][ind]
    lat=salid.variables["Lat_Stn_Model"][ind]
    return lon, lat


#####################
### Data Treatment###
#####################
### Get variable integrated in subbasin of a station - convert from mm/d to m^3/s - FOR PRECIPITATION AND EVAP
def get_varmask_stn(stname, chem_file, chem_grdc_rd, chem_grdc, chem_grid, namegr, variable): #ADD GRDC FILE !!!!! 
    Debug = False
    if Debug: print "start"
    # Get index
    index = stgrdcindex(stname,namegr) # pour grdc_rd
    
    if Debug: print "Indexes ok"
    # Get mask and time
    mask = getstn_grdc_rd(chem_grdc_rd, index)
    if mask is None:
        return None, None
    mask = ma.array(mask)
    dtime = importTIME(chem_file, "time_counter")

    if Debug: print "Grid area"
    # Get upstream area
    gridarea = importgridarea(chem_grid)
    
    if Debug: print "Get time serie"
    # Get time series of integrated value on upstream area
    M=ma.zeros((len(dtime)))
    i=0
    print "Dealing with subbasin integration : "
    while i<len(M):
        if i==(len(M)/4): print "25%"
        if i==(len(M)/2): print "50%"
        if i==(len(M)*3/4): print "75%"
        vari=importTIMEvalue(chem_file, variable, i)
        # mm/d = kg/m^2/d = m^3/1000/m^2/d = m/1000/s/86400 | * upstream en m^2
        M[i] = ma.sum(vari[:,:]*mask[:,:]*gridarea)/1000/86400 
        i=i+1
    return M, dtime
### Extract data from a station from GRDC output, rain/evap data, ORCHIDEE hydrographs (needs correspondance of lon/lat)
def extract_stn(stname, chem_file, filetype, namegr="", 
                chem_grid="", chem_grdc="", chem_grdc_rd=""): #Last part for rain/evap
    """
    Extract the hydrographs data for a specific station.
    stname: string, name of the station.
    var: array, array of hydrographs data
    filetype: string, values = "old","new","GRDCnew" // possible rajout ancien mais inutile
    """
    if filetype == "old": # for 0.5° simulation
        dtime1, lon1, lat1, hydro1 = importOLDSIM(chem_file)
        Ls=np.array(Lstold)
        i=np.where(Ls==stname)[0][0] # Give indexation Python
        nlon,nlat = lonlatij(lon1,lat1,LonMod[i],LatMod[i])
        hydrostn = hydro1[:,nlat,nlon]
        return hydrostn, dtime1
    
    elif filetype == "new" or filetype == "GRDCnew": # for GRDC output (FUXING Module)
        if filetype == "new" : dtime0, hydro0, index0 = importNEWSIM(chem_file)
        if filetype == "GRDCnew" : 
            dtime0, hydro0, index0 = importGRDCnew(chem_file)
            hydro0=ma.masked_where(hydro0<0, hydro0)
        ind=stoutputindex(stname, namegr, index0)
        if ind == None:
            return None, None
        hydrostn = hydro0[:,ind]
        return hydrostn, dtime0

    elif filetype == "rain" or filetype == "evap":
        outdata, dtime = get_varmask_stn(stname, chem_file, chem_grdc_rd, chem_grdc, chem_grid, namegr, filetype) 
        return outdata, dtime
    else:
        print error
        print "Filetype is not well defined"
### Extract list of data ! 
def extract_liststn(stname, Li, chem_GRDC, chem_grdcnew="", chem_grid="", chem_grdc_rd=""):
    """
    Construction of hydrographs for station stname
    from element of L - [[filedir (str dir of file), filetype (str "old" "new" "GRDCnew"), name_simu str simulation name]]
    stname: string, name of the station.
    chem_GRDC: string, path to GRDC original file (used in simulation).
    """
    Name=importGRDCname(chem_GRDC)
    li=len(Li)
    outlist=[]
    i=0
    while i<li:
        L0=Li[i]
        if L0[1]=="old":
            data, dtime = extract_stn(stname, L0[0], "old")
            outlist.append([L0[2], dtime, data])  
            
        ### mettre ensemble et mettre Name pour tous re grouper au moins
        elif L0[1]=="new" or L0[1]=="GRDCnew":
            data, dtime = extract_stn(stname, L0[0], L0[1], namegr=Name)
            outlist.append([L0[2], dtime, data])
            
        elif L0[1] == "rain" or L0[1] == "evap":
            data, dtime = extract_stn(stname, L0[0], L0[1], namegr=Name, chem_grid=chem_grid,
                                      chem_grdc=chem_grdcnew, chem_grdc_rd=chem_grdc_rd)
            outlist.append([L0[2], dtime, data])
            
        else:
            print error
            print "Not valid filetype"
            print L0[1]
            return None
        if data is None:
            print "Absent data in the list"
            return None
        i=i+1
    return outlist
### Extract time series from a list of data
def extract_timeseries(stname, L, chem_GRDC, y1, y2, chem_grid="", chem_grdc_rd=""):
    """
    Extract the time serie from a list of data
    """
    debug=False
    li=len(L)
    i=0
    isgrdc=False
    # Detect GRDC
    if debug: print "Detect GRDC"
    while i<li:
        if L[i][1]=="GRDCnew":          
            isgrdc=True
            grind=i # indice de GRDC
            break
        i=i+1
    if debug: print "Start"
    
    if isgrdc: #Cas ou GRDC donc oublier chiffre si absence de donné pour cohérence 
        outlist=extract_liststn(stname, L, chem_GRDC, L[grind][0], chem_grid, chem_grdc_rd)
    
        if outlist is None:
            return None
        H=[]    
        
        # Creation list recoupé
        if debug: print "Cut list"
        i=0
        while i<li:
            if debug: print i, L[i][2]
            A=outlist[i]
            tbeg=datebeg(A[1],y1)
            tend=finddate(A[1],31,12,y2)
            H.append([A[0],A[1][tbeg:(tend+1)],ma.array(A[2][tbeg:(tend+1)])])
            if debug: print A[1][tbeg]
            if debug: print A[1][tend]
            if debug: print len(H[i][2]) # Vérification même taille de liste
            i=i+1
        return H
    else:
        return None
### Extract annual cycle from a list of data
def extract_annualcycle(stname, L, chem_GRDC, y1, y2):
    """
    Extract the annual cycle from a list of data.
    """
    debug=False
    li=len(L)
    i=0
    isgrdc=False
    # Detect GRDC
    if debug: print "Detect GRDC"
    while i<li:
        if L[i][1]=="GRDCnew":
            isgrdc=True
            grind=i # indice de GRDC
            break
        i=i+1
    
    outlist=extract_liststn(stname, L, chem_GRDC)
    if outlist == None:
        return None,None,None,None
    
    if debug: print "Start"
    if isgrdc: #Cas ou GRDC donc oublier chiffre si absence de donné pour cohérence 
        H=[]    
        
        # Creation list recoupé
        if debug: print "Cut list"
        i=0
        while i<li:
            if debug: print i, L[i][2]
            A=outlist[i]
            tbeg=datebeg(A[1],y1)
            tend=finddate(A[1],31,12,y2)
            H.append([A[0],A[1][tbeg:(tend+1)],ma.array(A[2][tbeg:(tend+1)])])
            if debug: print A[1][tbeg]
            if debug: print A[1][tend]
            if debug: print len(H[i][2]) # Vérification même taille de liste
            i=i+1
        K=H
        
        
        # Masked where no data - étape pour presence grdc
        if debug: print "Mask data"
        i=0
        while i<li:
            if i != grind:
                H[i][2]=ma.masked_where(H[grind][2]<0,H[i][2])
            i=i+1
        H[grind][2]=ma.masked_where(H[grind][2]<0,H[grind][2])
        
        # Passer dans un tableau (sim,mois,année) H - M
        if debug: print "Month mean"
        dy = y2-y1+1
        M=ma.zeros([li,12,dy])
        y=0
        while y<dy:
            m=0
            while m<12:
                Sim=0
                while Sim<li:
                    M[Sim,m,y] = monthmean(H[Sim][2], H[Sim][1], m+1, y1+y)
                    Sim=Sim+1
                m=m+1
            y=y+1
               
        # Faire la moyenne sur toutes les années M - T
        if debug: print "Year mean"
        T = ma.mean(M, axis=2)       
        
        simname=[]
        i=0
        while i<li:
            simname.append(L[i][2])
            i=i+1
        
        return T, M, simname, K
        #Ajout list des noms
    else:
        print "NO GRDC"
        return None, None, None, None

#################
### PLOT DATA ###
#################
### TOOLS###
############
### Add map of GRDC station
# reprendre avec retrouver dans modèle location exact cf tools !!!!!! - GRDCnew 
def addcardgrdcnew(stname, chem_GRDC, basin):
    """
    Add the map with the indication of the station next to the graphic.
    k: int, index of the station
    fig: axis of the figure to plot the map
    """
    namegr = importGRDCname(chem_GRDC)
    i = stgrdcindex(stname, namegr)-1 # car index commence à 1

    lon = importvariable(chem_GRDC, "lon", 1)[i]
    lat = importvariable(chem_GRDC, "lat", 1)[i]
    ibas = np.where(Basins == basin)[0][0]
    Lbas= Basins[ibas]
    ax2 = plt.subplot2grid((3, 5), (2, 4),colspan=1)
    m = Basemap(projection="cyl", llcrnrlon=float(Lbas[1]), llcrnrlat=float(Lbas[2]),         \
                    urcrnrlon=float(Lbas[3]), urcrnrlat= float(Lbas[4]), resolution="i")
    m.arcgisimage(server='http://server.arcgisonline.com/ArcGIS', service = 'World_Physical_Map',epsg=3000,xpixels=300, dpi=300,verbose=True)
    m.drawcountries(linewidth=0.25)
    m.drawcoastlines(linewidth=0.25)
    m.drawrivers(linewidth=0.15,color="b")
    
    ax2.plot([lon],[lat],'o',markersize=2,color='r')
    return 
### Plot xticks timeseries
def xtickstimeMonth(y1, y2, ax):
    #T=np.arange(1,len(dtime)+1,1) #index nombre 
    dtime=[]
    y=y1
    while y<y2+1:
        m=1
        while m<13:
            dtime.append(date(y,m,1))
            m=m+1
        y=y+1
    
    Ti=[0]
    Tii=[dtime[0].year]
    k=1
    while Tii[k-1]<y2:
        Ti.append(datebeg(dtime,Tii[k-1]+1))
        Tii.append(dtime[Ti[k]].year)
        k=k+1
    Ti.append(len(dtime)-1)
    Tii.append(y2+1)

    plt.xticks(Ti, Tii, rotation='horizontal',fontsize=6)
    ax.xaxis.set_minor_locator(AutoMinorLocator(12))
    plt.tick_params(which='minor', length=2, color='grey')
    plt.tick_params(axis ='both', which='major', length=4)
    return 

### PLOT GRAPHS ###
###################
### Plot the annual cycle graphs for a station for a list of data
def plot_annualcyclestn(stname, L, chem_GRDC,y1,y2, dgraphs, basin): #style included
    """
    Plot the annual cycle between y1 and y2 for GRDC station stname
    between y1 and y2. stle define the style of the corresponding curves for L.
    chem_GRDC is the original GRDC observation file from which the simulation have been done.
    It's essential to find the index of the stations.
    When data is missing in GRDC, it's not taken in account in model output.
    stname: string, name of the station.
    L: array, output array (each output is define by [dir_output, filetype, simulation name])
    chem_GRDC: string, direction of GRDC file.
    y1: int, beginning year.
    y2: int, ending year.
    style: style for each output, define by : ([linestyle,marker,color])
    """
    debug = None
    print stname    

    X=np.arange(1,13,1)  
    # Get data
    if debug: print "Get data"
    T, M, simname, K = extract_annualcycle(stname, L, chem_GRDC, y1, y2)  
    if T is None:
        return None
    if type(T.mask) == np.ndarray:
        if not False in T.mask: 
            print "No data for the period"
            return None,None,None
    # Legend
    LEG=[]
    i=0
    while i<len(L):
        LEG.append(mlines.Line2D([], [], color=style[i][2], marker=style[i][1],label=simname[i],ls=style[i][0],ms=4))
        i=i+1

    if debug: print "Plot"
    LabMonths=["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","Jan"]
    fig=plt.figure(figsize=(4.5,2.5),dpi=250)
    #ax1=plt.gca()
    ax1 = plt.subplot2grid((1, 5), (0, 0), colspan=4)   
    i=0
    while i<len(simname):
        ax1.plot(X, T[i,:]/1000, color = style[i][2] , marker = style[i][1],ls=style[i][0], ms=2,lw=0.5)
        i=i+1
    
    plt.ylim( 0, np.max(T/1000)*1.1)
    ax1.set_ylabel('($10^3 m^3/s$)',fontsize=6,labelpad=3,rotation=90)
    plt.setp(ax1.get_yticklabels(), fontsize=4)
    
    plt.ylim( 0, np.max(T/1000)*1.1)
    
    ax1.set_xticks(X)
    ax1.set_xticklabels(LabMonths, fontsize=6)
    ax1.tick_params(axis='y', which='major',pad=0.1,labelsize=6)    
    
    addcardgrdcnew(stname, chem_GRDC, basin)
    
    legend=ax1.legend(bbox_to_anchor=(1.05, 0.6, 0.2, 0.4),handles=LEG,fontsize=4,title=r'Legend',loc = 1, edgecolor="none")
    Outnum = NumObsStn(chem_GRDC,[stname],y1,y2)
    a=np.sum(Outnum[0,:])
    txt=str(round(int(a),0))+"/"+str((y2-y1+1)*12)+" months"
    plt.text(1.045,0.45,"Avalaible data",transform=ax1.transAxes, fontsize=6)
    plt.text(1.045,0.4,txt, transform=ax1.transAxes, fontsize=6)
    plt.setp(legend.get_title(),fontsize=8)
    
    fig.subplots_adjust(left=0.1, right=0.99,bottom=0.1, top=0.93,wspace= 0.04)
    fig.suptitle(r'Annual Cycle '+stname, fontsize=8,y=0.985)#loc="left"
    fig.savefig(dgraphs+stname.replace(" ","-").replace("/","-")+"-Annual_cycle.jpg",dpi=350)
    plt.close()
    return T, M, K
def plotallstn_annualcycle(Lst, L, chem_GRDC, y1, y2, dgraphs, basin):
    i=0
    while i<len(Lst):
        print "####"
        print i+1,"/",len(Lst)
        plot_annualcyclestn(Lst[i], L, chem_GRDC,y1,y2, dgraphs, basin)
        i=i+1


### Plot the time serie for a station and a list of data
def plottimeserie(stname, L, chem_GRDC, y1, y2, dgraphs, basin, chem_grid="", chem_grdc_rd=""): #Style included
    """
    plot the time serie
    """
    debug = None
    print "####"
    print stname   
    if debug: print "Get data"
    # Get data
    out = extract_timeseries(stname, L, chem_GRDC, y1, y2, chem_grid, chem_grdc_rd)  
    
    if out is None:
        print "Error - closed"
        return None
    # LEGEND
    LEG=[]
    i=0
    while i<len(L):
        LEG.append(mlines.Line2D([], [], color=style[i][2], marker=style[i][1],label=L[i][2],ls=style[i][0],ms=4))
        i=i+1
        
    if debug: print "Plot"
    # PLOT
    fig=plt.figure(figsize=(4.5,2.5),dpi=250)
    ax1 = plt.subplot2grid((1, 5), (0, 0), colspan=4)   
    i=0
    #X=np.arange(0,len(out[i][1]))
    X=np.arange(0,(y2-y1+1)*12)
    
    while i<len(out):
        print L[i][2]
        out0=monthmeantot(out[i][2],out[i][1],y1,y2) #data dtime y1 y2
        print "Mean value : ",round(ma.mean(out0/1000),2)
        ax1.plot(X, out0/1000, color = style[i][2] , marker = style[i][1],ls=style[i][0], ms=1,lw=style[i][3], markevery = 10)
        i=i+1
    out00=[0]*len(X)
    ax1.plot(X, out00, color = "black" , ls="-", lw=0.2)
        
    # ytick    
    ax1.set_ylabel('($10^3 m^3/s$)',fontsize=6,labelpad=3,rotation=90)
    plt.setp(ax1.get_yticklabels(), fontsize=4)
    # xtick
    xtickstimeMonth(y1, y2 , ax1)
    # Map
    addcardgrdcnew(stname, chem_GRDC, basin)
    # Legend
    legend=ax1.legend(bbox_to_anchor=(1.05, 0.6, 0.2, 0.4),handles=LEG,fontsize=4,title=r'Legend',loc = 1, edgecolor="none")
    plt.setp(legend.get_title(),fontsize=8)
    # Finalize    
    fig.subplots_adjust(left=0.1, right=0.99,bottom=0.1, top=0.93,wspace= 0.04)
    fig.suptitle(r'Time series '+stname, fontsize=8,y=0.985)#loc="left"
    fig.savefig(dgraphs+stname.replace(" ","-").replace("/","-")+"-timeserie.jpg",dpi=350)
    plt.close()
    
    return 
### Plot the time series graph for a list of stations
def plotallstn_timeseries(Lst, L, chem_GRDC, y1, y2, dgraphs, basin, chem_grid="", chem_grdc_rd=""):
    i=0
    while i<len(Lst):
        print "####"
        print i+1,"/",len(Lst)
        plottimeserie(Lst[i], L, chem_GRDC, y1, y2, dgraphs, basin, chem_grid, chem_grdc_rd)
        i=i+1

#####################
### Help Analysis ###
#####################
def NumObsStn(chem_GRDC, ListStn, y1, y2):
    """
    Get the number of months with observation for a List of station ListStn
    between the year y1 and y2 (included).
    Output Outnum[stationindex (from ListStn), year (starting from y1)]
    chem_GRDC: string, direction of GRDC file with observations.
    ListStn: array(string), List of the name of the station we are interested in.
    y1: int, beginning year.
    y2: int, ending year.
    """
    if y2<y1:
        print "Error : y2 must be superior to y1"
        return None
    GRDC = NetCDFFile(chem_GRDC, 'r')
    timv=GRDC.variables['time']
    time=timv[:] #seconds since 1807-01-15 00:00:00
    dtimegr = num2date(time,timv.units)
    hydrog = GRDC.variables['hydrographs'][:,:] #(time,stations) Original Hydrographs
    Name = GRDC.variables['name'][:]
    
    OutNum=np.zeros((len(ListStn),(y2-y1+1)))
    i=0
    while i<len(ListStn):
        y=y1
        k=stgrdcindex(ListStn[i], Name) -1 
        while y<y2+1:
            a=0
            j0=datebeg(dtimegr,y)
            j=j0
            while j<j0+12:
                if hydrog[j,k] is not ma.masked:
                    a=a+1
                j=j+1
            OutNum[i,y-y1]=a
            y=y+1
        i=i+1
    return OutNum


def AvailableStn(chem_GRDC_rd, chem_GRDC, basin, AR=False, BR=False):
    # Get stations output list on Parana basin
    GRDC_rd = NetCDFFile(chem_GRDC_rd, 'r')
    ind=[]
    for l in GRDC_rd.variables.keys():
        if l[0:len(basin)]==basin:
            ind.append(GRDC_rd.variables[l].Index_of_GRDC_Station)
    # Get station output information
    GRDC = NetCDFFile(chem_GRDC, 'r')
    Name = GRDC.variables["name"][:]
    country = GRDC.variables["country"]
    lon = GRDC.variables["lon"][:]
    lat = GRDC.variables["lat"][:]
    area = GRDC.variables["area"][:]
    L=[]
    if AR==True: Lar=[]
    i=0
    if BR==True: Lbr=[]
    i=0
    while i<len(ind):
        index=ind[i]-1
        namest = convertgrdcname(Name[index])
        pays = convertgrdcname(country[index])
        L.append([namest,pays, round(area[index]/1000,2), lon[index], lat[index]])
        if AR==True:
            if pays=="AR": Lar.append([namest,pays, round(area[index]/1000,2), lon[index], lat[index]])
        if BR==True:
            if pays=="BR": Lbr.append([namest,pays, round(area[index]/1000,2), lon[index], lat[index]])
        # M[i,0]=namest ; M[i,1]=pays ; M[i,2]=lon[index]; M[i,3]=lat[index]
        i=i+1
    if AR==True :
        if BR==True :
            return L, Lar, Lbr
        else:
            return L, Lar
    elif BR==True : 
        return L, Lbr
    else:
        return L


#### Plot all stn available basin ####
# Time series
def plotallstn_timeseries_basin(L, chem_GRDC, y1, y2, dgraphs, basin, chem_grid="", chem_grdc_rd=""): #actualiser format
    print basin
    Lavst = AvailableStn(chem_grdc_rd, chem_GRDC, basin, AR=False, BR=False)
    Lst=[]
    i=0
    while i<len(Lavst):
        Lst.append(Lavst[i][0])
        i=i+1
    plotallstn_timeseries(Lst, L, chem_GRDC, y1, y2, dgraphs, basin, chem_grid, chem_grdc_rd)
    return

# Annual Cycle
def plotallstn_annualcycle_basin(L, chem_GRDC, y1, y2, dgraphs, basin, chem_grdc_rd=""):
    Lavst = AvailableStn(chem_grdc_rd, chem_GRDC, basin, AR=False, BR=False)
    Lst=[]
    i=0
    while i<len(Lavst):
        Lst.append(Lavst[i][0])
        i=i+1
    plotallstn_annualcycle(Lst, L, chem_GRDC, y1, y2, dgraphs, basin) # Pas besoin de chemgrid ni chemgrdcrd???
    return

    

### TO DO ###
#############
#*ajouter dépendance dir_GRDC et dir_GRDC_rd
#*Cas station pas disponible pour l'une des liste - message d'erreur mais faire le graphique sans
