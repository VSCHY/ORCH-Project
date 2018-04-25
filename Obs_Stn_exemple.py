#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Application of tools for observations on station's upstream area
Created on Wed Apr 25 14:22:20 2018

@author: Anthony Schrapffer
@institution: CIMA, Universidad de Buenos Aires (Argentina)
"""

"""
Convertir to monthly 
GPCC : 1º, precip(time, lat, lon) en mm/mn, time // (X.5º)
GPCP : 0.5º, precip(time, lat, lon) en mm/day, time // 
3CN : 0.5 º, rr(time, latitude, longitude) (mm), time, 
CRU : 0.5º, pre(time, lat, lon) mm/month
CPC_UNI : 0.5º, precip(time, level001, latitude, longitude) mm/day (meme si mensuel)
E2OFD : 0.25º, rain, a mettre sur le serveur (demander si possible) daily mm/d
"""

sys.path.append(os.getcwd())
import Obs_Stn_tools as Obs

################
#### Prueba #### 
################
#
### For one station ###
#
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
La=[chem_file, "GPCC", "precip", "time", k, True, chem_grid_alt]
# [direction, name, variablename, timename, k-coeff multiplication =1 daily =31 monthly, 3CN - altgrdc (True ou direction grdc rd), chemgrid si alt]
# Erreur creer par monthly, besoin de [31,28.25,31,30,31,30,31,31,30,31,30,31]
L=[La]

M0, M1= plotstn_obs_annualcycle("Corrientes", L, chem_grdc_rd, chem_grdc, chem_grid, dgraphs, y1, y2, style, basin)
##############################################################################
#
### Stations from a list ###
#
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

Lst=["Corrientes", "Posadas"]

plotallstn_obs_annualcycle(Lst L, chem_grdc_rd, chem_grdc, chem_grid, dgraphs, y1, y2, style, basin)

##############################################################################
#
### All station from a basin ###
#
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

plotallbasin_obs_annualcycle(L, chem_grdc_rd, chem_grdc, chem_grid, dgraphs, y1, y2, style, basin)
