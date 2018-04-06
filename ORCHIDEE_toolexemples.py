#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 07:10:37 2018

@author: Anthony Schrapffer
@institution: CIMA, Universidad de Buenos Aires (Argentina)
"""
import sys
sys.path.append("/home/anthony/TotiTools/ORCH-Project/")
import ORCHIDEE_tools as OR

### ORCHIDEE Valid Stations - La Plata Basin 0.5° ###
# List of name
L2=["Fecho Dos Morros","La Punilla","Andira","Barbosa Ferraz","Zanja Del Tigre",
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





### ANNUAL CYCLE ###
# Parameters
dgraphs="/home/direction/of/your/graph/"

# Data direction
dir_GRDC=dorig+"GRDC_Monthly_June14_LaPlataBasin.nc"
dir_E2OFD_GRDC=dorig+"SAM_E2OFD_orchidee_Dis_Stn_GRDC_1979_1990.nc"
dir_E2OFD=dorig+"SAM_E2OFD_orchidee_Dis_Stn_Model_1979_1990.nc"
dir_WFDEI=dorig+"SAM_WFDEI_GPCC_flood_orchidee_hydrographs_1979_1987.nc"

# Descriptive list for each simulation
# simulation direction, filetype, name of simulation
# filetype : GRDCnew, new, old, rain, evap
La = [dir_E2OFD_GRDC, "GRDCnew", "GRDC"]
Lb = [dir_E2OFD, "new", "E2OFD"]
Lc = [dir_WFDEI, "old", "WFDEI GPCC flood"]

# Creation of simulation list
L=[La, Lb, Lc]

# Plot the data
# stname: station string name (or first part) in GRDC
# L: list of simulation
# chem_GRDC: direction of GRDC original file
# y1, y2: first and last year (last year included)
# dgraphs: direction where the graph will be saved
# basin: name of the basin, must be in the graph
TT =OR.plot_annualcyclestn(stname, L, chem_GRDC, y1, y2, dgraphs, basin)





### TIME SERIE  ###
# Parameters
dgraphs="/home/direction/of/your/graph/"

# Data direction
dir_GRDC=dorig+"GRDC_Monthly_June14_LaPlataBasin.nc"
dir_grdc_rd="/dir/grdc_river_desc.nc"
dir_grid

dir_E2OFD_GRDC=dorig+"SAM_E2OFD_orchidee_Dis_Stn_GRDC_1979_1990.nc"
dir_E2OFD=dorig+"SAM_E2OFD_orchidee_Dis_Stn_Model_1979_1990.nc"
dir_WFDEI=dorig+"SAM_WFDEI_GPCC_flood_orchidee_hydrographs_1979_1987.nc"
dir_E2OFD_rain = dorig + "SAM_E2OFD_orchidee_rain_1979_1990.nc"
dir_E2OFD_evap = dorig + "SAM_E2OFD_orchidee_evap_1979_1990.nc"
dir_E2OFD_pr_ev = dorig + "SAM_E2OFD_orchidee_pr-ev_1979_1990.nc"


La = [dir_E2OFD_GRDC, "GRDCnew", "GRDC"]
Lb = [dir_E2OFD, "new", "E2OFD"]
Lc = [dir_WFDEI, "old", "WFDEI GPCC flood"]
Ld = [dir_E2OFD_rain, "rain", "E2OFD rain"]
Le = [dir_E2OFD_evap, "evap", "E2OFD evap"]
Lf = [dir_E2OFD_pr_ev, "rain", "E2OFD pr-ev"]

L1=[La, Lb, Lc, Lf]

# Plot time serie for list of simulation for a station
# stname: station string name (or first part) in GRDC
# L: list of simulation
# chem_GRDC: direction of GRDC original file
# y1, y2: first and last year (last year included)
# dgraphs: direction where the graph will be saved
# basin: name of the basin (must be in the list)
# chem_grid: direction of the cell area grid file corresponding to output
# chem_grdc_rd: direction of grdc_river_desc.nc file corresponding
OR.plottimeserie(stname, L, chem_GRDC, y1, y2, dgraphs, basin, chem_grid, chem_grdc_rd, basin)

# Plot all station from list Lst
Lst=[stname1, stname2,...]
OR.plotallstn_timeseries(Lst, L, chem_GRDC, y1, y2, dgraphs, basin, chem_grid, chem_grdc_rd)


#All station from a basin (Same parameter in the GRDC module and same original file)
OR.plotallstn_timeseries_basin(L, chem_GRDC, y1, y2, dgraphs, basin, chem_grid="", chem_grdc_rd=""): #actualiser format