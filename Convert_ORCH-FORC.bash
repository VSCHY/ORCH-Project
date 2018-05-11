#!/bin/bash
#
# Convert lon / lat for an ORCHIDEE Forcing
#
# Conf outil utilise (DIR)
PYTHON=""
#
# Directions
ORCHTools="/"
dir_Forc=$ORCHTools"FORCING_CONVERT_tools.py"
DIRIN="/bdd/ORCHIDEE_Forcing/BC/OOL/OL2/E2OFD/"
Input=$DIRIN"E2OFD_$YEAR.nc"
DIROUT="./"
Output="E2OFD_1979_convpart_"
datan="Rainf"
NC=".nc"
#
# Parameters
YEAR=1979
L=100
#
counter=0
while [$counter -le 28]; do
   echo $counter-$(($counter+2))
#
   $PYTHON $dir_Forc -o conversionlonlat_parallel -i $Input -w $DIROUT$Output$counter$NC -v $datan -l $L -p $counter
   $counter=$(($counter+1))
# 
   $PYTHON $dir_Forc -o conversionlonlat_parallel -i $Input -w $DIROUT$Output$counter$NC -v $datan -l $L -p $(($counter+1))
   $counter=$(($counter+1))
#
   $PYTHON $dir_Forc -o conversionlonlat_parallel -i $Input -w $DIROUT$Output$counter$NC -v $datan -l $L -p $(($counter+2))
   $counter=$(($counter+1))
done




