#!/bin/sh   

# change binning

KBIN=1.
KCENTER=1.

python test.py --kbin $KBIN --kcenter $KCENTER 

# change grid and LOS 
#(note that change in LOS is not compatible with Sancho dataset, here just to show a test case)

GRID=64
LOSx=0
LOSy=2
LOSz=1

python test.py --grid $GRID --los ${LOSx} ${LOSy} ${LOSz}

# How to put to false a parameter in the config file: 
#Ex: measurements with no interlacing

python test.py --interlacing ""