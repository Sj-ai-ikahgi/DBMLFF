import sys
import os
import argparse
import time
from copy import deepcopy
import numpy as np
from aux_dbmlff import *

isfile = open("./atom.config", 'r')
natoms, lattice, eindex, coords = Read_One_Config(isfile)
isfile.close

o2pId = np.arange(natoms)
forces = np.zeros([natoms,3], dtype=float)
velocity = np.zeros([natoms,3], dtype=float)

DBMLFF_MOVEMENT_headline = "%7d   1 %7d 1 \n" % (natoms, natoms)

for m in range(20):
  tcoords = deepcopy(coords)
  for i in range(len(tcoords)):
    tcoords[i][0] = tcoords[i][0] + m*(1.0/100.0/20.0)

  osfile = open("./shift_config/atom_" + str(m) + ".config", 'w')
  Write_MLFF_Init_Config(osfile, \
                      DBMLFF_MOVEMENT_headline, \
                      lattice, \
                      o2pId, \
                      eindex, \
                      tcoords, \
                      forces, \
                      velocity) 
  osfile.close()



