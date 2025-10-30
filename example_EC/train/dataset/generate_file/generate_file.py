import sys
import os
import numpy as np 
from aux_dbmlff import *


def Write_MOVEMENT_TYPE(molType):

  isfile = open("./reference_conf/atom_" + molType + ".config", 'r')
  natoms, lattice, eindex, coords = Read_One_Config(isfile)
  isfile.close()

  atomic_symmetry = MolAtomSysList[molType]

  osfile = open("./results/MOVEMENT.type_" + molType, 'w')
  osfile.write("%7d\n" % (natoms))
  for i in range(natoms):
    osfile.write("%5d   %3d   %12.9f %12.9f %12.9f    1  1  1\n" % \
      (eindex[i], atomic_symmetry[i], coords[i][0], coords[i][1], coords[i][2]))
  osfile.close()



def Write_Mol_Atom_Bond_Info(molType):

  eindex = MolEindexList[molType]
  atomic_sys = MolAtomSysList[molType]

  natoms = len(eindex)
  atomic_type_num = max(atomic_sys)


  isfile = open("./results/bond.molecule_" + molType + "_bt", 'r')
  token = isfile.readline().split(',')
  bond_num = int(token[0])
  bond_type_num = int(token[1])

  bond_type = np.zeros(bond_num, dtype=int)
  bond_pair = np.zeros([bond_num,2], dtype=int)
  for i in range(bond_num):
    token = isfile.readline().split()
    bond_type[i] = int(token[0])
    for j in range(2):
      bond_pair[i][j] = int(token[j+1])
  isfile.close()


  osfile = open("./results/mol_atom_bond_info.dat_" + molType, 'w')

  osfile.write("%7d %7d\n" % (natoms+bond_num, atomic_type_num+bond_type_num))
  for i in range(natoms):
    # print(i, eindex[i], zValenceList[eindex[i]])
    osfile.write("%3d %6.2f   1  %3d   1.000000\n" % \
     (atomic_sys[i], zValenceList[eindex[i]],i+1))

  for i in range(bond_num):
    osfile.write("%3d   0.0    2  %3d %3d   0.5000  0.5000\n" % \
      (atomic_type_num+bond_type[i], bond_pair[i][0], bond_pair[i][1]))

  osfile.close()


  

def Write_Mol_Info(molType):

  eindex = MolEindexList[molType]
  atomic_sys = MolAtomSysList[molType]

  natoms = len(eindex)
  atomic_type_num = max(atomic_sys)

  ofilename_1 = "./results/mol.info_" + molType
  ofilename_2 = "./results/mol_atom_bond_info.dat_" + molType

  osfile_1 = open(ofilename_1, 'w')
  osfile_1.write("%7d %7d \n" % (natoms, atomic_type_num))
  for i in range(natoms):
    osfile_1.write("%5d   %3d  %3d  %6.3f\n" % \
      (i+1, atomic_sys[i], eindex[i], atomicMassList[eindex[i]]))
  osfile_1.write("----------------------------\n")
  osfile_1.close()


  isfile = open("./results/bond.molecule_" + molType + "_bt", 'r')
  token = isfile.readline().split(',')
  bond_num = int(token[0])
  bond_type_num = int(token[1])

  bond_type = np.zeros(bond_num, dtype=int)
  bond_pair = np.zeros([bond_num,2], dtype=int)
  for i in range(bond_num):
    token = isfile.readline().split()
    bond_type[i] = int(token[0])
    for j in range(2):
      bond_pair[i][j] = int(token[j+1])
  isfile.close()


  osfile = open(ofilename_2, 'w')

  osfile.write("%7d %7d\n" % (natoms+bond_num, atomic_type_num+bond_type_num))
  for i in range(natoms):
    # print(i, eindex[i], zValenceList[eindex[i]])
    osfile.write("%3d %6.2f   1  %3d   1.000000  %3d\n" % \
     (atomic_sys[i], zValenceList[eindex[i]],i+1, eindex[i]))

  for i in range(bond_num):
    osfile.write("%3d   0.0    2  %3d %3d   0.5000  0.5000  0\n" % \
      (atomic_type_num+bond_type[i], bond_pair[i][0], bond_pair[i][1]))

  osfile.close()

  os.system("cat " + ofilename_2 + " >> " + ofilename_1)







molType = sys.argv[1]

print("   @@@ start generate file ")
print("   @@@ this process need    " +
      " 'atom_*.config', 'Read_One_Config()', " +
      " 'MolAtomSysList'")
Write_MOVEMENT_TYPE(molType)
print("   @@@ this process need    " +
      " 'MolEindexList', 'MolAtomSysList', " +
      " 'bond.molecule_*_bt'")
Write_Mol_Atom_Bond_Info(molType)
print("   @@@ this process need    " +
      " 'MolEindexList', 'MolAtomSysList', " +
      " 'atomicMassList', 'bond.molecule_*_bt' " +
      " 'zValenceList' \n")
Write_Mol_Info(molType)

osfile = open("./results/find_neighbore.in", 'w')
osfile.write("3\n")
osfile.close()
