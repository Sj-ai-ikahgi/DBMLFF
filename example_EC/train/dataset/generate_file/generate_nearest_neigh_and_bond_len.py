import sys
import numpy as np 
from aux_dbmlff import *




def Neigh_and_Bond_List(molType):

  neigh_cutoff_H_C = 1.3
  neigh_cutoff_X_X = 1.8

  isfile = open("./reference_conf/atom_" + molType + ".config", 'r')
  natoms, lattice, eindex, coords = Read_One_Config(isfile)
  isfile.close()

  neigh_num = np.zeros(natoms, dtype=int)
  neigh_list = [[] for i in range(natoms)]
  neigh_bond_len = [[] for i in range(natoms)]

  for i in range(natoms):
    for j in range(natoms):
      if i != j:
        dr, d = two_atom_distance(coords, lattice, i, j)
        if 1 in [eindex_list[i], eindex_list[j]]:
          if (d < neigh_cutoff_H_C):
            neigh_num[i] += 1
            neigh_list[i].append(j+1)
            neigh_bond_len[i].append(d)
        elif d < neigh_cutoff_X_X:
          neigh_num[i] += 1
          neigh_list[i].append(j+1)
          neigh_bond_len[i].append(d)

  return natoms, neigh_num, neigh_list, neigh_bond_len



def Write_Neigh_List(molType, natoms, neigh_num, neigh_list):

  print("   !!! create new file " + "./results/nearest_neigh_list_" + molType + ".dat")

  osfile = open("./results/nearest_neigh_list_" + molType + ".dat", 'w')
  osfile.write("   %7d\n" % (natoms))
  for i in range(natoms):
    osfile.write("%5d   %3d   " % (i+1, neigh_num[i]))
    for nid in range(neigh_num[i]):
      osfile.write("%3d " % (neigh_list[i][nid]))
    osfile.write("\n")
  osfile.close()



def Write_Bond_Len(molType, natoms, neigh_num, neigh_list, bond_len):

  bond_num = 0
  for i in range(natoms):
    bond_num += neigh_num[i]
  bond_num = round(bond_num*0.5)

  print("   !!! create new file " + "./results/bond.molecule." + molType + "_need_modify")

  osfile = open("./results/bond.molecule_" + molType + "_need_modify", 'w')
  osfile.write("   %7d\n" % (bond_num))
  for i in range(natoms):
    for nid in range(neigh_num[i]):
      j = neigh_list[i][nid]
      if i < j:
        osfile.write("%3d %3d   %7.3f  20.0 \n" % (i+1, j, bond_len[i][nid]))
  osfile.close()






print("   @@@ start generate nearest neighbor and those bond length")
print("   @@@ this process need    " + 
      " 'atom_*.config', 'MolEindexList', "+
      " 'Read_One_Config()', 'two_atom_distance()'")


molType = sys.argv[1]
eindex_list = MolEindexList[molType]



natoms, neigh_num, neigh_list, neigh_bond_len = Neigh_and_Bond_List(molType)
Write_Bond_Len(molType, natoms, neigh_num, neigh_list, neigh_bond_len)
Write_Neigh_List(molType, natoms, neigh_num, neigh_list)

print("   ### Please MODIFY the file ./results/bond.molecule." + molType)
print("       add bond type, bond type number and ")
print("       change the file name to ./results/bond.molecule." + molType + "_bt")

print("\n")

