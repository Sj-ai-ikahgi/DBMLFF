import sys
import numpy as np 

whome="/data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC"
sys.path.append(whome + "/dataset")

from aux_dbmlff import *


def Read_Dis_Matrix(molType,natoms):

  dis_matrix = np.zeros([natoms, natoms], dtype=int)

  ifilename = whome + "/dataset/generate_file/results/dis_matrix_" + molType + ".dat"
  isfile = open(ifilename, 'r')

  isfile.readline()
  for i in range(natoms):
    token = isfile.readline().split()
    for j in range(natoms):
      dis_matrix[i][j] = int(token[j+1])
  isfile.close()

  return dis_matrix


def Write_Neigh(molType, natoms, atom_dist, atomic_symmetry):

  bond_cutoff = 3.1

  for i in range(natoms):
    for j in range(i):
      atom_dist[i][j] = atom_dist[j][i]

  neigh_list = [[] for i in range(natoms)]
  for i in range(natoms):
    for j in range(natoms):
      if 0.1 < atom_dist[i][j] < bond_cutoff:
        neigh_list[i].append(j+1)

  num_type = max(atomic_symmetry)

  neigh_type_num = np.zeros([natoms, num_type], dtype=int)
  for i in range(natoms):
    for nid in range(len(neigh_list[i])):
      ntype = atomic_symmetry[neigh_list[i][nid]-1]
      neigh_type_num[i][ntype-1] += 1

  print("   !!! create new file " + "./results/IN.NEIGHBORE_" + molType)

  osfile = open("./IN.NEIGHBORE_all", 'w')

  osfile.write("%5d %5d \n" % (natoms, num_type))

  for cyc in range(2):
    for i in range(natoms):
      osfile.write("%5d " % (i+1))
      for j in range(num_type):
        osfile.write("%5d " % neigh_type_num[i][j])
      osfile.write("\n")
    for i in range(natoms):
      nId = 0
      for t in range(num_type):
        for j in range(neigh_type_num[i][t]):
          osfile.write("%5d %5d %5d %5d %5d\n" % \
           (i+1, t+1, t+1, neigh_list[i][nId], j+1))
          nId += 1
  osfile.close()



print("   @@@ start generate the neighbor list")
print("   @@@ this process need    " +
      " 'MolAtomSysList' 'dis_matrix_*.dat")

molType = "EC"

atomic_symmetry = MolAtomSysList[molType]
natoms = len(atomic_symmetry)

atom_dist = Read_Dis_Matrix(molType, natoms)
Write_Neigh(molType, natoms, atom_dist, atomic_symmetry)


print("\n")