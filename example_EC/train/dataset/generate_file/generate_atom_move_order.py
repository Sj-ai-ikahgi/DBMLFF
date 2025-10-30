import sys
import numpy as np 
from aux_dbmlff import *

def Read_Dis_Matrix(molType,natoms):

  dis_matrix = np.zeros([natoms, natoms], dtype=int)

  ifilename = "./results/dis_matrix_" + molType + ".dat"
  isfile = open(ifilename, 'r')

  isfile.readline()
  for i in range(natoms):
    token = isfile.readline().split()
    for j in range(natoms):
      dis_matrix[i][j] = int(token[j+1])
  isfile.close()

  return dis_matrix


def Write_Atom_Move_Order(molType, natoms, atom_dist, saId):

  visited = np.zeros(natoms, dtype=int)
  visited[saId - 1] = 1

  nid = saId - 1
  move_order = [nid+1] 

  loop_list = []
  for i in range(natoms):
    if (visited[i] == 0) and (atom_dist[nid][i] == 1):
      loop_list.append(i)

  while loop_list  != []:
    # print(np.array(loop_list)+1)
    # # BFS
    # nid = loop_list.pop(-1)
    # DFS
    nid = loop_list.pop(0)
    if visited[nid] == 0:
      move_order.append(nid+1)
    visited[nid] = 1

    for i in range(natoms):
      if (visited[i] == 0) and (atom_dist[nid][i] == 1):
        loop_list.append(i)

  # print(move_order)

  osfile = open("./results/move_atom_order.dat_" + molType, 'w')
  osfile.write("%7d\n" % (natoms))
  for i,aid in enumerate(move_order[:-1]):
    if i % 5 == 0:
      osfile.write("%5d," % (aid))
    else:
      osfile.write("%d," % (aid))
  osfile.write("%d" % move_order[-1])
  osfile.write("\n")
  osfile.close()




print("   @@@ start generate the atom move order")
print("   @@@ this process need    " +
      " 'MolAtomSysList' 'dis_matrix_*.dat")

molType = sys.argv[1]
saId = moveStartAtomId[molType]

atomic_symmetry = MolAtomSysList[molType]
natoms = len(atomic_symmetry)

atom_dist = Read_Dis_Matrix(molType, natoms)

Write_Atom_Move_Order(molType, natoms, atom_dist, saId)




