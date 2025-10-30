import sys
import numpy as np 
from aux_dbmlff import *


def Read_Nearest_Neigh(molType):

  isfile = open("./results/nearest_neigh_list_" + molType + ".dat", 'r')
  token = isfile.readline().split()
  natoms = int(token[0])

  nearest_neigh_num = np.zeros(natoms, dtype=int)
  nearest_neigh_list = np.zeros([natoms, natoms], dtype=int)

  for i in range(natoms):
    token = isfile.readline().split()
    aid = int(token[0]) - 1
    tnum = int(token[1])
    nearest_neigh_num[i] = tnum
    for j in range(tnum):
      nearest_neigh_list[i][j] = int(token[j+2]) - 1
  isfile.close()

  return natoms, nearest_neigh_num, nearest_neigh_list



def Shortest_Path(natoms, nearest_neigh_num, nearest_neigh_list, molType):

  dis_matrix = np.zeros([natoms, natoms], dtype=int)
  for i in range(natoms):
    for j in range(natoms):
      dis_matrix[i][j] = 10**2

  # BFS
  for cid in range(natoms):
    visited = np.zeros(natoms, dtype=int)
    visited[cid] = 1
    loop_list = []
    dis_list = []
    for i in range(nearest_neigh_num[cid]):
      loop_list.append(nearest_neigh_list[cid][i])
      dis_list.append(1)
      # visited[nearest_neigh_list[cid][i]] = 1
    while loop_list != []:
      nid = loop_list.pop(0)
      dis = dis_list.pop(0)
      visited[nid] = 1
      if dis_matrix[cid][nid] > dis:
        dis_matrix[cid][nid] = dis
      for i in range(nearest_neigh_num[nid]):
        if visited[nearest_neigh_list[nid][i]] == 0:
          loop_list.append(nearest_neigh_list[nid][i])
          dis_list.append(dis+1)

  for i in range(natoms):
    for j in range(i,natoms):
      if dis_matrix[i][j] != dis_matrix[j][i]:
        print("@@@ ERROR dis_matrix")

  print("   !!! create new file " + "./results/dis_matrix_" + molType + ".dat")

  osfile = open("./results/dis_matrix_" + molType + ".dat", 'w')
  osfile.write("    ")
  for i in range(natoms):
    osfile.write("%3d " % (i+1))
  osfile.write("\n")
  for i in range(natoms):
    osfile.write("%3d " % (i+1))
    for j in range(natoms):
      osfile.write("%3d " % dis_matrix[i][j])
    osfile.write("\n")
  osfile.close()







print("   @@@ start generate distance matrix (the shorest bond-path length)")
print("   @@@ this process need    " +
      " 'nearest_neigh_list_*.dat'")

molType = sys.argv[1]

natoms, nearest_neigh_num, nearest_neigh_list = Read_Nearest_Neigh(molType)
Shortest_Path(natoms, nearest_neigh_num, nearest_neigh_list, molType)


print("\n")



