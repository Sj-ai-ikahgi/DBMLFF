import os
from aux_dbmlff import *



def Get_Point_All(ifilename):

  isfile = open(ifilename, 'r')

  cpoint_list = []

  total_point = 0

  while True:
    token = isfile.readline()
    if token == "":
      break
    else:
      total_point += 1
      token = token.split()
      pnum = int(token[0])
      for p in range(pnum):
        token = isfile.readline().split()
        if pnum == 1 and total_point%2 == 0:
          continue
        tcoord = []
        for i in range(3):
          tcoord.append(float(token[i+2]))
        cpoint_list.append(tcoord)
      
  isfile.close()

  return cpoint_list


def Merge_Coord(natoms, eindex, coords, cpoint_list, num):

  oeindex = []
  ocoords = []

  for i in range(natoms):
    oeindex.append(eindex[i])
    ocoords.append(coords[i])

  for i in range(num):
    oeindex.append(3)
    ocoords.append(cpoint_list[i])

  return natoms+num, oeindex, ocoords





file_list = ["2A", "4A", "6A", "2A.2A", "2A.4A", "2A.6A"]

for fileId in file_list:

  isfile = open("../xatom0.config", 'r')
  natoms, lattice, eindex, coords = Read_One_Config(isfile)
  isfile.close()

  cpoint_list = Get_Point_All("../point." + fileId)
  num = len(cpoint_list)

  natoms, eindex, coords = Merge_Coord(natoms, eindex, coords, cpoint_list, num)

  ofilename = "charge_position_" + fileId + ".xyz"
  os.system("rm -rf " + ofilename)

  coords = Lamda2Cart(lattice, coords)
  Write_XYZ(1,ofilename,natoms,eindex,coords)



