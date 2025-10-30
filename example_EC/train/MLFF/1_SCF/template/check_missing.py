import os
import sys
import numpy as np 



natoms = int(sys.argv[1])
total_conf = int(sys.argv[2])
ssize = int(sys.argv[3])



pnum = int(total_conf/ssize)

# print(pnum)
spoint_list = []
epoint_list = []


for pId in range(pnum):
  spoint=int(pId*200+1)
  epoint=int(spoint + 200 - 1)
  spoint_list.append(spoint)
  epoint_list.append(epoint)

if pnum < total_conf/ssize:
  pnum += 1
  spoint = epoint_list[-1]+1
  epoint = total_conf
  spoint_list.append(spoint)
  epoint_list.append(epoint)


eindex = []
Etot_list = []
f_list = []

for pId in range(pnum):
  # print(pId, spoint_list[pId], epoint_list[pId])
  spoint = spoint_list[pId]
  epoint = epoint_list[pId]

  ifilehead = str(spoint) + "_" + str(epoint)
  for cId in range(spoint, epoint+1):
    ifilename = ifilehead + "/" + str(cId)

    # os.system("rm " + ifilename + "/OUT.GKK")
    # os.system("rm " + ifilename + "/OUT.RHO")
    # os.system("rm " + ifilename + "/OUT.VR")
    # os.system("rm " + ifilename + "/OUT.VR_hion")

    # E_tot = Get_Eng(ifilename + "/REPORT")
    # eindex, force = Get_Force(ifilename + "/OUT.FORCE", natoms)
    
    if (not os.path.exists(ifilename + "/OUT.FORCE")) :
      print(spoint, epoint, cId)