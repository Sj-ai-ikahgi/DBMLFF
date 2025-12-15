import sys
import numpy as np 


def Get_Eng(ifilename):

  E_tot = 0.0

  isfile = open(ifilename, 'r')
  while True:
    token = isfile.readline()
    if token == "":
      break
    else:
      if "E_tot(eV)    =" in token:
        token = token.split()
        E_tot = float(token[2])
        break
  isfile.close()

  return E_tot


def Get_Force(ifilename, natoms):

  force = np.zeros([natoms,3], dtype=float)
  eindex = np.zeros(natoms, dtype=int)

  isfile = open(ifilename, 'r')

  for i in range(natoms+4):
    token = isfile.readline()

  for i in range(natoms):
    token = isfile.readline().split()
    eindex[i] = int(token[0])
    for j in range(3):
      force[i][j] = float(token[j+1])

  isfile.close()

  return eindex, force


def Write_Eng_Force(ofilename, confNum, natoms, eindex, Etot_list, f_list):

  osfile = open(ofilename, 'w')

  for cId in range(confNum):
    osfile.write("%20.9f\n" % Etot_list[cId])
    f = f_list[cId]
    for i in range(natoms):
      osfile.write("%7d %20.9f %20.9f %20.9f\n" %\
        (eindex[i], f[i][0], f[i][1], f[i][2]))

  osfile.close()



print("### get_dft_eng_force.py start ###")


natoms = int(sys.argv[1])
total_conf = int(sys.argv[2])
ssize = int(sys.argv[3])


pnum = int(total_conf/ssize)

# print(pnum)
spoint_list = []
epoint_list = []


for pId in range(pnum):
  spoint=int(pId*ssize+1)
  epoint=int(spoint + ssize - 1)
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
confNum = 0

for pId in range(pnum):
# for pId in range(1):
  print(pId, spoint_list[pId], epoint_list[pId])
  spoint = spoint_list[pId]
  epoint = epoint_list[pId]

  # if (1 <= spoint <=  200) or (801<= spoint<=1000):
  #   continue

  ifilehead = str(spoint) + "_" + str(epoint)
  for cId in range(spoint, epoint+1):
    ifilename = ifilehead + "/" + str(cId)

    E_tot = Get_Eng(ifilename + "/REPORT")
    eindex, force = Get_Force(ifilename + "/OUT.FORCE", natoms)

    confNum += 1

    Etot_list.append(E_tot)
    f_list.append(force)

ofilename = "dft_eng_force.dat"
Write_Eng_Force(ofilename, confNum, natoms, eindex, Etot_list, f_list)

