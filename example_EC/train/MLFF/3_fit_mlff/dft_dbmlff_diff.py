import os
import numpy as np
from aux_dbmlff import *


def Get_Eng_Force(ifilename, natoms):

  confNum = 0

  Etot_list = []
  f_list = []
  eindex = np.zeros(natoms, dtype=int)

  isfile = open(ifilename, 'r')
  while True:
    token = isfile.readline()
    if token == "":
      break
    else:
      Etot_list.append(float(token.split()[0]))
      force = np.zeros([natoms,3], dtype=float)
      for i in range(natoms):
        token = isfile.readline().split()
        eindex[i] = int(token[0])
        for j in range(3):
          force[i][j] = float(token[j+1])
      f_list.append(force)
      confNum += 1
  isfile.close()

  Etot_list = np.array(Etot_list)
  f_list = np.array(f_list)

  return confNum, Etot_list, f_list


# get atom nubmer
isfile = open("../1_SCF/energy_shift/atom.config", 'r')
token = isfile.readline().split()
natoms = int(token[0])
isfile.close()

# get energy shift
isfile = open("../1_SCF/energy_shift/energy_shift.dat", 'r')
token = isfile.readline().split()
energy_shift = float(token[0])
isfile.close()

print("@@@ natoms %7d" % (natoms))
print("@@@ energy_shift %12.6f" % (energy_shift))


for it, T in enumerate(["400K", "700K", "1200K"]):

  # get large distrotion mol Id
  ld_Id = []
  ifilename = "../MLFF_dataset/check_conf/" + T + "_DFT_large_distortion.dat"
  isfile = open(ifilename, 'r')
  tokens = isfile.readlines()
  isfile.close()
  for l in tokens:
    ld_Id.append(int(l.split()[0]))


  ifilename = "../1_SCF/" + T + "/dft_eng_force.dat"
  dft_confNum, dft_etot_list, dft_f_list = Get_Eng_Force(ifilename, natoms)

  ifilename = "../2_dbmlff/" + T + "/OUT.DB"
  dbmlff_confNum, dbmlff_etot_list, dbmlff_f_list = Get_Eng_Force(ifilename, natoms)

  confNum = min([dft_confNum, dbmlff_confNum])

  print(it, T, confNum)
  
  diff_etot_list = []
  diff_f_list = []
  for s in range(confNum):
    diff_etot_list.append((dft_etot_list[s] - energy_shift) - dbmlff_etot_list[s])
    diff_f_list.append(dft_f_list[s] - dbmlff_f_list[s])


  odir = "./data_set/" + T
  if not os.path.isdir(odir):
    os.system("mkdir -p " +odir)
  osfile = open(odir + "/MOVEMENT", 'w')

  for cId in range(confNum):

    if cId + 1 in ld_Id:
      print(" @@@ jump over %s %7d" % (T, cId + 1))
      continue

    ifilename = "../MLFF_dataset/Conf_" + T + "/DFT_" + str(cId + 1) + ".config"
    isfile = open(ifilename, 'r')
    natoms, lattice, eindex, coords = Read_One_Config(isfile)
    isfile.close()

    oheadline = "%7d iteration \n" % natoms
    forces = diff_f_list[cId]
    velocity = []
    atomic_eng = diff_etot_list[cId]
    Write_One_MOVEMENT(osfile, oheadline, lattice, eindex, coords, forces, velocity, atomic_eng)

  osfile.close()


