import numpy as np
import matplotlib.pyplot as plt
from aux_dbmlff import *


def Get_MOVEMENT_ALL(ifilename):

  natoms = 0
  numConf = 0
  lattice = []
  eng_list = []
  coord_list = []
  f_list = []

  isfile = open(ifilename)

  while True:
    headline = New_MOVEMENT_Conf(isfile) 
    if headline == '': break

    ntimestep, Ep, natoms, lattice, eindex, \
    coords, forces, velocity, atomic_eng =  \
      Get_One_MOVEMENT_Conf(headline, isfile)

    numConf += 1
    eng_list.append(Ep)
    coord_list.append(coords)
    f_list.append(forces)

  isfile.close()
  
  eng_list = np.array(eng_list)
  coord_list = np.array(coord_list)
  f_list = np.array(f_list)

  return natoms, numConf, lattice, eng_list, coord_list, f_list


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


def Fcheck(natoms, numConf, lattice, eng_list, coord_list, f_list):

  ffc_list = []
  efc_list = []
  fcheck_list = []

  for cId in range(numConf - 2):

    etot1, etot2 = eng_list[cId+1], eng_list[cId+2]
    x1, x2 = coord_list[cId], coord_list[cId+1]
    f1, f2 = f_list[cId+1], f_list[cId+2]

    dx = np.zeros([natoms,3], dtype=float)
    for i in range(natoms):
      for j in range(3):
        dr = x1[i][j] - x2[i][j]
        if np.fabs(dr) > 0.5:
          dr = dr - np.sign(dr)
        dx[i][j] = dr*lattice[j][j]

    sf = (f1 + f2)*0.5
    ffc = 0.0
    for i in range(natoms):
      for j in range(3):
        ffc += dx[i][j]*sf[i][j]
    efc = etot1 - etot2

    if abs(efc) < 1.0e-6:
      efc = 1.0e-6
      print(" @@@ too small efc %7d %20.12f" % (cId+1, efc))
  
    ffc_list.append(ffc)
    efc_list.append(efc)

    fcheck_list.append(ffc/efc)
  
  return ffc_list, efc_list, fcheck_list


strategy = "md_test"

fig, ax = plt.subplots(4, 1, figsize=(15,6), constrained_layout=True)


for ei, epart in enumerate(["energy_force", "ML_FF", "DB", "POL"]):

  print(epart)

  ifilename = "./" + strategy + "/MOVEMENT"
  natoms, numConf, lattice, eng_list, coord_list, f_list = Get_MOVEMENT_ALL(ifilename)

  ifilename = "./" + strategy + "/OUT." + epart
  numConf_f, eng_list, f_list = Get_Eng_Force(ifilename, natoms)

  if numConf != numConf_f-1:
    print("@@@ numConf not same! %7d %7d" % (numConf, numConf_f))

  ffc_list, efc_list, fcheck_list = Fcheck(natoms, numConf, lattice, eng_list, coord_list, f_list)

  ax[ei].plot(fcheck_list, 'o', ms=1, label=epart)
    # ax[c].legend(title=elem_label[e] + "-" + elem_label[c],loc='lower right', ncol=1, fontsize=9)

  ax[ei].set_ylabel(epart, fontsize=12)
  ax[ei].set_ylim(-2,3)

ax[3].set_xlabel("steps", fontsize=16)
# ax[2][2].set_xlabel("dx percentage", fontsize=12)
# ax[0][0].set_ylabel("dx*F", fontsize=12)
# ax[1][0].set_ylabel("dE", fontsize=12)
# ax[2][0].set_ylabel("fcheck", fontsize=12)

# for r in range(3):
#   for c in range(5):
#     ax[r][c].legend(bbox_to_anchor=(0.2, 0.1),
#                             loc='upper left', borderaxespad=0.)

plt.suptitle("Fcheck " + strategy, fontsize=16)
plt.savefig("./results_fcheck/fcheck_" + strategy + ".png", dpi=600)