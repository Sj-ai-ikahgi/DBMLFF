
import argparse
from aux_file import *

from copy import deepcopy
import numpy as np
import sys
import os
import time




class sMSD:


  def __init__(self,infilename,ofilehead):

    self.dt_num = 50
    self.infilename = infilename
    self.ofilehead = ofilehead
    self.all_coords = []
    self.all_atype = []

    self.total_conf_num = 0
    self.half_conf_num = 0
    self.natoms = 0

    self.dt_list = []
    self.timeStep = 0.5

    self.stime = time.time()


  def Read_Data(self):

    isfile = open(self.infilename,'r')
    cid = 0
    New_Conf = New_DUMP_Conf(isfile)
    while New_Conf:
      ntimestep, natoms, lattice, atype, unwrapped_coords, forces = \
        Get_One_DUMP_Conf(isfile)
      New_Conf = New_DUMP_Conf(isfile)

      if cid%1000 == 0:
        print('-'*10 + str(ntimestep) + " start " + '-'*10)

      cid += 1
      if cid == 1:
        self.charge = np.zeros(natoms,dtype=float)
        for i in range(natoms):
          if atype[i] == 0:
            atype[i] = 1
            self.charge[i] = 1
          elif atype[i] == 1:
            atype[i] = 2
            self.charge[i] = -1
        self.all_atype = deepcopy(atype)

      # if cid < 10000 or cid > 35000 : continue
      if cid < 3000 : continue

      self.all_coords.append(unwrapped_coords)
    isfile.close()

    if len(self.all_coords) < 1:
      print("There is no coordiante!")
      sys.exit()

    self.total_conf_num = np.shape(self.all_coords)[0]
    self.half_conf_num = int(self.total_conf_num*0.5)

    atype_list, each_atype_num = np.unique(atype, return_counts=True)
    self.atype_num = len(atype_list)
    self.natoms = np.shape(self.all_coords)[1]

    print(self.total_conf_num)
    print(np.shape(self.all_coords))


  def Get_Time_List(self):

    elog = np.log(self.half_conf_num)
    point = np.linspace(1,elog,self.dt_num)

    point = pow(np.e,point)
    for i,p in enumerate(point):
      tx = int(p)
      if len(self.dt_list) != 0:
        if tx != self.dt_list[-1]:
          self.dt_list.append(tx)
      else:
        self.dt_list.append(tx)


  def Calculate(self):

    self.lamda_app = np.zeros([self.dt_num], dtype=float)
    self.self_cation = np.zeros([self.dt_num], dtype=float)
    self.self_anion = np.zeros([self.dt_num], dtype=float)
    self.dist_cation = np.zeros([self.dt_num], dtype=float)
    self.dist_anion = np.zeros([self.dt_num], dtype=float)
    self.dist_ca = np.zeros([self.dt_num], dtype=float)

    cnum = np.zeros([self.dt_num,6],dtype=int)

    for idt,dt in enumerate(self.dt_list):
      print(idt, dt)
      tnum = 0
      # for sid in range(0,self.half_conf_num,100):
      for sid in range(0,self.total_conf_num-dt,100):
        eid = sid + dt
        scoord = self.all_coords[sid]
        ecoord = self.all_coords[eid]

        for i in range(self.natoms):
          for j in range(self.natoms):

            dr2 = 0.0
            for k in range(3):
              dx1 = ecoord[i][k] - scoord[i][k]
              dx2 = ecoord[j][k] - scoord[j][k]
              dr2 += dx1*dx2

            zizj = self.charge[i]*self.charge[j]
            dr2 = dr2*zizj
            self.lamda_app[idt] += dr2
            cnum[idt][0] += 1

            if (i == j) and (self.all_atype[i] == 1):
              self.self_cation[idt] += dr2
              cnum[idt][1] += 1
              
            elif (i == j) and (self.all_atype[i] == 2):
              self.self_anion[idt] += dr2
              cnum[idt][2] += 1

            elif (i != j) and (self.all_atype[i] == 1) and (self.all_atype[j] == 1):
              self.dist_cation[idt] += dr2
              cnum[idt][3] += 1

            elif (i != j) and (self.all_atype[i] == 2) and (self.all_atype[j] == 2):
              self.dist_anion[idt] += dr2
              cnum[idt][4] += 1

            elif (self.all_atype[i] != self.all_atype[j]):
              self.dist_ca[idt] += dr2
              cnum[idt][5] += 1
        tnum += 1

    print(" @@@ confNum %7d " % (tnum))

    for idt in range(len(self.dt_list)):
      if ( cnum[idt][0] != sum(cnum[idt][1:6]) ):
        print(idt, cnum[idt,:])
      # self.lamda_app[idt] /= cnum[idt][0]
      # self.self_cation[idt] /= cnum[idt][1]
      # self.self_anion[idt] /= cnum[idt][2]
      # self.dist_cation[idt] /= cnum[idt][3]
      # self.dist_anion[idt] /= cnum[idt][4]
      # self.dist_ca[idt] /= cnum[idt][5]
      self.lamda_app[idt] /= tnum
      self.self_cation[idt] /= tnum
      self.self_anion[idt] /= tnum
      self.dist_cation[idt] /= tnum
      self.dist_anion[idt] /= tnum
      self.dist_ca[idt] /= tnum

  def Write(self):

    osfile = open(self.ofilehead + "_compose_msd.dat", 'w')

    s = "t lamda_app self_cation self_anion dist_cation dist_anion dist_ca\n"
    osfile.write(s)

    for idt, dt in enumerate(self.dt_list):
      s = "%12.6f %16.9f %16.9f %16.9f %16.9f %16.9f %16.9f\n"% \
                                    (dt*self.timeStep, \
                                     self.lamda_app[idt], \
                                     self.self_cation[idt], \
                                     self.self_anion[idt], \
                                     self.dist_cation[idt], \
                                     self.dist_anion[idt], \
                                     self.dist_ca[idt])
      osfile.write(s)
    osfile.close()

    self.etime = time.time() - self.stime
    print("Total using time %12.6f" % (self.etime))



if __name__ == "__main__":

  sysId = sys.argv[1]
  idir = sys.argv[2]

  infilename = idir + "/" + sysId + "_lammps_unwrapped.dump"
  ofilehead  = idir + "/" + sysId

  print("Start MSD alpha")
  print(infilename)

  msd = sMSD(infilename, ofilehead)
  msd.Read_Data()
  msd.Get_Time_List()
  msd.Calculate()
  msd.Write()


