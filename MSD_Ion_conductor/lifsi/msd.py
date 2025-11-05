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
        for i in range(natoms):
          if atype[i] == 0:
            atype[i] = 1
          elif atype[i] == 1:
            atype[i] = 2
        self.all_atype = deepcopy(atype)

      if cid < 3000: continue

      self.all_coords.append(unwrapped_coords)
    isfile.close()

    if len(self.all_coords) < 1:
      print("There is no coordiante!")
      sys.exit()

    self.total_conf_num = np.shape(self.all_coords)[0]
    self.half_conf_num = int(self.total_conf_num*0.5)

    atype_list, self.each_atype_num = np.unique(atype, return_counts=True)
    self.atype_num = len(atype_list)
    self.natoms = np.shape(self.all_coords)[1]

    print(self.total_conf_num)
    print(np.shape(self.all_coords))

  def Get_Time_List(self):
    point = np.linspace(1,self.total_conf_num*0.5, self.dt_num)

    for i,p in enumerate(point):
      tx = int(p)
      if len(self.dt_list) != 0:
        if tx != self.dt_list[-1]:
          self.dt_list.append(tx)
      else:
        self.dt_list.append(tx)

  def Calculate(self):
    self.tmsd = np.zeros([self.atype_num,len(self.dt_list)], dtype=float)
    self.pmsd = np.zeros([self.natoms, len(self.dt_list)], dtype=float)

    tnum = np.zeros([self.dt_num,self.atype_num], dtype=int)
    pnum = np.zeros([self.dt_num,self.natoms], dtype=int)

    for idt,dt in enumerate(self.dt_list):
      print(idt, dt)
      for sid in range(0,self.total_conf_num-dt,1):
        eid = sid + dt
        scoord = self.all_coords[sid]
        ecoord = self.all_coords[eid]
        for i in range(self.natoms):
          dr = 0.0
          for j in range(3):
            dx = ecoord[i][j] - scoord[i][j]
            dr += dx*dx
          t = self.all_atype[i] - 1
          self.tmsd[t][idt] += dr
          self.pmsd[i][idt] += dr
          tnum[idt][t] += 1
          pnum[idt][i] += 1
    
    for idt in range(self.dt_num):
      for t in range(self.atype_num):
        self.tmsd[t][idt] /= tnum[idt][t]
      for i in range(self.natoms):
        self.pmsd[i][idt] /= pnum[idt][i]

  def Write(self):
    LiNum = self.each_atype_num[0]

    osfile = open(self.ofilehead + "_Li_msd.dat", 'w')
    s = "t  Li "
    for ti in range(LiNum):
      s += " Li_%d " % (ti+1)
    s += "\n"
    osfile.write(s)

    for idt, dt in enumerate(self.dt_list):
      s = "%12.6f %16.9f "%(dt*self.timeStep, self.tmsd[0][idt])
      for i in range(LiNum):
        s += "%16.9f "% self.pmsd[i][idt]
      s += "\n"
      osfile.write(s)
    osfile.close()

    # 修改点1：将PF6改为FSI (文件名)
    osfile = open(self.ofilehead + "_FSI_msd.dat", 'w')
    # 修改点2：将PF6改为FSI (表头)
    s = "t  FSI "
    FSI_Num = self.each_atype_num[1]  # 修改点3：变量名重命名
    for ti in range(FSI_Num):  # 修改点4：使用新变量名
      s += " FSI_%d " % (ti+1)  # 修改点5：原子标签改为FSI
    s += "\n"
    osfile.write(s)

    for idt, dt in enumerate(self.dt_list):
      s = "%12.6f %16.9f "%(dt*self.timeStep, self.tmsd[1][idt])
      for i in range(LiNum, LiNum+FSI_Num):  # 修改点6：使用新变量名
        s += "%16.9f "% self.pmsd[i][idt]
      s += "\n"
      osfile.write(s)
    osfile.close()

    self.etime = time.time() - self.stime
    print("Total using time %12.6f" % (self.etime))

if __name__ == "__main__":
  sysId = sys.argv[1]
  idir = sys.argv[2]

  infilename = idir + "/" + sysId + "_lammps_unwrapped.dump"
  ofilehead  = idir + "/" + sysId

  print("Start MSD")
  print(infilename)

  msd = sMSD(infilename, ofilehead)
  msd.Read_Data()
  msd.Get_Time_List()
  msd.Calculate()
  msd.Write()