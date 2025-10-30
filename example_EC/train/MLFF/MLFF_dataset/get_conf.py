import sys
import os
import argparse
import time
from copy import deepcopy
import numpy as np
from aux_dbmlff import *




class Configuration:


  def __init__(self,
               sysName,
               infile_format,
               infilename,
               idir,
               odir,
               sconfId = 1,
               econfId = 1e12,
               confDelt = 1,
               molNum_list = [],
               molSize_list = []):

    self.sysName = sysName
    self.infile_format = infile_format
    self.infilename = infilename
    self.idir = idir
    self.odir = odir
    self.sconfId = sconfId
    self.econfId = econfId
    self.confDelt = confDelt

    self.molNum_list = molNum_list
    self.molSize_list = molSize_list

    self.isfile = '' 

    # the number of configuration write out
    self.confId = 0
    self.headline = []
    #$ this for output sorted *.config and *.pwmat    
    # o2pId[oid] = pid
    self.o2pId = []
    #$ this for molId
    # p2oId[pid] = oid
    self.p2oId = []
    self.molId = []
    
    self.init_data()

    self.stime = time.time()

    # for check
    if not os.path.isdir("./check_conf"):
      os.system("mkdir ./check_conf")

    self.lammps_dump_fname = "./check_conf/" + sysName + "_" + \
      infile_format + "_MLFF_dataset.dump"
    os.system("rm -f " + self.lammps_dump_fname)

    self.large_distrotion_fname = "./check_conf/" + sysName + "_" + \
      infile_format + "_large_distortion.dat"
    os.system("rm -f " + self.large_distrotion_fname)
    self.osfile_large_distrotion = open(self.large_distrotion_fname, 'w')

    self.bond_length_fname = "./check_conf/" + sysName + "_" + \
      infile_format + "_bond_length.dat"
    os.system("rm -f " + self.bond_length_fname)
    self.osfile_bond_length = open(self.bond_length_fname, 'w')

    pre_coords = []
    unwrapped_coords = []

    self.total_DMCNum = 0


  def init_data(self):

    self.ntimestep = 0.0
    self.Ep = 0.0
 
    self.natoms = 0
    self.lattice = []
    self.eindex = []
    self.coords = []
    self.forces = []
    self.velocity = []


  def Get_atomId(self):

    self.o2pId = np.repeat(-1,self.natoms)
    self.p2oId = np.repeat(-1,self.natoms)
    self.molId = np.zeros([self.natoms], dtype=int)

    # set the mapping of atom id between the ORIGIN and PRESENT configuration
    if self.infile_format == "DBMLFF":
      print("\n------ process DBMLFF format file ------\n")
      self.o2pId = np.arange(0,self.natoms)
      self.p2oId = np.arange(0,self.natoms)
    else:
      print("\n------ process PWMAT format file ------\n")

      tisfile = open(self.idir + "/ORIGIN.INDEX", 'r')

      tnum = 0
      while True:
        token = tisfile.readline()
        if token == "":
          break
        else:
          token = token.split()
          oid = int(token[0]) - 1
          pid = int(token[4]) - 1
          self.o2pId[oid] = pid
          self.p2oId[pid] = oid
      tisfile.close()

    # set the molecular id of each atom
    i = 0
    mid = 0
    for mtype in range(len(self.molNum_list)):
      for molNum in range(self.molNum_list[mtype]):
        mid += 1
        for aid in range(self.molSize_list[mtype]):
          pid = self.o2pId[i]
          self.molId[pid] = mid
          i += 1

    # result check
    for i in range(self.natoms):
      if self.p2oId[i] == -1:
        print("haha read p2oId wrong")
        sys.exit()



  def Get_Mol(self,mid_list):

    teindex = []
    tcoords = []
    tforces = []
    tvelocity = []

    natoms = 0
    for m in range(len(mid_list)):
      mid = mid_list[m]
      for i in range(self.natoms):
        pid = self.o2pId[i]
        if self.molId[pid] == mid:
          natoms += 1
          teindex.append(self.eindex[pid])
          tcoords.append(self.coords[pid])
          tforces.append(self.forces[pid])
          tvelocity.append(self.velocity[pid])

    return natoms, teindex, tcoords, tforces, tvelocity


  def Get_ref_bond_length(self):

    isfile = open("./parameters/bond.molecule", "r")

    token = isfile.readline().split(',')
    self.bondNum = int(token[0])
    self.atomPair = np.zeros([self.bondNum,2], dtype=int)
    self.bondLength = np.zeros(self.bondNum, dtype=float)

    for ib in range(self.bondNum):
      token = isfile.readline().split()
      self.atomPair[ib][0] = int(token[1])
      self.atomPair[ib][1] = int(token[2])
      self.bondLength[ib] = float(token[3])

    isfile.close()


  def Get_bond_length(self,lattice,coords):

    for ib in range(self.bondNum):
      i = self.atomPair[ib][0] - 1
      j = self.atomPair[ib][1] - 1
      dr = np.zeros(3, dtype=float)
      for k in range(3):
        dr[k] = (coords[i][k] - coords[j][k])*lattice[k][k]
      dist = np.sqrt(np.dot(dr,dr))

      if abs(dist - self.bondLength[ib]) / self.bondLength[ib] > 0.2:
        print(" @@@ large distrotion mol %7d %7d %12.6f %12.6f" % \
          (self.total_DMCNum, ib, dist, self.bondLength[ib]))
        self.osfile_large_distrotion.write(" %7d \n" % self.total_DMCNum)

      self.osfile_bond_length.write(" %12.6f " % (dist))

    self.osfile_bond_length.write("\n")
    

  def Get_Move_Atom_Order(self):

    isfile = open("./parameters/move_atom_order.dat", "r")

    token = isfile.readline().split()
    tnatoms = int(token[0])
    self.mao = []
    token = isfile.readline().split(',')
    for i in range(tnatoms):
      self.mao.append(int(token[i]))

    isfile.close()

    print("@@@ move atom order is: ")
    print("   ", self.mao)


  def Move_Mol(self, icoords):

    mol_type_list = [0]
    mol_atom_list = self.molSize_list

    mol_num = len(mol_type_list)

    ocoords = []

    for mid in range(mol_num):
      sp = sum(mol_atom_list[0:mid])
      ep = sum(mol_atom_list[0:(mid+1)])
      tcoords = icoords[sp:ep]

      for ii in range(mol_atom_list[mid]-1):
        caId = self.mao[ii] - 1
        i = self.mao[ii+1] - 1
        for j in range(3):
          dr = tcoords[i][j] - tcoords[caId][j]
          if np.fabs(dr) > 0.5:
            tcoords[i][j] = tcoords[i][j] - np.sign(dr)
      
      for i in range(mol_atom_list[mid]):
        ocoords.append(tcoords[i])

    icoords = ocoords



  def Write_Mol(self):
    ofilehead =  self.odir + "/DFT_"

    for m in range(self.molNum_list[0]):

      mid_list = [m+1]

      self.total_DMCNum += 1

      natoms, teindex, tcoords, tforces, tvelocity = self.Get_Mol(mid_list)

      self.Move_Mol(tcoords)
      lattice, tcoords = Move_Enlarge_Box(self.lattice[0][0], self.lattice, natoms, tcoords)

      DBMLFF_MOVEMENT_headline = "%7d     1 %3d   1 \n" % (natoms, natoms)

      to2pId = np.arange(natoms)


      self.Get_bond_length(self.lattice, tcoords)


      osfile = open(ofilehead + str(self.total_DMCNum) + ".config", 'w')
      Write_MLFF_Init_Config(osfile, \
                          DBMLFF_MOVEMENT_headline, \
                          lattice, \
                          to2pId, \
                          teindex, \
                          tcoords, \
                          tforces, \
                          tvelocity) 
      osfile.close()

      elem_list, eindex_list, emass_list, natype, atype = eindex2atype(teindex)
      coords = Lamda2Cart(lattice, tcoords)
      Write_LAMMPS_Dump(self.total_DMCNum, self.lammps_dump_fname, \
        natoms, lattice, atype, coords)



  def Process_One_MOVEMENT_File(self):

    self.Get_Move_Atom_Order()

    self.isfile = open(self.idir + "/" + self.infilename, 'r')

    if not os.path.isdir(self.odir):
      os.system("mkdir -p " + self.odir)
  
    cand_list = np.arange(self.sconfId,self.econfId,self.confDelt)
    # print("cand_list ", cand_list)

    # tconfId is the total read in configuration;
    # while self.confId is the total configuration will further proceed
    tconfId = 0
    while True:
      tconfId += 1
      # whether there is new movement
      self.headline = New_MOVEMENT_Conf(self.isfile) 
      if self.headline == '': break

      # clear the configuration data and load a new configuration data
      self.init_data()
      self.ntimestep, self.Ep, self.natoms, self.lattice, self.eindex, \
      self.coords, self.forces, self.velocity, _ =  \
        Get_One_MOVEMENT_Conf(self.headline, self.isfile)


      # read over the configuration
      if tconfId < self.sconfId:
        if tconfId % 1000 == 0:
          print('-'*10 + str(tconfId) + ' - ' + \
                str(self.ntimestep) + " read over " + '-'*10)
        continue
      elif tconfId > self.econfId: break

      # we will deal with those configuration that we wanted
      if tconfId not in cand_list: continue

      # given the information of progress
      if (tconfId-1) % 1000 == 0:
        print("-"*10 + " start output " + str(tconfId) + "-"*10)

      self.confId += 1


      # 1) set the atomic mapping between DFT and DBMLFF,
      # when we using a atom.config of DBMLFF as the input configuration of PWMAT,
      # the atomic order in the MOVEMENT of PWMAT will change with repect to atom.config.
      # 2) set the molId of each atom
      if self.confId == 1:
        self.Get_atomId()
        self.Get_ref_bond_length()


      # ----------- pick up some molecular -----------

      self.Write_Mol()

    self.isfile.close()
    self.osfile_large_distrotion.close()
    self.osfile_bond_length.close()
    self.etime = time.time() - self.stime

    print("Total output configuration number is %9d" % (self.confId))
    print("Total using time %12.6f\n\n" % (self.etime))


if __name__ == "__main__":

  narg = len(sys.argv)


  sysName = sys.argv[1]
  infile_format = "DFT"
  infilename    = "MOVEMENT"
  idir = sys.argv[2]
  odir = sys.argv[3]

  molNum = int(sys.argv[4])
  molSize = int(sys.argv[5])
  confDelt = int(sys.argv[6])

  print("sId           ", sysName)
  print("infile format ", infile_format)
  print("infile dir    ", idir)
  print("outfile dir   ", odir)

  sconfId = 1
  econfId = 40000

  molNum_list = [molNum]

  molSize_list = [molSize]

  conf = Configuration(sysName,
                      infile_format,
                      infilename,
                      idir = idir,
                      odir = odir,
                      sconfId = sconfId,
                      econfId = econfId,
                      confDelt = confDelt,
                      molNum_list = molNum_list,
                      molSize_list = molSize_list)

  conf.Process_One_MOVEMENT_File()


