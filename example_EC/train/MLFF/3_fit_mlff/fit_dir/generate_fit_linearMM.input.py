import os
import sys
import numpy as np


def collectAllSourceFiles():
    '''
    搜索工作文件夹，得到所有MOVEMENT文件的路径，并将之存储在 sourceFileList中
    '''
    sourceFileList = []
    for path,dirList,fileList in os.walk("./"):
        if "MOVEMENT" in fileList:
            sourceFileList.append(os.path.abspath(path))
    return sourceFileList


def Write_fit_linearMM_input(natype, m_neigh, natoms, sysNum, sysName_list, weight_list):

  osfile = open("./fit_linearMM.input", 'w')
  osfile.write("%7d, %7d, %7d,      1, 0, 1 \
    ! natype, m_neigh, natoms, bond, vdw, pair\n" % (natype, m_neigh, natoms))
  for i in range(1,natype+1):
    osfile.write("%-7d, 0.0, 0.0            ! itype, rad_atom,wp_atom\n" % (i))
  osfile.write("0.0, 2.0, 0.01, 0.00001  ! E_weight ,Etot_weight, F_weight, delta\n")
  osfile.write("4.0  ! dwidth for energy_group \n")

  osfile.write("%-7d     ! sysNum\n" % (sysNum))
  for i in range(sysNum):
    osfile.write("%-5.3f  ! %s\n" % (weight_list[i], sysName_list[i]))


natype =  int(sys.argv[1])
m_neigh = int(sys.argv[2])
natoms =  int(sys.argv[3])

weight_map = {"400K":1.0, "700K": 1.0, "1200K":1.0, \
              "400K_s": 1.0, "600K_s":0.2, "1900K_s":0.02}


sourceFileList = collectAllSourceFiles()
sysNum = len(sourceFileList)

sysName_list = []
weight_list = []
for s in sourceFileList:
  sysName = s.split("/")[-1]
  sysName_list.append(sysName)
  weight_list.append(weight_map[sysName])


Write_fit_linearMM_input(natype, m_neigh, natoms, sysNum, sysName_list, weight_list)








