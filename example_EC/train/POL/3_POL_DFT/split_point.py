import sys
import os
import numpy as np 



def Get_Point_All(ifilename):

  isfile = open(ifilename, 'r')

  pnum_list = []
  cpoint_list = []

  while True:
    token = isfile.readline()
    if token == "":
      break
    else:
      token = token.split()
      pnum = int(token[0])
      pnum_list.append(pnum)
      for p in range(pnum):
        token = isfile.readline()
        cpoint_list.append(token)
      
  isfile.close()

  return pnum_list, cpoint_list


def Write_Point(ofilename, pnum_list, cpoint_list, prange):

  osfile = open(ofilename, 'w')

  # print(prange[1]-prange[0])
  for i in range(prange[0], prange[1]):
    osfile.write("%3d %5d \n" % (pnum_list[i], i+1))
    for p in range(pnum_list[i]):
      osfile.write(cpoint_list[i+p])



# ----- split the file start

odir = "./split_point_result"
if not os.path.isdir(odir):
  os.system("mkdir -p " + odir)
else:
  os.system("rm " + odir + "/*")


# missing point will no more then 5
sysId_list = ["2A", "4A", "6A", "2A.2A", "2A.4A", "2A.6A"]
split_num_list = [3, 3, 3, 4, 4, 5]
confNum_list = []

for si, sysId in enumerate(sysId_list):

  split_num = split_num_list[si]

  pnum_list, cpoint_list = Get_Point_All("../2_gen_points/point." + sysId)

  pnum = len(pnum_list)

  # print(pnum, split_num)

  ssize = int(pnum/split_num)
  confNum_list.append(ssize)

  prange = np.zeros([2,split_num], dtype=int)

  for i in range(split_num):
    prange[0,i] = ssize*i
    prange[1,i] = ssize*(i+1)

  prange[1,-1] = pnum

  for i in range(split_num):
    ofilename = odir + "/point." + sysId + "." + str(i+1)
    Write_Point(ofilename, pnum_list, cpoint_list, prange[:,i])

# ----- split the file end


# ----- set the variable in qsub_*.sh start

print("### Subsitute the variable in qsub.sh and merge_polarization.sh \
 with the following value! ### \n")

print(" sys_list=(", end="")
for sysId in sysId_list:
  print("%s " % (sysId), end="")
print(")")

print(" split_num=(", end="")
for split_num in split_num_list:
  print("%d " % (split_num), end="")
print(")")

print(" confNum_list=(", end="")
for confNum in confNum_list:
  print("%d " % (confNum), end="")
print(")\n")

os.system("cp -p qsub_dft_temp.sh qsub_dft.sh")
os.system("cp -p qsub_pxyz_temp.sh qsub_pxyz.sh")
os.system("cp -p merge_polarization_temp.sh merge_polarization.sh")

sys_list = "("
for sysId in sysId_list:
  sys_list += "%s " % (sysId)
sys_list += ")"
os.system('sed -i \"s!@sys_list!' + sys_list + '!g\" qsub_dft.sh')
os.system('sed -i \"s!@sys_list!' + sys_list + '!g\" qsub_pxyz.sh')
os.system('sed -i \"s!@sys_list!' + sys_list + '!g\" merge_polarization.sh')

osplit_num_list = "("
for sn in split_num_list:
  osplit_num_list += "%s " % sn
osplit_num_list += ")"
os.system('sed -i \"s!@split_num_list!' + osplit_num_list + '!g\" qsub_dft.sh')
os.system('sed -i \"s!@split_num_list!' + osplit_num_list + '!g\" qsub_pxyz.sh')
os.system('sed -i \"s!@split_num_list!' + osplit_num_list + '!g\" merge_polarization.sh')

oconfNum_list = "("
for confNum in confNum_list:
  oconfNum_list += "%d " % (confNum)
oconfNum_list += ")"
os.system('sed -i \"s!@confNum_list!' + oconfNum_list + '!g\" qsub_dft.sh')
os.system('sed -i \"s!@confNum_list!' + oconfNum_list + '!g\" qsub_pxyz.sh')

# ----- set the variable in qsub_*.sh end


  # sys_list=@sys_list
  # split_num_list=@split_num_list
  # confNum_list=@confNum_list

