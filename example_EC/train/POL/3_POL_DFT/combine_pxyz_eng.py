import numpy as np 
import sys


def Read_All_Pxyz(ifilename):

  confNum = 0
  oneConfLine = 0

  isfile = open(ifilename, 'r')

  all_lines = []

  while True:
    token = isfile.readline()
    if token == "":
      break
    else:
      if "---------------" in token:
        confNum += 1
        lines = []
        lines.append(token)
        for i in range(2):
          token = isfile.readline()
          lines.append(token)
        token = isfile.readline()
        ln = int(token.split()[0])
        lines.append(token)
        for i in range(ln):
          token = isfile.readline()
          lines.append(token)
        oneConfLine = 1+2+1+ ln
        all_lines.append(lines)

  isfile.close()

  return confNum, all_lines, oneConfLine       


def Read_All_Eng(ifilename):

  confNum = 0 

  isfile = open(ifilename, 'r')

  all_lines = []

  while True:
    token = isfile.readline()
    if token == "":
      break
    else:
      if "E_tot(eV)" in token:
        confNum += 1
        all_lines.append(token)
    
  isfile.close()

  return confNum, all_lines


def Write_All(ofilename, confNum, all_pxyz,\
  oneConfLine, all_eng):

  osfile = open(ofilename, 'a')

  for cId in range(confNum):
    for ln in range(oneConfLine):
      osfile.write(all_pxyz[cId][ln])
    osfile.write(all_eng[cId])

  osfile.close()


ifilename = sys.argv[1]
ofilename = sys.argv[2]

confNum_pxyz, all_pxyz, oneConfLine = Read_All_Pxyz(ifilename + "_pxyz/polarization.out")
confNum_eng, all_eng = Read_All_Eng(ifilename + "_dft/polarization.out")


if confNum_pxyz != confNum_eng:
  print("@@@  confNum_pxyz not equal confNum_eng!")
  print("@@@  confNum_pxyz %7d" % confNum_pxyz)
  print("@@@  confNum_eng  %7d" % confNum_eng)
  sys.exit()

print(ifilename, confNum_pxyz)

Write_All(ofilename, confNum_pxyz, all_pxyz, oneConfLine, all_eng)



# ofilename = "./temp.out"
# osfile = open(ofilename, 'w')
# for l in all_lines:
#   osfile.write(l)
# osfile.close()

# print(ifilename, confNum)

