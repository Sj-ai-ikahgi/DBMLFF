import numpy as np 
import sys

def Read_All_Pol(ifilename):

  confNum = 0 

  isfile = open(ifilename, 'r')

  all_lines = []
  lines = []

  while True:
    token = isfile.readline()
    if token == "":
      break
    elif "E_tot(eV) " in token:
      confNum += 1
      lines.append(token)
      for l in lines:
        all_lines.append(l)
    else:
      if "---------------" in token:
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
    
  isfile.close()

  return confNum, all_lines       


ifilename = sys.argv[1]

confNum, all_lines = Read_All_Pol(ifilename)

ofilename = "./temp.out"
osfile = open(ofilename, 'w')
for l in all_lines:
  osfile.write(l)
osfile.close()


print(ifilename, confNum)


  # rm -f polarization.out

  # pnum_list=(3 3 3 4 4)

  # i=0
  # for cut in "2A" #"4A" "6A" "2A.2A" "2A.4A"
  # do
  #   pnum=${pnum_list[$i]}
  #   i=$((${i} + 1))
  #   for p in `seq 1 ${pnum}`
  #   do
  #     echo $cut $p
  #     cat ${cut}_${p}/polarization.out >> polarization.out
  #   done
  # done

  # mv polarization.out ../fitting

