
import sys
import os
import numpy as np
import matplotlib.pyplot as plt



def Get_Dbmlff(ifilename,natoms):
  isfile = open(ifilename)

  numConf = 0
  db_eng_list = []
  db_f_list = []
  while True:
    token = isfile.readline()
    if token == "":
      break
    else:
      db_eng_list.append(float(token.split()[0]))
      f = []
      for i in range(natoms):
        token = isfile.readline().split()
        f.append([float(token[j]) for j in range(1,4)])
      db_f_list.append(f)
      numConf += 1
  isfile.close()

  db_eng_list = np.array(db_eng_list)
  db_f_list = np.array(db_f_list)
  return numConf, db_eng_list, db_f_list




isfile = open("./atom.config", 'r')
token = isfile.readline().split()
natoms = int(token[0])
isfile.close()


n123_list = [70,80,90,100,120,140]
n123Num = len(n123_list)
delx = 1.0/100.0

db_eng_list = []
db_force_list = []

for n123 in n123_list:
  ifilename = str(n123) + "/OUT.DB"
  numConf, eng, force = Get_Dbmlff(ifilename,natoms)
  db_eng_list.append(eng)
  db_force_list.append(force)

db_eng_list = np.array(db_eng_list)
db_force_list = np.array(db_force_list)

print(np.shape(db_eng_list), np.shape(db_force_list))




markers = ['s', 'o', '*', '<', '>', 'd']
axis_dir = ['x', 'y', 'z']
color_list = [(31/255,119/255,180/255),\
              (255/255,127/255,14/255),\
              (255/255,127/255,255/255),\
              (255/255,0/255,255/255)]


fig, ax = plt.subplots(1, 1, figsize=(6,4), constrained_layout=True)
for n in range(n123Num):

  y = db_eng_list[n]
  x = [i*delx for i in range(len(y))]
  
  ax.plot(x,y, marker=markers[0], mfc="None", \
    linestyle="-", ms=3, label=str(n123_list[n]))
  
ax.set_ylabel("$E_{DB}$", fontsize=12)
ax.set_xlabel("$displacement\ (\AA)$", fontsize=12)

#ax[0].legend(loc='upper right', ncol=3, fontsize=9)
ax.legend(ncol=6, fontsize=6)
#ax[0].set_title("first point value is " + str(db_eng[0]), fontsize=9)
#ax[1].legend(loc='upper right', ncol=3, fontsize=9)
# ax.legend(ncol=6, fontsize=6)
#ax[1].set_title("first point value is " + str(pol_eng[0]), fontsize=9)

plt.suptitle("DB Energy")
    
plt.savefig("./data_out/egg_box_energy.png", dpi=300)
plt.close()