
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

slabel = "rho_r4_c"
cnum = int(sys.argv[1])

dr = [ [] for i in range(2)]
rho = [ [[] for i in range(cnum)] for j in range(2) ]

file_list = ["plot.funcr", "plot.funcr_new"]

# ------- read in DB of DBMLFF

for i,ifilename in enumerate(file_list):
  isfile = open("./" + ifilename, 'r')
  lines = isfile.readlines()
  isfile.close()

  for line in lines[:]:
    token = line.split()
    dr[i].append(float(token[0])*0.529177208)
    for c in range(cnum):
      # rho[i][c].append(float(token[1+c]))
      rho[i][c].append(float(token[1+c])*dr[i][-1]**4)

markers = ['s', 'o', '*', '<', '>', 'd']
marker_size = 0.5 
MD_label = ["old", "new"]


for c in range(cnum):
  print(c)
  fig, ax = plt.subplots(1, 1, figsize=(6,4), constrained_layout=True)

  for i in range(2):
    ax.plot(dr[i],rho[i][c], marker=markers[i], mfc="None", \
     linestyle="-", ms=marker_size, label="center_" + str(c+1) + "-" + MD_label[i])

  ax.plot([4.5,4.5],[-0.001,0.001], 'r--',lw=0.5)
  ax.plot([0.3,0.3],[-0.001,0.001], 'r--',lw=0.5)

  ax.set_xlabel("$r\ \AA $", fontsize=12)
  #ax.set_xlim(-0.5,5.0)
  # ax.set_ylim(-0.001,0.001)
  ax.legend(loc='upper right', ncol=1, fontsize=12)

  ax.set_ylabel("$\\rho (r)*r^4$", fontsize=12)
  # ax.set_ylabel("$\\rho (r)$", fontsize=12)


  plt.savefig("./data_out/"+ slabel + "_" + str(c+1) + ".png", dpi=300)


