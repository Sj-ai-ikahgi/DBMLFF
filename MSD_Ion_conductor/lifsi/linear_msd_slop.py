import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def func(x, a, b):
  return a*x + b

# 修改点1：将PF6改为FSI
MD_label = ["Li", "FSI"]  # 将"PF6"改为"FSI"
out_label = ["$Li^{+}$", "$FSI^{-}$"]  # 将"$PF6^{-}$"改为"$FSI^{-}$"

color_list = [(31/255,119/255,180/255),\
              (255/255,127/255,14/255),\
              (255/255,127/255,255/255)]

markers = ['s', 'o', '*', '<', '>', 'd']

molTypeNum = len(MD_label)

colNum = 3
msd = [ [[] for i in range(colNum)] 
        for j in range(molTypeNum)]
dt_list = [[] for i in range(molTypeNum)]


system = "20LiFSI-24DME-60TTE"

ofname = "./data_out/" + system

for ml in range(molTypeNum):

  # 修改点2：文件名中的PF6改为FSI
  ifilename = "data_out/" + system + "_DBMLFF_" + MD_label[ml] + "_msd.dat"

  isfile = open(ifilename, 'r')
  lines = isfile.readlines()
  isfile.close()
  
  for line in lines[1:]:
    token = line.split()
    dt_list[ml].append(float(token[0])*100*0.001)
    for col in range(1,colNum+1):
      msd[ml][col-1].append(float(token[col]))

dt_list = np.array(dt_list)
msd = np.array(msd)


fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(7,3.5), constrained_layout=True)

D_list = []


for mtId in range(molTypeNum):

  x = dt_list[mtId]

  col = mtId
  # for cId in range(0,colNum):
  for cId in range(1):
    y = msd[mtId][cId]

    sp = 5
    ep = sp + 50

    xx = x[sp:ep]
    yy = y[sp:ep]

    popt, pcov = curve_fit(func, xx, yy)
    yy = func(xx, *popt)
    D = popt[0]/6.0 * 100

    ax[col].plot(x, y, marker=markers[mtId], mfc="None", \
        linestyle="-", ms=5, label="$D= %3.2f *10^{-10}\ m^2/s$" % (D))

    xx = x[sp:ep]
    yy = y[sp:ep]

    popt, pcov = curve_fit(func, xx, yy)
    yy = func(xx, *popt)
    D = popt[0]/6.0 * 100

    ax[col].plot(xx, yy, 'r--')

    print(" %5s D=%12.7f (*10^-10 m^2/s) %12.6f  " % (MD_label[mtId], D, np.log(D)))
    D_list.append(D)


  ax[col].tick_params(which='both',top=True,right = True)
  ax[col].legend(loc='upper left', ncol=2, fontsize=9)

  col += 1

# 修改点3：标题中的PF6改为FSI
ax[0].set_title("$Li^{+}$")
ax[1].set_title("$FSI^{-}$")  # 将"$PF6^{-}$"改为"$FSI^{-}$"

ax[0].set_xlabel("$t\ (ps)$")
ax[1].set_xlabel("$t\ (ps)$")
ax[0].set_ylabel("$MSD\ ({\AA}^2/s)$")

plt.suptitle( system + " @298K")
plt.savefig(ofname + "_msd_linear.png", dpi=600)

# --- Ionic conductivity

c_e = 1.6021766208      # C
kB = 8.617343*10**-5    # eV/K
T = 298                 # K
Vol = pow(27.815314,3.0)   # A**3

def Ionic_Conductivity(c_e,kB,T,V,N,q,D):
  return c_e/(kB*T*V)*(N*q**2*D)

# 修改点4：变量名和输出信息中的PF6改为FSI
D_Li = D_list[0]
N_Li = 20  # 假设有20个Li离子

kai_Li = Ionic_Conductivity(c_e,kB,T,Vol,N_Li,1,D_Li)
kai_Li *= 10  # 转换为mS/cm (S/m * 10 = mS/cm)
print(" @@@ Li+ conductivity     %12.6f (mS/cm)" % (kai_Li))

D_FSI = D_list[1]
N_FSI = 20  # 假设有20个FSI离子

kai_FSI = Ionic_Conductivity(c_e,kB,T,Vol,N_FSI,1,D_FSI)
kai_FSI *= 10  # 转换为mS/cm
print(" @@@ FSI- conductivity     %12.6f (mS/cm)" % (kai_FSI))

# 总离子电导率是两者的平均值
kai_total = (kai_Li + kai_FSI) / 2
print(" @@@ Total ionic conductivity %12.6f (mS/cm)" % (kai_total))