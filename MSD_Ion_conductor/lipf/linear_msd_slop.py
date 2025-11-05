import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# 单位转换系数
# m -> A: 10^-10
# s -> ps: 10^-12
# A²/ps -> m²/s: 10^-8

def func(x, a, b):
    return a*x + b

MD_label = ["Li", "PF6"]
out_label = ["$Li^{+}$", "$PF6^{-}$"]

color_list = [(31/255,119/255,180/255),
              (255/255,127/255,14/255),
              (255/255,127/255,255/255)]

markers = ['s', 'o', '*', '<', '>', 'd']

molTypeNum = len(MD_label)

colNum = 3
msd = [ [[] for i in range(colNum)] 
        for j in range(molTypeNum)]
dt_list = [[] for i in range(molTypeNum)]


system = "80EC-10LiPF6-60DMC"

ofname = "./data_out/" + system

for ml in range(molTypeNum):
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
    
    for cId in range(1):
        y = msd[mtId][cId]

        sp = 1
        ep = sp + 50

        xx = x[sp:ep]
        yy = y[sp:ep]

        popt, pcov = curve_fit(func, xx, yy)
        yy = func(xx, *popt)
        D = popt[0]/6.0 * 100
        D_rounded = round(D, 3)  # 保留三位小数


        ax[col].plot(x, y, marker=markers[mtId], mfc="None", 
            linestyle="-", ms=5, label="$D= %3.3f \\times 10^{-10}\ m^2/s$" % (D_rounded))

        xx = x[sp:ep]
        yy = y[sp:ep]

        popt, pcov = curve_fit(func, xx, yy)
        yy = func(xx, *popt)
        D = popt[0]/6.0 * 100
        D_rounded = round(D, 3)  # 保留三位小数

        ax[col].plot(xx, yy, 'r--')

        print(" %5s D=%10.3f (*10^-10 m^2/s) %12.6f  " % (MD_label[mtId], D_rounded, np.log(D)))
        D_list.append(D_rounded)


    ax[col].tick_params(which='both',top=True,right = True)
    ax[col].legend(loc='upper left', ncol=2, fontsize=9)
    col += 1


ax[0].set_title("$Li^{+}$")
ax[1].set_title("$PF6^{-}$")

ax[0].set_xlabel("$t\ (ps)$")
ax[1].set_xlabel("$t\ (ps)$")
ax[0].set_ylabel("$MSD\ ({\AA}^2/s)$")

plt.suptitle( system + " @298K")
plt.savefig(ofname + "_msd_linear.png", dpi=600)


# 离子电导率计算参数
c_e = 1.6021766208      # C
kB = 8.617343*10**-5    # eV/K
T = 298                 # K
Vol = pow(26.589,3.0)   # A³

def Ionic_Conductivity(c_e,kB,T,V,N,q,D):
    return c_e/(kB*T*V)*(N*q**2*D)


# 计算并输出离子电导率，保留三位小数
D = D_list[0]
N = 10
kai = Ionic_Conductivity(c_e,kB,T,Vol,N,1,D)
kai *= 10*10  # 单位转换
print(" @@@ Li_total     %10.3f (mS/cm)" % (round(kai, 3)))


D = D_list[1]
N = 10
kai = Ionic_Conductivity(c_e,kB,T,Vol,N,1,D)
kai *= 10*10  # 单位转换
print(" @@@ PF6_total    %10.3f (mS/cm)" % (round(kai, 3)))
