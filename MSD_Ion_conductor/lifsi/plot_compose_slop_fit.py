
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# m -> A: 10-10
# s -> ps: 10-12
# A**2/ps -> m**2/s: 10-8

def func(x, a, b):
  return a*x + b


MD_label = ["haha"]
out_label = [""]

color_list = [(31/255,119/255,180/255),\
              (255/255,127/255,14/255),\
              (255/255,127/255,255/255)]

markers = ['s', 'o', '*', '<', '>', 'd']

molTypeNum = len(MD_label)

colNum = 6
msd = [ [[] for i in range(colNum)] 
        for j in range(molTypeNum)]
dt_list = [[] for i in range(molTypeNum)]


label_list = ['lamda_app', 'self_cation', 'self_anion', 'dist_cation', 'dist_anion', 'dist_ca']
system = "20LiFSI-24DME-60TTE"

ofname = "./data_out/" + system

for ml in range(molTypeNum):

  ifilename = "data_out/" + system + "_DBMLFF_compose_msd.dat"

  isfile = open(ifilename, 'r')
  lines = isfile.readlines()
  isfile.close()
  
  for line in lines[1:35]:
    token = line.split()
    dt_list[ml].append(float(token[0])*0.1)
    for col in range(1,colNum+1):
      msd[ml][col-1].append(float(token[col]))

msd = np.array(msd)

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(7,3.5), constrained_layout=True)

D_list = []
s_list = []
colNum = 3
for i in range(colNum):
  xx = dt_list[0]
  yy = msd[0][i]
  tplot = ax[0].plot(xx, yy, marker=markers[0], mfc="None", \
    linestyle="-", ms=5, lw=0.5, label=label_list[i])

  xx = xx[-8:]
  yy = msd[0][i][-8:]
  popt, pcov = curve_fit(func, xx, yy)

  x = np.linspace(xx[0],xx[-1], 100)
  y = func(x,*popt)
  ax[0].plot(x,y,'--', lw=2, color=tplot[0].get_color())

  s_list.append(popt[0])

slop_alpha = s_list[0]/(s_list[1] + s_list[2])
print(" @@@ slop alpha %12.6f " % ( slop_alpha ))
ax[0].text(0.05,0.75, " slop_alpha = %6.3f" % slop_alpha, transform=ax[0].transAxes)

# xx = dt_list[0][:35]
# yy = msd[0][0][:35]/(msd[0][1][:35] + msd[0][2][:35])

xx = dt_list[0][:]
yy = msd[0][0][:]/(msd[0][1][:] + msd[0][2][:])

print(" @@@ min(alpha) % 12.6f %12.6f" % (min(yy), np.mean(yy[-6:-1])))

ax[1].plot(xx, yy, marker=markers[0], mfc="None", \
  linestyle="-", ms=5)


# x = []
# y = []

# sp = 18
# for ep in range(sp+3, len(msd[0][0])):
#   xx = dt_list[0][sp:ep]
#   yy = msd[0][0][sp:ep]
#   popt, pcov = curve_fit(func, xx, yy)
#   slop_0 = popt[0]
#   # print(" @@@ 0 slopt %12.6f " % (popt[0]))

#   xx = dt_list[0][sp:ep]
#   yy = msd[0][1][sp:ep]
#   popt, pcov = curve_fit(func, xx, yy)
#   slop_1 = popt[0]
#   # print(" @@@ 1 slopt %12.6f " % (popt[0]))

#   xx = dt_list[0][sp:ep]
#   yy = msd[0][2][sp:ep]
#   popt, pcov = curve_fit(func, xx, yy)
#   slop_2 = popt[0]
#   # print(" @@@ 2 slopt %12.6f " % (popt[0]))

#   # print(" @@@ alpha %12.6f" % (slop_0/(slop_1 + slop_2)))

#   x.append(dt_list[0][ep])
#   y.append(slop_0/(slop_1+slop_2))

# ax[1].plot(x, y, 'r--', marker=markers[0], mfc="None", \
#   linestyle="-", ms=5, label='slop')

# print(" @@@ min(alpha), slop % 12.6f" % (min(y)))



ax[0].tick_params(which='both',top=True,right = True)
ax[0].legend(loc='upper left', ncol=2, fontsize=8)

# ax.set_title("apparent")

ax[0].set_xlabel("t (ps)")
ax[0].set_ylabel("$\sum_i \sum_j(r_i[t] - r_i[0])(r_j[t]-r_j[0]) $")
ax[1].set_xlabel("t (ps)")
ax[1].set_ylabel("$\\alpha (t)$")

# ax.set_xlim(0,50)
# ax.set_ylim(-0.3,1)


#plt.legend(bbox_to_anchor=(0.5,0.95),bbox_transform=ax.transAxes, ncol=2)
#plt.legend(title='DBMLFF DFT', loc='upper left', ncol=2, fontsize=9)
#plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)
plt.suptitle( system + " @298")
plt.savefig(ofname + "_compose.png", dpi=600)



# isfile = open("./data_out/D_list.dat", 'w')
# for D in D_list:
#   isfile.write("%12.7f %12.6f\n" % (D, np.log(D)))
# isfile.close()

# # --- mean and std

# D_list = np.array(D_list)
# logD_list = np.log(D_list)

# x = D_list[0:10]
# mean_x = np.mean(x)
# std_x = np.std(x)
# print(" @@@ D_Li mean %12.6f " % (mean_x))
# print(" @@@ D_Li std  %12.6f " % (std_x))

# x = logD_list[0:10]
# mean_x = np.mean(x)
# std_x = np.std(x)
# print(" @@@ logD_Li mean %12.6f " % (mean_x))
# print(" @@@ logD_Li std  %12.6f " % (std_x))


# x = D_list[10:20]
# print(len(x))
# mean_x = np.mean(x)
# std_x = np.std(x)
# print(" @@@ D_PF6 mean %12.6f " % (mean_x))
# print(" @@@ D_PF6 std  %12.6f " % (std_x))

# x = logD_list[10:20]
# mean_x = np.mean(x)
# std_x = np.std(x)
# print(" @@@ logD_PF6 mean %12.6f " % (mean_x))
# print(" @@@ logD_PF6 std  %12.6f " % (std_x))




# # --- Ionic conductivity

def Ionic_Conductivity(c_e,kB,T,V,N,q,D):
  return c_e/(kB*T*V)*(N*q**2*D)


c_e = 1.6021766208      # C
kB = 8.617343*10**-5    # eV/K
T = 300                 # K
Vol = pow(27.815314,3.0)   # A**3


coeff = c_e/(kB*T*Vol)
# print(" @@@ coeff %15.6f " % (coeff) )


# for di,D in enumerate([D_list[0]]):
#   kai = Ionic_Conductivity(c_e,kB,T,Vol,10,1,D)
#   # one 10 is from exponential,
#   # the other is from S/m -> mS/cm
#   kai *= 10*10
#   print(" @@@ Li_total     %12.6f (mS/cm)" % (kai))

# kai = 0.0
# for di,D in enumerate(D_list[1:11]):
#   kai += Ionic_Conductivity(c_e,kB,T,Vol,1,1,D)
# kai *= 10*10
# print(" @@@ Li_seperate  %12.6f (mS/cm)" % (kai))


# for di,D in enumerate([D_list[11]]):
#   kai = Ionic_Conductivity(c_e,kB,T,Vol,10,1,D)
#   kai *= 10*10
#   print(" @@@ PF6_total    %12.6f (mS/cm)" % (kai))

# kai = 0.0
# for di,D in enumerate(D_list[12:22]):
#   kai += Ionic_Conductivity(c_e,kB,T,Vol,1,1,D)
# kai *= 10*10
# print(" @@@ PF6_seperate %12.6f (mS/cm)" % (kai))


# for di,D in enumerate([4.03]):
#   kai_1 = Ionic_Conductivity(c_e,kB,T,Vol,10,1,D)
# kai_1 *= 10*10
# print(" @@@ Li+ %12.6f (mS/cm)" % (kai_1))

# for di,D in enumerate([5.06]):
#   kai_2 = Ionic_Conductivity(c_e,kB,T,Vol,10,1,D)
# kai_2 *= 10*10
# print(" @@@ PF6- %12.6f (mS/cm)" % (kai_2))

# # print(" @@@ total %12.6f (mS/cm)" % ((kai_1 + kai_2)*0.5))



