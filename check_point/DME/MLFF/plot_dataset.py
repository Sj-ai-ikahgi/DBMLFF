import os
import sys
from os.path import join

import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt

import numpy as np

fig, ax = plt.subplots(1, 2, figsize=(8, 4.5), constrained_layout=True)


ifile = open("./fread_dfeat/energyL.pred.tot", 'r')
lines = ifile.readlines()
ifile.close

energy = [[],[]]
for line in lines:
    token = line.split()
    energy[0].append(float(token[0]))
    energy[1].append(float(token[1]))

energy = np.array(energy)
diff = energy[0] - energy[1]
eng_rmse = np.sqrt(np.dot(diff,diff)/len(diff))

ifile = open("./fread_dfeat/forceL.pred.all", 'r')
lines = ifile.readlines()
ifile.close

force = [[],[]]
for line in lines:
    token = line.split()
    force[0].append(float(token[0]))
    force[1].append(float(token[1]))

force = np.array(force)
diff = force[0] - force[1]
force_rmse = np.sqrt(np.dot(diff,diff)/len(diff))

# The difference eng_rmse with the AEM_Etot
# is caused by weight_system in calc_lin_forceMM.f 
ax[0].plot(energy[0], energy[1], 'o', ms=3, label="energy rmse = %.3e" % eng_rmse)
ax[1].plot(force[0], force[1], '<', ms=3, label="force rmse = %.3e" % force_rmse)


spoint_0 = np.min(energy)
epoint_0 = np.max(energy)
delt = epoint_0 - spoint_0
spoint_0 -= delt*0.1
epoint_0 += delt*0.1

x = [spoint_0, epoint_0]
y = [spoint_0, epoint_0]
ax[0].plot(x, y, 'r--')

spoint_1 = np.min(force)
epoint_1 = np.max(force)
delt = epoint_1 - spoint_1
spoint_1 -= delt*0.1
epoint_1 += delt*0.1
x = [spoint_1, epoint_1]
y = [spoint_1, epoint_1]
ax[1].plot(x, y, 'r--')


ax[0].set_xlim([spoint_0, epoint_0])
ax[0].set_ylim([spoint_0, epoint_0])
ax[1].set_xlim([spoint_1, epoint_1])
ax[1].set_ylim([spoint_1, epoint_1])


ax[0].set_xlabel("DFT")
ax[1].set_xlabel("DFT")
ax[0].set_ylabel("DBMLFF")

ax[0].legend(loc=(0.05,0.85), fontsize=9)
ax[1].legend(loc=(0.05,0.85), fontsize=9)

ax[0].set_title("energy")
ax[1].set_title("force")
plt.suptitle("MLFF fitting results")

fig.savefig("./eng_force_dataset.png", dpi=300)



