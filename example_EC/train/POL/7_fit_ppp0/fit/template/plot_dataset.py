import os
import sys
from os.path import join

import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt

import numpy as np


system = sys.argv[1]

fig, ax = plt.subplots(1, 1, figsize=(5, 4), constrained_layout=True)


ifile = open("./temp_folder/fitting/fread_dfeat/energyL.pred.tot", 'r')
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


# ifile = open("./temp_folder/fitting/fread_dfeat/forceL.pred.all", 'r')
# lines = ifile.readlines()
# ifile.close

# force = [[],[]]
# for line in lines:
#     token = line.split()
#     force[0].append(float(token[0]))
#     force[1].append(float(token[1]))

# force = np.array(force)
# diff = force[0] - force[1]
# force_rmse = np.sqrt(np.dot(diff,diff)/len(diff))

ax.plot(energy[0], energy[1], 'o', ms=3, label="energy rmse = %.3e" % eng_rmse)
# ax[1].plot(force[0], force[1], '<', ms=3, label="force rmse = %.3e" % force_rmse)

# ax[0].text(-0.3,0.6, "energy rmse = %.3e" % eng_rmse)
# ax[1].text(-2.5,2.5, "force rmse = %.3e" % force_rmse)

ax.set_xlabel("DBMLFF")
ax.set_ylabel("predict")

ax.legend()
# ax[0].text(-3928, -3928, "energy rmse = %.3e" % eng_rmse)
# ax[1].text(-20.5,20, "force rmse = %.3e" % force_rmse)


ax.set_title("energy " + system)
# ax[1].set_title("force")
# plt.suptitle("DMC MLFF")

fig.savefig("./eng_force_dataset.png", dpi=300)



