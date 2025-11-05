import os
from os.path import join

import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt

import numpy as np

md_type = 'dbmlff'
max_time = 60000

etot_range_length = 50
ep_range_length = 50
t_range_length = 200
avet_range_length =50

figsize = (16, 8)

# --------------------------------------------------

with open(join('./', 'MDSTEPS')) as f:
    data = f.read().splitlines()

temp = {
    'time' : [], 
    'Etot' : [], 
    'Ep' : [], 
    'Ek' : [], 
    'T' : [], 
    'aveT' : [], 
    'Fcheck' : [], 
}

for i in data:
    split = i.split()
    temp['time'].append(float(split[1]))
    temp['Etot'].append(float(split[3]))
    temp['Ep'].append(float(split[4]))
    temp['Ek'].append(float(split[5]))
    temp['T'].append(float(split[7]))
    temp['aveT'].append(float(split[9]))
    if md_type == 'dbmlff':
        temp['Fcheck'].append(float(split[13]))
    else:
        temp['Fcheck'].append(float(split[-1]))

data = temp

fig, ax = plt.subplots(5, 1, figsize=figsize, constrained_layout=True)

ax[0].plot(data['time'], data['Fcheck'], 'o', ms=1)
ax[1].plot(data['time'], data['Etot'])
ax[2].plot(data['time'], data['Ep'])
ax[3].plot(data['time'], data['T'])
ax[4].plot(data['time'], data['aveT'])

ax[0].set_xlabel('Time (fs)')
ax[0].set_ylabel('Fcheck')
ax[1].set_xlabel('Time (fs)')
ax[1].set_ylabel('Etot (eV)')
ax[2].set_xlabel('Time (fs)')
ax[2].set_ylabel('Ep (eV)')
ax[3].set_xlabel('Time (fs)')
ax[3].set_ylabel('T (K)')
ax[4].set_xlabel('Time (fs)')
ax[4].set_ylabel('aveTemp (K)')

stime = 0
max_time = data['time'][-1] 
# max_time = 20000

ax[0].set_xlim(stime, max_time)
ax[0].set_ylim(-1, 3)
ax[1].set_xlim(stime, max_time)
ax[1].set_ylim(np.array(data['Etot']).mean()-etot_range_length//2, np.array(data['Etot']).mean()+etot_range_length//2)
ax[2].set_xlim(stime, max_time)
ax[2].set_ylim(np.array(data['Ep']).mean()-ep_range_length//2, np.array(data['Ep']).mean()+ep_range_length//2)
ax[3].set_xlim(stime, max_time)
ax[3].set_ylim(np.array(data['T']).mean()-t_range_length//2, np.array(data['T']).mean()+t_range_length//2)
ax[4].set_xlim(stime, max_time)
ax[4].set_ylim(np.array(data['aveT']).mean()-avet_range_length//2, np.array(data['aveT']).mean()+avet_range_length//2)

fig.savefig('./basic_info.png', dpi=300)




