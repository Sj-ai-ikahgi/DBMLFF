import os
from os.path import join

import shutil

import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt

import numpy as np

# pymatgen.core.structure.Structure
def get_starting_index_of_element(structure, element):
    for i in range(len(structure)):
        if structure[i].specie.name == element:
            return i
    return 0


# bash
def cp_r(src, des):
    if os.path.isdir(src):
        shutil.copytree(join(src), join(des))
    elif os.path.isfile(src):
        shutil.copy(join(src), join(des))
    else:
        raise ValueError('%s represents neither a directory nor a file' % (src))


def rm_fr(path):
    if not os.path.exists(path):
        pass
    else:
        if os.path.isdir(path):
            shutil.rmtree(join(path))
        elif os.path.isfile(path):
            os.remove(join(path))
        else:
            raise ValueError('%s represents neither a directory nor a file' % (path))

# matplotlib
def plot_45_degree_line(x, y, title, xlabel='x', ylabel='y', fontsize=20, check_outlier=False, threshold=None, relative_threshold=False, highlight_outlier=False):
    fig, ax = plt.subplots(1, 1, figsize=(10, 10), tight_layout=True)
    
    x_info = np.array([[''] for i in range(len(x))])
    if len(x.shape) > 1:
        x_info = x[:,1:]
        x = x[:,0].astype(np.float64)
    
    y_info = np.array([[''] for i in range(len(y))])
    if len(y.shape) > 1:
        y_info = y[:,1:]
        y = y[:,0].astype(np.float64)
    
    max_data = max(x.max(), y.max())
    min_data = min(x.min(), y.min())
    
    max_lim = max_data + (max_data - min_data) * 0.1
    min_lim = min_data - (max_data - min_data) * 0.1
    
    ax.set_xlim(min_lim, max_lim)
    ax.set_ylim(min_lim, max_lim)
    
    ax.set_title(title, fontsize=fontsize)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)
    
    ax.plot(x, y, 'bo')
    ax.plot([min_lim, max_lim], [min_lim, max_lim], 'r-')
    
    if check_outlier and (threshold != None):
        if relative_threshold:
            threshold = threshold * abs(x-y).mean()
        x_outlier = []
        y_outlier = []
        for i in range(len(x)):
            if abs(x[i]-y[i]) > threshold:
                x_outlier.append(x[i])
                y_outlier.append(y[i])
                print('%23.16E ' % (x[i]), end='')
                for j in x_info[i]:
                    print('%s ' % (j), end='')
                print()
                print('%23.16E ' % (y[i]), end='')
                for j in y_info[i]:
                    print('%s ' % (j), end='')
                print()
        
        if highlight_outlier:
            ax.plot(x_outlier, y_outlier, 'ro')
            ax.plot([min_lim, max_lim-threshold], [min_lim+threshold, max_lim], 'r--')
            ax.plot([min_lim+threshold, max_lim], [min_lim, max_lim-threshold], 'r--')
    
    if title == '':
        fig.savefig('45.png')
    else:
        fig.savefig('45_%s.png' % (title))
