import os
import sys
from os.path import join

import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt

import numpy as np

fig, ax = plt.subplots(1, 2, figsize=(8, 4.5), constrained_layout=True)

# 读取能量数据
ifile = open("./fread_dfeat/energyL.pred.tot", 'r')
lines = ifile.readlines()
ifile.close()

energy = [[],[]]
for line in lines:
    token = line.split()
    energy[0].append(float(token[0]))  # DFT Energy
    energy[1].append(float(token[1]))  # ML-FF Energy

energy = np.array(energy)
diff = energy[0] - energy[1]

# 计算 RMSE, MAE 和 R² (能量)
eng_rmse = np.sqrt(np.dot(diff, diff) / len(diff))
eng_mae = np.mean(np.abs(diff))
eng_r2 = 1 - (np.sum(diff**2) / np.sum((energy[0] - np.mean(energy[0]))**2))

# 读取力数据
ifile = open("./fread_dfeat/forceL.pred.all", 'r')
lines = ifile.readlines()
ifile.close()

force = [[],[]]
for line in lines:
    token = line.split()
    force[0].append(float(token[0]))  # DFT Force
    force[1].append(float(token[1]))  # ML-FF Force

force = np.array(force)
diff = force[0] - force[1]

# 计算 RMSE, MAE 和 R² (力)
force_rmse = np.sqrt(np.dot(diff, diff) / len(diff))
force_mae = np.mean(np.abs(diff))
force_r2 = 1 - (np.sum(diff**2) / np.sum((force[0] - np.mean(force[0]))**2))

# 绘制能量对比图
ax[0].plot(energy[0], energy[1], 'o', ms=3, label="RMSE = %.3f eV" % eng_rmse)
ax[0].text(0.05, 0.85, "MAE = %.3f eV\nR² = %.3f" % (eng_mae, eng_r2), transform=ax[0].transAxes, fontsize=9)

# 绘制力对比图
ax[1].plot(force[0], force[1], '<', ms=3, label="RMSE = %.3f eV/Å" % force_rmse)
ax[1].text(0.05, 0.85, "MAE = %.3f eV/Å\nR² = %.3f" % (force_mae, force_r2), transform=ax[1].transAxes, fontsize=9)

# 绘制1:1参考线
spoint_0 = np.min(energy)
epoint_0 = np.max(energy)
delt = epoint_0 - spoint_0
spoint_0 -= delt * 0.1
epoint_0 += delt * 0.1
x = [spoint_0, epoint_0]
y = [spoint_0, epoint_0]
ax[0].plot(x, y, 'r--')

spoint_1 = np.min(force)
epoint_1 = np.max(force)
delt = epoint_1 - spoint_1
spoint_1 -= delt * 0.1
epoint_1 += delt * 0.1
x = [spoint_1, epoint_1]
y = [spoint_1, epoint_1]
ax[1].plot(x, y, 'r--')

# 设置坐标轴范围
ax[0].set_xlim([spoint_0, epoint_0])
ax[0].set_ylim([spoint_0, epoint_0])
ax[1].set_xlim([spoint_1, epoint_1])
ax[1].set_ylim([spoint_1, epoint_1])

# 设置坐标轴标签
ax[0].set_xlabel("DFT Etot (eV)")
ax[1].set_xlabel("DFT Atomic Force (eV/Å)")
ax[0].set_ylabel("ML-FF Etot (eV)")
ax[1].set_ylabel("ML-FF Atomic Force (eV/Å)")

# 添加图例
ax[0].legend(loc=(0.05, 0.75), fontsize=9)
ax[1].legend(loc=(0.05, 0.75), fontsize=9)

# 设置标题
plt.suptitle("MLFF Fitting Results")

# 保存图像
fig.savefig("./eng_force_dataset.png", dpi=300)
