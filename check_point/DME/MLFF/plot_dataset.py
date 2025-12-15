import os
import sys
from os.path import join

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# 读取能量数据
with open("./energyL.pred.tot", 'r') as ifile:
    lines = ifile.readlines()

energy = [[], []]
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
with open("./forceL.pred.all", 'r') as ifile:
    lines = ifile.readlines()

force = [[], []]
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

# 创建图形
fig, ax = plt.subplots(1, 2, figsize=(12, 6), constrained_layout=True)

# 设置刻度方向为朝内，并调整字体大小和刻度线宽度
for axis in ax:
    axis.tick_params(axis='both', which='major', direction='in', labelsize=12, width=2)
    axis.tick_params(axis='both', which='minor', direction='in', width=2)

# 绘制能量对比图
ax[0].scatter(energy[0], energy[1], c='#5392CE', s=10, alpha=0.7, label="")
ax[0].text(0.05, 0.75, "RMSE = %.3f eV\nMAE = %.3f eV\nR² = %.3f" % (eng_rmse, eng_mae, eng_r2), transform=ax[0].transAxes, fontsize=16, bbox=None)

# 绘制力对比图
ax[1].scatter(force[0], force[1], c='#5392CE', s=10, alpha=0.7, label="")
ax[1].text(0.05, 0.75, "RMSE = %.3f eV/Å\nMAE = %.3f eV/Å\nR² = %.3f" % (force_rmse, force_mae, force_r2), transform=ax[1].transAxes, fontsize=16, bbox=None)

# 绘制1:1参考线
spoint_0 = np.min(energy)
epoint_0 = np.max(energy)
delt = epoint_0 - spoint_0
spoint_0 -= delt * 0.1
epoint_0 += delt * 0.1
x = [spoint_0, epoint_0]
y = [spoint_0, epoint_0]
ax[0].plot(x, y, 'r--', linewidth=2)

spoint_1 = np.min(force)
epoint_1 = np.max(force)
delt = epoint_1 - spoint_1
spoint_1 -= delt * 0.1
epoint_1 += delt * 0.1
x = [spoint_1, epoint_1]
y = [spoint_1, epoint_1]
ax[1].plot(x, y, 'r--', linewidth=2)

# 设置坐标轴范围
ax[0].set_xlim([spoint_0, epoint_0])
ax[0].set_ylim([spoint_0, epoint_0])
ax[1].set_xlim([spoint_1, epoint_1])
ax[1].set_ylim([spoint_1, epoint_1])

# 设置坐标轴标签
ax[0].set_xlabel("DFT Etot (eV)", fontsize=16)
ax[1].set_xlabel("DFT Atomic Force (eV/Å)", fontsize=16)
ax[0].set_ylabel("ML-FF Etot (eV)", fontsize=16)
ax[1].set_ylabel("ML-FF Atomic Force (eV/Å)", fontsize=16)

# 加粗坐标轴线条
for axis in ax:
    for spine in axis.spines.values():
        spine.set_linewidth(2)

# 移除网格线
# 注释掉或删除以下两行来移除网格线
# ax[0].grid(True, linestyle='--', alpha=0.5)
# ax[1].grid(True, linestyle='--', alpha=0.5)

# 移除图例
ax[0].legend().remove()
ax[1].legend().remove()

# 设置标题
plt.suptitle("MLFF Fitting Results", fontsize=20, fontweight='bold')

# 保存图像
fig.savefig("./eng_force_dataset.png", dpi=600, bbox_inches='tight')