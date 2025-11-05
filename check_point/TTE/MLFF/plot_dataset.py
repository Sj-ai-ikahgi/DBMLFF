import os
import sys
from os.path import join

import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import numpy as np

def calculate_metrics(true, pred):
    diff = true - pred
    rmse = np.sqrt(np.mean(diff**2))
    mae = np.mean(np.abs(diff))
    r2 = 1 - (np.sum(diff**2) / np.sum((true - np.mean(true))**2))
    return rmse, mae, r2

# 初始化画布
fig, ax = plt.subplots(1, 2, figsize=(8, 4.5), constrained_layout=True)

# 读取能量数据
with open("./fread_dfeat/energyL.pred.tot", 'r') as ifile:
    lines = ifile.readlines()

energy = [[], []]
for line in lines:
    token = line.split()
    energy[0].append(float(token[0]))  # DFT Energy
    energy[1].append(float(token[1]))  # ML-FF Energy

energy = np.array(energy)
eng_rmse, eng_mae, eng_r2 = calculate_metrics(energy[0], energy[1])

# 读取力数据
with open("./fread_dfeat/forceL.pred.all", 'r') as ifile:
    lines = ifile.readlines()

force = [[], []]
for line in lines:
    token = line.split()
    force[0].append(float(token[0]))  # DFT Force
    force[1].append(float(token[1]))  # ML-FF Force

force = np.array(force)
force_rmse, force_mae, force_r2 = calculate_metrics(force[0], force[1])

# 绘制能量对比图
ax[0].plot(energy[0], energy[1], 'o', ms=3, label=f"RMSE = {eng_rmse:.3f} eV")
ax[0].text(0.05, 0.90, f"MAE = {eng_mae:.3f} eV\nR² = {eng_r2:.3f}", transform=ax[0].transAxes, fontsize=9)

# 绘制力对比图
ax[1].plot(force[0], force[1], '<', ms=3, label=f"RMSE = {force_rmse:.3f} eV/Å")
ax[1].text(0.05, 0.90, f"MAE = {force_mae:.3f} eV/Å\nR² = {force_r2:.3f}", transform=ax[1].transAxes, fontsize=9)

# 绘制1:1参考线
for i, data in enumerate([energy, force]):
    spoint = np.min(data)
    epoint = np.max(data)
    delt = epoint - spoint
    spoint -= delt * 0.1
    epoint += delt * 0.1
    x = [spoint, epoint]
    y = [spoint, epoint]
    ax[i].plot(x, y, 'r--')

    # 设置坐标轴范围
    ax[i].set_xlim([spoint, epoint])
    ax[i].set_ylim([spoint, epoint])

# 设置坐标轴标签
ax[0].set_xlabel("DFT Etot (eV)")
ax[1].set_xlabel("DFT Atomic Force (eV/Å)")
ax[0].set_ylabel("ML-FF Etot (eV)")
ax[1].set_ylabel("ML-FF Atomic Force (eV/Å)")

# 添加图例
ax[0].legend(loc='upper left', fontsize=9)
ax[1].legend(loc='upper left', fontsize=9)

# 设置标题
plt.suptitle("MLFF Fitting Results")

# 保存图像
output_file = "./eng_force_dataset.png"
fig.savefig(output_file, dpi=300)
print(f"Figure saved to {output_file}")