import os
import sys
from os.path import join

import matplotlib as mpl
mpl.use('Agg')  # 非交互式后端，适合服务器环境
from matplotlib import pyplot as plt

import numpy as np

# -------------------------- 核心修改：从目标文件读取数据 --------------------------
# 读取能量数据文件（格式：第1列=DFT能量（已反转），第2列=DB能量（已反转））
energy_data = np.loadtxt("energyL.pred.tot")
dft_energy_inverted = energy_data[:, 0]  # 提取DFT能量列
db_energy_inverted = energy_data[:, 1]   # 提取DB能量列

# 读取力数据文件（格式：第1列=DFT力，第2列=DB力）
force_data = np.loadtxt("forceL.pred.all")
dft_force = force_data[:, 0]  # 提取DFT力列
db_force = force_data[:, 1]   # 提取DB力列
# --------------------------------------------------------------------------------

# 初始化画布（2个子图，保持原尺寸和布局）
fig, ax = plt.subplots(1, 2, figsize=(8, 4.5), constrained_layout=True)

# -------------------------- 能量误差计算（逻辑不变） --------------------------
energy_diff_inverted = dft_energy_inverted - db_energy_inverted
eng_rmse_inverted = np.sqrt(np.dot(energy_diff_inverted, energy_diff_inverted) / len(energy_diff_inverted))  # RMSE
eng_mae_inverted = np.mean(np.abs(energy_diff_inverted))  # MAE
ss_res_eng = np.sum(energy_diff_inverted ** 2)  # 残差平方和
ss_tot_eng = np.sum((dft_energy_inverted - np.mean(dft_energy_inverted)) ** 2)  # 总平方和
eng_r2_inverted = 1 - (ss_res_eng / ss_tot_eng) if ss_tot_eng != 0 else 0  # 避免除零错误

# -------------------------- 力误差计算（逻辑不变） --------------------------
force_diff = dft_force - db_force
force_rmse = np.sqrt(np.dot(force_diff, force_diff) / len(force_diff))  # RMSE
force_mae = np.mean(np.abs(force_diff))  # MAE
ss_res = np.sum(force_diff ** 2)  # 残差平方和
ss_tot = np.sum((dft_force - np.mean(dft_force)) ** 2)  # 总平方和
force_r2 = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0  # 避免除零错误

# -------------------------- 画图逻辑（完全不变） --------------------------
# 绘制能量对比图（含误差指标）
ax[0].plot(dft_energy_inverted, db_energy_inverted, 'o', ms=3, 
           label=f"RMSE = {eng_rmse_inverted:.3f} eV\nMAE = {eng_mae_inverted:.3f} eV\nR² = {eng_r2_inverted:.3f}")
# 绘制力图（含误差指标）
ax[1].plot(dft_force, db_force, '<', ms=3, 
           label=f"RMSE = {force_rmse:.3f} eV/Å\nMAE = {force_mae:.3f} eV/Å\nR² = {force_r2:.3f}")

# 能量图1:1参考线（自适应数据范围）
spoint_0 = np.min(dft_energy_inverted)
epoint_0 = np.max(dft_energy_inverted)
delt = epoint_0 - spoint_0
spoint_0 -= delt * 0.1  # 扩展10%范围，避免点贴边
epoint_0 += delt * 0.1
ax[0].plot([spoint_0, epoint_0], [spoint_0, epoint_0], 'r--')

# 力图1:1参考线（自适应数据范围）
spoint_1 = np.min(dft_force)
epoint_1 = np.max(dft_force)
delt = epoint_1 - spoint_1
spoint_1 -= delt * 0.1
epoint_1 += delt * 0.1
ax[1].plot([spoint_1, epoint_1], [spoint_1, epoint_1], 'r--')

# 禁用能量图坐标偏移量（避免科学计数法干扰）
ax[0].ticklabel_format(useOffset=False)

# 设置坐标轴范围（与参考线范围一致）
ax[0].set_xlim([spoint_0, epoint_0])
ax[0].set_ylim([spoint_0, epoint_0])
ax[1].set_xlim([spoint_1, epoint_1])
ax[1].set_ylim([spoint_1, epoint_1])

# 设置标签、标题和图例
ax[0].set_xlabel("DFT Energy (eV)")
ax[1].set_xlabel("DFT Atomic Force (eV/Å)")
ax[0].set_ylabel("DB Energy (eV)")
ax[1].set_ylabel("DB Atomic Force (eV/Å)")

ax[0].legend(loc=(0.05, 0.7), fontsize=8)  # 图例位置适配多行文本
ax[1].legend(loc=(0.05, 0.7), fontsize=8)

ax[0].set_title("Energy")
ax[1].set_title("Force")
plt.suptitle("results")  # 总标题

# 保存图像（高分辨率，适合论文/报告）
fig.savefig("eng_force_dataset.png", dpi=600, bbox_inches='tight')  # bbox_inches避免标签被截断