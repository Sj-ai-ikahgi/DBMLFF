import os
import sys
from os.path import join
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# 设置全局绘图参数
plt.rcParams.update({
    'font.family': 'Arial',  # 优先使用Arial或Helvetica
    'font.size': 12,         # 正文字体大小基准
    'axes.labelsize': 14,    # 坐标轴标签字体大小
    'axes.titlesize': 16,    # 子图标题字体大小
    'xtick.labelsize': 12,   # x轴刻度字体
    'ytick.labelsize': 12,   # y轴刻度字体
    'lines.linewidth': 1.5,  # 线条粗细
    'axes.linewidth': 1.5,   # 坐标轴线宽
    'savefig.dpi': 600,      # 输出分辨率
    'savefig.bbox': 'tight', # 自动裁剪白边
    'legend.frameon': False  # 图例无边框
})

# 读取能量数据
def read_data(filepath):
    data = [[], []]
    with open(filepath, 'r') as f:
        for line in f:
            tokens = line.strip().split()
            data[0].append(float(tokens[0]))  # DFT数据
            data[1].append(float(tokens[1]))  # ML预测数据
    return np.array(data)

energy = read_data("./fread_dfeat/energyL.pred.tot")
force = read_data("./fread_dfeat/forceL.pred.all")

# 计算统计指标
def calc_metrics(y_true, y_pred):
    diff = y_true - y_pred
    rmse = np.sqrt(np.mean(diff**2))
    mae = np.mean(np.abs(diff))
    r2 = 1 - (np.sum(diff**2) / np.sum((y_true - np.mean(y_true))**2))
    return rmse, mae, r2

eng_rmse, eng_mae, eng_r2 = calc_metrics(energy[0], energy[1])
force_rmse, force_mae, force_r2 = calc_metrics(force[0], force[1])

# 创建画布
fig, axs = plt.subplots(1, 2, figsize=(10, 4.5), gridspec_kw={'wspace':0.25})

# 能量对比图
ax = axs[0]
ax.scatter(energy[0], energy[1], s=20, alpha=0.6, 
          c='#2c7bb6', edgecolor='none', label='Data points')
ax.plot([energy.min(), energy.max()], [energy.min(), energy.max()], 
       '--', color='#d7191c', label='Ideal fit')
ax.text(0.05, 0.85, 
       f'RMSE = {eng_rmse:.3f} eV\nMAE = {eng_mae:.3f} eV\n$R^2$ = {eng_r2:.3f}', 
       transform=ax.transAxes, fontsize=12)
ax.set(xlabel='DFT Energy (eV)', ylabel='MLFF Energy (eV)')

# 力对比图
ax = axs[1]
ax.scatter(force[0], force[1], s=20, alpha=0.6, 
          c='#2c7bb6', edgecolor='none')
ax.plot([force.min(), force.max()], [force.min(), force.max()], 
       '--', color='#d7191c')
ax.text(0.05, 0.85, 
       f'RMSE = {force_rmse:.3f} eV/Å\nMAE = {force_mae:.3f} eV/Å\n$R^2$ = {force_r2:.3f}', 
       transform=ax.transAxes, fontsize=12)
ax.set(xlabel='DFT Atomic Force (eV/Å)', ylabel='MLFF Atomic Force (eV/Å)')

# 统一坐标轴范围
for ax in axs:
    lim_min = min(ax.get_xlim()[0], ax.get_ylim()[0])
    lim_max = max(ax.get_xlim()[1], ax.get_ylim()[1])
    padding = (lim_max - lim_min)*0.05
    ax.set_xlim(lim_min - padding, lim_max + padding)
    ax.set_ylim(lim_min - padding, lim_max + padding)
    ax.set_aspect('equal')  # 确保1:1比例
    
# 添加统一图例
fig.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), 
          ncol=2, frameon=False)

# 保存图像
fig.savefig("MLFF_validation.png")
fig.savefig("MLFF_validation.pdf")  # 矢量图格式更佳