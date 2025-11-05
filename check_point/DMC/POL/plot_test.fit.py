import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# 读取数据
dft = []
predict = []

with open("test.fit", 'r') as isfile:
    for line in isfile:
        token = line.split()
        dft.append(float(token[0]))
        predict.append(float(token[1]))

dft = np.array(dft)
predict = np.array(predict)

# 计算 RMSE, MAE 和 R²
diff = dft - predict
rmse = np.sqrt(np.mean(diff**2))
mae = np.mean(np.abs(diff))
ss_total = np.sum((dft - np.mean(dft))**2)  # 总平方和
ss_residual = np.sum(diff**2)               # 残差平方和
r2 = 1 - (ss_residual / ss_total)

# 创建图形
fig, ax = plt.subplots(figsize=(5, 5), constrained_layout=True)

# 设置刻度方向为朝内，并调整字体大小和刻度线宽度
ax.tick_params(axis='both', which='major', direction='in', labelsize=12, width=2)
ax.tick_params(axis='both', which='minor', direction='in', width=2)

# 绘制散点图
ax.scatter(dft, predict, c='#5392CE', s=10, alpha=0.7)

# 添加 1:1 参考线
spoint = np.min(dft)
epoint = np.max(predict)
delt = epoint - spoint
spoint -= delt * 0.1
epoint += delt * 0.1
ax.plot([spoint, epoint], [spoint, epoint], 'r--', linewidth=2)

# 设置坐标轴范围
ax.set_xlim(spoint, epoint)
ax.set_ylim(spoint, epoint)

# 设置坐标轴标签
ax.set_xlabel("DFT", fontsize=16)
ax.set_ylabel("Predict", fontsize=16)

# 加粗坐标轴线条
for spine in ax.spines.values():
    spine.set_linewidth(2)

# 移除图例
ax.legend().remove()

# 直接添加统计信息文本而不使用文本框边框
textstr = f"RMSE = {rmse:.3f} eV\nMAE = {mae:.3f} eV\nR² = {r2:.3f}"
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=16, verticalalignment='top')

# 设置标题
ax.set_title("EC Polarization Fitting Result", fontsize=20, fontweight='bold')

# 保存图像
plt.savefig('plot_test.fit.png', dpi=600, bbox_inches='tight')