import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import r2_score, mean_absolute_error

# 读取数据
dft = []
predict = []

with open("E_dft_dbmlff.txt", 'r') as isfile:
    for line in isfile:
        tokens = line.split()
        if len(tokens) >= 2:
            dft.append(float(tokens[0]))
            predict.append(float(tokens[1]))

# 转换为numpy数组
dft = np.array(dft)
predict = np.array(predict)

# 计算统计指标
rmse = np.sqrt(np.mean((dft - predict)**2))
mae = mean_absolute_error(dft, predict)
r2 = r2_score(dft, predict)

print(f"RMSE: {rmse:.4f}")
print(f"MAE: {mae:.4f}")
print(f"R²: {r2:.4f}")

# 绘图
fig, ax = plt.subplots(figsize=(6, 6), constrained_layout=True)

# 散点图
ax.scatter(dft, predict, s=20, alpha=0.7, label=f'Data points (n={len(dft)})')

# 理想拟合线
spoint = min(np.min(dft), np.min(predict))
epoint = max(np.max(dft), np.max(predict))
delt = epoint - spoint
spoint -= delt * 0.1
epoint += delt * 0.1

ax.plot([spoint, epoint], [spoint, epoint], 'r--', linewidth=2, label='Ideal fit')

# 设置坐标轴范围和标签
ax.set_xlim(spoint, epoint)
ax.set_ylim(spoint, epoint)
ax.set_xlabel("DFT Energy", fontsize=12)
ax.set_ylabel("DBMLFF Energy", fontsize=12)

# 添加统计信息文本框
stats_text = f'RMSE = {rmse:.4f}\nMAE = {mae:.4f}\nR² = {r2:.4f}'
ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

ax.set_title("DFT vs DBMLFF Interaction Energy Comparison")
plt.legend(loc="lower right")
plt.grid(True, alpha=0.3)

plt.savefig('E_dft_dbmlff_comparison.png', dpi=300, bbox_inches='tight')
plt.show()