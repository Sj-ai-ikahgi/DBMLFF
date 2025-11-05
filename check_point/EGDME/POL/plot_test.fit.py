import matplotlib.pyplot as plt
import numpy as np

dft = []
predict = []

# 读取数据
with open("test.fit", 'r') as isfile:
    for line in isfile:
        token = line.split()
        dft.append(float(token[0]))
        predict.append(float(token[1]))

dft = np.array(dft)
predict = np.array(predict)

# 计算 RMSE
rmse = np.sqrt(np.mean((dft - predict)**2))

# 计算 MAE
mae = np.mean(np.abs(dft - predict))

# 计算 R²
ss_total = np.sum((dft - np.mean(dft))**2)  # 总平方和
ss_residual = np.sum((dft - predict)**2)    # 残差平方和
r2 = 1 - (ss_residual / ss_total)

# 绘制图表
fig, ax = plt.subplots(figsize=(5, 5), constrained_layout=True)

# 绘制散点图
ax.plot(dft, predict, "o", ms=3, label="RMSE = %.3f" % rmse)

# 添加 1:1 参考线
spoint = np.min(dft)
epoint = np.max(predict)
delt = epoint - spoint
spoint -= delt * 0.1
epoint += delt * 0.1
ax.plot([spoint, epoint], [spoint, epoint], 'r--')

# 设置坐标轴范围
ax.set_xlim(spoint, epoint)
ax.set_ylim(spoint, epoint)

# 设置坐标轴标签
ax.set_xlabel("DFT", fontsize=12)
ax.set_ylabel("Predict", fontsize=12)

# 添加文本框显示 MAE 和 R²
ax.text(0.05, 0.85, f"MAE = {mae:.3f}\nR² = {r2:.3f}", transform=ax.transAxes, fontsize=10)

# 添加图例
plt.legend(loc="upper left")

# 设置标题
ax.set_title("EA Polarization Fitting Result")

# 保存图像
plt.savefig('plot_test.fit.png')
