import matplotlib.pyplot as plt
import numpy as np


dft = []
predict = []

isfile = open("test.fit", 'r')
while True:
  token = isfile.readline()
  if token == "":
    break
  else:
    token = token.split()
    dft.append(float(token[0]))
    predict.append(float(token[1]))
isfile.close()

rmse = 0.0
for i in range(len(dft)):
  # if dft[i] < -0.6:
  #   print(i, dft[i])
  rmse += (dft[i] - predict[i])*(dft[i] - predict[i])

rmse = np.sqrt(rmse/len(dft))


fig, ax = plt.subplots(figsize=(5,5), constrained_layout=True)
  
ax.plot(dft, predict, "o", ms=3, label = "rmse = %6.3f" % rmse)

spoint = np.min(dft)
epoint = np.max(predict)
delt = epoint - spoint
spoint -= delt*0.1
epoint += delt*0.1

ax.plot([spoint,epoint], [spoint,epoint], 'r--')

ax.set_xlim(spoint, epoint)
ax.set_ylim(spoint, epoint)

ax.set_xlabel("DFT",fontsize=12)
ax.set_ylabel("predict",fontsize=12)

plt.legend(loc="upper left")

ax.set_title("DMC polarization fitting result")
plt.savefig('plot_test.fit.png')



