  #!/usr/bin/env python3
  # -*- coding: utf-8 -*-

import numpy as np

#f_npfile='data_scaler.npy'
f_npfile='Wij.npy'
dsnp = np.load(f_npfile, allow_pickle=True)
#print('test',f_npfile,dsnp.dtype,dsnp.shape)
#print('shape',dsnp.shape)
#print('dim',dsnp.ndim)
print('size',dsnp.size)
n=dsnp.size
count = 0
while count <= n-1:
  a=dsnp[count]
  m1=a.shape[0]
  m2=a.shape[1]
  print('m12=',m1,m2)
  for i in range(0,m1):
    for j in range(0,m2):
      print(i,j,a[i,j])
#  print(a)
  count = count+1

