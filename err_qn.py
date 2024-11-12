# Run on Python3: python3 err_qn.py

""" Plotting the absolute difference between q_n(x) computed
in quadruple and double precisions """

import linecache
import numpy as np
import qn_cfs # or alternatively qn
import matplotlib.pyplot as plt

n = int(linecache.getline("qnx.dat",1))

x = np.linspace(n-1.,n+1.,20001)

y = []
for i in x:
    y.append(qn_cfs.schlf(n,i)) # or alternatively qn.schlf(n,i)

xf = []
yf = []

for i in range(2,20003):
    l = linecache.getline("qnx.dat",i).split()
    xf.append(np.float128(l[0]))
    yf.append(np.float128(l[1]))

dy = []

for i in range(0,len(x)):
    dy.append(abs(y[i]-yf[i]) / (n*1.e-16)) # or (n*1.e-14) for qn.schlf(n,i)

z = np.full((20001,),1.0)

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(x,dy)

ax.set_xlabel('x')
ax.set_ylabel(r'$\delta q_n(x)/(n\times 10^{-16})$') # or (n\times 10^{-14})
fig.suptitle(f'n = {n}')
plt.show()
#fig.savefig(f"error_q_{n}.pdf")
