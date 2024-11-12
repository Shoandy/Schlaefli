# Run on Python3: python3 plot_qn.py

""" Plotting q_n(x) """

import numpy as np
import qn_cfs # or alternatively qn
import matplotlib.pyplot as plt

n = int(input("Enter n: "))
x = np.linspace(n-1.,n+1.,2001)

y = []
for i in x:
    y.append(qn_cfs.schlf(n,i)) # or alternatively qn.schlf(n,i)

plt.plot(x,y)
plt.xlabel('x')
plt.ylabel(r'$q_n(x)$')
plt.title(f'n = {n}')
plt.show()
#plt.savefig(f"q_{n}.pdf")
