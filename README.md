# Schlaefli

Computation of the Schläfli function

FILES DESCRIPTION: 

qn.f90 Fortran program computes qn(x) for x in [n-1,n+1] and n = 4,5,...,10000 
with quadruple precision.

qn.py Python program computes qn(x) for x in [n-1,n+1] and n = 4,5,...,10000 
with double precision.

cheb_cfs.f90 Fortran program computes Chebyshev coefficients for n = 4,5,...,10000
and saves them in cfs.dat (in the first column is the corresponding n).

qn_cfs.py Python program computes qn(x) by reading Chebyshev coefficients from cfs.dat

plot_qn.py plots qn(x) for given n = 4,5,...,10000 and x in [n-1,n+1]. 

err_qn.py plots the absolute difference between q_n(x) computed in quadruple and double
precisions for each n separately and x in [n-1,n+1] (import modules: qn_cfs or qn and used data from qnx.dat computed in qnx.f90)

Remark: One can take n > 10000 with the cost of precision and computation time.
