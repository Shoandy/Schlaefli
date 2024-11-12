# Run on Python3: python3 qn_cfs.py

""" Compute function q_n(x) defined as schlf(n, x) """

import numpy as np
import time
from itertools import islice

""" INITIALISATION AND FUNCTIONS DEFINITION: """

N = 24 # Number of Chebyshev coefficients
""" N = 54 gives minimal Chebyshev truncation error for n = 4 and 5,
which is of order |a_N| which is just above machine accuracy ~ 1.9e-34 for
quadruple precision. """

# Data: key = dimension of Euclidean space, value = Chebyshev coefficients a_i:
dic = {}

def qn(y, a):
    """ Compute Chebyshev expansion of q_n(y) """
    r = [0] * (N+2) # Initialize array
    r[N+1] = 0
    r[N] = 0
    for i in range(N-1, 0, -1):
        r[i] = 2*y*r[i+1] - r[i+2] + a[i]
    q = y*r[1] - r[2] + a[0]/2
    return q

""" COMPUTATION OF q_n(x): """

def schlf(n, x):
    """ q_n(x) function, related to the Schläfli function f_n(x) as follows:
    f_n(x) = 2**n*sqrt(n)*(n/2)!*(x-n+1)**((n-1)/2)/(pi**(n/2)*n!*n!)*q_n(x)
    for even n, and
    f_n(x) = sqrt(n)*(x-n+1)**((n-1)/2)/(pi**((n-1)/2)*n!*((n-1)/2)!)*q_n(x)
    for odd n.
    Here n is positive integer (Euclidean space dimension) and x is real
    from the interval [n-1, n+1] """
    _n = int(n)
    if _n != n:
        raise ValueError("n must be positive integer")
    elif n < 1:
        raise ValueError("n must be positive integer")
    elif n == 1:
        return 1

    _x = float(x)
    if _x != x:
        raise ValueError("x must be real")
    elif x < n-1. or x > n+1.:
        raise ValueError(f"x must lie in [{n-1}, {n+1}]")
    elif x == n - 1.:
        return 1

    y = x - n # Variable change: -1 <= y <= 1

    if n == 2: # q_2(x) is known:
        q = np.arccos(1/(2+y)) / np.sqrt(2+2*y)
    elif n == 3: # q_3(x) is known:
        q = 2*np.sqrt(3)*(np.arccos(1/(3 + y))-np.pi/3) / (1+y)
    else: # q_n(x) with n > 4:
        if n in dic.keys(): # Check if Chebyshev coefficients were computed
            a = dic[n] # Read Chebyshev coefficients from dictionary
        else: # Read Chebyshev coefficients:
            f = open("cfs.dat", "r")
            a = []
            lines = list(islice(f,(n-4)*N,(n-3)*N))
            for line in lines:
                a.append(np.float128(line.split()[1]))
                f.close()
        dic[n] = a # Record Chebyshev coefficients
        q = qn(y, a)
    # Correct digits for double precision maximal roundoff error ~ n*1.e-16:
    sf = int(16-np.floor(np.log10(n)))
    return round(q, sf)

def main():
    n = int(input("Dimension n = "))
    x = float(input("x = "))
    start = time.process_time()
    qnx = schlf(n,x)
    end = time.process_time()
    print(f"q_{n}({x}) = {qnx}")
    print(f"Comutation time = {end-start:0.5f}sec")

if __name__ == "__main__":
    main()
