# Run on Python3: python3 qn.py

""" Compute function q_n(x) defined as schlf(n, x) """

import numpy as np
import time

""" INITIALISATION AND FUNCTIONS DEFINITION: """

N = 24 # Number of Chebyshev coefficients
""" N = 24 gives minimal Chebyshev truncation error for n = 4 and 5, which is
of order |a_N|, which is just above the machine accuracy ~ 2.2e-16 """

# Chebyshev expansion matrix:
A = np.fromfunction(lambda i, j: np.cos(np.pi/N*i*(j+1/2)), (N, N), dtype = int)

# Data: key = dimension of Euclidean space, value = Chebyshev coefficients a_i:
dic = {}

# Functions:
def gq_2(y):
    """ Function G_4(y)*Q_2(y) """
    q = 18 * np.arccos(1/(2+y)) / ((4+y) * np.sqrt(1+y) * np.sqrt((4+y)**2-1))
    return q

def gq_3(y):
    """ Function G_5(y)*Q_3(y) """
    q = 96 * np.sqrt(5) * (np.arccos(1/(3+y)) - np.pi/3) \
      / ((1+y) * (5+y) * np.sqrt((5+y)**2-1))
    return q

def g(n, y):
    """ Function G_n(y) """
    g = (n-1)**2 * np.sqrt(n*(n-2)) / ((y+n) * np.sqrt((y+n)**2-1))
    return g

def cheb_f(f):
    """ Compute Chebyshev coefficients a_i of a function f(y) """
    dg = np.fromfunction(lambda j: f(A[1,j]), (N, ), dtype = int)
    a = 2/N*np.dot(A,dg)
    return a

def dca(a, c):
    """ Compute Chebyshev coefficients d_i via the recurrence relation """
    d = [0] * N # Initialize array
    d[0] = np.dot(c,a) - 1/2 * c[0] * a[0]
    d[1] = (np.dot(a[:-1], c[1:]) + np.dot(c[:-1], a[1:])) / 2
    for i in range(2, N, 1):
        d[i] = (np.dot(a[:-i], c[i:]) + np.dot(c[:-i], a[i:]) \
             + np.dot(c[1:i], a[i-1:0:-1])) / 2
    return d

def cheb_a(d, n):
    """ Compute Chebyshev coefficients a_i of q_n(y) """
    a = [0] * N # Initialize array
    a[N-1] = d[N-1] / (2*N+n-3)
    a[N-2] = (d[N-2]-4*(N-1)*a[N-1]) / (2*N+n-5)
    for i in range(N-3, -1, -1):
        a[i] = (d[i]-d[i+2] + (n-2*(i+1)-3)*a[i+2] - 4*(i + 1)*a[i+1]) \
             / (2*(i+1)+n-3)
    return a

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
    elif n % 2 == 0: # Even n > 2
        if n in dic.keys(): # Check if Chebyshev coefficients were computed
            a = dic[n] # Read Chebyshev coefficients from dictionary
        else: # Compute Chebyshev coefficients:
            d = cheb_f(gq_2)
            a = cheb_a(d, 4)
            dic[4] = a
        if n > 4:
            if n in dic.keys(): # Check if Chebyshev coefficients were computed
                a = dic[n] # Read Chebyshev coefficients
            else: # Compute Chebyshev coefficients:
                for k in range(4, n, 2):
                    c = cheb_f(lambda y: g(k+2, y))
                    d = dca(a, c)
                    a = cheb_a(d, k+2)
                    dic[k+2] = a
        dic[n] = a # Record Chebyshev coefficients
        q = qn(y, a)
    else: # Odd n > 3
        if n in dic.keys(): # Check if Chebyshev coefficients were computed
            a = dic[n] # Read Chebyshev coefficients
        else: # Compute Chebyshev coefficients:
            d = cheb_f(gq_3)
            a = cheb_a(d, 5)
            dic[5] = a
        if n > 5:
            if n in dic.keys(): # Check if Chebyshev coefficients were computed
                a = dic[n] # Read Chebyshev coefficients
            else: # Compute Chebyshev coefficients:
                for k in range(5, n, 2):
                    c = cheb_f(lambda y: g(k+2, y))
                    d = dca(a, c)
                    a = cheb_a(d, k + 2)
                    dic[k+2] = a
        dic[n] = a # Record Chebyshev coefficients
        q = qn(y, a)
    # Correct digits for double precision maximal roundoff error ~ n*1.e-14:
    sf = int(16-2-np.floor(np.log10(n)))
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
