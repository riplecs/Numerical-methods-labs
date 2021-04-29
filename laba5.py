# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 22:22:06 2021

@author: RIPLECS
"""


import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.interpolate import CubicSpline

n=12
func = lambda x: 2*(np.sin(x/2))**2
X=np.linspace(-np.pi, np.pi, n)
Y=func(X)

df=pd.DataFrame({'x':[i for i in X], 'f(x)':[l for l in Y]})
print(df.set_index('x'))

x=np.linspace(-np.pi, np.pi)
plt.plot(x, func(x))
plt.scatter(X, Y)


#df = pd.DataFrame({'x':[-1.5, -1, -0.5, 0, 0.5, 1, 1.5], 'f(x)':[-0.7, 0, 0.7, 1, 0.7, 0, -0.7]})
#df = pd.DataFrame({'x':[1, 2, 3, 4, 5], 'f(x)':[5, 3, 2,]})

a = [i for i in df['f(x)'][:-1]]
y = [i for i in df['f(x)']]
h = [df['x'][i]-df['x'][i-1] for i in range(1, len(df))]

#h_1 = [i*h[i] for i in range(0, 4)] 
#h_2 = [(i-1)*h_1[i] for i in range(0, 4)]



def Thomas_method(mat, v):
    n=len(mat)
    p, q = [0]*n, [0]*n
    p[0] = -mat[0][1]/mat[0][0]
    q[0] = -(-v[0]/mat[0][0])
    for i in range(1, n-1):
        p[i] = mat[i][i+1]/(-mat[i][i]-mat[i+1][i]*p[i-1])
        q[i] = (mat[i][i-1]*q[i-1] - v[i])/(-mat[i][i] - mat[i][i-1]*p[i-1])
    x = [0]*n
    x[n-1] = (mat[-1][n-2]*q[n-2]-v[-1])/(-mat[-1][-1]-mat[-1][n-2]*p[n-2])
    for i in reversed(range(n-1)):
        x[i]=p[i]*x[i+1]+q[i]
    return x

def cubic_spline(f, r):
    A=np.eye((len(df)-1))
    n=len(A)
    '''
    for i in range(1, n-1):
       A[i][i] = 2*(r[i] + r[(i+1)])
    for i in range(1, n-1):
       A[i][i-1] = r[i]
    for i in range(1, n-1):
       A[i][i+1] = r[i]
    '''
    for i in range(2, n):
        A[i][i]=2*(h[i-1]+h[i])
        A[i][i-1]=h[i-1]
        A[i-1][i]=h[i]
    w = [3*((f[i]-f[i-1])/r[i]-(f[i-1]-f[i-2])/r[i-1]) for i in range(2, n)]
    w.append(0.0)
    w.insert(0, 0.0)
    c = Thomas_method(A, w)
    c[0]=0.0
    c[n-1]=0.0
    d = [(c[i+1] - c[i])/(3*r[i]) for i in range(n-1)]
    d.append(-c[n-1]/(3*r[n-1]))
    b = [(f[i+1]-f[i])/r[i] - (c[i+1]+2*c[i])*r[i]/3 for i in range(n-1)]
    b.append((f[n-1]-f[n-2])/r[n-1] - 2*c[n-1]*r[n-1]/3)
    return b, c, d

bb, cc, dd = cubic_spline(y, h)
cs=CubicSpline(df['x'], y)
print(cubic_spline(y, h))
'''  

print(a[0])
print(bb[0])
print(cc[0])
print(dd[0])

 
print(cs.c.item(3, 0))  
print(cs.c.item(2, 0))   
print(cs.c.item(1, 0))   
print(cs.c.item(0, 0))   
'''

def spline(dot):
    for i in range(len(df['x'])):
        if df['x'][i]>dot:
            x0=df['x'][i-1]
            j=i-1
            break
    f = lambda x: a[j]+bb[j]*(x-x0)+cc[j]*(x-x0)**2 + dd[j]*(x-x0)**3
    return f(dot)



X2=np.linspace(-np.pi+0.1, np.pi-0.4, n+1)
print(X2)
Y2=[spline(i) for i in X2]
plt.plot(X2, Y2)
X3=np.linspace(-np.pi, np.pi, 5*n)

ress=[]
xx=[]
for j in range(1, len(df['x'])):
    res=[abs(func(i)-spline(i)) for i in X3 if i<df['x'][j] and i>df['x'][j-1]]
    ress.append(max(res))
    for i in res:
        if i==max(res):
            xx.append(i)

        
plt.show()
x=np.linspace(-np.pi, np.pi, n-1)
plt.plot(x[:-1], ress[:-1])
plt.show()

'''
A = np.array([[2, 1, 0, 0],
             [2, 3, -1, 0],
             [0, 1, -1, 3],
             [0, 0, 1, -1]])

b = np.array([4, 9, 12, -4])

def Laugrange(x, y, k):
    res = lambda z: 1
    for i in range(len(x)):
        if i!=k:
            res = lambda z: res*(z-x[i])/(x[k]-x[i])
    return lambda z: res*y[k]

def Laugrange(x, y, k):
    #np.prod([(x-X[j])/(X[i]-X[j]) for j in range(len(X)) if i != j]) * Y[i]
    return lambda z: np.prod([(z-x[i])/(x[k]-x[i]) for i in range(len(x)) if i!=k])*y[k]
        
def total(x, y):
    summ=[Laugrange(x, y, i) for i in range(len(x))] 
    return lambda x: np.sum([s(x) for s in summ])

F=total(X, Y)
x_range = np.linspace(X[0], X[-1], 100)
plt.plot(X, Y, 'ro')
plt.plot(x_range, map(F, x_range))
plt.xlabel(r'$x$')
plt.ylabel(r'$F(x)$')
plt.title('Lagrange polynomial interpolation')
plt.grid(True)
plt.show()
'''