# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 21:42:47 2021

@author: RIPLECS
"""
import numpy as np

def join(A, w):
    D=np.zeros((len(A), len(A)+1))
    for i in range(len(A)):
        D[i]=np.append(A[i], w[i])
    return D

def change(mat, p, q):
    for i in range(len(mat)):
        v=mat[i][p].copy()
        mat[i][p]=mat[i][q]
        mat[i][q]=v
    return mat

def check_diagonal(mat):
    k=0
    for i in range(len(mat)):
        res=0
        for j in range(len(mat)):
            res=res+abs(mat[i][j])
        if abs(mat[i][i])>=res-abs(mat[i][i]): k=k+1
    if k==len(mat): print('\nМатриця А - з діагональною перевагою.')
    else: print('ERROR')

def split(mat):
    v=[]
    for i in range(len(mat)):
        v.append(mat[i][-1])
    D=np.zeros((len(mat), len(mat)))
    for i in range(len(mat)):
        D[i]=mat[i][:-1]
    return D, v

def iter_Zeidel(mat, v, x0, e):
    true = True
    iter=0
    while true:
        iter=iter+1
        x1 = np.copy(x0)
        for i in range(len(mat)):
            summ=0
            for j in range(0, i):
                summ=summ+mat[i][j]*x1[j]
            x1[i]=v[i]-summ
            summ=0
            for j in range(i+1, len(mat)):
                summ=summ+mat[i][j]*x0[j]
            x1[i]=x1[i]-summ
            x1[i]=x1[i]/mat[i][i]
        true = max(abs(x1[i]-x0[i]) for i in range(len(mat)))>e
        x0 = x1
        r=abs(v-np.dot(mat, x0))
        print(f'{iter}) x = ', x0)
        gap=(len(str(iter))+2)*' '
        print(f'{gap}r = |A-b*x| = ', r)
    print('\nВідповідь: x = ', x0)
    print('\n           r = |A-b*x| = ', r)
    return x0

AA=np.array([[3.81, 0.25, 1.28, 1.75],
            [2.25, 1.32, 5.58, 0.49],
            [5.31, 7.28, 0.98, 1.04],
            [10.39, 2.45, 3.35, 2.28]])

b=np.array([4.21, 8.97, 2.38, 12.98])

print(AA, ' = A\n')
print(b, ' = b\n')

A=join(AA, b)
print('Матриця до змін:\n\n', A)    
A=change(A, 1, 2)   
print('\nПоміняємо місцями 2й і 3й стовпці:\n\n', A)
A=change(A, 0, 3)
print('\nПоміняємо місцями 1й і 4й стовпці:\n\n', A)

for i in range(len(A)+1):
    A[0][i]=3*A[0][i]-A[-1][i]*1.1
    
print('\nДомножимо 1шу строку на 3 і віднімемо від неї останню, помножену на 1.1:\n\n', A)
for i in range(len(A)+1):
    A[2][i]=A[2][i]-A[0][i]
    
print('\nВід 3ї строки віднімемо 1шу:\n\n', A)

check_diagonal(A)  
        
    
M=split(A)[0]
b1=split(A)[1]

print('\nОстаточні результати: \n\n', M, ' = A\n\n b = ', b1)

print('\nA^T*A = \n', np.dot(np.transpose(AA), AA))
print('\nA^T*b = \n', np.dot(np.transpose(AA), b))

X = [0.001, 0.001, 0.001, 0.001]           
eps=0.0000001

#print('\n')
#print(iter_Zeidel(M, b1, X, eps))
print('\n')
iter_Zeidel(np.dot(np.transpose(AA), AA), np.dot(np.transpose(AA), b), X, eps)

