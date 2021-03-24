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

def make_C(mat):
    c=np.zeros((len(mat), len(mat)))
    for i in range(len(c)):
        for j in range(len(c)):
            if i==j: continue
            else: c[i][j]=-mat[i][j]/mat[i][i]
    return c

def check_simetric(mat):
    for i in range(len(mat)):
        for j in range(len(mat)):
            if mat[i][j]==mat[j][i]:
                continue
            else: 
                return False
                break
    return True

def iter_Zeidel(mat, v, x0, e):
    true = True
    check=check_simetric(mat)
    iter=0
    while true:
        iter=iter+1
        x1 = np.copy(x0)
        for i in range(len(mat)):
            summ=0
            for j in range(0, i):
                if check is True: summ=summ-mat[i][j]*x1[j]
                else:summ=summ+mat[i][j]*x1[j]
            x1[i]=v[i]+summ
            summ=0
            for j in range(i+1, len(mat)):
                if check is True: summ=summ-mat[i][j]*x0[j]
                else:summ=summ+mat[i][j]*x0[j]
            x1[i]=x1[i]+summ
            if check is True: x1[i]=x1[i]/mat[i][i]
        true = max(abs(x1[i]-x0[i]) for i in range(len(mat)))>=e
        x0 = x1
        if check is True: r=abs(v-np.dot(mat, x0))
        else: r=abs(b1-np.dot(M, x0))
        print(f'{iter}) x = ', x0[::-1])
        gap=(len(str(iter))+2)*' '
        print(f'{gap}r = |A-b*x| = ', r)
    res=x0 if check is True else x0[::-1]
    print('\nВідповідь: x = ', res)
    print('\n           r = |A-b*x| = ', r)
    return res

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

C=make_C(M)
b2=[]
for i in range(len(M)):
    b2.append(b1[i]/M[i][i])


X = [0.001, 0.001, 0.001, 0.001]           
eps=0.0001
print('\nУтворимо матрицю С:\n')
print(C)
print('\nТоді вектор b:\n')
print(b2)
#for i in range(len(C)):
#    res=0
#    for j in range(len(C)):
#        res=res+C[i][j]
#    print(res)
#print('\n')

print('\nМетод Зейделя для матриці С:\n')
iter_Zeidel(C, b2, X, eps)

print('\nДомножимо рівняння Ax=b на А транспоновану і отримаємо наступні матриці:')
A=np.dot(np.transpose(AA), AA)
b=np.dot(np.transpose(AA), b)
print('\nA^T*A = \n', A)
print('\nA^T*b = \n', b)
print('\nМетод Зейделя для матриці A^T*A:\n')
iter_Zeidel(A, b, X, eps)
print("\nРозв'язок за допомогою програмного пакету Python:\n")
print(np.linalg.solve(AA, np.array([4.21, 8.97, 2.38, 12.98])))
