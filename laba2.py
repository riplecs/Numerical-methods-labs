# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 16:03:36 2021

@author: RIPLECS
"""
import numpy as np

def Gauss(AA):
    A = AA
    print(AA, ' = A')
    B = []
    while len(A) != 0:
        m = A[0][0]
        for i in range(len(A)):
            for j in range(len(A)):
                if abs(A[i][j]) > abs(m): 
                    m = A[i][j]
                if A[i][j] == m: 
                    p, q = i, j
        print(f'\na_main = {m}, p = {p}, q = {q}\n')
        for i in range(len(A[0])): 
            A[p][i] /= m
        print('Нормуємо головний рядок:\n')
        print(A, ' = A ')
        v = A[p].copy()
        A[p] = A[0] 
        A[0] = v
        for i in range(len(A)):
            v = A[i][q].copy()
            A[i][q] = A[i][0]
            A[i][0] = v
        M = len(A)*[0]
        for k in range(len(A)):
            if k == 0: 
                continue
            M[k] = -A[k][0]/A[0][0]
        for i in range(len(A)):
            for j in range(len(A[0])):
                if i == 0: 
                    continue
                A[i][j] += A[0][j]*M[i]
        B.insert(0, A[0])
        print('\nПересуваємо a_main у верхній лівий куток і утворюємо під ним нулі:\n')
        print(A, ' = A ')
        idx = np.where(A == A[0][0])
        A = np.delete(A, idx[0][0], axis = 0)
        A = np.delete(A, idx[1][0], axis = 1)
    k = 0
    A = np.zeros((len(AA), len(AA) + 1))
    for i in reversed(B):
        x = list(i[:(len(i) - 1)])
        y = list(i[(len(i) - 1):])
        A[k] = np.array(list((len(A[0]) - len(i))*[0]) + x + y).reshape(1, -1)
        k += 1
    print('\nРезультат:\n')
    print(A, ' = A')
    print(X)
    res = []
    n = len(A)
    for i in reversed(range(n)):
        result = A[i][-1]
        k = 0
        for j in reversed((range(n - i))):
            r = A[i][n - j]
            if r != A[i][-1]:
                result -= r*res[k]
                k += 1
        res.insert(0, result)
    print('\nx = ', res)
    return A, res
    
def checking(mas, x, v):
    return mas.dot(np.transpose((np.matrix(x)))) - np.transpose((np.matrix(v)))

def join(A, w):
    D = np.zeros((len(A), len(A) + 1))
    for i in range(len(A)):
        D[i] = np.append(A[i], w[i])
    return D

M = np.array([[3.81, 0.25, 1.28, 1.75],
              [2.25, 1.32, 5.58, 0.49],
              [5.31, 7.28, 0.98, 1.04],
              [10.39, 2.45, 3.35, 2.28]])

b = [4.21, 8.97, 2.38, 12.98] 


A = join(M, b)
r = Gauss(A)
print("\nРозв'язок за допомогою програмного пакету Python:\n")
print(np.linalg.solve(M, b))
print('\nAx-b = \n', checking(M, r[1], b))


