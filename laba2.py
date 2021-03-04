# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 16:03:36 2021

@author: RIPLECS
"""
import numpy as np

def Gauss(AA):
    A=AA
    print(AA, ' = A')
    B=[]; X=[]
    x=['x1', 'x2', 'x3', 'x4', 'x5']
    while len(A)!=0:
        m=A[0][0]
        for i in range(len(A)):
            for j in range(len(A)):
                if abs(A[i][j])>abs(m): m=A[i][j]
                if A[i][j]==m: p=i; q=j
        print(f'\na_main = {m}, p={p}, q={q}\n')
        for i in range(len(A[0])): A[p][i]=A[p][i]/m
        print('Нормуємо головний рядок:\n')
        print(A, ' = A ')
        v=A[p].copy()
        A[p]=A[0] 
        A[0]=v
        for i in range(len(A)):
            v=A[i][q].copy()
            A[i][q]=A[i][0]
            A[i][0]=v
        r=x[q]
        x[q]=x[0]
        x[0]=r
        X.append(x)
        x=x[1:]
        M=len(A)*[0]
        for k in range(len(A)):
            if k==0: continue
            M[k]=-A[k][0]/A[0][0]
        for i in range(len(A)):
            for j in range(len(A[0])):
                if i==0: continue
                A[i][j]=A[i][j]+A[0][j]*M[i]
        B.insert(0, A[0])
        print('\nПересуваємо a_main у верхній лівий куток і утворюємо під ним нулі:\n')
        print(A, ' = A ')
        idx=np.where(A==A[0][0])
        A=np.delete(A, idx[0][0], axis=0)
        A=np.delete(A, idx[1][0], axis=1)
    k=0
    A=np.zeros((len(AA), len(AA)+1))
    for i in reversed(B):
        x=list(i[:(len(i)-1)])
        y=list(i[(len(i)-1):])
        A[k]=np.array(list((len(A[0])-len(i))*[0])+x+y).reshape(1, -1)
        k=k+1
    print('\nРезультат:\n')
    print(A, ' = A')
    print(X)
    res=[]
    n=len(A)
    for i in reversed(range(n)):
        result=A[i][-1]
        k=0
        for j in reversed((range(n-i))):
            r=A[i][n-j]
            if r!=A[i][-1]:
                result=result-r*res[k]
                k=k+1
        res.insert(0, result)
    B=res[::-1]
    mas=zip(reversed(X), B)
    xs=sorted(mas, key=lambda tup: tup[0])
    B=[b[1] for b in xs]
    print('res = ', B)
    return A, B
    
B=np.array([[1.0, 5.0, 3.0, -4.0, 20.0],
            [3.0, 1.0, -2.0, 0.0, 9.0],
            [5.0, -7.0, 0.0, 10.0, -9.0],
            [0.0, 3.0, -5.0, 0.0, 1.0]])

Gauss(B)

