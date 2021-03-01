# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 16:03:36 2021

@author: RIPLECS
"""
import numpy as np

def Gauss(AA):
    A=AA
    print('A = ', AA)
    strs=[]
    while len(A)!=1:
        m=A[0][0]
        for i in range(len(A)):
            for j in range(len(A)):
                if abs(A[i][j])>abs(m):  m=A[i][j]
        for i in range(len(A)):
            for j in range(len(A)):
                if A[i][j]==m: p=i; q=j
        M=len(A)*[0]
        for k in range(len(A)):
            if k==p: continue
            M[k]=-A[k][q]/m
        for i in range(len(A)):
            for j in range(len(A[0])):
                if i==p: continue
                A[i][j]=A[i][j]+A[p][j]*M[i]
        strs.append(A[p])
        idx=np.where(A==A[p][q])
        A=np.delete(A, idx[0][0], axis=0)
        A=np.delete(A, idx[1][0], axis=1)
        print('A = ', A)
    strs.append(A[0])
    k=0
    B=np.zeros((len(AA), len(AA)+1))
    for i in reversed(strs):
        x=list(i[:(len(i)-1)])
        y=list(i[(len(i)-1):])
        B[k]=np.array(x+list((len(B[0])-len(i))*[0])+y).reshape(1, -1)
        k=k+1
    print('B = ', B)
    x1=B[0][-1]/B[0][0]
    x2=(B[1][-1]-x1*B[1][0])/B[1][1]
    x3=(B[2][-1]-x1*B[2][0]-x2*B[2][1])/B[2][2]
    x4=(B[3][-1]-x1*B[3][0]-x2*B[3][1]-x3*B[3][2])/B[3][3]
    print('x1 = ', x1)
    print('x2 = ', x2)
    print('x3 = ', x3)
    print('x4 = ', x4)
    return A
    
B=np.array([[1.0, 5.0, 3.0, -4.0, 20.0],
            [3.0, 1.0, -2.0, 0.0, 9.0],
            [5.0, -7.0, 0.0, 10.0, -9.0],
            [0.0, 3.0, -5.0, 0.0, 1.0]])

Gauss(B)

