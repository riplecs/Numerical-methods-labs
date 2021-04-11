# -*- coding: utf-8 -*-
"""
Created on Sun Mar 28 22:48:24 2021

@author: RIPLECS
"""

import numpy as np
import math

A=np.array([[7.41, 1.13, 0.93, 1.22],
            [1.13, 8.31, 1.30, 0.16],
            [0.93, 1.30, 5.42, 2.10],
            [1.22, 0.16, 2.10, 11.10]])

def max_up_diag(mat):
    triangle=[]
    for i in range(len(mat)):
        for j in range(len(mat)):
            if i>j:
                triangle.append(mat[i][j])
    return max(triangle, key=abs)


def vectors(mat, k):
    return np.matrix([mat[i][k] for i in range(len(mat))])
eps=0.00001

def Jakobi(M):
    v=np.eye(len(M))
    M0=M.copy()
    main=max_up_diag(M0)
    while abs(main)>eps:
        Sd, Snd, SA = 0, 0, 0
        for i in range(len(M0)):
            for j in range(len(M0)):
                SA=SA+M0[i][j]**2
        print('SA = ', SA)
        for i in range(len(M0)):
            Sd=Sd+M0[i][i]
        print('Sd = ', Sd)
        for i in range(len(M0)):
            for j in range(len(M0)):
                if i!=j:
                    Snd=Snd+M0[i][j]**2
        print('Snd = ', Snd)
        main=max_up_diag(M0)
        for k in range(len(M0)):
            for l in range(len(M0)):
                if k>l:
                    if M0[k][l]==main:
                        i, j = k, l
                        break
        nu=2*main/(M0[i][i]-M0[j][j])
        H=np.eye(len(M0))
        c=np.sqrt(0.5*(1+1/np.sqrt(1+nu**2)))
        s=np.sign(nu)*np.sqrt(0.5*(1-1/np.sqrt(1+nu**2)))
        H[i][i]=H[j][j]=c
        H[i][j]=-s
        H[j][i]=s
        print(H, ' = T')
        v=np.dot(v, H)
        M1=((np.transpose(H)).dot(M0)).dot(H)
        main=max_up_diag(M1)
        M0=M1
    lambdas=[]
    for i in range(len(M0)):
        lambdas.append(M0[i][i])
    for i in range(len(lambdas)):
        print(f'λ_{i} = ', lambdas[i])
    print('Власні вектори:')
    for i in range(len(v)):
        print(np.transpose(vectors(v, i)), f' = v_{i+1}')


Jakobi(A)
print('By python: \n', np.linalg.eigh(A, UPLO='U'))