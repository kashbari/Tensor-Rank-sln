# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 23:11:06 2017
@author: Kash
"""

import numpy as np
import sympy as sp
import math
from itertools import *
from sympy import Matrix

def nCr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)
    
def a(x,y):
    a = np.zeros( (int(nCr(x,y+1)), int(nCr(x,y)) * x) )
    pdict = {s:i for i,s in enumerate(combinations(range(x),y))}
    pp1dict = {s:i for i,s in enumerate(combinations(range(x),y+1))}

    for S in combinations(range(x),y):
        for j in range(x):
            if j in S: continue
            Sp = tuple(sorted(S + (j,)))
            ix = pp1dict[Sp]
            jx = pdict[S] * x + j
            sgn = (-1) ** Sp.index(j)
            a[ix,jx] = sgn
    return a    
    
    
# COMPUTE STRUCTURE TENSOR OF SL3 LIE ALGEBRA
H1 = np.array([[1,0,0],[0,-1,0],[0,0,0]])
H2 = np.array([[0,0,0],[0,1,0],[0,0,-1]])
E12 = np.array([[0,1,0],[0,0,0],[0,0,0]])
E13 = np.array([[0,0,1],[0,0,0],[0,0,0]])
E23 = np.array([[0,0,0],[0,0,1],[0,0,0]])
E21 = E12.T
E31 = E13.T
E32 = E23.T

def c(x,y):
    c = np.dot(x,y) - np.dot(y,x)
    return c

#print(c(H1,H2)) # = 0
#print(c(H1,E12)) # = 2E12
#print(c(H1,E13)) # = E13
#print(c(H1,E23)) # = -E23
#print(c(H1,E21)) # = -2E21
#print(c(H1,E31)) # = -E31
#print(c(H1,E32)) # = E32

T = np.zeros((8,8,8))
T[0,2,2] = 2
T[0,3,3] = 1
T[0,4,4] = -1
T[0,5,5] = -2
T[0,6,6] = -1
T[0,7,7] = 1

T[2,0,2] = -2
T[3,0,3] = -1
T[4,0,4] = 1
T[5,0,5] = 2
T[6,0,6] = 1
T[7,0,7] = -1

#print(c(H2,E12)) # = -E12
#print(c(H2,E13)) # = E13
#print(c(H2,E23)) # = 2E23
#print(c(H2,E21)) # = E23
#print(c(H2,E31)) # = -E31
#print(c(H2,E32)) # = -2E32

T[1,2,2] = -1
T[1,3,3] = 1
T[1,4,4] = 2
T[1,5,5] = 1
T[1,6,6] = -1
T[1,7,7] = -2

T[2,1,2] = 1
T[3,1,3] = -1
T[4,1,4] = -2
T[5,1,5] = -1
T[6,1,6] = 1
T[7,1,7] = 2

#print(c(E21,E12)) # = -H1
#print(c(E23,E13)) # = 0
#print(c(E21,E13)) # = E23
#print(c(E31,E13)) # = -H1 - H2
#print(c(E32,E23)) # = -H2
#print(c(E13,E21)) # = -E23
#print(c(E31,E12)) # = E32
#print(c(E12,E23)) # = E13

T[5,2,0] = -1
T[6,3,0] = -1
T[6,3,1] = -1
T[2,5,0] = 1
T[7,4,1] = -1
T[4,7,1] = 1
T[3,6,0] = 1
T[3,6,1] = 1

T[5,3,4] = 1 #E21 E13
T[3,5,4] = -1
T[6,2,7] = 1 #E31 E12
T[2,6,7] = -1
T[2,4,3] = 1 #E12 E23
T[4,2,3] = -1
T[7,5,6] = 1 #E32 E21
T[5,7,6] = -1
T[3,7,2] = 1 #E13 E32
T[7,3,2] = -1
T[4,6,5] = 1 #E23 E31
T[6,4,5] = -1
print('spaghetti')
''' 
Basis:
0 = H1
1 = H2
2 = E12
3 = E13
4 = E23
5 = E21
6 = E31
7 = E32    
'''

# COMPUTE THE RANK OF THE KOSZUL FLATTENING
n=8
p=1
 
d = int(nCr(n,p))

# S: B* \to A \otimes C
Tf = T.reshape(n,n*n)
print(Tf.shape)
print(np.linalg.matrix_rank(Tf))
# Id_Skew \otimes S: L^p A \otimes B^* \to L^p A \otimes A \otimes C
K = np.kron(np.eye(d), Tf)
print(K.shape)
#Projection L^p \otimes A to L^{p+1}A    
aa = a(int(n),int(p))
#Kronecker product with C
P = np.kron(aa.T,np.eye(n))
print(P.shape)
TAp = np.dot(K,P)
print(TAp.shape)
print(np.linalg.matrix_rank(TAp))    
print('spaghetti2')
#RESTRICTION TO SUBSPACE -------------------- WIP
def res(T,k,n):
    R = np.random.rand(k,n)
    Q = np.zeros((n,k,n))
    for j in range(k):
        s=[]
        for i in range(n):
            s.append(R[j,i]*T[:,:,i])
        s1 = np.zeros(n)
        for l in range(n):
            s1 = s1+ s[l]
        s1 = np.asarray(s1)
        Q[:,j,:]= s1   
    return Q
    
k=3
p=(k-1)/2
d = int(nCr(k,p))    
S = res(T,k,n)
# S: B* \to A \otimes C
Sf = S.reshape(n,k*n)
# Id_Skew \otimes S: L^p A \otimes B^* \to L^p A \otimes A \otimes C
SK = np.kron(np.eye(d), Sf)
#Projection L^p \otimes A to L^{p+1}A    
aa = a(int(k),int(p))
#Kronecker product with C
SP = np.kron(aa.T,np.eye(n))
SAp = np.dot(SK,SP)
print(SAp.shape)
print(np.linalg.matrix_rank(SAp))    


'''
NOTES:
   sl3 structure tensor (8x8x8)
   p=3 (448,560) -> 440     dim Ker = 8
   p=2 (224,448) -> 223     dim Ker = 1
   p=1 (64,224) -> 64       FULL RANK!!!
   Restrict A (dim n) to subspace A' (dim k)
   k=7 p=3 (280,280) -> 273    dim Ker = 7
   k=5 p=2 (80,80) -> 76        dim Ker = 4
   k=3 p=1 (24,24) -> 24        FULL RANK!!!
    
'''

#def null(a, rtol=1e-5):
#    u, s, v = np.linalg.svd(a)
#    rank = (s > rtol*s[0]).sum()
#    return rank, v[rank:].T.copy()
#    
#M = Matrix(TAp)
#KerTAp = M.nullspace()
