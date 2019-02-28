# -*- coding: utf-8 -*-
"""
Created on Sat Apr 14 21:56:01 2018

@author: kashb
"""

from datetime import datetime
import matplotlib.pylab as plt #for spy sparse matry
import sys #for exiting script early
import numpy as np
#import sympy as sp
import scipy as scipy
from scipy import sparse
from scipy.sparse.linalg import eigs
from numpy import linalg as LA
import math
from itertools import *
#from sympy import Matrix
#from scipy.sparse.linalg import inv


#// division is better for integers
def nCr(n,r):
    f = math.factorial
    return f(n) // f(r) // f(n-r)
    
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
    
    
''' COMPUTE STRUCTURE TENSOR OF SL4 LIE ALGEBRA '''
n=4
m=(n**2)-1

def ind(x):
    a = int(np.floor(x/n))
    b= x % n
    c = [a,b]
    return c



def iind(x):
    if x == 0:
        c = -1
    else:
        c = n*x[0] +x[1]
    return c



''' Computes Lie bracket E = [m,n] a = [i,j] '''
def LB(E,a):
    c1 = 0
    c2 = 0
    if a[0] ==E[1]:
        c1 = [E[0],a[1]]
        if a[1] == E[0]:
            c2 = [a[0],E[1]]
        else:
            c2 = 0
    else:
        c1 = 0
        if a[1] == E[0]:
            c2 = [a[0],E[1]]
        else:
            c2 = 0
    c = [c1,c2]
    return c

def LB1(E,a):
    [c1,c2] = LB(E,a)
    if c1 == 0 and c2 == 0:
        d = print(0)    
    elif c1 == 0 and c2 != 0:
        d = print("-",c2, sep="")
    elif c2 == 0 and c1 != 0:
        d = print(c1)
    else:
        d = print(c1,c2,sep="-")
    return d


def LBD(E,a):
    if a[0] == a[1]:
        c = LB(E,a)
        a1 = [a[0]+1,a[1]+1]
        d = LB(E,a1)
        c = c+ [d[1],d[0]]
    else:
        c=LB(E,a)
    return c



''' Computes Lie bracket E = [m,n] a = [a1,a2,...] = [i,j,k,l,...], returns + - + - '''
def replace(q,v,k):
    q1 = q
    q1[k] = v[0]
    q1[k+1] = v[1]
    return q1

def LieBracketTensor(E,a):
    b = []
    l = len(a)
    l2 = int(l/2)
    for k in range(0,l2):
        a1 = [a[2*k],a[2*k+1]]
        L = LBD(E,a1)
        for i in range(0,len(L)):
            if L[i] != 0:
                q = replace(a,L[i],2*k)
                b = b + [q]
            else:
                b = b + [0]
    return b

def ComputeLieTensor(x,y):
    x1 = ind(x)
    y1 = ind(y)
    z1 = LieBracketTensor(x1,y1)
    z = [0]*len(z1)
    for i in range(len(z1)):
        z[i] = iind(z1[i])
    return z

''' Actually construct structure tensor '''
print('Calculating Tensor')
Tensor = np.zeros((m,m,m))
for i in range(m):
    for j in range(m):
        z = ComputeLieTensor(i,j)
        for k in range(len(z)):
            if z[k] != -1 and z[k] < m:
                Tensor[i,j,z[k]] += (-1)**k
                
print('     Tensor memory: ', Tensor.data.nbytes, 'bytes')


# COMPUTE THE RANK OF THE KOSZUL FLATTENING
p = 7
print('p = ',p,',',' m = ', m)

d = int(nCr(m,p))
# S: B* \to A \otimes C
#FlattenedTensor = Tensor.reshape(m,m*m)
print('Calculating FlattenedTensor')
FlattenedTensor = sparse.csr_matrix(Tensor.reshape(m,m*m), dtype=np.int32)
del Tensor
print('     FlattenedTensor shape: ', FlattenedTensor.shape)
print('     FlattenedTensor memory: ', FlattenedTensor.data.nbytes, 'bytes')

#plt.spy(FlattenedTensor, markersize = 2, aspect = 'auto')
#plt.show()

# Id_Skew \otimes S: L^p A \otimes B^* \to L^p A \otimes A \otimes C
print('Calculating K')
K = sparse.kron(np.eye(d), FlattenedTensor, 'csr')
print('     K shape: ', K.shape)
print('     K memory: ', K.data.nbytes, 'bytes')
del FlattenedTensor #delete objects from memory

#Projection L^p \otimes A to L^{p+1}A    
aa = a(int(m),int(p))
aa = sparse.csr_matrix(aa, dtype=np.int32)

#Kronecker product with C
print('Calculating P')
P = sparse.kron(aa.T,np.eye(m),'csr')
print('     P shape: ', P.shape)
print('     P memory: ', P.data.nbytes, 'bytes')
del aa #delete objects from memory

print('Calculating TAp')
TAp = sparse.csr_matrix(K @ P, dtype=np.int32)

del P, K #delete objects from memory

print('     TAp shape: ', TAp.shape)
print('     TAp memory: ', TAp.data.nbytes, 'bytes')
#print('TAp condition number: ', LA.cond(TAp))

#plt.spy(TAp, markersize = 2)
#plt.show()
sys.exit()

"""Various attempts to compute the rank"""
def rank(A, eps=1e-12):
    u, s, vh = LA.svd(A)
    return len([x for x in s if abs(x) > eps])

#import scipy.linalg.interpolative as sli
#import scipy.sparse.linalg.interpolative as spl
def TApMatrixVectorProduct(v):
    return TAp@v


startTime = datetime.now()
print(scipy.sparse.linalg.svds(TAp, k=6, which='SM'))

#print(rank(TAp))

#print(sli.estimate_rank(spl.LinearOperator((TAp.shape[0], TAp.shape[1]),matvec = TApMatrixVectorProduct), 1e-8))
print(datetime.now() - startTime)
"""===================================="""
#sys.exit()

"""
print('Now the annoyingly long wait...')
startTime = datetime.now()
#GIMME RANK OF SPARSE MATRIX print(np.linalg.matrix_rank(TAp))
vals, vecs = eigs(TAp, k=1, which = 'SM', tol = 0.00000001)
print(vals)
print(datetime.now() - startTime)
"""





 