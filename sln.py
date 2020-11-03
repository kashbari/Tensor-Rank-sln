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

def e(i,j,n):
    return matrix(QQ,n,n,{(i,j):1})

def basis(i,n):
    assert i < n**2-1
    i,j = i//n, i%n
    if i==j:
        return e(i,i,n) - e(i+1,i+1,n)
    else:
        return e(i,j,n)

# e_ij
def Tsln(n):
    T = [{} for i in range(n**2-1)]
    B = matrix(QQ,[ basis(i,n).list() for i in range(n**2-1)])
    for i in range(n**2-1):
        a = basis(i,n)
        for j in range(n**2-1):
            b = basis(j,n)
            c = a*b - b*a
            c = B.solve_left(vector(QQ,c.list()))
            for k in c.nonzero_positions():
                T[i][(j,k)] = c[k]
    T = [matrix(QQ,n**2-1,n**2-1,m) for m in T]

    return T
print('Koszul Flattenings for sln')
n = input('n = ')
#print('n = %i' %n)
m = n**2-1
p = input('p = ')
#print('p = %i' %p)
d = int(nCr(m,p))
#print('Computing sln')


T = Tsln(n)
T = np.array([np.array(T[i]) for i in range(0,m)])

#print('T in np.array')

# S: B* \to A \otimes C
Tf = T.reshape(m,m*m)
#print(Tf.shape)
#print(np.linalg.matrix_rank(Tf))
# Id_Skew \otimes S: L^p A \otimes B^* \to L^p A \otimes A \otimes C
K = np.kron(np.eye(d), Tf)
#print(K.shape)
#Projection L^p \otimes A to L^{p+1}A
aa = a(int(m),int(p))
#Kronecker product with C
P = np.kron(aa.T,np.eye(m))
#print(P.shape)
TAp = np.dot(K,P)
print('TAp has dimensions and is of rank:')
print(TAp.shape)
print(np.linalg.matrix_rank(TAp))
#print('spaghetti2')



#RESTRICTION TO SUBSPACE
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


REST = raw_input('Restrict to subspace? (y/n) ')
if REST == 'y':
	k = input('Dimension of subspace = ')
	print('k = %i' %k)
	p=(k-1)/2
	d = int(nCr(k,p))
	S = res(T,k,m)
	# S: B* \to A \otimes C
	Sf = S.reshape(m,k*m)
	# Id_Skew \otimes S: L^p A \otimes B^* \to L^p A \otimes A \otimes C
	SK = np.kron(np.eye(d), Sf)
	#Projection L^p \otimes A to L^{p+1}A
	aa = a(int(k),int(p))
	#Kronecker product with C
	SP = np.kron(aa.T,np.eye(m))
	SAp = np.dot(SK,SP)
	print('RESTRICTION to %i dim subspace' %k)
	print(SAp.shape)
	print(np.linalg.matrix_rank(SAp))


print('DONE!')

