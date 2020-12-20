# Code for constructing irreducible representations
# of semisimple lie algebras

# Distinguish primitive positive root vectors x_1, ..., x_r. These determine
# corresponding h_i and y_i so that [x_i,y_i] = h_i, [h_i, x_i] = 2x_i and
# [h_i,y_i] = -2y_i.

# Representations will be specified by giving the action of these three lists of
# elements. Explicitly, if V is a representation of L of dimension n, V will be
# described by [ [ X_1, ..., X_r], [Y_1, ..., Y_r], [H_1, ..., H_r] ], where
# the X_i, Y_i and H_i are n x n matrices.

# In addition to the representations in the above form, routines of this file
# also require the cartan matrix of L, with rows and columns ordered to
# correspond with the x_i.

# A_n = sl_n+1
# B_n = so_2n+1
# C_n = sp_2n
# D_n = so_2n

from itertools import *
import numpy as np
import collections


class ssLieAlgebra():

	def __init__(self,SLA):
		if SLA[0] not in {'A','B','C','D','E','F','G'}:
			raise NameError("Not a valid Simple Lie Algebra")
		try:
			self.rank = int(SLA[1:])
		except ValueError:
			raise NameError("Not a valid Simple Lie Algebra")
		self.name = SLA
		self.type = SLA[0]
	
	def cartan(self):
		rk = self.rank
		t = self.type
		C= np.zeros((rk+1,rk+1))
		for i in range(rk+1):
			C[i,i] = 2
			if i < (rk-1):
				C[i,i+1] = -1
				C[i+1,i] = -1
		if t == 'A':
			C[rk,rk-1] = -1
			C[rk-1,rk] = -1
		elif t == 'B':
			C[rk,rk-1] = -2
			C[rk-1,rk] = -1
		elif t == 'C':
			C[rk,rk-1] = -1
			C[rk-1,rk] = -2
		elif t == 'D':
			C[rk,rk-2] = -1
			C[rk-2,rk] = -1
		self.cartan = C
		return self.cartan

	def reps(self):
		R = [ [], [], [] ]
		self.reps = R
		self.X = self.reps[0]
		self.Y = self.reps[1]
		self.H = self.reps[2]
		return reps
	

	def structure_tensor(self,lie_alg):
		pass
			

class Representations():

	def __init__(self,sslie_alg,dim):
		self.dim = dim
		self.lie_alg = ssLieAlgebra(sslie_alg)
		
	
	def tensor_product(self,a):
		pass

	def dual(self):
		return [[ -x.transpose() for x in s] for s in self]


# Tensor product of modules
def module_product(a,b):
    return [[x.tensor_product(identity_matrix(QQ,y.dimensions()[0],sparse=True)) +
            identity_matrix(QQ,x.dimensions()[0],sparse=True).tensor_product(y)
        for x,y in zip(s1,s2)] for s1,s2 in zip(a,b)]

# Dual of module
def module_dual(a):
    return [[-x.transpose() for x in s] for s in a]

# Symmetrization of module (S^k A)
def module_sym(a,k):
    pass

# Skew symmetrization of module (Lambda^k A)
def module_skew(a,k):
    pass

# B is a basis of column vectors for the desired submodule
# being a submodule is not checked
def submodule_from_basis(a,B):
    return [[restrict_map(m,B,B) for m in s] for s in a]

# B: n x k matrix
# C: m x l matrix
# m: m x n matrix mapping the column space of B to the column space of C
# returns the l x k matrix representing m in the bases B and C
def restrict_map(m,C,B):
    I = C.pivot_rows()
    return C[I,:].solve_right((m*B)[I,:])

# a : a representation of L
# C : the cartan matrix of the semisimple part of L
#
# computes a dictionary summarizing info of this representation. In particular,
# the keys of the returned dictionary are
# 'M': a dictionary mapping weights to a distinguished basis of the corresponding weight space
# 'C': The Cartan matrix of L
# 'a': the representation
# 'tot': a distinguished total ordering of the weights compatible with P; for
#     convenience weights are given along with the multiplicity of the weight space
# 'wtg': directed graph with a vertex for each weight and an edge for each
#     raising operator. The edges are labelled with the raising operator between the
#     corresponding weight spaces in the distinguished bases of those spaces
#
def weight_decomposition(a,C):
    tot = simultaneous_eigenspace(a[2])

    totkey = lambda wt: tuple(C.solve_right(vector(QQ,wt)).list())
    tot.sort(key=lambda p: totkey(p[1]))

    M = {wt:V for V,wt in tot}
    tot = [(wt,V.dimensions()[1]) for V,wt in tot]

    wtg = DiGraph()
    for p in tot:
        wt, mult = p
        B = M[wt]
        for i,x in enumerate(a[0]):
            wt2 = tuple((C[:,i] + vector(wt).column()).list())
            if wt2 in M:
                wtg.add_edge(p,(wt2,M[wt2].dimensions()[1]),restrict_map(x,M[wt2],B))

    return { 'M' : M, 'tot' : tot, 'C': C, 'a': a, 'wtg': wtg }



# ms : list of matrices over QQ which commute and are diagonalizable
#
# computes the simultaneous eigenspaces of ms
def simultaneous_eigenspace(ms):
    n = ms[0].dimensions()[0]
    spaces = [(identity_matrix(n,sparse=True) ,())]
    for m in ms:
        nspaces = []
        for B,wt in spaces:
            for wt0,W in restrict_map(m,B,B).eigenspaces_right():
                W = W.basis_matrix().transpose().sparse_matrix()
                nspaces.append((B*W,wt+(wt0,)))
        spaces = nspaces
    return spaces




def weight_to_reps(wt):
    pass


